use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};
use std::mem;
use std::path::Path;

use byteorder::{ByteOrder, LittleEndian};
use flate2::read::GzDecoder;

use flate2::write::GzEncoder;
use flate2::Compression;

use sprs::CsMat;

pub type MatValT = f32;

fn get_reserved_spaces(
    num_bit_vecs: usize,
    num_rows: usize,
    mut file: GzDecoder<BufReader<File>>,
) -> Result<Vec<usize>, Box<dyn Error>> {
    let mut bit_vec_lengths: Vec<usize> = vec![0; num_rows + 1];
    let mut bit_vec = vec![0; num_bit_vecs];
    let mut running_sum = 0;

    for i in 0..num_rows {
        file.read_exact(&mut bit_vec[..])?;
        let mut num_ones = 0;
        for bits in bit_vec.iter() {
            num_ones += bits.count_ones() as usize;
        }

        running_sum += num_ones;
        bit_vec_lengths[i + 1] = running_sum;

        // no seek command yet
        // copied from https://github.com/rust-lang/rust/issues/53294#issue-349837288
        io::copy(
            &mut file
                .by_ref()
                .take((num_ones * mem::size_of::<MatValT>()) as u64),
            &mut io::sink(),
        )?;
        //file.seek(SeekFrom::Current((num_ones * 4) as i64))?;
    }

    Ok(bit_vec_lengths)
}

// reads the EDS format single cell matrix from the given path
pub fn reader(
    file_path: &Path,
    num_rows: usize,
    num_cols: usize,
) -> Result<CsMat<MatValT>, Box<dyn Error>> {
    // reading the matrix
    let file_handle = File::open(file_path.to_owned())?;
    let buffered = BufReader::new(file_handle);
    let file = GzDecoder::new(buffered);

    let num_bit_vecs: usize = (num_cols + 7) / 8;
    let bit_vector_lengths = get_reserved_spaces(num_bit_vecs, num_rows, file)?;

    let total_nnz = bit_vector_lengths[num_rows];
    let mut data: Vec<MatValT> = vec![0.0; total_nnz];
    let mut indices: Vec<usize> = vec![0; total_nnz];

    let file_handle = File::open(file_path)?;
    let buffered = BufReader::new(file_handle);
    let mut file = GzDecoder::new(buffered);

    let mut global_pointer = 0;
    let mut bit_vec = vec![0; num_bit_vecs];
    for i in 0..num_rows {
        file.read_exact(&mut bit_vec[..])?;
        let num_ones = bit_vector_lengths[i + 1] - bit_vector_lengths[i];

        let mut one_validator = 0;
        for (j, flag) in bit_vec.iter().enumerate() {
            if *flag != 0 {
                for (i, bit_id) in format!("{:8b}", flag).chars().enumerate() {
                    if let '1' = bit_id {
                        let offset = i + (8 * j);
                        indices[global_pointer + one_validator] = offset;

                        one_validator += 1;
                    };
                }
            }
        }
        assert_eq!(num_ones, one_validator);

        let mut expression: Vec<u8> = vec![0; mem::size_of::<MatValT>() * (num_ones as usize)];
        let mut float_buffer: Vec<MatValT> = vec![0.0; num_ones as usize];
        file.read_exact(&mut expression[..])?;

        // NOTE: if we change MatValT, double check below line
        LittleEndian::read_f32_into(&expression, &mut float_buffer);

        for value in float_buffer {
            data[global_pointer] = value as MatValT;
            global_pointer += 1;
        }
    }

    assert_eq!(global_pointer, total_nnz);
    let matrix = CsMat::new((num_rows, num_cols), bit_vector_lengths, indices, data);
    Ok(matrix)
}

// writes the EDS format single cell matrix into the given path
pub fn writer(file_path: &Path, matrix: &CsMat<MatValT>) -> Result<(), Box<dyn Error>> {
    let file_handle = File::create(file_path)?;
    let buffered = BufWriter::new(file_handle);
    let mut file = GzEncoder::new(buffered, Compression::default());

    let num_bit_vecs: usize = (matrix.cols() + 7) / 8;
    let mut bit_vecs: Vec<u8> = vec![0; num_bit_vecs];

    for row_vec in matrix.outer_iterator() {
        let mut positions = Vec::new();
        let mut values = Vec::new();

        for (col_ind, &val) in row_vec.iter() {
            positions.push(col_ind);
            values.push(val as MatValT);
        }

        // clearing old bit vector
        bit_vecs.iter_mut().for_each(|x| *x = 0);

        // refilling bit vector
        for pos in positions {
            let i = pos / 8;
            let j = pos % 8;

            bit_vecs[i] |= 128u8 >> j;
        }

        let mut bin_exp: Vec<u8> = vec![0_u8; values.len() * 4];

        // NOTE: if we change MatValT, double check below line
        LittleEndian::write_f32_into(&values, &mut bin_exp);
        file.write_all(&bit_vecs)?;
        file.write_all(&bin_exp)?;
    }

    Ok(())
}

// writes the EDS format single cell matrix into the given path
pub fn as_bytes(matrix_row: &[MatValT], num_cols: usize) -> Result<Vec<u8>, Box<dyn Error>> {
    let num_bit_vecs: usize = (num_cols + 7) / 8;
    let mut bit_vecs: Vec<u8> = vec![0_u8; num_bit_vecs];
    let mut values = Vec::new();

    let error = f32::EPSILON;
    for (col_ind, &val) in matrix_row.iter().enumerate() {
        let val = val as MatValT;

        if (val - 0.0).abs() < error {
            continue;
        }
        values.push(val);

        let i = col_ind / 8;
        let j = col_ind % 8;
        bit_vecs[i] |= 128u8 >> j;
    }

    let mut bin_exp: Vec<u8> = vec![0_u8; values.len() * 4];
    if std::any::TypeId::of::<MatValT>() == std::any::TypeId::of::<f32>() {
        LittleEndian::write_f32_into(&values, &mut bin_exp);
        bit_vecs.append(&mut bin_exp);
    } else {
        unreachable!("EDS not supported for non float objects");
    }

    Ok(bit_vecs)
}

#[cfg(test)]
mod tests {
    use crate::single_cell::single_cell_experiment::SingleCellExperiment;
    use sprs::CsMat;

    fn get_test_matrix_f32() -> CsMat<f32> {
        CsMat::new(
            (3, 3),
            vec![0, 2, 4, 5],
            vec![0, 1, 0, 2, 2],
            vec![1.2, 2.3, 3.4, 4.5, 5.6],
        )
    }

    fn get_test_sce_f32_data() -> (CsMat<f32>, Vec<String>, Vec<String>) {
        let a = get_test_matrix_f32();
        let b: Vec<String> = vec!["1", "2", "3"]
            .into_iter()
            .map(|x| x.to_string())
            .collect();
        let c: Vec<String> = vec!["4", "5", "6"]
            .into_iter()
            .map(|x| x.to_string())
            .collect();

        (a, b, c)
    }

    fn get_temp_file(ext: String) -> std::path::PathBuf {
        let mut dir = std::env::temp_dir();
        dir.push(format!("foo{}", ext));
        dir
    }

    #[test]
    fn test_eds_f32() {
        let (a, b, c) = get_test_sce_f32_data();
        let sce = match SingleCellExperiment::from_csr(a, b.clone(), c.clone()) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        let file = get_temp_file(".eds.gz".to_owned());
        let fname = file.to_str().unwrap();
        match sce.to_eds(fname) {
            Ok(_) => (),
            Err(y) => panic!("ERROR: {}", y),
        };
        println!("{:?}", fname);

        let sce_eds: SingleCellExperiment<f32> = match SingleCellExperiment::from_eds(fname, b, c) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        assert_eq!(sce, sce_eds);
        std::fs::remove_file(fname).expect("can't remove temp file");
    }
}
