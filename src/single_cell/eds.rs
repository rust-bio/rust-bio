use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};
use std::mem;
use std::path::Path;

use byteorder::{ByteOrder, LittleEndian};
use flate2::read::GzDecoder;
use math::round;

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

    let num_bit_vecs: usize = round::ceil(num_cols as f64 / 8.0, 0) as usize;
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

    let num_bit_vecs: usize = round::ceil(matrix.cols() as f64 / 8.0, 0) as usize;
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
            let i = round::floor(pos as f64 / 8.0, 0) as usize;
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
