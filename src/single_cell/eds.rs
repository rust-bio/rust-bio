use std::error::Error;
use std::fs::{canonicalize, File};
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

use byteorder::{ByteOrder, LittleEndian};
use flate2::read::GzDecoder;
use math::round;

use flate2::write::GzEncoder;
use flate2::Compression;

use crate::single_cell::scmatrix::{MatValT, ScMatrix};
use sprs::CsMatBase;

fn get_file_names(path: &str) -> Result<(PathBuf, PathBuf, PathBuf), Box<dyn Error>> {
    // futuristic customization of alevin output name
    let eds_type = 0_14;
    let quants_mat: PathBuf;
    let mut quants_mat_rows: PathBuf;
    let mut quants_mat_cols: PathBuf;
    {
        // assigning names of the file to read
        match eds_type {
            0_14 => {
                quants_mat = canonicalize(path)?;

                quants_mat_rows = quants_mat.clone();
                quants_mat_rows.set_file_name("quants_mat_rows.txt");

                quants_mat_cols = quants_mat.clone();
                quants_mat_cols.set_file_name("quants_mat_cols.txt");
            }
            _ => unreachable!(),
        };
    }

    Ok((quants_mat, quants_mat_rows, quants_mat_cols))
}

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
            &mut file.by_ref().take((num_ones * 4) as u64),
            &mut io::sink(),
        )?;
        //file.seek(SeekFrom::Current((num_ones * 4) as i64))?;
    }

    Ok(bit_vec_lengths)
}

// reads the EDS format single cell matrix from the given path
pub fn reader(input_path: &str) -> Result<ScMatrix, Box<dyn Error>> {
    // extracting the path of the files
    let (quants_mat, quants_mat_rows, quants_mat_cols) = get_file_names(input_path)?;

    let mut row_names: Vec<String> = Vec::with_capacity(1_000);
    {
        // reading rows file
        let file_handle = File::open(quants_mat_rows)?;
        let buffered = BufReader::new(file_handle);
        for line in buffered.lines() {
            row_names.push(line?);
        }
    }

    let mut column_names: Vec<String> = Vec::with_capacity(1_000);
    {
        // reading columns file
        let file_handle = File::open(quants_mat_cols)?;
        let buffered = BufReader::new(file_handle);
        for line in buffered.lines() {
            column_names.push(line?);
        }
    }

    let matrix = {
        // reading the matrix
        let num_rows = row_names.len();
        let num_columns = column_names.len();

        let file_handle = File::open(quants_mat.clone())?;
        let buffered = BufReader::new(file_handle);
        let file = GzDecoder::new(buffered);

        let num_bit_vecs: usize = round::ceil(num_columns as f64 / 8.0, 0) as usize;
        let bit_vector_lengths = get_reserved_spaces(num_bit_vecs, num_rows, file)?;

        let total_nnz = bit_vector_lengths[num_rows];
        let mut data: Vec<MatValT> = vec![0.0; total_nnz];
        let mut indices: Vec<usize> = vec![0; total_nnz];

        let file_handle = File::open(quants_mat)?;
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
                        match bit_id {
                            '1' => {
                                let offset = i + (8 * j);
                                indices[global_pointer + one_validator] = offset;

                                one_validator += 1;
                            }
                            _ => (),
                        };
                    }
                }
            }
            assert_eq!(num_ones, one_validator);

            let mut expression: Vec<u8> = vec![0; 4 * (num_ones as usize)];
            let mut float_buffer: Vec<f32> = vec![0.0_f32; num_ones as usize];
            file.read_exact(&mut expression[..])?;
            LittleEndian::read_f32_into(&expression, &mut float_buffer);

            for i in 0..float_buffer.len() {
                data[global_pointer] = float_buffer[i] as MatValT;
                global_pointer += 1;
            }
        }

        assert_eq!(global_pointer, total_nnz);
        CsMatBase::new((num_rows, num_columns), bit_vector_lengths, indices, data)
    };

    Ok(ScMatrix::new(matrix, row_names, column_names))
}

// writes the EDS format single cell matrix into the given path
pub fn writer(matrix: ScMatrix, path_str: &str) -> Result<(), Box<dyn Error>> {
    let quants_file_handle = File::create(path_str)?;
    let (_, quants_mat_rows, quants_mat_cols) = get_file_names(path_str)?;

    {
        // writing row file
        let file_handle = File::create(quants_mat_rows)?;
        let mut file = BufWriter::new(file_handle);

        let row_names = matrix.row_names();
        for name in row_names {
            write!(&mut file, "{}\n", name)?;
        }
    }

    {
        // writing columns file
        let file_handle = File::create(quants_mat_cols)?;
        let mut file = BufWriter::new(file_handle);

        let column_names = matrix.column_names();
        for name in column_names {
            write!(&mut file, "{}\n", name)?;
        }
    }

    {
        // writing matrix
        let buffered = BufWriter::new(quants_file_handle);
        let mut file = GzEncoder::new(buffered, Compression::default());

        let num_bit_vecs: usize = round::ceil(matrix.num_columns() as f64 / 8.0, 0) as usize;
        let mut bit_vecs: Vec<u8> = vec![0; num_bit_vecs];

        for row_vec in matrix.data().outer_iterator() {
            let mut positions = Vec::new();
            let mut values = Vec::new();

            for (col_ind, &val) in row_vec.iter() {
                positions.push(col_ind);
                values.push(val as f32);
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
            LittleEndian::write_f32_into(&values, &mut bin_exp);
            file.write_all(&bit_vecs)?;
            file.write_all(&bin_exp)?;
        }
    }

    Ok(())
}
