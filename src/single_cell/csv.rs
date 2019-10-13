use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::error::Error;
use std::fs::{canonicalize, File};
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::PathBuf;

use crate::single_cell::scmatrix::{MatValT, ScMatrix};
use sprs::CsMatBase;

fn get_file_names(path: &str) -> Result<(PathBuf, PathBuf, PathBuf), Box<dyn Error>> {
    // futuristic customization of alevin output name
    let mtx_type = "alevin";

    let quants_mat: PathBuf;
    let mut quants_mat_rows: PathBuf;
    let mut quants_mat_cols: PathBuf;
    {
        // assigning names of the file to read
        match mtx_type {
            "alevin" => {
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

// reads the CSV format single cell matrix from the given path
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
        let ballpark: usize = ((num_rows * num_columns) as f64 / 8.0) as usize;

        let file_handle = File::open(quants_mat.clone())?;
        let file = GzDecoder::new(file_handle);
        let buffered = BufReader::new(file);

        let mut ind_ptr = vec![0; num_rows + 1];
        let mut data = Vec::with_capacity(ballpark);
        let mut indices = Vec::with_capacity(ballpark);

        let mut total_entries = 0;
        for (row_id, line) in buffered.lines().enumerate() {
            let record = line?;
            let values: Vec<MatValT> = record.split(",").flat_map(str::parse::<MatValT>).collect();

            total_entries += values.len();
            ind_ptr[row_id + 1] = total_entries;

            assert_eq!(values.len(), num_columns);
            for (column_id, val) in values.into_iter().enumerate() {
                indices.push(column_id);
                data.push(val as MatValT);
            }
        }

        CsMatBase::new((num_rows, num_columns), ind_ptr, indices, data)
    };

    Ok(ScMatrix::new(matrix, row_names, column_names))
}

// writes the CSV format single cell matrix into the given path
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

        let num_columns = matrix.num_columns();
        for row_vec in matrix.data().outer_iterator() {
            let mut it = row_vec.iter();
            let mut column_id_validator = 0;

            match it.next() {
                Some((0, &val)) => {
                    column_id_validator += 1;
                    write!(&mut file, "{}", val)?;
                }
                Some((idx, &val)) => {
                    for _ in 0..idx {
                        write!(&mut file, "{},", 0 as MatValT)?;
                    }

                    column_id_validator += idx + 1;
                    write!(&mut file, "{}", val)?;
                }
                None => {
                    for _ in 1..num_columns {
                        write!(&mut file, "{},", 0 as MatValT)?;
                    }

                    write!(&mut file, "{}", 0 as MatValT)?;
                    continue;
                }
            };

            while let Some((column_ind, &val)) = it.next() {
                while column_id_validator < column_ind {
                    write!(&mut file, ",{}", 0 as MatValT)?;
                    column_id_validator += 1;
                }
                write!(&mut file, ",{}", val)?;
                column_id_validator += 1;
            }

            while column_id_validator < num_columns {
                write!(&mut file, ",{}", 0 as MatValT)?;
                column_id_validator += 1;
            }
        } // end row iterator
    }

    Ok(())
}
