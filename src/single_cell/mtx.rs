use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, BufWriter, Write};

use sprs::CsMatBase;
use std::error::Error;
use std::fs::{canonicalize, File};
use std::path::PathBuf;

use crate::single_cell::scmatrix::{MatIdxT, MatValT, ScMatrix};

fn get_file_names(path: &str) -> Result<(PathBuf, PathBuf, PathBuf), Box<dyn Error>> {
    // futuristic customization of alevin output name
    let mtx_type = "tenx";

    let quants_mat: PathBuf;
    let mut quants_mat_rows: PathBuf;
    let mut quants_mat_cols: PathBuf;
    {
        // assigning names of the file to read
        match mtx_type {
            "tenx" => {
                quants_mat = canonicalize(path)?;

                quants_mat_rows = quants_mat.clone();
                quants_mat_rows.set_file_name("barcodes.tsv");

                quants_mat_cols = quants_mat.clone();
                quants_mat_cols.set_file_name("features.tsv");
            }
            _ => unreachable!(),
        };
    }

    Ok((quants_mat, quants_mat_rows, quants_mat_cols))
}

fn triplets_to_sparse_matrix(
    triplets: Vec<BTreeMap<usize, MatValT>>,
    num_columns: usize,
) -> Result<CsMatBase<MatValT, MatIdxT, Vec<MatIdxT>, Vec<MatIdxT>, Vec<MatValT>>, Box<dyn Error>> {
    let num_rows = triplets.len();
    let ballpark: usize = ((num_rows * num_columns) as f64 / 8.0) as usize;

    let mut ind_ptr = vec![0; num_rows + 1];
    let mut data = Vec::with_capacity(ballpark);
    let mut indices = Vec::with_capacity(ballpark);

    let mut total_entries = 0;
    for (row_id, column_data) in triplets.into_iter().enumerate() {
        total_entries += column_data.len();
        ind_ptr[row_id + 1] = total_entries;

        for (column_id, val) in column_data {
            indices.push(column_id);
            data.push(val as MatValT);
        }
    }

    Ok(CsMatBase::new(
        (num_rows, num_columns),
        ind_ptr,
        indices,
        data,
    ))
}

// reads the MTX format single cell matrix from the given path
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

    // Todo: make this customizable
    let row_indexed = true;
    let (first_index, second_index) = match row_indexed {
        true => (0, 1),
        false => (1, 0),
    };

    let matrix = {
        // reading the matrix
        let num_rows = row_names.len();
        let num_columns = column_names.len();

        let file_handle = File::open(quants_mat)?;
        let file = GzDecoder::new(file_handle);
        let buffered = BufReader::new(file);

        let mut header = true;
        let mut triplets: Vec<BTreeMap<usize, MatValT>> = vec![BTreeMap::new(); num_rows];

        for line in buffered.lines() {
            let record = line?;
            if record.chars().nth(0).unwrap() == '%' {
                continue;
            }

            let vals: Vec<&str> = record.split_whitespace().collect();

            let row_id = vals[first_index]
                .parse::<usize>()
                .expect("can't convert row id");

            let columns_id = vals[second_index]
                .parse::<usize>()
                .expect("can't convert column id");

            let value = vals[2].parse::<MatValT>().expect("can't convert value");

            if header {
                header = false;

                assert!(num_rows == row_id);
                assert!(num_columns == columns_id);
                continue;
            }

            triplets[row_id - 1].insert(columns_id - 1, value);
        }

        triplets_to_sparse_matrix(triplets, num_columns)?
    };

    Ok(ScMatrix::new(matrix, row_names, column_names))
}

// writes the MTX format single cell matrix into the given path
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

        let header = "%%MatrixMarket\tmatrix\tcoordinate\treal\tgeneral".to_string();
        write!(
            &mut file,
            "{}\n{}\t{}\t{}\n",
            header,
            matrix.num_rows(),
            matrix.num_columns(),
            matrix.nnz()
        )?;

        for (row_ind, row_vec) in matrix.data().outer_iterator().enumerate() {
            for (column_ind, &val) in row_vec.iter() {
                write!(&mut file, "{}\t{}\t{}\n", row_ind + 1, column_ind + 1, val)?;
            }
        }
    }

    Ok(())
}
