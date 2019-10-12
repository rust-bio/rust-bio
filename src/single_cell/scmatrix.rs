use sprs::CsMatBase;
use std::error::Error;

use crate::single_cell::eds;

// currently fixing the generic to f64
#[derive(Debug)]
pub struct ScMatrix {
    data: CsMatBase<f64, usize, Vec<usize>, Vec<usize>, Vec<f64>>,
    row_names: Vec<String>,
    column_names: Vec<String>,
}

#[derive(Debug)]
pub enum ScMatType {
    EDS,
    MTX,
    H5,
    CSV,
    Dummy(String),
}

impl ScMatrix {
    // returns number of rows in the matrix
    pub fn num_rows(&self) -> usize {
        return self.row_names.len();
    }

    // returns number of columns in the matrix
    pub fn num_columns(&self) -> usize {
        return self.column_names.len();
    }

    // returns the reference to the row names
    pub fn row_names(&self) -> &Vec<String> {
        return &self.row_names;
    }

    // returns the reference to the column names
    pub fn column_names(&self) -> &Vec<String> {
        return &self.column_names;
    }

    // returns the reference to the sparse matrix data
    pub fn data(&self) -> &sprs::CsMat<f64> {
        return &self.data;
    }

    pub fn new(
        data: sprs::CsMat<f64>,
        row_names: Vec<String>,
        column_names: Vec<String>,
    ) -> ScMatrix {
        ScMatrix {
            data: data,
            row_names: row_names,
            column_names: column_names,
        }
    }

    // single cell matrix reader based on the file format
    pub fn reader(mat_type: ScMatType, input_path: &str) -> Result<ScMatrix, Box<dyn Error>> {
        match mat_type {
            ScMatType::EDS => eds::reader(input_path),
            _ => unimplemented!(),
        }
    }

    // single cell matrix writer based on the file format
    pub fn writer(
        mat_type: ScMatType,
        matrix: ScMatrix,
        output_path: &str,
    ) -> Result<(), Box<dyn Error>> {
        match mat_type {
            ScMatType::EDS => eds::writer(matrix, output_path),
            _ => unimplemented!(),
        }
    }
}
