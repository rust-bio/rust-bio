use sprs::CsMat;
use std::error::Error;
use std::path::Path;

use crate::single_cell::{csv, eds, mtx};

#[derive(Debug, PartialEq)]
pub struct SingleCellExperiment<T> {
    counts: CsMat<T>,
    rows: Vec<String>,
    cols: Vec<String>,
}

impl<T> SingleCellExperiment<T> {
    pub fn cols(&self) -> usize {
        self.counts.cols()
    }

    pub fn rows(&self) -> usize {
        self.counts.rows()
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.counts.rows(), self.counts.cols())
    }

    pub fn nnz(&self) -> usize {
        self.counts.nnz()
    }

    pub fn row_names(&self) -> &Vec<String> {
        &self.rows
    }

    pub fn get_row_name(&self, index: usize) -> Option<&String> {
        if self.rows() <= index {
            return None;
        }
        Some(&self.rows[index])
    }

    pub fn col_names(&self) -> &Vec<String> {
        &self.cols
    }

    pub fn counts(&self) -> &CsMat<T> {
        &self.counts
    }

    pub fn transpose_into(self) -> SingleCellExperiment<T> {
        let rows = self.row_names().to_owned();
        let cols = self.col_names().to_owned();
        SingleCellExperiment {
            counts: self.counts.transpose_into(),
            rows,
            cols,
        }
    }

    pub fn from_csr(
        counts: CsMat<T>,
        rows: Vec<String>,
        cols: Vec<String>,
    ) -> Result<SingleCellExperiment<T>, String> {
        if rows.len() != counts.rows() {
            return Err("Number of rows in the matrix doesn't match the row names".to_owned());
        }

        if cols.len() != counts.cols() {
            return Err(
                "Number of columns in the matrix doesn't match the column names".to_owned(),
            );
        }

        Ok(SingleCellExperiment { counts, rows, cols })
    }

    pub fn from_mtx(
        file_path: &str,
        rows: Vec<String>,
        cols: Vec<String>,
    ) -> Result<SingleCellExperiment<T>, Box<dyn Error>>
    where
        T: Clone + num_traits::Num + num_traits::NumCast,
    {
        let file = Path::new(file_path);
        let counts: CsMat<T> = mtx::reader(file)?;

        Ok(SingleCellExperiment { counts, rows, cols })
    }

    pub fn from_csv(
        file_path: &str,
        rows: Vec<String>,
        cols: Vec<String>,
    ) -> Result<SingleCellExperiment<T>, Box<dyn Error>>
    where
        T: std::str::FromStr + num::Num + Clone,
    {
        let file = Path::new(file_path);
        let counts: CsMat<T> = csv::reader(file, rows.len(), cols.len())?;

        Ok(SingleCellExperiment { counts, rows, cols })
    }

    pub fn to_mtx(&self, file_path: &str) -> Result<(), Box<dyn Error>>
    where
        T: std::fmt::Display + Copy + sprs::num_kinds::PrimitiveKind,
    {
        let file = Path::new(file_path);
        mtx::writer(file, self.counts())
    }

    pub fn to_csv(&self, file_path: &str) -> Result<(), Box<dyn Error>>
    where
        T: Copy + num::traits::Zero + std::fmt::Display,
    {
        let file = Path::new(file_path);
        csv::writer(file, self.counts())
    }
}

impl SingleCellExperiment<f32> {
    pub fn from_eds(
        file_path: &str,
        rows: Vec<String>,
        cols: Vec<String>,
    ) -> Result<SingleCellExperiment<f32>, Box<dyn Error>> {
        let file = Path::new(file_path);
        let counts: CsMat<f32> = eds::reader(file, rows.len(), cols.len())?;

        Ok(SingleCellExperiment { counts, rows, cols })
    }

    pub fn to_eds(&self, file_path: &str) -> Result<(), Box<dyn Error>> {
        let file = Path::new(file_path);
        eds::writer(file, self.counts())
    }
}

#[cfg(test)]
mod tests {
    use super::SingleCellExperiment;
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

    #[test]
    fn test_single_cell_experiment() {
        let (a, b, c) = get_test_sce_f32_data();
        let sce = match SingleCellExperiment::from_csr(a.clone(), b.clone(), c.clone()) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        assert_eq!(sce.cols(), 3);
        assert_eq!(sce.rows(), 3);

        assert_eq!(sce.row_names().to_owned(), b);
        assert_eq!(sce.col_names().to_owned(), c);
        assert_eq!(sce.shape(), (3, 3));
        assert_eq!(sce.nnz(), 5);
        assert_eq!(sce.counts().clone(), get_test_matrix_f32());

        let sce_t = match SingleCellExperiment::from_csr(a.transpose_into(), b.clone(), c.clone()) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };
        assert_eq!(sce.transpose_into(), sce_t);
    }
}
