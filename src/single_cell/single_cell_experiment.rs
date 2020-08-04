extern crate byteorder;
extern crate csv as ext_csv;
extern crate flate2;
extern crate math;
extern crate num;
extern crate num_traits;
extern crate sprs;

pub mod csv;
pub mod eds;
pub mod mtx;

use sprs::CsMat;
use std::error::Error;
use std::path::Path;

#[derive(Debug, PartialEq)]
pub struct SingleCellExperiment<T> {
    counts: CsMat<T>,
    rows: Vec<String>,
    cols: Vec<String>,
}

impl<T> SingleCellExperiment<T> {
    pub fn cols(&self) -> usize {
        return self.counts.cols();
    }

    pub fn rows(&self) -> usize {
        return self.counts.rows();
    }

    pub fn shape(&self) -> (usize, usize) {
        return (self.counts.rows(), self.counts.cols());
    }

    pub fn nnz(&self) -> usize {
        return self.counts.nnz();
    }

    pub fn row_names(&self) -> &Vec<String> {
        return &self.rows;
    }

    pub fn col_names(&self) -> &Vec<String> {
        return &self.cols;
    }

    pub fn counts(&self) -> &CsMat<T> {
        return &self.counts;
    }

    pub fn transpose_into(self) -> SingleCellExperiment<T> {
        let row_names = self.row_names().to_owned();
        let col_names = self.col_names().to_owned();
        SingleCellExperiment {
            counts: self.counts.transpose_into(),
            rows: col_names,
            cols: row_names,
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

        Ok(SingleCellExperiment {
            counts: counts,
            rows: rows,
            cols: cols,
        })
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
        let counts_matrix: CsMat<T> = mtx::reader(file)?;

        Ok(SingleCellExperiment {
            counts: counts_matrix,
            rows: rows,
            cols: cols,
        })
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
        let counts_matrix: CsMat<T> = csv::reader(file, rows.len(), cols.len())?;

        Ok(SingleCellExperiment {
            counts: counts_matrix,
            rows: rows,
            cols: cols,
        })
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
        let counts_matrix: CsMat<f32> = eds::reader(file, rows.len(), cols.len())?;

        Ok(SingleCellExperiment {
            counts: counts_matrix,
            rows: rows,
            cols: cols,
        })
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

    fn get_test_matrix_usize() -> CsMat<usize> {
        CsMat::new(
            (3, 3),
            vec![0, 2, 4, 5],
            vec![0, 1, 0, 2, 2],
            vec![1, 2, 3, 4, 5],
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

    fn get_test_sce_usize_data() -> (CsMat<usize>, Vec<String>, Vec<String>) {
        let a = get_test_matrix_usize();
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

        let sce_t = match SingleCellExperiment::from_csr(a.transpose_into(), c.clone(), b.clone()) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };
        assert_eq!(sce.transpose_into(), sce_t);
    }

    fn get_temp_file(ext: String) -> std::path::PathBuf {
        let mut dir = std::env::temp_dir();
        dir.push(format!("foo{}", ext));
        dir
    }

    #[test]
    fn test_csv_f32() {
        let (a, b, c) = get_test_sce_f32_data();
        let sce = match SingleCellExperiment::from_csr(a, b.clone(), c.clone()) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        let file = get_temp_file(".csv.gz".to_owned());
        let fname = file.to_str().unwrap();
        match sce.to_csv(fname) {
            Ok(_) => (),
            Err(y) => panic!("ERROR: {}", y),
        };
        println!("{:?}", fname);

        let sce_csv: SingleCellExperiment<f32> = match SingleCellExperiment::from_csv(fname, b, c) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        assert_eq!(sce, sce_csv);
        std::fs::remove_file(fname).expect("can't remove temp file");
    }

    #[test]
    fn test_csv_usize() {
        let (a, b, c) = get_test_sce_usize_data();
        let sce = match SingleCellExperiment::from_csr(a, b.clone(), c.clone()) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        let file = get_temp_file(".csv.usize.gz".to_owned());
        let fname = file.to_str().unwrap();
        match sce.to_csv(fname) {
            Ok(_) => (),
            Err(y) => panic!("ERROR: {}", y),
        };
        println!("{:?}", fname);

        let sce_csv: SingleCellExperiment<usize> = match SingleCellExperiment::from_csv(fname, b, c)
        {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        assert_eq!(sce, sce_csv);
        std::fs::remove_file(fname).expect("can't remove temp file");
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

    #[test]
    fn test_mtx_f32() {
        let (a, b, c) = get_test_sce_f32_data();
        let sce = match SingleCellExperiment::from_csr(a, b.clone(), c.clone()) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        let file = get_temp_file(".mtx".to_owned());
        let fname = file.to_str().unwrap();
        match sce.to_mtx(fname) {
            Ok(_) => (),
            Err(y) => panic!("ERROR: {}", y),
        };
        println!("{:?}", fname);

        let sce_mtx: SingleCellExperiment<f32> = match SingleCellExperiment::from_mtx(fname, b, c) {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        assert_eq!(sce, sce_mtx);
        std::fs::remove_file(fname).expect("can't remove temp file");
    }

    #[test]
    fn test_mtx_usize() {
        let (a, b, c) = get_test_sce_usize_data();
        let sce = match SingleCellExperiment::from_csr(a, b.clone(), c.clone()) {
            Ok(x) => x,
            Err(_) => unreachable!(),
        };

        let file = get_temp_file("_usize.mtx".to_owned());
        let fname = file.to_str().unwrap();
        match sce.to_mtx(fname) {
            Ok(_) => (),
            Err(y) => panic!("ERROR: {}", y),
        };
        println!("{:?}", fname);

        let sce_mtx: SingleCellExperiment<usize> = match SingleCellExperiment::from_mtx(fname, b, c)
        {
            Ok(x) => x,
            Err(y) => panic!("ERROR: {}", y),
        };

        assert_eq!(sce, sce_mtx);
        std::fs::remove_file(fname).expect("can't remove temp file");
    }
}
