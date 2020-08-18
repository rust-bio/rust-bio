extern crate csv as ext_csv;

use std::error::Error;
use std::path::Path;

use sprs::{CsMat, TriMat};

// reads the CSV format single cell matrix from the given path
pub fn reader<MatValT>(
    file_path: &Path,
    num_rows: usize,
    num_cols: usize,
) -> Result<CsMat<MatValT>, Box<dyn Error>>
where
    MatValT: std::str::FromStr + num::Num + Clone,
{
    let mut tri_matrix = TriMat::new((num_rows, num_cols));
    let mut rdr = ext_csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(file_path)?;

    for (row_id, line) in rdr.records().enumerate() {
        let record = line?;
        let values: Vec<MatValT> = record.into_iter().flat_map(str::parse::<MatValT>).collect();
        assert_eq!(values.len(), num_cols);

        for (column_id, val) in values.into_iter().enumerate() {
            if val == MatValT::zero() {
                continue;
            }
            tri_matrix.add_triplet(row_id, column_id, val);
        }
    }

    Ok(tri_matrix.to_csr())
}

// writes the CSV format single cell matrix into the given path
pub fn writer<MatValT>(file_path: &Path, matrix: &CsMat<MatValT>) -> Result<(), Box<dyn Error>>
where
    MatValT: Copy + num::traits::Zero + std::fmt::Display,
{
    // writing matrix
    let mut wtr = ext_csv::WriterBuilder::new().from_path(file_path)?;

    let num_columns = matrix.cols();
    let zero: MatValT = MatValT::zero();
    let mut columns: Vec<MatValT> = vec![zero; num_columns];

    for row_vec in matrix.outer_iterator() {
        columns.iter_mut().for_each(|x| *x = zero);

        let mut it = row_vec.iter();
        while let Some((col_idx, &val)) = it.next() {
            columns[col_idx] = val;
        }

        let record: Vec<String> = columns.iter().map(|x| x.to_string()).collect();
        wtr.write_record(&record)?;
    } // end row iterator

    wtr.flush()?;
    Ok(())
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
}
