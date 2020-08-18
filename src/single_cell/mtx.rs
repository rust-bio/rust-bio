use sprs::CsMat;
use std::error::Error;
use std::path::Path;

// reads the MTX format single cell matrix from the given path
pub fn reader<MatValT>(file_path: &Path) -> Result<CsMat<MatValT>, Box<dyn Error>>
where
    MatValT: Clone + num_traits::Num + num_traits::NumCast,
{
    let matrix = sprs::io::read_matrix_market(file_path)?;
    Ok(matrix.to_csr())
}

// writes the MTX format single cell matrix from the given path
pub fn writer<MatValT>(file_path: &Path, matrix: &CsMat<MatValT>) -> Result<(), Box<dyn Error>>
where
    MatValT: std::fmt::Display + Copy + sprs::num_kinds::PrimitiveKind,
{
    let file_path = file_path.to_str().expect("can't extract file path");
    sprs::io::write_matrix_market(file_path, matrix)?;
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
