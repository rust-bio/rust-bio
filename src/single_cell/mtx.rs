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
