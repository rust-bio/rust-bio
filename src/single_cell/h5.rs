use crate::single_cell::scmatrix::ScMatrix;
use hdf5;
use std::error::Error;

// writes the MTX format single cell matrix into the given path
pub fn writer(matrix: ScMatrix, path_str: &str) -> Result<(), Box<dyn Error>> {
    let file = hdf5::File::open(path_str, "w").expect("can't create output file");

    let group = file
        .create_group("matrix")
        .expect("can't create group in h5");

    {
        let shape = group
            .new_dataset::<u64>()
            .create("shape", 2)
            .expect("can't write shape in h5");

        shape
            .write(&[matrix.num_rows(), matrix.num_columns()])
            .expect("error writing shape");
    }

    // can't figure out string writing into h5
    //{
    //    let row_names = group
    //        .new_dataset::<u8>()
    //        .create("barcodes", matrix.num_rows())
    //        .expect("can't write barcodes in h5");

    //    let names_as_bytes:Vec<Vec<u8>> = matrix.row_names()
    //        .iter()
    //        .map(|x| x.as_bytes())
    //        .collect();

    //    row_names
    //        .write_raw(&names_as_bytes)
    //        .expect("error writing row names");
    //}

    //{
    //    let column_names = group
    //        .new_dataset::<u8>()
    //        .create("features", matrix.num_columns())
    //        .expect("can't write features in h5");

    //    column_names
    //        .write_raw(&matrix.column_names().as_bytes())
    //        .expect("error writing column names");
    //}

    {
        let indptr = group
            .new_dataset::<usize>()
            .gzip(6)
            .create("indptr", matrix.num_rows())
            .expect("can't write indptr in h5");

        indptr
            .write_raw(&matrix.data().indptr())
            .expect("error writing indptr");
    } // end writing indptr

    {
        let data = group
            .new_dataset::<f64>()
            .gzip(6)
            .create("data", matrix.nnz())
            .expect("can't write data in h5");

        data.write_raw(&matrix.data().data())
            .expect("can't write data");
    } // end writing data

    {
        let indices = group
            .new_dataset::<u32>()
            .gzip(6)
            .create("indices", matrix.nnz())
            .expect("can't write positions in h5");

        indices
            .write_raw(&matrix.data().indices())
            .expect("can't write indices");
    } // end writing indices

    Ok(())
}
