use crate::single_cell::single_cell_experiment::SingleCellExperiment;

impl<'a, T> IntoIterator for &'a SingleCellExperiment<T> {
    type Item = SingleCellExperimentRow<'a, T>;
    type IntoIter = SingleCellExperimentIntoIterator<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        SingleCellExperimentIntoIterator {
            sce: self,
            row_id: 0,
            row_it: self.counts().outer_iterator(),
        }
    }
}

pub struct SingleCellExperimentRow<'a, T> {
    row_id: usize,
    row_counts: sprs::CsVecBase<&'a [usize], &'a [T]>,
    sce: &'a SingleCellExperiment<T>,
}

impl<'a, T> SingleCellExperimentRow<'a, T> {
    pub fn name(&'a self) -> &'a String {
        &self.sce.row_names()[self.row_id]
    }

    pub fn id(&'a self) -> usize {
        self.row_id
    }
}

impl<'a, T> IntoIterator for &'a SingleCellExperimentRow<'a, T> {
    type Item = SingleCellExperimentEntry<'a, T>;
    type IntoIter = SingleCellExperimentIntoRow<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        SingleCellExperimentIntoRow {
            sce: self.sce,
            row_it: self.row_counts.iter(),
        }
    }
}

pub struct SingleCellExperimentEntry<'a, T> {
    sce: &'a SingleCellExperiment<T>,
    count: &'a T,
    col_id: usize,
}

impl<'a, T> SingleCellExperimentEntry<'a, T> {
    pub fn id(&self) -> usize {
        self.col_id
    }

    pub fn name(&self) -> &'a String {
        &self.sce.col_names()[self.col_id]
    }

    pub fn count(&self) -> &T {
        self.count
    }
}

pub struct SingleCellExperimentIntoRow<'a, T> {
    sce: &'a SingleCellExperiment<T>,
    row_it: sprs::vec::VectorIterator<'a, T, usize>,
}

impl<'a, T> Iterator for SingleCellExperimentIntoRow<'a, T> {
    type Item = SingleCellExperimentEntry<'a, T>;

    fn next(&mut self) -> Option<Self::Item> {
        let nval = self.row_it.next();
        match nval {
            Some(x) => Some(SingleCellExperimentEntry {
                sce: self.sce,
                count: x.1,
                col_id: x.0,
            }),
            None => None,
        }
    }
}

pub struct SingleCellExperimentIntoIterator<'a, T> {
    row_id: usize,
    sce: &'a SingleCellExperiment<T>,
    row_it: sprs::OuterIterator<'a, T, usize>,
}

impl<'a, T> Iterator for SingleCellExperimentIntoIterator<'a, T> {
    type Item = SingleCellExperimentRow<'a, T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.row_id >= self.sce.rows() {
            return None;
        }

        let row_id = self.row_id;
        let row_counts = self.row_it.next().expect("can't get rowdata");

        self.row_id += 1;
        Some(SingleCellExperimentRow {
            row_id,
            row_counts,
            sce: self.sce,
        })
    }
}
