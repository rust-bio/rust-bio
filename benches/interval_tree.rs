#![feature(test)]

extern crate test;

use test::Bencher;

use bio::data_structures::interval_tree::*;
use bio::utils::Interval;
use std::cmp::{max, min};
use std::ops::Range;

#[bench]
fn bench_interval_few_large_queries(b: &mut Bencher) {
    // insert 100_000 intervals of size 10
    // do 1000 queries, each resulting in 1000 matches
    b.iter(|| test_insert_query(10, &(100_000..200_000), 1000, 105_000..106_000));
}

#[bench]
fn bench_interval_many_small_queries(b: &mut Bencher) {
    // insert 100_000 intervals of size 10
    // do 100_000 queries, each resulting in at most 10 matches
    b.iter(|| test_insert_query(10, &(100_000..200_000), 10, 99_995..199_995));
}

fn test_insert_query(
    insert_size: i64,
    insert_bounds: &Range<i64>,
    query_size: i64,
    query_bounds: Range<i64>,
) {
    let mut tree: IntervalTree<i64, Range<i64>> = IntervalTree::new();

    for i in insert_bounds.clone() {
        tree.insert(i..i + insert_size, i..i + insert_size);
    }
    for i in query_bounds {
        let lower_bound = i;
        let upper_bound = i + query_size;
        let smallest_start = max(lower_bound - insert_size + 1, insert_bounds.start);
        let largest_start = min(upper_bound, insert_bounds.end);
        let mut expected_intersections = vec![];
        for j in smallest_start..largest_start {
            expected_intersections.push(j..j + insert_size);
        }
        assert_intersections(&tree, lower_bound..upper_bound, &expected_intersections);
    }
}

fn assert_intersections(
    tree: &IntervalTree<i64, Range<i64>>,
    target: Range<i64>,
    expected_results: &[Range<i64>],
) {
    let mut actual_entries: Vec<_> = tree.find(target).collect();
    actual_entries.sort_by(|x1, x2| x1.data().start.cmp(&x2.data().start));
    let mut expected_entries: Vec<_> = expected_results
        .iter()
        .map(|x| (x.clone(), Interval::from(x.clone())))
        .collect();
    expected_entries.sort_by(|x1, x2| x1.0.start.cmp(&x2.0.start));
    assert_eq!(actual_entries.len(), expected_entries.len());
    for (actual, expected) in actual_entries.iter().zip(expected_entries.iter()) {
        assert_eq!(actual.interval(), &expected.1);
        assert_eq!(actual.data(), &expected.0);
    }
}
