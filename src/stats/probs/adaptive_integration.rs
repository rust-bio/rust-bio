// Copyright 2021 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::collections::HashMap;
use std::convert::Into;
use std::hash::Hash;
use std::{
    fmt::Debug,
    ops::{Add, Div, Mul, Sub},
};

use bio::stats::probs::LogProb;
use itertools::Itertools;
use itertools_num::linspace;
use ordered_float::NotNan;

/// Integrate over an interval of type T with a given density function while trying to minimize
/// the number of grid points evaluated and still hit the maximum likelihood point.
/// This is achieved via a binary search over the grid points.
/// The assumption is that the density is unimodal. If that is not the case,
/// the binary search will not find the maximum and the integral can miss features.
pub fn ln_integrate_exp<T, F>(
    mut density: F,
    min_point: T,
    max_point: T,
    max_resolution: T,
) -> LogProb
where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Div<Output = T>
        + Div<NotNan<f64>, Output = T>
        + Mul<Output = T>
        + Into<f64>
        + From<f64>
        + Ord
        + Debug
        + Hash,
    F: FnMut(T) -> LogProb,
    f64: From<T>,
{
    let mut probs = HashMap::new();

    let mut grid_point = |point, probs: &mut HashMap<_, _>| {
        probs.insert(point, density(point));
        point
    };
    let middle_grid_point = |left: T, right: T| (right + left) / NotNan::new(2.0).unwrap();
    // METHOD:
    // Step 1: perform binary search for maximum likelihood point
    // Remember all points.
    let mut left = grid_point(min_point, &mut probs);
    let mut right = grid_point(max_point, &mut probs);
    let mut first_middle = None;
    let mut middle = None;

    while (((right - left) >= max_resolution) && left < right) || middle.is_none() {
        middle = Some(grid_point(middle_grid_point(left, right), &mut probs));

        if first_middle.is_none() {
            first_middle = middle;
        }

        let left_prob = probs.get(&left).unwrap();
        let right_prob = probs.get(&right).unwrap();

        if left_prob > right_prob {
            // investigate left window more closely
            right = middle.unwrap();
        } else {
            // investigate right window more closely
            left = middle.unwrap();
        }
    }
    // METHOD: add additional grid point in the initially abandoned arm
    if middle < first_middle {
        grid_point(
            middle_grid_point(first_middle.unwrap(), max_point),
            &mut probs,
        );
    } else {
        grid_point(
            middle_grid_point(min_point, first_middle.unwrap()),
            &mut probs,
        );
    }
    // METHOD additionally investigate small interval around the optimum
    for point in linspace(
        cmp::max(
            middle.unwrap() - (max_resolution.into() * 3.0).into(),
            min_point,
        )
        .into(),
        middle.unwrap().into(),
        4,
    )
    .take(3)
    .chain(
        linspace(
            middle.unwrap().into(),
            cmp::min(
                middle.unwrap() + (max_resolution.into() * 3.0).into(),
                max_point,
            )
            .into(),
            4,
        )
        .skip(1),
    ) {
        grid_point(point.into(), &mut probs);
    }

    let sorted_grid_points: Vec<f64> = probs.keys().sorted().map(|point| (*point).into()).collect();

    // METHOD:
    // Step 2: integrate over grid points visited during the binary search.
    LogProb::ln_trapezoidal_integrate_grid_exp::<f64, _>(
        |_, g| *probs.get(&T::from(g)).unwrap(),
        &sorted_grid_points,
    )
}
