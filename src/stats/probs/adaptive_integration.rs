// Copyright 2021-2022 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp;
use std::collections::HashMap;
use std::convert::{Into, TryFrom};
use std::hash::Hash;
use std::{
    fmt::Debug,
    ops::{Add, Div, Mul, Sub},
};

use crate::stats::probs::LogProb;
use itertools::Itertools;
use itertools_num::linspace;
use ordered_float::NotNan;

/// Integrate over an interval of type T with a given density function while trying to minimize
/// the number of grid points evaluated and still hit the maximum likelihood point.
/// This is achieved via a binary search over the grid points.
/// The assumption is that the density is unimodal. If that is not the case,
/// the binary search will not find the maximum and the integral can miss features.
///
/// # Example
///
/// ```rust
/// use approx::abs_diff_eq;
/// use bio::stats::probs::adaptive_integration::ln_integrate_exp;
/// use bio::stats::probs::{LogProb, Prob};
/// use ordered_float::NotNan;
/// use statrs::distribution::{Continuous, Normal};
/// use statrs::statistics::Distribution;
///
/// let ndist = Normal::new(0.0, 1.0).unwrap();
///
/// let integral = ln_integrate_exp(
///     |x| LogProb::from(Prob(ndist.pdf(*x))),
///     NotNan::new(-1.0).unwrap(),
///     NotNan::new(1.0).unwrap(),
///     NotNan::new(0.01).unwrap(),
/// );
/// abs_diff_eq!(integral.exp(), 0.682, epsilon = 0.01);
/// ```
pub fn ln_integrate_exp<T, F, E>(
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
        + TryFrom<f64, Error = E>
        + Ord
        + Debug
        + Hash,
    E: Debug,
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
    // After that loop, we are guaranteed that middle.is_some().
    let middle = middle.unwrap();
    let first_middle = first_middle.unwrap();
    // METHOD: add additional grid point in the initially abandoned arm
    if middle < first_middle {
        grid_point(middle_grid_point(first_middle, max_point), &mut probs);
    } else {
        grid_point(middle_grid_point(min_point, first_middle), &mut probs);
    }
    // METHOD additionally investigate small interval around the optimum
    for point in linspace(
        cmp::max(
            T::try_from(middle.into() - max_resolution.into() * 3.0).unwrap(),
            min_point,
        )
        .into(),
        middle.into(),
        4,
    )
    .take(3)
    .chain(
        linspace(
            middle.into(),
            cmp::min(
                T::try_from(middle.into() + max_resolution.into() * 3.0).unwrap(),
                max_point,
            )
            .into(),
            4,
        )
        .skip(1),
    ) {
        grid_point(T::try_from(point).unwrap(), &mut probs);
    }

    let sorted_grid_points: Vec<f64> = probs.keys().sorted().map(|point| (*point).into()).collect();

    // METHOD:
    // Step 2: integrate over grid points visited during the binary search.
    LogProb::ln_trapezoidal_integrate_grid_exp::<f64, _>(
        |_, g| *probs.get(&T::try_from(g).unwrap()).unwrap(),
        &sorted_grid_points,
    )
}
