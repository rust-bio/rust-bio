// Copyright 2018 Manuel Holtgrewe, Berlin Institute of Health.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! An implementation of Hidden Markov Models in Rust.
//!
//! ## Examples
//!
//! ### Discrete Emission Distribution
//!
//! We construct the example from Borodovsky & Ekisheva (2006), pp. 80 (also see
//! [these slides](http://cecas.clemson.edu/~ahoover/ece854/refs/Gonze-ViterbiAlgorithm.pdf).
//!
//! ```rust
//! use approx::assert_relative_eq;
//! use bio::stats::hmm::discrete_emission::Model as DiscreteEmissionHMM;
//! use bio::stats::hmm::viterbi;
//! use bio::stats::Prob;
//! use ndarray::array;
//!
//! let transition = array![[0.5, 0.5], [0.4, 0.6]];
//! let observation = array![[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]];
//! let initial = array![0.5, 0.5];
//!
//! let hmm = DiscreteEmissionHMM::with_float(&transition, &observation, &initial)
//!     .expect("Dimensions should be consistent");
//! let (path, log_prob) = viterbi(&hmm, &vec![2, 2, 1, 0, 1, 3, 2, 0, 0]);
//! let prob = Prob::from(log_prob);
//! assert_relative_eq!(4.25e-8_f64, *prob, epsilon = 1e-9_f64);
//! ```
//!
//! ### Continuous (Gaussian) Emission Distribution
//!
//! ```rust
//! use approx::assert_relative_eq;
//! use bio::stats::hmm::univariate_continuous_emission::GaussianModel as GaussianHMM;
//! use bio::stats::hmm::viterbi;
//! use bio::stats::Prob;
//! use ndarray::array;
//! use statrs::distribution::Normal;
//!
//! let transition = array![[0.5, 0.5], [0.4, 0.6]];
//! let observation = vec![
//!     Normal::new(0.0, 1.0).unwrap(),
//!     Normal::new(2.0, 1.0).unwrap(),
//! ];
//! let initial = array![0.5, 0.5];
//!
//! let hmm = GaussianHMM::with_float(&transition, observation, &initial)
//!     .expect("Dimensions should be consistent");
//! let (path, log_prob) = viterbi(
//!     &hmm,
//!     &vec![-0.1, 0.1, -0.2, 0.5, 0.8, 1.1, 1.2, 1.5, 0.5, 0.2],
//! );
//! let prob = Prob::from(log_prob);
//! assert_relative_eq!(2.64e-8_f64, *prob, epsilon = 1e-9_f64);
//! ```
//!
//! ## Numeric Stability
//!
//! The implementation uses log-scale probabilities for numeric stability.
//!
//! ## Limitations
//!
//! Currently, only discrete and single-variate Gaussian continuous HMMs are implemented.
//! Also, only dense transition matrices are supported.
//!
//! ## References
//!
//! - Rabiner, Lawrence R. "A tutorial on hidden Markov models and selected applications
//!   in speech recognition." Proceedings of the IEEE 77, no. 2 (1989): 257-286.

pub mod errors;

use std::cmp::Ordering;

use ndarray::prelude::*;
use num_traits::Zero;
use ordered_float::OrderedFloat;
use statrs::distribution::Continuous;

pub use self::errors::{Error, Result};

use super::LogProb;

custom_derive! {
    /// A newtype for HMM states.
    #[derive(
        NewtypeFrom,
        NewtypeDeref,
        PartialEq,
        Copy,
        Clone,
        Debug
    )]
    // #[derive(Serialize, Deserialize)]
    pub struct State(pub usize);
}

/// Iterate over the states of a `Model`.
pub struct StateIter {
    nxt: usize,
    max: usize,
}

impl StateIter {
    /// Constructor.
    pub fn new(num_states: usize) -> Self {
        Self {
            nxt: 0,
            max: num_states,
        }
    }
}

impl Iterator for StateIter {
    type Item = State;

    fn next(&mut self) -> Option<State> {
        if self.nxt < self.max {
            let cur = self.nxt;
            self.nxt += 1;
            Some(State(cur))
        } else {
            None
        }
    }
}

/// Transition between two states in a `Model`.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct StateTransition {
    /// Source of the transition.
    pub src: State,
    /// Destination of the transition.
    pub dst: State,
}

impl StateTransition {
    /// Constructor.
    pub fn new(src: State, dst: State) -> Self {
        Self { src, dst }
    }
}

/// Iterate over all state transitions of a `Model`.
pub struct StateTransitionIter {
    nxt_a: usize,
    nxt_b: usize,
    max: usize,
}

impl StateTransitionIter {
    /// Constructor.
    pub fn new(num_states: usize) -> Self {
        Self {
            nxt_a: 0,
            nxt_b: 0,
            max: num_states,
        }
    }
}

impl Iterator for StateTransitionIter {
    type Item = StateTransition;

    fn next(&mut self) -> Option<StateTransition> {
        let cur_b = self.nxt_b;
        let cur_a = self.nxt_a;
        if self.nxt_b < self.max {
            self.nxt_b += 1;
            Some(StateTransition::new(State(cur_a), State(cur_b)))
        } else if self.nxt_a < self.max {
            self.nxt_b = 0;
            self.nxt_a += 1;
            Some(StateTransition::new(State(cur_a), State(cur_b)))
        } else {
            None
        }
    }
}

/// A trait for Hidden Markov Models (HMM) with generic `Observation` type.
///
/// Rabiner (1989) defines a Hidden Markov Model λ as the tiple (*A*, *B*, π) of transition matrix
/// *A*, emission probabilities *B*, and initial state distribution π.  This has been generalized
/// in `Model` such that you implement `transition_prob()`, `observation_prob()`, and
/// `initial_prob()` (and the other methods; implementation of `transition_prob_idx()` can
/// optionally be implemented and your implementation of `transition_prob()` can then panic).
///
/// The inference algorithm implementations `viterbi()`, `forward()`, and `backward()` will work
/// with any implementation.
///
/// Consequently, this allows for the implementation of HMMs with both discrete and continuous
/// emission distributions.
pub trait Model<Observation> {
    /// The number of states in the model.
    fn num_states(&self) -> usize;

    /// Return iterator over the states of an HMM.
    fn states(&self) -> StateIter;

    /// Returns an iterator of all transitions.
    fn transitions(&self) -> StateTransitionIter;

    /// Transition probability between two states `from` and `to`.
    fn transition_prob(&self, from: State, to: State) -> LogProb;

    /// Transition probability between two states `from` and `to` for observation with index
    /// `_to_idx` (index of `to`).
    ///
    /// This feature comes in handy in several applications of HMMs to biological sequences.
    /// One prominent one is how XHMM by Fromer et al. (2014) uses the distance between target
    /// regions for adjusting the transition probabilities.
    ///
    /// The default implementation return the result of the position-independent
    /// `transition_prob()`.
    fn transition_prob_idx(&self, from: State, to: State, _to_idx: usize) -> LogProb {
        self.transition_prob(from, to)
    }

    /// Initial probability given the HMM `state`.
    fn initial_prob(&self, state: State) -> LogProb;

    /// Probability for the given observation in the given state.
    fn observation_prob(&self, state: State, observation: &Observation) -> LogProb;
}

/// Compute the probability Viterbi matrix and the pointers to the origin.
fn viterbi_matrices<O, M: Model<O>>(
    hmm: &M,
    observations: &[O],
) -> (Array2<LogProb>, Array2<usize>) {
    // The matrix with probabilities.
    let mut vals = Array2::<LogProb>::zeros((observations.len(), hmm.num_states()));
    // For each cell in `vals`, a pointer to the row in the previous column (for the traceback).
    let mut from = Array2::<usize>::zeros((observations.len(), hmm.num_states()));

    // Compute matrix.
    for (i, o) in observations.iter().enumerate() {
        if i == 0 {
            // Initial column.
            for s in hmm.states() {
                vals[[0, *s]] = hmm.initial_prob(s) + hmm.observation_prob(s, o);
                from[[0, *s]] = *s;
            }
        } else {
            // Subsequent columns.
            for j in hmm.states() {
                let x = vals
                    .index_axis(Axis(0), i - 1)
                    .iter()
                    .enumerate()
                    .map(|(a, p)| (State(a), p))
                    .max_by(|(a, &x), (b, &y)| {
                        if x.is_zero() && y.is_zero() {
                            Ordering::Equal
                        } else if x.is_zero() {
                            Ordering::Less
                        } else if y.is_zero() {
                            Ordering::Greater
                        } else {
                            (x + hmm.transition_prob_idx(*a, j, i))
                                .partial_cmp(&(y + hmm.transition_prob_idx(*b, j, i)))
                                .unwrap()
                        }
                    })
                    .map(|(x, y)| (x, *y))
                    .unwrap();
                vals[[i, *j]] =
                    x.1 + hmm.transition_prob_idx(x.0, j, i) + hmm.observation_prob(j, o);
                from[[i, *j]] = *x.0;
            }
        }
    }

    (vals, from)
}

fn viterbi_traceback(vals: Array2<LogProb>, from: Array2<usize>) -> (Vec<State>, LogProb) {
    // Traceback through matrix.
    let n = vals.len_of(Axis(0));
    let mut result: Vec<State> = Vec::new();
    let mut curr = 0;
    let mut res_prob = LogProb::ln_zero();
    for (i, col) in vals.axis_iter(Axis(0)).rev().enumerate() {
        if i == 0 {
            let tmp = col
                .iter()
                .enumerate()
                .max_by_key(|&(_, item)| OrderedFloat(**item))
                .unwrap();
            curr = tmp.0;
            res_prob = *tmp.1;
        } else {
            curr = from[[n - i, curr]];
        }
        result.push(State(curr));
    }
    result.reverse();

    (result, res_prob)
}

/// Execute Viterbi algorithm on the given slice of `Observation` values to get the maximum a
/// posteriori (MAP) probability.
///
/// ## Arguments
///
/// - `hmm` - the `Model` to run the Viterbi algorithm on
/// - `observations` - a slice of observation values to use in the algorithm
///
/// ## Result
///
/// The resulting pair *(s, p)* is the `Vec<State>` of most probable states given `hmm`
/// and `observations` as well as the probability (as `LogProb`) of path `s`.
///
/// ## Type Parameters
///
/// - `O` - the observation type
/// - `M` - type `Model` type
pub fn viterbi<O, M: Model<O>>(hmm: &M, observations: &[O]) -> (Vec<State>, LogProb) {
    let (vals, from) = viterbi_matrices(hmm, observations);
    viterbi_traceback(vals, from)
}

/// Execute the forward algorithm and return the forward probabilites as `LogProb` values
/// and the resulting forward probability.
///
/// ## Arguments
///
/// - `hmm` - the `Model` to run the forward algorithm on
/// - `observations` - a slice of observation values to use in the algorithm
///
/// ## Result
///
/// The resulting pair (*P*, *p*) is the forward probability table (`P[[s, o]]` is the entry
/// for state `s` and observation `o`) and the overall probability for `observations` (as
/// `LogProb`).
///
/// ## Type Parameters
///
/// - `O` - the observation type
/// - `M` - type `Model` type
pub fn forward<O, M: Model<O>>(hmm: &M, observations: &[O]) -> (Array2<LogProb>, LogProb) {
    // The matrix with probabilities.
    let mut vals = Array2::<LogProb>::zeros((observations.len(), hmm.num_states()));

    // Compute matrix.
    for (i, o) in observations.iter().enumerate() {
        if i == 0 {
            // Initial column.
            for s in hmm.states() {
                vals[[0, *s]] = hmm.initial_prob(s) + hmm.observation_prob(s, o);
            }
        } else {
            // Subsequent columns.
            for j in hmm.states() {
                let xs = hmm
                    .states()
                    .map(|k| {
                        vals[[i - 1, *k]]
                            + hmm.transition_prob_idx(k, j, i)
                            + hmm.observation_prob(j, o)
                    })
                    .collect::<Vec<LogProb>>();
                vals[[i, *j]] = LogProb::ln_sum_exp(&xs);
            }
        }
    }

    // Compute final probability.
    let prob = LogProb::ln_sum_exp(vals.row(observations.len() - 1).to_slice().unwrap());

    (vals, prob)
}

/// Execute the backward algorithm and return the backward probabilities as `LogProb` values
/// and the resulting backward probability.
///
/// ## Arguments
///
/// - `hmm` - the `Model` to run the backward algorithm on
/// - `observations` - a slice of observation values to use in the algorithm
///
/// ## Result
///
/// The resulting pair (*P*, *p*) is the backward probability table (`P[[s, o]]` is the entry
/// for state `s` and observation `o`) and the overall probability for `observations` (as
/// `LogProb`).
///
/// ## Type Parameters
///
/// - `O` - the observation type
/// - `M` - type `Model` type
pub fn backward<O, M: Model<O>>(hmm: &M, observations: &[O]) -> (Array2<LogProb>, LogProb) {
    // The matrix with probabilities.
    let mut vals = Array2::<LogProb>::zeros((observations.len(), hmm.num_states()));

    // Compute matrix.
    let n = observations.len();
    for (i, o) in observations.iter().rev().enumerate() {
        if i == 0 {
            for j in hmm.states() {
                let maybe_initial = if i == observations.len() - 1 {
                    hmm.initial_prob(j)
                } else {
                    LogProb::ln_one()
                };
                vals[[0, *j]] = LogProb::ln_one() + hmm.observation_prob(j, o) + maybe_initial;
            }
        } else {
            // Previous columns.
            for j in hmm.states() {
                let maybe_initial = if i == observations.len() - 1 {
                    hmm.initial_prob(j)
                } else {
                    LogProb::ln_one()
                };
                let xs = hmm
                    .states()
                    .map(|k| {
                        vals[[i - 1, *k]]
                            + hmm.transition_prob_idx(j, k, n - i)
                            + hmm.observation_prob(j, o)
                            + maybe_initial
                    })
                    .collect::<Vec<LogProb>>();
                vals[[i, *j]] = LogProb::ln_sum_exp(&xs);
            }
        }
    }

    // Compute final probability.
    let prob = LogProb::ln_sum_exp(vals.row(observations.len() - 1).to_slice().unwrap());

    (vals, prob)
}

/// Implementation of Hidden Markov Model with emission values from discrete distributions.
pub mod discrete_emission {
    use super::super::{LogProb, Prob};
    use super::*;

    /// Implementation of a `hmm::Model` with emission values from discrete distributions.
    ///
    /// Log-scale probabilities are used for numeric stability.
    ///
    /// In Rabiner's tutorial, a discrete emission value HMM has `N` states and `M` output symbols.
    /// The state transition matrix with dimensions `NxN` is `A`, the observation probability
    /// distribution is the matrix `B` with dimensions `NxM` and the initial state distribution `pi`
    /// has length `N`.
    #[derive(Debug, PartialEq)]
    pub struct Model {
        /// The state transition matrix (size `NxN`), `A` in Rabiner's tutorial.
        transition: Array2<LogProb>,

        /// The observation symbol probability distribution (size `NxM`), `B` in Rabiner's tutorial.
        observation: Array2<LogProb>,

        /// The initial state distribution (size `N`), `pi` in Rabiner's tutorial.
        initial: Array1<LogProb>,
    }

    impl Model {
        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors already in log-probability space.
        pub fn new(
            transition: Array2<LogProb>,
            observation: Array2<LogProb>,
            initial: Array1<LogProb>,
        ) -> Result<Self> {
            let (an0, an1) = transition.dim();
            let (bn, bm) = observation.dim();
            let pin = initial.dim();

            if an0 != an1 || an0 != bn || an0 != pin {
                Err(Error::InvalidDimension {
                    an0,
                    an1,
                    bn,
                    bm,
                    pin,
                })
            } else {
                Ok(Self {
                    transition,
                    observation,
                    initial,
                })
            }
        }

        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors already as `Prob` values.
        pub fn with_prob(
            transition: &Array2<Prob>,
            observation: &Array2<Prob>,
            initial: &Array1<Prob>,
        ) -> Result<Self> {
            Self::new(
                transition.map(|x| LogProb::from(*x)),
                observation.map(|x| LogProb::from(*x)),
                initial.map(|x| LogProb::from(*x)),
            )
        }

        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors with probabilities as `f64` values.
        pub fn with_float(
            transition: &Array2<f64>,
            observation: &Array2<f64>,
            initial: &Array1<f64>,
        ) -> Result<Self> {
            Self::new(
                transition.map(|x| LogProb::from(Prob(*x))),
                observation.map(|x| LogProb::from(Prob(*x))),
                initial.map(|x| LogProb::from(Prob(*x))),
            )
        }
    }

    impl super::Model<usize> for Model {
        fn num_states(&self) -> usize {
            self.transition.dim().0
        }

        fn states(&self) -> StateIter {
            StateIter {
                nxt: 0,
                max: self.num_states(),
            }
        }

        fn transitions(&self) -> StateTransitionIter {
            StateTransitionIter {
                nxt_a: 0,
                nxt_b: 0,
                max: self.num_states(),
            }
        }

        fn transition_prob(&self, from: State, to: State) -> LogProb {
            self.transition[[*from, *to]]
        }

        fn initial_prob(&self, state: State) -> LogProb {
            self.initial[[*state]]
        }

        fn observation_prob(&self, state: State, observation: &usize) -> LogProb {
            self.observation[[*state, *observation]]
        }
    }
}

/// Implementation of Hidden Markov Models with emission values from univariate continuous
/// distributions.
pub mod univariate_continuous_emission {
    use super::super::{LogProb, Prob};
    use super::*;

    /// Implementation of a `hmm::Model` with emission values from univariate continuous distributions.
    ///
    /// Log-scale probabilities are used for numeric stability.
    pub struct Model<Dist: Continuous<f64, f64>> {
        /// The state transition matrix (size `NxN`), `A` in Rabiner's tutorial.
        transition: Array2<LogProb>,

        /// The emission probability distributions.
        observation: Vec<Dist>,

        /// The initial state distribution (size `N`), `pi` in Rabiner's tutorial.
        initial: Array1<LogProb>,
    }

    impl<Dist: Continuous<f64, f64>> Model<Dist> {
        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors already in log-probability space.
        pub fn new(
            transition: Array2<LogProb>,
            observation: Vec<Dist>,
            initial: Array1<LogProb>,
        ) -> Result<Self> {
            let (an0, an1) = transition.dim();
            let bn = observation.len();
            let pin = initial.dim();

            if an0 != an1 || an0 != bn || an0 != pin {
                Err(Error::InvalidDimension {
                    an0,
                    an1,
                    bn,
                    bm: bn,
                    pin,
                })
            } else {
                Ok(Self {
                    transition,
                    observation,
                    initial,
                })
            }
        }

        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors already as `Prob` values.
        pub fn with_prob(
            transition: &Array2<Prob>,
            observation: Vec<Dist>,
            initial: &Array1<Prob>,
        ) -> Result<Self> {
            Self::new(
                transition.map(|x| LogProb::from(*x)),
                observation,
                initial.map(|x| LogProb::from(*x)),
            )
        }

        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors with probabilities as `f64` values.
        pub fn with_float(
            transition: &Array2<f64>,
            observation: Vec<Dist>,
            initial: &Array1<f64>,
        ) -> Result<Self> {
            Self::new(
                transition.map(|x| LogProb::from(Prob(*x))),
                observation,
                initial.map(|x| LogProb::from(Prob(*x))),
            )
        }
    }

    impl<Dist: Continuous<f64, f64>> super::Model<f64> for Model<Dist> {
        fn num_states(&self) -> usize {
            self.transition.dim().0
        }

        fn states(&self) -> StateIter {
            StateIter {
                nxt: 0,
                max: self.num_states(),
            }
        }

        fn transitions(&self) -> StateTransitionIter {
            StateTransitionIter {
                nxt_a: 0,
                nxt_b: 0,
                max: self.num_states(),
            }
        }

        fn transition_prob(&self, from: State, to: State) -> LogProb {
            self.transition[[*from, *to]]
        }

        fn initial_prob(&self, state: State) -> LogProb {
            self.initial[[*state]]
        }

        fn observation_prob(&self, state: State, observation: &f64) -> LogProb {
            LogProb::from(Prob::from(self.observation[*state].pdf(*observation)))
        }
    }

    /// Shortcut for HMM with emission values from a Gaussian distribution.
    pub type GaussianModel = Model<statrs::distribution::Normal>;
}

#[cfg(test)]
mod tests {
    use super::super::Prob;
    use ndarray::array;
    use statrs::distribution::Normal;

    use super::discrete_emission::Model as DiscreteEmissionHMM;
    use super::univariate_continuous_emission::GaussianModel as GaussianHMM;
    use super::*;

    #[test]
    fn test_discrete_viterbi_toy_example() {
        // We construct the toy example from Borodovsky & Ekisheva (2006), pp. 80.
        //
        // http://cecas.clemson.edu/~ahoover/ece854/refs/Gonze-ViterbiAlgorithm.pdf
        //
        // States: 0=High GC content, 1=Low GC content
        // Symbols: 0=A, 1=C, 2=G, 3=T
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = array![[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]];
        let initial = array![0.5, 0.5];

        let hmm = DiscreteEmissionHMM::with_float(&transition, &observation, &initial)
            .expect("Dimensions should be consistent");
        let (path, log_prob) = viterbi(&hmm, &vec![2, 2, 1, 0, 1, 3, 2, 0, 0]);
        let prob = Prob::from(log_prob);

        let expected = vec![0, 0, 0, 1, 1, 1, 1, 1, 1]
            .iter()
            .map(|i| State(*i))
            .collect::<Vec<State>>();
        assert_eq!(expected, path);
        assert_relative_eq!(4.25e-8_f64, *prob, epsilon = 1e-9_f64);
    }

    #[test]
    fn test_discrete_forward_toy_example() {
        // Same toy example as above.
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = array![[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]];
        let initial = array![0.5, 0.5];

        let hmm = DiscreteEmissionHMM::with_float(&transition, &observation, &initial)
            .expect("Dimensions should be consistent");
        let log_prob = forward(&hmm, &vec![2, 2, 1, 0]).1;
        let prob = Prob::from(log_prob);

        assert_relative_eq!(0.0038432_f64, *prob, epsilon = 0.0001);
    }

    #[test]
    fn test_discrete_backward_toy_example() {
        // Same toy example as above.
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = array![[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]];
        let initial = array![0.5, 0.5];

        let hmm = DiscreteEmissionHMM::with_float(&transition, &observation, &initial)
            .expect("Dimensions should be consistent");
        let log_prob = backward(&hmm, &vec![2, 2, 1, 0]).1;
        let prob = Prob::from(log_prob);

        assert_relative_eq!(0.0038432_f64, *prob, epsilon = 0.0001);
    }

    #[test]
    fn test_discrete_forward_equals_backward_toy_example() {
        // Same toy example as above.
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = array![[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]];
        let initial = array![0.5, 0.5];
        let hmm = DiscreteEmissionHMM::with_float(&transition, &observation, &initial)
            .expect("Dimensions should be consistent");

        for len in 1..10 {
            let mut seq: Vec<usize> = vec![0; len];
            while seq.iter().sum::<usize>() != len {
                for i in 0..len {
                    if seq[i] == 0 {
                        seq[i] = 1;
                        break;
                    } else {
                        seq[i] = 0;
                    }
                }

                let prob_fwd = *Prob::from(forward(&hmm, &seq).1);
                let prob_bck = *Prob::from(backward(&hmm, &seq).1);
                assert_relative_eq!(prob_fwd, prob_bck, epsilon = 0.00001);
            }
        }
    }

    #[test]
    fn test_gaussian_viterbi_simple_example() {
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = vec![
            Normal::new(0.0, 1.0).unwrap(),
            Normal::new(2.0, 1.0).unwrap(),
        ];
        let initial = array![0.5, 0.5];

        let hmm = GaussianHMM::with_float(&transition, observation, &initial)
            .expect("Dimensions should be consistent");
        let (path, log_prob) = viterbi(
            &hmm,
            &vec![-0.1, 0.1, -0.2, 0.5, 0.8, 1.1, 1.2, 1.5, 0.5, 0.2],
        );
        let prob = Prob::from(log_prob);

        let expected = vec![0, 0, 0, 0, 0, 1, 1, 1, 0, 0]
            .iter()
            .map(|i| State(*i))
            .collect::<Vec<State>>();
        assert_eq!(expected, path);
        assert_relative_eq!(2.64e-8_f64, *prob, epsilon = 1e-9_f64);
    }

    #[test]
    fn test_gaussian_forward_simple_example() {
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = vec![
            Normal::new(0.0, 1.0).unwrap(),
            Normal::new(2.0, 1.0).unwrap(),
        ];
        let initial = array![0.5, 0.5];

        let hmm = GaussianHMM::with_float(&transition, observation, &initial)
            .expect("Dimensions should be consistent");
        let log_prob = forward(&hmm, &vec![0.1, 1.5, 1.8, 2.2, 0.5]).1;
        let prob = Prob::from(log_prob);

        assert_relative_eq!(7.820e-4_f64, *prob, epsilon = 1e-5_f64);
    }

    #[test]
    fn test_gaussian_backward_simple_example() {
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = vec![
            Normal::new(0.0, 1.0).unwrap(),
            Normal::new(2.0, 1.0).unwrap(),
        ];
        let initial = array![0.5, 0.5];

        let hmm = GaussianHMM::with_float(&transition, observation, &initial)
            .expect("Dimensions should be consistent");
        let log_prob = backward(&hmm, &vec![0.1, 1.5, 1.8, 2.2, 0.5]).1;
        let prob = Prob::from(log_prob);

        assert_relative_eq!(7.820e-4_f64, *prob, epsilon = 1e-5_f64);
    }

    #[test]
    fn test_gaussian_forward_equals_backward_simple_example() {
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = vec![
            Normal::new(0.0, 1.0).unwrap(),
            Normal::new(2.0, 1.0).unwrap(),
        ];
        let initial = array![0.5, 0.5];
        let hmm = GaussianHMM::with_float(&transition, observation, &initial)
            .expect("Dimensions should be consistent");

        let seqs = vec![vec![0.1, 0.5, 1.0, 1.5, 1.8, 2.1]];
        for seq in &seqs {
            let prob_fwd = *Prob::from(forward(&hmm, &seq).1);
            let prob_bck = *Prob::from(backward(&hmm, &seq).1);
            assert_relative_eq!(prob_fwd, prob_bck, epsilon = 0.00001);
        }
    }
}
