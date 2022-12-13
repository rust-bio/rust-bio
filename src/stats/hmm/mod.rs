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
//! ### Trainning a Discrete Emission Model with Baum-Welch algorithm
//!
//! We construct the example from Jason Eisner lecture which can be followed along
//! with his spreadsheet ([link](http://www.cs.jhu.edu/~jason/papers/#eisner-2002-tnlp)).
//! Take a look at tests in source file.
//!
//! ## Numeric Stability
//!
//! The implementation uses log-scale probabilities for numeric stability.
//!
//! ## Limitations
//!
//! Currently, only discrete and single-variate Gaussian continuous HMMs are implemented.
//! Also, only dense transition matrices and trainning of discrete models are supported.
//!
//! ## References
//!
//! - Rabiner, Lawrence R. "A tutorial on hidden Markov models and selected applications
//!   in speech recognition." Proceedings of the IEEE 77, no. 2 (1989): 257-286.
//! - Eisner, Jason "An interactive spreadsheet for teaching the forward-backward algorithm.
//!   in speech recognition." In ACL Workshop on Teaching NLP and CL (2002).
pub mod errors;

use std::cmp::Ordering;

use ndarray::prelude::*;
use num_traits::Zero;
use ordered_float::OrderedFloat;
use statrs::distribution::Continuous;

pub use self::errors::{Error, Result};

use super::LogProb;

use std::cmp::Eq;
use std::fmt::Debug;
use std::hash::Hash;

use super::Prob;
use std::cell::RefCell;
use std::collections::BTreeMap;

custom_derive! {
    /// A newtype for HMM states.
    // #[derive(
    //     NewtypeFrom,
    //     NewtypeDeref,
    //     Default,
    //     Copy,
    //     Clone,
    //     Eq,
    //     PartialEq,
    //     Ord,
    //     PartialOrd,
    //     Hash,
    //     Debug,
    // )]
    // #[derive(Serialize, Deserialize)]
    #[derive(
        NewtypeFrom,
        NewtypeDeref,
        Default,
        Copy,
        Clone,
        Eq,
        PartialEq,
        Ord,
        PartialOrd,
        Hash,
        Debug,
    )]
    #[derive(serde::Serialize, serde::Deserialize)]
    pub struct State(pub usize);
}

/// Iterate over the states of a `Model`.
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
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
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
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
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
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
#[allow(unconditional_recursion)]
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

    /// End probability given the HMM `state`.
    fn end_prob(&self, _state: State) -> LogProb {
        self.end_prob(_state)
    }

    fn has_end_state(&self) -> bool {
        false
    }
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
                        // + maybe_initial
                    })
                    .collect::<Vec<LogProb>>();
                vals[[i, *j]] = LogProb::ln_sum_exp(&xs);
            }
        }
    }

    // Compute final probability.

    let prob_vec_final = hmm
        .states()
        .map(|k| vals[[observations.len() - 1, *k]] + hmm.end_prob(k))
        .collect::<Vec<LogProb>>();

    let prob = LogProb::ln_sum_exp(&prob_vec_final);
    // let prob = LogProb::ln_sum_exp(vals.row(observations.len() - 1).to_slice().unwrap());

    (vals, prob)
}

pub fn backward<O, M: Model<O>>(hmm: &M, observations: &[O]) -> (Array2<LogProb>, LogProb) {
    // The matrix with probabilities.
    let mut vals = Array2::<LogProb>::zeros((observations.len(), hmm.num_states()));
    let mut prob_vec_final = vec![];

    // Compute matrix.
    let n = observations.len();
    for (i, o) in observations.iter().rev().enumerate() {
        if i == 0 {
            for j in hmm.states() {
                vals[[0, *j]] = hmm.end_prob(j);
            }

            for j in hmm.states() {
                let xs = hmm
                    .states()
                    .map(|k| {
                        vals[[i, *k]]
                            + hmm.transition_prob_idx(j, k, n - i)
                            + hmm.observation_prob(k, o)
                    })
                    .collect::<Vec<LogProb>>();
                if observations.len() > 1 {
                    vals[[i + 1, *j]] = LogProb::ln_sum_exp(&xs);
                } else {
                    prob_vec_final = hmm
                        .states()
                        .map(|k| vals[[i, *k]] + hmm.initial_prob(k) + hmm.observation_prob(k, o))
                        .collect::<Vec<LogProb>>();
                }
            }
        } else if i == (observations.len() - 1) {
            prob_vec_final = hmm
                .states()
                .map(|k| vals[[i, *k]] + hmm.initial_prob(k) + hmm.observation_prob(k, o))
                .collect::<Vec<LogProb>>();
        } else {
            // Previous columns.
            for j in hmm.states() {
                let xs = hmm
                    .states()
                    .map(|k| {
                        vals[[i, *k]]
                            + hmm.transition_prob_idx(j, k, n - i)
                            + hmm.observation_prob(k, o)
                    })
                    .collect::<Vec<LogProb>>();
                vals[[i + 1, *j]] = LogProb::ln_sum_exp(&xs);
            }
        }
    }

    let prob = LogProb::ln_sum_exp(&prob_vec_final);

    (vals, prob)
}

/// Execute **one step** of Baum-Welch algorithm to find the maximum likelihood estimate of the parameters of a HMM given a set of observed
/// feature vector and return the estimated initial state distribution (*π**), estimated transition matrix (*A**),
///  estimated emission probabilities matrix (*B**) and end probabilities vector (if the model has declared an end state beforehand).
/// This function doesn't update the HMM parameters in the model and has been implemented for Discrete Emissions Models only.
/// It return values as `LogProb`.
///
/// ## Arguments
///
/// - `hmm` - the `Model` to run the baum-welch or expected maximization algorithm on. It has to be a Discrete Model with a `Trainable` trait implemented.
/// - `observations` - a slice of observation values to use in the algorithm
///
/// ## Result
///
/// The resulting tuple (*π**, *A**, B*, E*) is the estimated initial probability table (`P[s]`),
/// the estimated transitions probability table (`P[[s, o]]` is the entry
/// for state `s1` and other state `s2`), the estimated emission probability table (`P[[s, s]]` is the entry
/// for state `s` and observation class `o`) and if we specify an end probability when building the model,
/// E* is the estimated end probabilities. Otherwise, E* is a vector with size equal to initial probability
/// and all values set to LogProb(1.0).  The values in all outputs are shown as `LogProb`.
///
/// ## Type Parameters
///
/// - `O` - the observation type (only discrete emissions)
/// - `M` - type `Model` type
pub fn baum_welch<O: Debug + Eq + Hash + Ord, M: Model<O>>(
    hmm: &M,
    observations: &[O],
) -> (
    Array1<LogProb>,
    Array2<LogProb>,
    Array2<LogProb>,
    Array1<LogProb>,
) {
    // Execute forward algorithm to calculate the alpha probabilities for each time-step.
    // Ignore P(x)
    let (prob_table_f, _) = forward(hmm, observations);
    // Execute backward algorithm to calculate the beta probabilities for each time-step.
    let (prob_table_b_cor, _) = backward(hmm, observations);

    // The time-step in forward matrix goes from t=0 to t=N. Coversely, the backward matrix
    // entries are in reverse order, i.e., from t=N to t=0.
    // The following code puts the backward matrix in the same order as forward matrix.
    // This is usefull when computing the alfa times beta probabilities.
    let n = observations.len() - 1;
    let mut prob_table_b = Array2::<LogProb>::zeros((observations.len(), hmm.num_states()));
    for (j, el) in prob_table_b_cor.axis_iter(Axis(0)).enumerate() {
        prob_table_b.row_mut(n - j).assign(&el);
    }

    // Product of the forward probability and backward probability which is translated
    // into a sum in log space. Results in gamma array which is an array of size = len(Observation) x N_states representing
    // the probability of being in state t at time j for this observations.
    let alpha_betas = &prob_table_f + &prob_table_b;

    // Log-likelihood of the observation given the model.
    // This probability must be the same as probabilities calculated by forward or backward functions
    let probx = LogProb::ln_sum_exp(alpha_betas.row(observations.len() - 1).to_slice().unwrap());

    // Dictionary-like vec_hashs_prob_obs stores b_{i}^{*}(v_{k}) which is the expected number of times the output observations
    // have been equal to observation v_{k} while in state i over the expected
    // total number of times in state i.
    let mut vec_hashs_prob_obs: Vec<BTreeMap<&O, LogProb>> = Vec::new();

    let mut distinct_obs = 0;

    for h in hmm.states() {
        let mut probs_observations: BTreeMap<&O, LogProb> = BTreeMap::new();

        // First, for the numerator of b_{i}^{^}(v_{k}) we sum gamma for all time steps t in which the observation o_{t}
        // is the symbol v_{k}
        for (t, o) in observations.iter().enumerate() {
            // For example, in einsner example, it calculates the probabilities P(->state, observation), for example p(->C,1) or p(->C,2) or p(->C,3)
            let p = probs_observations.entry(o).or_insert_with(LogProb::ln_zero);
            *p = (*p).ln_add_exp(alpha_betas[[t, *h]] - probx);
        }
        distinct_obs = probs_observations.len();
        vec_hashs_prob_obs.push(probs_observations);
    }

    let mut vals_xi =
        Array2::<LogProb>::zeros((observations.len(), hmm.num_states() * hmm.num_states()));

    // Calculate the arc probabilities.
    // It's the xi equation which it's entries represent the probability of being in state i and j,
    // in times t and t+1, respectively, given the observed sequence and parameters.
    // The vals_xi matrix is of size (time-steps, N_states^2) to store for each observation (or time-step)
    // the probability of being at state i and j in times t and t+1. The element vals_xi[t, k] indicates
    // the probability in time t and k is the index encoding of the pair of states (i,j).
    for (t, o) in observations.iter().enumerate() {
        if t == 0 {
            continue;
        } else {
            for (idxstate, j) in hmm.states().enumerate() {
                let numerador = hmm
                    .states()
                    .map(|i| {
                        prob_table_f[[t - 1, *j]]
                            + hmm.transition_prob_idx(j, i, t)
                            + prob_table_b[[t, *i]]
                            + hmm.observation_prob(i, o)
                            - probx
                    })
                    .collect::<Vec<LogProb>>();

                vals_xi
                    .slice_mut(s![
                        t,
                        idxstate * hmm.num_states()..(idxstate + 1) * hmm.num_states()
                    ])
                    .assign(&Array1::from(numerador));
            }
        }
    }

    let mut sum_p_states = Vec::new();
    for k in hmm.states() {
        let sum_prob_states =
            LogProb::ln_sum_exp(&alpha_betas.column(*k).map(|x| x - probx).to_vec());
        sum_p_states.push(sum_prob_states);
    }

    let mut observations_hat = Array2::<LogProb>::zeros((hmm.num_states(), distinct_obs));
    let mut transitions_hat = Array2::<LogProb>::zeros((hmm.num_states(), hmm.num_states()));

    for (idxstate, i) in hmm.states().enumerate() {
        let gamma_i = LogProb::ln_sum_exp(
            &alpha_betas
                .column(*i)
                .to_vec()
                .iter()
                .map(|x| *x - probx)
                .collect::<Vec<LogProb>>(),
        );
        let end_i = if hmm.has_end_state() {
            LogProb::ln_zero()
        } else {
            alpha_betas[[observations.len() - 1, *i]] - probx
        };
        let q = vals_xi.slice(s![
            ..,
            idxstate * hmm.num_states()..(idxstate + 1) * hmm.num_states()
        ]);
        for k in hmm.states() {
            let sa = LogProb::ln_sum_exp(&q.column(*k).to_vec());
            transitions_hat[[*i, *k]] = sa - gamma_i.ln_sub_exp(end_i);
        }

        let mut ind_probs = vec![];
        for v in vec_hashs_prob_obs[*i].values() {
            ind_probs.push(*v - gamma_i);
        }
        observations_hat
            .row_mut(*i)
            .assign(&Array1::from(ind_probs));
    }

    let pi_hat = Array1::from(
        alpha_betas
            .row(0)
            .to_vec()
            .iter()
            .map(|x| *x - probx)
            .collect::<Vec<LogProb>>(),
    );

    let end_hat = if hmm.has_end_state() {
        Array1::from(
            sum_p_states
                .iter()
                .zip(alpha_betas.row(observations.len() - 1).to_vec().iter())
                .map(|(sg, g)| (*g - probx) - *sg)
                .collect::<Vec<LogProb>>(),
        )
    } else {
        Array1::from(
            pi_hat
                .iter()
                .map(|_x| LogProb::from(Prob(1.0)))
                .collect::<Vec<LogProb>>(),
        )
    };

    (pi_hat, transitions_hat, observations_hat, end_hat)
}

/// A trait for trainning Hidden Markov Models (HMM) with generic `Observation` type using Baum-Welch algorithm.
pub trait Trainable<Observation> {
    ///  Iterative procedure to train the model using Baum-Welch algorithm given the training sequences.
    ///
    /// As arguments, a set of sequences (observations) and two optional argumets: maximum number of iterations (`n_iter`) and tolerance (`tol`).
    /// The baum-welch iterative training procedure will stop either if it reaches the tolerance of the relative log-likelihood augmentation (default `1e-6`) or
    /// exceed the maximum number of iterations (default `500`).
    fn train_baum_welch(
        &self,
        observations: &[Vec<Observation>],
        n_iter: Option<usize>,
        tol: Option<f64>,
    ) -> (
        Array1<LogProb>,
        Array2<LogProb>,
        Array2<LogProb>,
        Array1<LogProb>,
    );

    /// This feature comes in handy in Bam-Welch algorithm when doing an update of HMM parameters.
    ///
    /// After receiving the estimated parameters found after trainning, this method updates the values in the
    /// HMM model.
    fn update_matrices(
        &self,
        transition_hat: Array2<LogProb>,
        observation_hat: Array2<LogProb>,
        initial_hat: Array1<LogProb>,
        end_hat: Array1<LogProb>,
    );
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
    #[derive(Default, Clone, PartialEq, Debug)]
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

        fn end_prob(&self, _state: State) -> LogProb {
            LogProb::ln_one()
        }
    }
}

/// Implementation of Hidden Markov Model with emission values from discrete distributions and an optional explicity end state.
/// This module also implements the `Trainable` trait allowing to be trainned by Baum-Welch algorithm.
pub mod discrete_emission_opt_end {
    use super::super::{LogProb, Prob};
    use super::*;

    /// Implementation of a `hmm::Model` with emission values from discrete distributions and an optional declared end state.
    ///
    /// Log-scale probabilities are used for numeric stability.
    ///
    /// In Rabiner's tutorial, a discrete emission value HMM has `N` states and `M` output symbols.
    /// The state transition matrix with dimensions `NxN` is `A`, the observation probability
    /// distribution is the matrix `B` with dimensions `NxM` and the initial state distribution `pi`
    /// has length `N`. We also included a silent end state `ε` with vector length `N` that do not emit symbols for
    /// modelling the end of sequences. It's optional to supply the end probabilities at the creation of the model.
    /// If this happens, we'll create a dummy end state to simulate as if the end state has not been included.
    #[derive(Default, Clone, PartialEq, Debug)]
    pub struct Model {
        /// The state transition matrix (size `NxN`), `A` in Rabiner's tutorial.
        transition: RefCell<Array2<LogProb>>,

        /// The observation symbol probability distribution (size `NxM`), `B` in Rabiner's tutorial.
        observation: RefCell<Array2<LogProb>>,

        /// The initial state distribution (size `N`), `pi` in Rabiner's tutorial.
        initial: RefCell<Array1<LogProb>>,

        end: RefCell<Array1<LogProb>>,

        has_end_state: bool,
    }

    impl Model {
        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors already in log-probability space.
        pub fn new(
            transition: RefCell<Array2<LogProb>>,
            observation: RefCell<Array2<LogProb>>,
            initial: RefCell<Array1<LogProb>>,
            end: RefCell<Array1<LogProb>>,
            has_end_state: bool,
        ) -> Result<Self> {
            let (an0, an1) = transition.borrow().dim();
            let (bn, bm) = observation.borrow().dim();
            let pin = initial.borrow().dim();

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
                    end,
                    has_end_state,
                })
            }
        }

        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors already as `Prob` values.
        pub fn with_prob(
            transition: &Array2<Prob>,
            observation: &Array2<Prob>,
            initial: &Array1<Prob>,
            end: Option<&Array1<Prob>>,
        ) -> Result<Self> {
            let initial_dim = initial.dim();
            let end_possible = (0..initial_dim)
                .map(|_x| Prob(1.0))
                .collect::<Array1<Prob>>();

            let has_end_state = end.is_some();

            let end_un = end.unwrap_or(&end_possible);

            Self::new(
                RefCell::new(transition.map(|x| LogProb::from(*x))),
                RefCell::new(observation.map(|x| LogProb::from(*x))),
                RefCell::new(initial.map(|x| LogProb::from(*x))),
                RefCell::new(end_un.map(|x| LogProb::from(*x))),
                has_end_state,
            )
        }

        /// Construct new Hidden MarkovModel with the given transition, observation, and initial
        /// state matrices and vectors with probabilities as `f64` values.
        pub fn with_float(
            transition: &Array2<f64>,
            observation: &Array2<f64>,
            initial: &Array1<f64>,
            end: Option<&Array1<f64>>,
        ) -> Result<Self> {
            let initial_dim = initial.dim();
            let end_possible = (0..initial_dim).map(|_x| 1.0).collect::<Array1<f64>>();
            let end_un = end.unwrap_or(&end_possible);

            let has_end_state = end.is_some();

            Self::new(
                RefCell::new(transition.map(|x| LogProb::from(Prob(*x)))),
                RefCell::new(observation.map(|x| LogProb::from(Prob(*x)))),
                RefCell::new(initial.map(|x| LogProb::from(Prob(*x)))),
                RefCell::new(end_un.map(|x| LogProb::from(Prob(*x)))),
                has_end_state,
            )
        }
    }

    impl super::Model<usize> for Model {
        fn num_states(&self) -> usize {
            let transition = self.transition.borrow();
            transition.dim().0
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
            let transition_prob = self.transition.borrow();
            transition_prob[[*from, *to]]
        }

        fn initial_prob(&self, state: State) -> LogProb {
            let initial_prob = self.initial.borrow();
            initial_prob[[*state]]
        }

        fn end_prob(&self, _state: State) -> LogProb {
            let end_prob = self.end.borrow();
            end_prob[[*_state]]
        }

        fn observation_prob(&self, state: State, observation: &usize) -> LogProb {
            let observation_prob = self.observation.borrow();
            observation_prob[[*state, *observation]]
        }

        fn has_end_state(&self) -> bool {
            self.has_end_state
        }
    }

    impl super::Trainable<usize> for Model {
        fn update_matrices(
            &self,
            transition_hat: Array2<LogProb>,
            observation_hat: Array2<LogProb>,
            initial_hat: Array1<LogProb>,
            end_hat: Array1<LogProb>,
        ) {
            let mut end_prob = self.end.borrow_mut();
            *end_prob = end_hat;

            let mut transition_prob = self.transition.borrow_mut();
            *transition_prob = transition_hat;

            let mut observation_prob = self.observation.borrow_mut();
            *observation_prob = observation_hat;

            let mut initial_prob = self.initial.borrow_mut();
            *initial_prob = initial_hat;
        }

        fn train_baum_welch(
            &self,
            observations: &[Vec<usize>],
            n_iter: Option<usize>,
            tol: Option<f64>,
        ) -> (
            Array1<LogProb>,
            Array2<LogProb>,
            Array2<LogProb>,
            Array1<LogProb>,
        ) {
            let hmm = self;

            let tol = tol.unwrap_or(1e-6_f64);
            let n_iter = n_iter.unwrap_or(500);

            let (pi_hat, transitions_hat, observations_hat, end_hat) =
                super::baum_welch(hmm, &observations[0]);
            let (
                mut pi_hat_ref,
                mut transitions_hat_ref,
                mut observations_hat_ref,
                mut end_hat_ref,
            ) = (pi_hat, transitions_hat, observations_hat, end_hat);
            let (_, mut prob_fwd_new) = super::forward(hmm, &observations[0]);

            // initialize var llh to store the log likelihood of trainned model
            let mut llh: LogProb = LogProb::ln_one();

            // Get the number of observations in first "sample"
            let mut obs_n = observations[0].len() as f64;

            // Normalize the previous log likelihood computed by number of observations
            let mut nllh_o = LogProb::from(Prob((*prob_fwd_new / obs_n).exp()));

            println!(
                "Iter num {:?} - LLH: {:?} - Normalized LLH: {:?}",
                0, prob_fwd_new, nllh_o
            );

            for iteration in 0..(n_iter - 1) {
                for obs in observations.iter() {
                    let (pi_hat, transitions_hat, observations_hat, end_hat) = baum_welch(hmm, obs);

                    pi_hat_ref = pi_hat.to_owned();
                    transitions_hat_ref = transitions_hat.to_owned();
                    observations_hat_ref = observations_hat.to_owned();
                    end_hat_ref = end_hat.to_owned();

                    hmm.update_matrices(transitions_hat, observations_hat, pi_hat, end_hat);

                    let (_, prob_fwd_new_ii) = super::forward(hmm, obs);

                    llh = prob_fwd_new_ii;

                    // Get the number of observations
                    obs_n = obs.len() as f64;
                }

                // Normalize the new log likelihood computed in this iteration by number of observations
                let nllh = LogProb::from(Prob((*llh / obs_n).exp()));

                println!(
                    "Iter num {:?} - LLH: {:?} - Normalized LLH: {:?}",
                    iteration + 1,
                    llh,
                    nllh
                );

                if nllh_o >= nllh {
                    prob_fwd_new = llh;
                    // Normalize the previous log likelihood computed by number of observations
                    nllh_o = LogProb::from(Prob((*prob_fwd_new / obs_n).exp()));
                    // Skip to next iteration
                    continue;
                }

                if nllh.ln_sub_exp(nllh_o) < LogProb::from(Prob(tol)) {
                    // Stop trainning if the difference between the new log like and the old is less than a threshold
                    break;
                } else {
                    // Otherwise, set the old log like as the new log like
                    prob_fwd_new = llh;

                    // Normalize the previous log likelihood computed by number of observations
                    nllh_o = LogProb::from(Prob((*prob_fwd_new / obs_n).exp()));
                }
            }
            (
                pi_hat_ref,
                transitions_hat_ref,
                observations_hat_ref,
                end_hat_ref,
            )
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
    #[derive(Default, Clone, PartialEq, Debug)]
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

        fn end_prob(&self, _state: State) -> LogProb {
            LogProb::ln_one()
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
    use super::discrete_emission_opt_end::Model as DiscreteEmissionHMMoptEND;
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
        let (path, log_prob) = viterbi(&hmm, &[2, 2, 1, 0, 1, 3, 2, 0, 0]);
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
        let log_prob = forward(&hmm, &[2, 2, 1, 0]).1;
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
        let (_, log_prob) = backward(&hmm, &[2, 2, 1, 0]);
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
        let (path, log_prob) = viterbi(&hmm, &[-0.1, 0.1, -0.2, 0.5, 0.8, 1.1, 1.2, 1.5, 0.5, 0.2]);
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
        let log_prob = forward(&hmm, &[0.1, 1.5, 1.8, 2.2, 0.5]).1;
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
        let log_prob = backward(&hmm, &[0.1, 1.5, 1.8, 2.2, 0.5]).1;
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

    #[test]
    fn test_recriate_discrete_backward_toy_example() {
        // Same toy example as above.
        let transition = array![[0.5, 0.5], [0.4, 0.6]];
        let observation = array![[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]];
        let initial = array![0.5, 0.5];

        let hmm = DiscreteEmissionHMMoptEND::with_float(&transition, &observation, &initial, None)
            .expect("Dimensions should be consistent");

        let (_, log_prob) = backward(&hmm, &[2, 2, 1, 0]);
        let prob = Prob::from(log_prob);

        assert_relative_eq!(0.0038432_f64, *prob, epsilon = 0.0001);
    }

    #[test]
    fn test_discrete_with_end_backward_toy_example() {
        // We construct Jason Eisner's dairy ice-cream consumption example.
        //
        // http://www.cs.jhu.edu/~jason/papers/#eisner-2002-tnlp
        //
        // States: 0=Hot day, 1=Cold day
        // Symbols: 0=one ice creams, 1=two ice creams, 2=three ice creams
        let transition = array![[0.8, 0.1], [0.1, 0.8]];
        let observation = array![[0.7, 0.2, 0.1], [0.1, 0.2, 0.7]];
        let initial = array![0.5, 0.5];
        let end = array![0.1, 0.1];

        let ices = vec![
            1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 0, 1, 2, 2,
            1, 2, 1, 1,
        ];

        let hmm =
            DiscreteEmissionHMMoptEND::with_float(&transition, &observation, &initial, Some(&end))
                .expect("Dimensions should be consistent");

        let (_, log_backward_probok) = backward(&hmm, &ices);
        let prob = Prob::from(log_backward_probok);

        assert_relative_eq!(0.912e-18_f64, *prob, epsilon = 0.1e-20_f64);
    }

    #[test]
    fn test_baum_welch_one_iter_example() {
        // Same Jason Eisner's ice cream example as above with little bias
        // towards initial day being a hot day (p=0.7)
        let transition = array![[0.8, 0.1], [0.1, 0.8]];
        let observation = array![[0.7, 0.2, 0.1], [0.1, 0.2, 0.7]];

        let initial = array![0.3, 0.7];
        let end = array![0.1, 0.1];

        let ices = vec![
            1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 0, 1, 2, 2,
            1, 2, 1, 1,
        ];

        let hmm =
            DiscreteEmissionHMMoptEND::with_float(&transition, &observation, &initial, Some(&end))
                .expect("Dimensions should be consistent");

        // Apply one iteration of baum-welch algorithm and return estimated parameters without
        // updating the model in loco.
        let (pi_hat, transitions_hat, observations_hat, end_hat) = baum_welch(&hmm, &ices);

        // Transform estimated log prob values into prob
        let pi_hat_vec = pi_hat.iter().map(|x| Prob::from(*x)).collect::<Vec<Prob>>();
        let transitions_hat_vec = transitions_hat
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let observations_hat_vec = observations_hat
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let end_hat_vec = end_hat
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();

        // Compare with the example given in Jason Eisner's spreadsheet (http://www.cs.jhu.edu/~jason/papers/#eisner-2002-tnlp)
        let pi_hat_vec_ori = vec![0.0597, 0.9403]
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let transitions_hat_vec_ori = vec![0.8797, 0.1049, 0.0921, 0.8658]
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let observations_hat_vec_ori = vec![0.6765, 0.2188, 0.1047, 0.0584, 0.4251, 0.5165]
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let end_hat_vec_ori = vec![0.0153, 0.0423]
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();

        assert!(pi_hat_vec
            .iter()
            .zip(pi_hat_vec_ori.iter())
            .all(|(a, b)| (*b - Prob::from(0.001) <= *a) && (*a <= *b + Prob::from(0.001))));

        assert!(transitions_hat_vec
            .iter()
            .zip(transitions_hat_vec_ori.iter())
            .all(|(a, b)| (*b - Prob::from(0.001) <= *a) && (*a <= *b + Prob::from(0.001))));

        assert!(observations_hat_vec
            .iter()
            .zip(observations_hat_vec_ori.iter())
            .all(|(a, b)| (*b - Prob::from(0.01) <= *a) && (*a <= *b + Prob::from(0.01))));

        assert!(end_hat_vec
            .iter()
            .zip(end_hat_vec_ori.iter())
            .all(|(a, b)| (*b - Prob::from(0.01) <= *a) && (*a <= *b + Prob::from(0.01))));
    }

    #[test]
    fn test_baum_welch_train_example() {
        // Same Jason Eisner's ice cream example as above with little bias
        // towards initial day being a hot day (p=0.7)
        let transition = array![[0.8, 0.1], [0.1, 0.8]];
        let observation = array![[0.7, 0.2, 0.1], [0.1, 0.2, 0.7]];

        let initial = array![0.3, 0.7];
        let end = array![0.1, 0.1];

        let observations = [vec![
            1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 0, 1, 2, 2,
            1, 2, 1, 1,
        ]];

        let hmm =
            DiscreteEmissionHMMoptEND::with_float(&transition, &observation, &initial, Some(&end))
                .expect("Dimensions should be consistent");

        // Run 10 iterations of Baum-Welch algorithm
        let (pi_hat, transitions_hat, observations_hat, end_hat) =
            hmm.train_baum_welch(&observations, Some(10), None);

        // Transform the obtained LogProb values into Prob.
        let pi_hat_vec = pi_hat.iter().map(|x| Prob::from(*x)).collect::<Vec<Prob>>();
        let transitions_hat_vec = transitions_hat
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let observations_hat_vec = observations_hat
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let end_hat_vec = end_hat
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();

        // Example based on Jason Eisner spreadsheet (http://www.cs.jhu.edu/~jason/papers/#eisner-2002-tnlp)
        let pi_hat_vec_ori = vec![0.0, 1.0]
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let transitions_hat_vec_ori = vec![0.9337, 0.0663, 0.0718, 0.865]
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let observations_hat_vec_ori = vec![0.6407, 0.1481, 0.2112, 1.5e-4, 0.5341, 0.4657]
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();
        let end_hat_vec_ori = vec![0.0, 0.0632]
            .iter()
            .map(|x| Prob::from(*x))
            .collect::<Vec<Prob>>();

        assert!(pi_hat_vec
            .iter()
            .zip(pi_hat_vec_ori.iter())
            .all(|(a, b)| (*b - Prob::from(0.001) <= *a) && (*a <= *b + Prob::from(0.001))));

        assert!(transitions_hat_vec
            .iter()
            .zip(transitions_hat_vec_ori.iter())
            .all(|(a, b)| (*b - Prob::from(0.001) <= *a) && (*a <= *b + Prob::from(0.001))));

        assert!(observations_hat_vec
            .iter()
            .zip(observations_hat_vec_ori.iter())
            .all(|(a, b)| (*b - Prob::from(0.01) <= *a) && (*a <= *b + Prob::from(0.01))));

        assert!(end_hat_vec
            .iter()
            .zip(end_hat_vec_ori.iter())
            .all(|(a, b)| (*b - Prob::from(0.01) <= *a) && (*a <= *b + Prob::from(0.01))));
    }
}
