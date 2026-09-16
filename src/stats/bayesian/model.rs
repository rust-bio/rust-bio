// Copyright 2019 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A trait system for Bayesian statistical modelling.
//!
//! A [`Model`] is assembled from three pieces that you implement for your own
//! event and data types: a [`Likelihood`] (`Pr(data | event)`), a [`Prior`]
//! (`Pr(event)`) and a [`Posterior`], which turns base events into the
//! (possibly aggregated) events you actually want a posterior probability
//! for. When there is no aggregation to perform, `Posterior::compute` can
//! simply forward to the provided `joint_prob` closure, as shown below.
//!
//! # Example
//!
//! Here, we infer the bias (`Pr(heads)`) of a coin from a series of flips,
//! using a uniform prior over a small discretized universe of candidate
//! biases and a binomial likelihood.
//!
//! ```
//! use approx::assert_relative_eq;
//! use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
//! use bio::stats::LogProb;
//! use ordered_float::NotNan;
//!
//! // The observed data: how many heads and tails we saw.
//! struct Flips {
//!     heads: u32,
//!     tails: u32,
//! }
//!
//! // A hypothesis about the coin, i.e. Pr(heads). We use `NotNan` so that
//! // biases can be hashed and compared, as required by `Model`.
//! type Bias = NotNan<f64>;
//!
//! struct BinomialLikelihood;
//!
//! impl Likelihood for BinomialLikelihood {
//!     type Event = Bias;
//!     type Data = Flips;
//!
//!     fn compute(&self, event: &Bias, data: &Flips, _payload: &mut ()) -> LogProb {
//!         let p = **event;
//!         LogProb(f64::from(data.heads) * p.ln() + f64::from(data.tails) * (1. - p).ln())
//!     }
//! }
//!
//! // A uniform prior over `n` candidate biases.
//! struct UniformPrior {
//!     n: usize,
//! }
//!
//! impl Prior for UniformPrior {
//!     type Event = Bias;
//!
//!     fn compute(&self, _event: &Bias) -> LogProb {
//!         LogProb((1. / self.n as f64).ln())
//!     }
//! }
//!
//! // With nothing to aggregate, the posterior of an event is just its joint probability.
//! struct DirectPosterior;
//!
//! impl Posterior for DirectPosterior {
//!     type Event = Bias;
//!     type BaseEvent = Bias;
//!     type Data = Flips;
//!
//!     fn compute<F: FnMut(&Bias, &Flips) -> LogProb>(
//!         &self,
//!         event: &Bias,
//!         data: &Flips,
//!         joint_prob: &mut F,
//!     ) -> LogProb {
//!         joint_prob(event, data)
//!     }
//! }
//!
//! // Candidate biases from 0.1 to 0.9 in steps of 0.1.
//! let universe: Vec<Bias> = (1..10).map(|i| NotNan::new(i as f64 / 10.).unwrap()).collect();
//! let model = Model::new(BinomialLikelihood, UniformPrior { n: universe.len() }, DirectPosterior);
//!
//! // Out of 10 flips, 8 came up heads.
//! let data = Flips { heads: 8, tails: 2 };
//! let instance = model.compute(universe, &data);
//!
//! // The candidate bias closest to the observed frequency (0.8) should be most probable.
//! let map = instance.maximum_posterior().unwrap();
//! assert_relative_eq!(**map, 0.8);
//! ```

use std::cmp::Eq;
use std::collections::HashMap;
use std::hash::Hash;
use std::marker::PhantomData;

use itertools::Itertools;
use ordered_float::NotNan;

use crate::stats::LogProb;

pub type JointProbUniverse<Event> = HashMap<Event, LogProb>;

/// Likelihood model.
pub trait Likelihood<Payload = ()> {
    type Event;
    type Data;

    /// Compute likelihood of event given the data. Optionally, the passed payload can be used
    /// to e.g., cache intermediate results. One payload corresponds to one model instance.
    fn compute(&self, event: &Self::Event, data: &Self::Data, payload: &mut Payload) -> LogProb;
}

/// Prior model.
pub trait Prior {
    type Event;

    fn compute(&self, event: &Self::Event) -> LogProb;
}

/// Posterior model.
pub trait Posterior {
    type Event;
    type BaseEvent;
    type Data;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb;
}

/// Bayesian model, consisting of a prior, a posterior and a likelihood model.
/// Thereby, `Payload` is a custom payload of the model instance.
/// This can be used to define custom caching mechanisms. See
/// [here](https://github.com/varlociraptor/varlociraptor/blob/694e994547e8f523e5b0013fdf951b694f3870fa/src/model/modes/generic.rs#L200)
/// for an example.
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct Model<L, Pr, Po, Payload = ()>
where
    L: Likelihood<Payload>,
    Pr: Prior,
    Po: Posterior,
    Payload: Default,
{
    likelihood: L,
    prior: Pr,
    posterior: Po,
    payload: PhantomData<Payload>,
}

impl<Event, PosteriorEvent, Data, L, Pr, Po, Payload> Model<L, Pr, Po, Payload>
where
    Payload: Default,
    Event: Hash + Eq + Clone,
    PosteriorEvent: Hash + Eq + Clone,
    L: Likelihood<Payload, Event = Event, Data = Data>,
    Pr: Prior<Event = Event>,
    Po: Posterior<BaseEvent = Event, Event = PosteriorEvent, Data = Data>,
{
    /// Create new instance.
    pub fn new(likelihood: L, prior: Pr, posterior: Po) -> Self {
        Model {
            likelihood,
            prior,
            posterior,
            payload: PhantomData,
        }
    }

    pub fn likelihood(&self) -> &L {
        &self.likelihood
    }

    pub fn likelihood_mut(&mut self) -> &mut L {
        &mut self.likelihood
    }

    pub fn prior(&self) -> &Pr {
        &self.prior
    }

    pub fn prior_mut(&mut self) -> &mut Pr {
        &mut self.prior
    }

    pub fn posterior(&self) -> &Po {
        &self.posterior
    }

    pub fn posterior_mut(&mut self) -> &mut Po {
        &mut self.posterior
    }

    /// Calculate joint probability, i.e. `Pr(event) * Pr(data | event)`.
    fn joint_prob(&self, event: &Event, data: &Data, payload: &mut Payload) -> LogProb {
        self.prior.compute(event) + self.likelihood.compute(event, data, payload)
    }

    /// Compute model for a given universe of events.
    ///
    /// # Complexity
    ///
    /// Calls `Posterior::compute` once per event in `universe`; each such call may in turn
    /// invoke the `joint_prob` closure (and hence `Likelihood::compute` and `Prior::compute`)
    /// an arbitrary number of times, depending on the `Posterior` implementation. Every
    /// distinct base event that is evaluated is cached, so `joint_prob` is never recomputed
    /// for the same base event. Requires `O(n + b)` additional space, where `n = universe.len()`
    /// and `b` is the number of distinct base events visited.
    pub fn compute<U: IntoIterator<Item = PosteriorEvent>>(
        &self,
        universe: U,
        data: &Data,
    ) -> ModelInstance<Event, PosteriorEvent> {
        let mut joint_probs = HashMap::new();
        let mut payload = Payload::default();
        let (posterior_probs, marginal) = {
            let mut joint_prob = |event: &Event, data: &Data| {
                let p = self.joint_prob(event, data, &mut payload);
                joint_probs.insert(event.clone(), p);
                p
            };

            let posterior_probs: HashMap<PosteriorEvent, LogProb> = universe
                .into_iter()
                .map(|event| {
                    let p = self.posterior.compute(&event, data, &mut joint_prob);
                    (event, p)
                })
                .collect();
            let marginal = LogProb::ln_sum_exp(&posterior_probs.values().cloned().collect_vec());

            (posterior_probs, marginal)
        };

        ModelInstance {
            joint_probs,
            posterior_probs,
            marginal,
        }
    }

    /// Compute model via the exploration of the marginal distribution of the data.
    ///
    /// # Complexity
    ///
    /// Same characteristics as [`Model::compute`], except that the events explored (and thus
    /// the number of calls to `joint_prob`) are driven by the given `Marginal` implementation
    /// rather than an explicit universe.
    pub fn compute_from_marginal<M>(
        &self,
        marginal: &M,
        data: &Data,
    ) -> ModelInstance<Event, PosteriorEvent>
    where
        M: Marginal<Data = Data, Event = PosteriorEvent, BaseEvent = Event>,
    {
        let mut joint_probs = HashMap::new();
        let mut posterior_probs = HashMap::new();
        let mut payload = Payload::default();
        let marginal = {
            let mut joint_prob = |event: &Event, data: &Data| {
                let p = self.joint_prob(event, data, &mut payload);
                joint_probs.insert(event.clone(), p);
                p
            };

            let mut joint_prob_posterior = |event: &PosteriorEvent, data: &Data| {
                let p = self.posterior.compute(event, data, &mut joint_prob);
                posterior_probs.insert(event.clone(), p);
                p
            };

            marginal.compute(data, &mut joint_prob_posterior)
        };

        ModelInstance {
            joint_probs,
            posterior_probs,
            marginal,
        }
    }
}

/// A trait for the exploration of the marginal distribution of the data.
pub trait Marginal {
    type Event;
    type BaseEvent;
    type Data;

    fn compute<F: FnMut(&Self::Event, &Self::Data) -> LogProb>(
        &self,
        data: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb;
}

/// Instance of a model for given data and event universe.
/// From the instance, posterior, marginal and MAP can be computed.
#[derive(Default, Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct ModelInstance<Event, PosteriorEvent>
where
    Event: Hash + Eq,
    PosteriorEvent: Hash + Eq,
{
    joint_probs: HashMap<Event, LogProb>,
    posterior_probs: HashMap<PosteriorEvent, LogProb>,
    marginal: LogProb,
}

impl<Event, PosteriorEvent> ModelInstance<Event, PosteriorEvent>
where
    Event: Hash + Eq,
    PosteriorEvent: Hash + Eq,
{
    /// Posterior probability of given event.
    ///
    /// Runs in `O(1)` average time, a single lookup in the underlying `HashMap`.
    pub fn posterior(&self, event: &PosteriorEvent) -> Option<LogProb> {
        self.posterior_probs.get(event).map(|p| p - self.marginal)
    }

    /// Marginal probability.
    ///
    /// Runs in `O(1)` time; the marginal is computed once by [`Model::compute`] and simply
    /// returned here.
    pub fn marginal(&self) -> LogProb {
        self.marginal
    }

    /// Maximum a posteriori estimate.
    ///
    /// Since dividing by the (fixed) marginal does not change which event maximizes the
    /// probability, this compares joint probabilities directly rather than posteriors.
    /// Runs in `O(b)` time, where `b` is the number of distinct base events cached by
    /// [`Model::compute`].
    pub fn maximum_posterior(&self) -> Option<&Event> {
        self.joint_probs
            .iter()
            .max_by_key(|(_, prob)| NotNan::new(***prob).unwrap())
            .map(|(event, _)| event)
    }

    /// Event posteriors sorted in descending order.
    ///
    /// Runs in `O(b log b)` time and requires `O(b)` additional space, where `b` is the number
    /// of distinct base events cached by [`Model::compute`], due to sorting.
    pub fn event_posteriors(&self) -> impl Iterator<Item = (&Event, LogProb)> {
        self.joint_probs
            .iter()
            .map(|(event, prob)| (event, prob - self.marginal))
            .sorted_by_key(|(_, prob)| -NotNan::new(**prob).unwrap())
    }
}

impl<PosteriorEvent> ModelInstance<NotNan<f64>, PosteriorEvent>
where
    PosteriorEvent: Hash + Eq,
{
    pub fn expected_value(&self) -> NotNan<f64> {
        self.joint_probs
            .iter()
            .map(|(event, prob)| *event * NotNan::new(**prob).unwrap())
            .fold(NotNan::default(), |s, e| s + e)
    }
}

mod tests {}
