use std::cmp::Ord;
use std::collections::BTreeMap;

use itertools::Itertools;
use ordered_float::NotNan;

use stats::LogProb;

pub type JointProbUniverse<Event> = BTreeMap<Event, LogProb>;

/// Likelihood model.
pub trait Likelihood {
    type Event;
    type Data;

    fn compute(&self, event: &Self::Event, data: &Self::Data) -> LogProb;
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

pub struct Model<L, Pr, Po>
where
    L: Likelihood,
    Pr: Prior,
    Po: Posterior,
{
    likelihood: L,
    prior: Pr,
    posterior: Po,
}

impl<Event, PosteriorEvent, Data, L, Pr, Po> Model<L, Pr, Po>
where
    Event: Ord + Clone,
    PosteriorEvent: Ord + Clone,
    L: Likelihood<Event = Event, Data = Data>,
    Pr: Prior<Event = Event>,
    Po: Posterior<BaseEvent = Event, Event = PosteriorEvent, Data = Data>,
{
    pub fn new(likelihood: L, prior: Pr, posterior: Po) -> Self {
        Model {
            likelihood,
            prior,
            posterior,
        }
    }

    pub fn compute<U: IntoIterator<Item=PosteriorEvent>>(
        &self,
        universe: U,
        data: &Data,
    ) -> ModelInstance<Event, PosteriorEvent> {
        let mut joint_probs = BTreeMap::new();
        let (posterior_probs, marginal) = {
            let mut joint_prob = |event: &Event, data: &Data| {
                let p = self.prior.compute(event) + self.likelihood.compute(event, data);
                joint_probs.insert(event.clone(), p);
                p
            };

            let posterior_probs: BTreeMap<PosteriorEvent, LogProb> = universe
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
}

/// Instance of a model for given data and event universe.
/// From the instance, posterior, marginal and MAP can be computed.
pub struct ModelInstance<Event, PosteriorEvent>
where
    Event: Ord,
    PosteriorEvent: Ord,
{
    joint_probs: BTreeMap<Event, LogProb>,
    posterior_probs: BTreeMap<PosteriorEvent, LogProb>,
    marginal: LogProb,
}

impl<Event, PosteriorEvent> ModelInstance<Event, PosteriorEvent>
where
    Event: Ord,
    PosteriorEvent: Ord,
{
    /// Posterior probability of given event.
    pub fn posterior(&self, event: &PosteriorEvent) -> Option<LogProb> {
        self.posterior_probs.get(event).map(|p| p - self.marginal)
    }

    /// Marginal probability.
    pub fn marginal(&self) -> LogProb {
        self.marginal
    }

    /// Maximum a posteriori estimate.
    pub fn maximum_posterior(&self) -> Option<&Event> {
        self.joint_probs
            .iter()
            .max_by_key(|(_, prob)| NotNan::new(***prob).unwrap())
            .map(|(event, _)| event)
    }
}

impl<PosteriorEvent> ModelInstance<NotNan<f64>, PosteriorEvent>
where
    PosteriorEvent: Ord,
{
    pub fn expected_value(&self) -> NotNan<f64> {
        self.joint_probs
            .iter()
            .map(|(event, prob)| *event * NotNan::new(**prob).unwrap())
            .fold(NotNan::default(), |s, e| s + e)
    }
}

mod tests {}
