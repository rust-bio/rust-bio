use std::cmp::Ord;
use std::collections::BTreeMap;

use itertools::Itertools;
use ordered_float::NotNan;

use stats::LogProb;

/// Likelihood model.
pub trait Likelihood {
    type Event;
    type Data;

    fn get(&self, event: &Self::Event, data: &Self::Data) -> LogProb;
}

/// Prior model.
pub trait Prior {
    type Event;

    fn get(&self, event: &Self::Event) -> LogProb;
}

/// A bayesian model.
pub trait Model<Event, Data>
where
    Event: Ord + Clone,
{
    /// Joint probability (prior * likelihood).
    fn joint_prob(&self, event: &Event, data: &Data) -> LogProb;

    /// Compute an instance of the model for the given event universe and data.
    fn compute(&self, universe: &[Event], data: &Data) -> ModelInstance<Event> {
        let joint_probs: BTreeMap<Event, LogProb> = universe
            .into_iter()
            .cloned()
            .map(|event| {
                let p = self.joint_prob(&event, data);
                (event, p)
            })
            .collect();
        let marginal = LogProb::ln_sum_exp(&joint_probs.values().cloned().collect_vec());
        ModelInstance {
            joint_probs,
            marginal,
        }
    }
}

/// Instance of a model for given data and event universe.
/// From the instance, posterior, marginal and MAP can be computed.
pub struct ModelInstance<Event>
where
    Event: Ord,
{
    joint_probs: BTreeMap<Event, LogProb>,
    marginal: LogProb,
}

impl<Event> ModelInstance<Event>
where
    Event: Ord,
{
    /// Posterior probability of given event.
    pub fn posterior(&self, event: &Event) -> Option<LogProb> {
        self.joint_probs.get(event).map(|p| p - self.marginal)
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

mod tests {}
