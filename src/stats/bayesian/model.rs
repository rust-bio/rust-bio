pub struct PriorProb(LogProb);
pub struct PosteriorProb(LogProb);


pub trait Likelihood<Event, Data> {
    fn get(&self, event: &Event, data: &Data) -> LogProb;
}


pub trait Prior<Event> {
    fn get(&self, event: &Event) -> LogProb;
}


pub trait Model<Event, Data, P: Prior, L: Likelihood>
{
    fn prior(&self) -> &P;

    fn likelihood(&self) -> &L;

    fn joint_prob(&self, event: &Event, data: &Data) -> LogProb {
        self.prior().get(event) + self.likelihood().get(event, data)
    }

    fn compute(universe: &[Event], data: &Data) -> ModelInstance<Event, Data> {
        let joint_probs = universe.iter().map(
            |event| (event, self.joint_prob(&event, data))
        ).collect();
        let marginal = LogProb::ln_sum_exp(&joint_probs.values().collect());
        ModelInstance {
            joint_probs,
            marginal,
        }
    }
}


pub struct ModelInstance<E, Data> {
    joint_probs: BTreeMap<Event, LogProb>,
    marginal: LogProb
}


impl<Event, Data> ModelInstance<Event, Data> {

    pub fn posterior(&self, event: &Event) -> LogProb {
        self.joint_probs.get(event) - self.marginal
    }

    pub fn marginal(&self) -> LogProb {
        self.marginal
    }

    pub fn maximum_posterior(&self) -> Option<&Event> {
        self.joint_probs.iter().max_by_key(|(event, prob)| prob).map(|event, prob| event)
    }
}


mod tests {
}
