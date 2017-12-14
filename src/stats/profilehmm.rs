// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Profile Hidden Markov Model. Can be used for construct profile of multiple sequence alignments

use stats::LogProb;


pub struct ProfileHMM {
    /// number of possible observations
    pub observation_count: usize,
    /// probability that x is the initial state
    pub initial_states_prob: Vec<LogProb>,
    /// probability of transiting from state x to y
    pub state_transitions: Vec<Vec<LogProb>>,
    /// probability of emiting observation y at state x
    pub emission_matrix: Vec<Vec<LogProb>>
}

impl ProfileHMM {
    pub fn new() -> Self {
        ProfileHMM {
            observation_count: 0,
            initial_states_prob: Vec::new(),
            state_transitions: Vec::new(),
            emission_matrix: Vec::new()
        }
    }

    pub fn viterbi(&self, observations : Vec<usize>)->Vec<usize> {
        //Viterbi Reader's Guide:
        //The term "Winning" transition means the transition selected to be taken by the engine
        //The term scoring is to multiply a transition by it's correspondance i the Emission Matrix
        //path_finder[k][i]=MAX(all states:l) {Score[l][i-1] X weight (l,k,i-1)} i:Column, k:Node Serial
        let mut path_finder:Vec<Vec<LogProb>>=Vec::new(); //path_finder is the memory used for dynamic programming

        let mut result_states:Vec<usize>=Vec::new(); //Output Vector
        let mut previous_state : Vec<Vec<usize>> = Vec::new();
        let trans_count=self.state_transitions[0].len();; //Number of transition states
        //Loop handling the first column
        let mut path_segment:Vec<LogProb>=Vec::new(); //In order to fill a Vector or Vectors we need a temp vector to push
        previous_state.push(Vec::new());
        for i in 0..trans_count {
                //The "winning" segment is certainly the ith
            path_segment.push(self.initial_states_prob[i] + self.emission_matrix[i][observations[0]]);
        }
        path_finder.push(path_segment); //<= Like here

        //Loop handling the other columns
        //For each state
        for i in 1..observations.len() {
            //The parameter is pushed here and updated later to ensure the equality of parameters and request for observed states
            path_segment=Vec::new(); //renewing path_segment
            previous_state.push(Vec::new());
            for j in 0..trans_count { //each transition-*
                path_segment.push(LogProb::ln_zero());
                previous_state[i].push(0);
                for k in 0..trans_count { //*-is compared to the other transitions
                    //value of the transition=the probability of the previous value X transition probability to this value
                    let tmp= path_finder[i-1][k] + self.state_transitions[k][j];
                    if tmp > path_segment[j] {
                        previous_state[i][j] = k; //Registering Path For result_states
                        path_segment[j]=tmp; //Registering Path Value for further path finding
                    }
                }
                path_segment[j] = path_segment[j] + self.emission_matrix[j][observations[i]]; //Scoring the "winning" transition
            }
            path_finder.push(path_segment); //2D Code fill process
        }
        let mut max:LogProb=LogProb::ln_zero();
        result_states = vec![0; observations.len()]; //pushing a value to keep editing it
        for i in 0..trans_count {//Finding Maximum path by iteratively searching the final entry in the path_finder
            if path_finder[observations.len()-1][i] >= max {
                result_states[observations.len()-1]=i;
                max=path_finder[observations.len()-1][i];
            }
        }
        for i in 2..observations.len()+1 {
            result_states[observations.len()-i]= previous_state[observations.len()-i+1][result_states[observations.len()-i+1]];
        }
        result_states
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use stats::{Prob, LogProb};
}