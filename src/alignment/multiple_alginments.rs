


use stats::{profilehmm, LogProb, Prob};
use alphabets::Alphabet;

struct MultipleAlignment {
    /// alphabet used
    alphabet : Alphabet,
    /// length of alignments
    length : usize,
    /// aligned sequences
    alignments: Vec<Vec<u8>>,
    /// consensus sequence
    consensus_sequence: Vec<u8>,
    /// underlining ProfileHMM
    profile_hmm : profilehmm::ProfileHMM,
}

impl MultipleAlignment {
    pub fn new(used_alphabet : Alphabet) -> Self {
        MultipleAlignment {
            alphabet: used_alphabet,
            length: 0,
            alignments: Vec::new(),
            consensus_sequence: Vec::new(),
            profile_hmm: profilehmm::ProfileHMM::new()
        }
    }

    /// Construct multiple alignments based on Hughey R, Krogh A (1996)
    /// "Hidden Markov models for sequence analysis: extension and analysis
    /// of the basic method". CABIOS. 12 (2): 95–107. doi:10.1093/bioinformatics/12.2.95
    pub fn construct(&mut self, sequences : Vec<Vec<u8>>) {
        // estimate HMM of the given data

        // get alignments using viterbi

    }

    /// Construct profile HMM if the alignment is already known
    /// based on Byung-Jun Yoon () "Hidden Markov Models and their Applications in Biological
    /// Sequence Analysis". Curr Genomics. 10(6): 402–415. doi:  10.2174/138920209789177575
    pub fn construct_profile_hmm(&self) -> profilehmm::ProfileHMM {
        // construct base frequency
        let mut base_frequency : Vec<Vec<LogProb>> = Vec::new();
        let mut symbols : Vec<u8> = Vec::new();
        let column_count : usize = self.alphabet.len();
        for symbol in &self.alphabet.symbols {
            symbols.push(symbol as u8);
        }
        for index in 0..self.length {
            let mut column_probability : Vec<LogProb> = vec![LogProb::ln_zero(); column_count];
            let mut column_observation_count : Vec<usize> = vec![0; column_count + 1];
            for alignment in &self.alignments {
                let symbol = alignment.get(index).unwrap();
                if self.alphabet.symbols.contains(*symbol as usize) {
                    column_observation_count[symbols.iter().position(|&s| s == *symbol).unwrap()] += 1;
                } else {
                    column_observation_count[column_count] += 1;
                }
            }
            for symbol_index in 0..self.alphabet.len() {
                column_probability[symbol_index] = LogProb::from(Prob((column_observation_count[symbol_index] as f64) / (column_count as f64)));
            }
            base_frequency.push(column_probability);
        }
        // construct ungapped
        let mut profile_hmm = profilehmm::ProfileHMM::new();
        let state_count = self.consensus_sequence.len() * 3 + 1; // match - delete - insert
        profile_hmm.observation_count = self.alphabet.len();
        profile_hmm.initial_states_prob = vec![LogProb::ln_zero(); state_count];
        profile_hmm.state_transitions = vec![vec![LogProb::ln_zero(); state_count]; state_count];
        let half = LogProb::from(Prob(0.5f64));
        let mut consensus_sequence_index : usize = 0;
        let mut current_symbol = self.consensus_sequence[0];
        for index in 0..self.length {
            let probability = base_frequency[index][symbols.iter().position(|&s| s == current_symbol).unwrap()];
            if probability >= half {
                profile_hmm.emission_matrix.push(base_frequency[index].clone());
                if consensus_sequence_index > 0 {
                    profile_hmm.state_transitions[consensus_sequence_index - 1][consensus_sequence_index] = probability;
                } else {
                    profile_hmm.initial_states_prob[consensus_sequence_index] = probability;
                }
                consensus_sequence_index += 1;
                if consensus_sequence_index == self.consensus_sequence.len() {break;}
                current_symbol = self.consensus_sequence[consensus_sequence_index];
            }
        }
        assert!(consensus_sequence_index == self.consensus_sequence.len());
        // construct profile (delete and insert)

        profile_hmm
    }

    /// Get score for sequence in the alignment
    pub fn search(sequence : Vec<u8>) -> LogProb {
        LogProb::from(Prob(1.0))
    }



    pub fn add_sequence(&mut self, sequence : Vec<u8>) {
        // get original sequences
        let mut sequences : Vec<Vec<u8>> = Vec::new();
        {
            let alignments = &self.alignments;
            for alignment in alignments {
                let mut ungapped = alignment.clone();
                ungapped.retain(|&symbol| self.alphabet.symbols.contains(symbol as usize));
                if ungapped.len() > 0 {
                    sequences.push(ungapped);
                }
            }
        }
        sequences.push(sequence);
        // construct again
        self.construct(sequences);
    }
}