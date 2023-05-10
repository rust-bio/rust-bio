//use bio_types::genome::Length;

use crate::utils::Text;
//use crate::alphabets;


/// In place implementation of scan over a slice.
pub fn transcription(dna_str : Text) -> Text {
    //let alphab = alphabets::dna::alphabet();
    let mut rna_str : Text = vec![];
    for i in &dna_str {
        let nuc_dna = i;
        match nuc_dna{
            84=>rna_str.push(85),
            116=>rna_str.push(117),
            _=>rna_str.push(*nuc_dna),
        }
    }
    rna_str
}




#[cfg(test)]
mod tests {
    use super::*;
    use std::ops::Deref;

    fn print_sequence<Item: Deref<Target = u8>, T: IntoIterator<Item = Item>>(sequence: T) {
        for c in sequence {
            println!("{}", *c);
        }
    }

    #[test]
    fn test_1() {
        let seq_ex = b"ACGT";
        let seq_vec = &vec![b'A',b'C'];
        print_sequence(seq_ex.iter().step_by(1));
        print_sequence(seq_vec);
    }

    #[test]
    fn test_2(){
        let seq_vec : Text= vec![b'A',b'C',b'T'];
        let new_arn = transcription(seq_vec);
        assert_eq!(new_arn,vec![b'A',b'C',b'U']);
    }
}
