//use bio_types::genome::Length;

use crate::utils::Text;
//use crate::alphabets;



/// # Example
/// ```
/// use bio::blast;
/// let dna_string = vec![b'C',b'T',b'A'];
/// let rna_string = blast::transcription(dna_string);
/// assert!(rna_string == vec![b'C',b'U',b'A']);
///  
/// let prot_string = blast::traduction(rna_string);
/// assert!(prot_string == vec![b'L']);
/// ```

/// Transcription of a dna (Text :type) into the rna version
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

/// traduction of a rna (Text :type) into the aminoacid strand
pub fn traduction(rna_str : Text) -> Text{
    let mut codo : Text = vec![];
    let mut protein_str : Text = vec![];
    for i in &rna_str {
        if codo.len() < 3 {
            codo.push(*i);
        }
        if codo.len() == 3{
            match codo[0]{
                b'C'=> match codo[1]{
                    b'A'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'H'),
                        _=>protein_str.push(b'Q'),
                    },
                    b'U'=>protein_str.push(b'L'),
                    b'C'=>protein_str.push(b'P'),
                    _=>protein_str.push(b'R'), //G
                },
                b'G'=> match codo[1]{
                    b'A'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'D'),
                        _=>protein_str.push(b'E'),
                    },
                    b'U'=>protein_str.push(b'V'),
                    b'C'=>protein_str.push(b'A'),
                    _=>protein_str.push(b'G'), //G
                },
                b'A'=> match codo[1]{
                    b'A'=> match codo[2]{
                        b'G' | b'A' =>protein_str.push(b'K'),
                        _=>protein_str.push(b'N'),
                    },
                    b'G'=> match codo[2]{
                        b'G' | b'A' =>protein_str.push(b'R'),
                        _=>protein_str.push(b'S'),
                    },
                    b'U'=>match codo[2]{
                        b'G' =>protein_str.push(b'M'),
                        _=>protein_str.push(b'I'),
                    },
                    _=>protein_str.push(b'P'), //C
                },
                b'U'=> match codo[1]{
                    b'U'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'F'),
                        _=>protein_str.push(b'L'),
                    },
                    b'A'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'Y'),
                        _=>protein_str.push(b'O'),
                    },
                    b'G'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'C'),
                        b'G'=>protein_str.push(b'W'),
                        _=>protein_str.push(b'O'),
                    },
                    _=>protein_str.push(b'S'), //C
                },
                _=>protein_str.push(b'X'),
            }

            codo.clear();
        }
    }
    
    protein_str
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
    fn printing_strand() {
        let seq_ex = b"ACGT";
        let seq_vec = &vec![b'A',b'C'];
        print_sequence(seq_ex.iter().step_by(1));
        print_sequence(seq_vec);
    }

    #[test]
    fn transcript_simple_dna(){
        let seq_vec : Text= vec![b'A',b'C',b'T'];
        let new_arn = transcription(seq_vec);
        assert_eq!(new_arn,vec![b'A',b'C',b'U']);
    }

    #[test]
    fn traduction_one_codon_protein(){
        let seq_arn : Text = vec![b'U',b'U',b'U'];
        let new_prot : Text = traduction(seq_arn);
        assert_eq!(new_prot,vec![b'F']);
    }
}
