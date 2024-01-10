use bio::alignment::{poa::*, TextSlice, pairwise::Scoring};

fn main() {
    let seqs = vec![
            "AACCCCTGTCAAAAA".to_string(),
            "TGT".to_string(),
                ];
    let mut seqs_bytes = vec![];
    for seq in seqs.iter() {
        seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
    }

    let mut scoring = Scoring::new(-12, -6, |a: u8, b: u8| if a == b { 3 } else { -4 });
    scoring.xclip_prefix = 0;
    scoring.yclip_prefix = 0;
    scoring.xclip_suffix = 0;
    scoring.yclip_suffix = 0;
    let mut aligner = Aligner::new(scoring, &seqs_bytes[0]);
    for seq in seqs_bytes.iter().skip(1) {
        aligner.custom(seq);
    }
    println!("here1");
    let test_ops = aligner.alignment();
    println!("here2");
    for op in test_ops.operations {
        match op {
            AlignmentOperation::Match(None) => {
                println!("Match none");
            }
            AlignmentOperation::Match(Some((_, p))) => {
                println!("Match some");
            }
            AlignmentOperation::Ins(None) => {
                println!("Ins none");
            }
            AlignmentOperation::Ins(Some(_)) => {
                println!("Ins some");
            }
            AlignmentOperation::Del(_) => {
                println!("Del");
            }
            AlignmentOperation::Xclip(p) => {
                println!("X clip");
            }
            AlignmentOperation::Yclip(p) => {
                println!("Y clip");
                
            }
        }
    }
    //let consensus = aligner.consensus();
}