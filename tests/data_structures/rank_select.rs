use proptest::prelude::*;
use std::collections::HashMap;
use std::hash::Hash;

// Naive Rank datastructure.
//
// Each index in the `ranks` vector is the index in the sequence.
// At each of the indexes the overall rank of all values seen up to
// that point is stored in a HashMap.
#[derive(Default)]
struct Rank<T: Hash + Eq + Clone> {
    ranks: Vec<HashMap<T, usize>>,
    len: usize,
}

// Naive Select data structure.
//
// All values in the sequence are keys in the HashMap.
// Each value has an associated vector where the indexes are the ranks
// and the value at each index is the index in the sequence at which the
// rank occured.
#[derive(Default)]
struct Select<T: Hash + Eq + Clone> {
    select: HashMap<T, Vec<usize>>,
    len: usize,
}

impl<T: Hash + Eq + Clone> Select<T> {
    fn push(&mut self, val: T) {
        self.select
            .entry(val)
            .or_insert_with(Vec::new)
            .push(self.len);
        self.len += 1;
    }

    fn select(&self, rank: usize, val: &T) -> Option<u64> {
        if rank > 0 {
            if let Some(ranks) = self.select.get(val) {
                if let Some(idx) = ranks.get(rank - 1) {
                    return Some(*idx as u64);
                }
            }
        }
        return None;
    }
}

impl<T: Hash + Eq + Clone> Rank<T> {
    fn push(&mut self, val: T) {
        if let Some(prev) = self.ranks.last() {
            let mut nxt: HashMap<_, _> = prev.clone();
            *nxt.entry(val).or_insert(0) += 1;
            self.ranks.push(nxt.clone());
        } else {
            let mut nxt = HashMap::new();
            nxt.insert(val, 1);
            self.ranks.push(nxt);
        }
        self.len += 1;
    }

    fn get_rank(&self, idx: usize, val: T) -> Option<u64> {
        if idx < self.len {
            let rank = self
                .ranks
                .get(idx)
                .map(|x| *x.get(&val).unwrap_or(&0))
                .unwrap_or(0);
            Some(rank as u64)
        } else {
            None
        }
    }
}

proptest! {
    #[test]
    fn rank_same_as_naive(input: Vec<bool>, k in 1usize..1000usize) {
        use bio::data_structures::rank_select::RankSelect;
        use bv::BitVec;
        let mut bv = BitVec::new();
        let mut naive = Rank::default();
        for &i in &input {
            naive.push(i);
            bv.push(i)
        }
        let rs = RankSelect::new(bv, k);
        for idx in 0..=input.len() {
            assert_eq!(naive.get_rank(idx, false), rs.rank_0(idx as u64));
            assert_eq!(naive.get_rank(idx, true), rs.rank_1(idx as u64));
        }
    }

    #[test]
    fn select_same_as_naive(input: Vec<bool>, k in 1usize..1000usize) {
        use bio::data_structures::rank_select::RankSelect;
        use bv::BitVec;
        let mut bv = BitVec::new();
        let mut naive = Select::default();
        for &i in &input {
            naive.push(i);
            bv.push(i)
        }
        let rs = RankSelect::new(bv, k);
        for idx in 0..=input.len() {
            assert_eq!(naive.select(idx, &false), rs.select_0(idx as u64));
            assert_eq!(naive.select(idx, &true), rs.select_1(idx as u64));
        }
    }
}
