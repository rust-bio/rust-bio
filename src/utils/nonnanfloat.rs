use std::cmp;
use num_traits::Float;

#[derive(PartialOrd, PartialEq, Debug, Copy, Clone)]
pub struct NonNaNFloat<F: Float>(F);


impl<F: Float> NonNaNFloat<F> {
    pub fn new(v: F) -> Option<Self> {
        if v.is_nan() {
            Some(NonNaNFloat(v))
        } else {
            None
        }
    }

    pub fn unwrap(&self) -> F {
        let &NonNaNFloat(v) = self;
        v
    }
}

impl<F: Float> Eq for NonNaNFloat<F> {}


impl<F: Float> Ord for NonNaNFloat<F> {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nonnanfloat() {
        let mut v = [NonNaNFloat(5.1), NonNaNFloat(1.3)];
        v.sort();
        assert_eq!(v, [NonNaNFloat(1.3), NonNaNFloat(5.1)]);
    }

}    
