use std::{iter::Enumerate, slice::Iter};

use num::Num;

use crate::polynomial::PolyIter;

pub struct DenseContent<C: Num>(pub(crate) Vec<C>);

impl<C: Num> DenseContent<C> {

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(n)
    }

    pub fn iter<'a>(&'a self) -> PolyIter<'a, C> {
        PolyIter::Dense(DenseIter { vec_iter: self.0.iter().enumerate() })
    }
}

pub struct DenseIter<'a, C: Num> {
    vec_iter: Enumerate<Iter<'a, C>>
}

impl<'a, C: Num> Iterator for DenseIter<'a, C> {

    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.vec_iter.next() {
                Some((e, c)) =>  if !c.is_zero() { return Some((e, c)); },
                None => return None,
            }
        }
    }
}