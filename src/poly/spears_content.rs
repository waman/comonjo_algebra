use std::collections::{btree_map::Iter, BTreeMap};
use num::Num;

use crate::polynomial::PolyIter;

pub(crate) struct SpearsContent<C: Num>(pub(crate) BTreeMap<usize, C>);

impl<C: Num> SpearsContent<C> {

    pub fn degree(&self) -> usize {
        *self.0.keys().max().unwrap()
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(&n)
    }

    pub fn iter<'a>(&'a self) -> PolyIter<'a, C> {
        PolyIter::Spears(SpearsIter { map_iter: self.0.iter() })
    }
}

pub(crate) struct SpearsIter<'a, C> {
    map_iter: Iter<'a, usize, C>

}

impl<'a, C: Num> Iterator for SpearsIter<'a, C> {

    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        self.map_iter.next().map(|(e, c)| (*e, c))
    }
}