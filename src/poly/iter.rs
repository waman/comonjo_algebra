use std::{collections::btree_map, iter::Enumerate, slice::Iter, vec};

use crate::algebra::Semiring;

use super::{CoeffsIter, Polynomial};

pub enum IntoNonzeroCoeffsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(Enumerate<vec::IntoIter<C>>),
    Sparse(btree_map::IntoIter<usize, C>),
}

impl<C> Iterator for IntoNonzeroCoeffsIter<C> where C: Semiring {

    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoNonzeroCoeffsIter::Zero() => None,
            IntoNonzeroCoeffsIter::Constant(value) => 
                match value.take() {
                    Some(c) => Some((0, c)),
                    None => None,
                },
            IntoNonzeroCoeffsIter::Dense(vec_iter) => 
                loop {
                    match vec_iter.next() {
                        Some((e, c)) =>  if !c.is_zero() { return Some((e, c)); },
                        None => return None,
                    }
                },
            IntoNonzeroCoeffsIter::Sparse(map_iter) => map_iter.next(),
        }
    }
}

pub enum NonzeroCoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<&'a C>),
    Dense(Enumerate<Iter<'a, C>>),
    Sparse(btree_map::Iter<'a, usize, C>),
}

impl<'a, C> Iterator for NonzeroCoeffsIter<'a, C> where C: Semiring {

    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            NonzeroCoeffsIter::Zero() => None,
            NonzeroCoeffsIter::Constant(value) => 
                match value.take() {
                    Some(c) => Some((0, c)),
                    None => None,
                },
            NonzeroCoeffsIter::Dense(vec_iter) => 
                loop {
                    match vec_iter.next() {
                        Some((e, c)) =>  if !c.is_zero() { return Some((e, c)); },
                        None => return None,
                    }
                },
            NonzeroCoeffsIter::Sparse(map_iter) => map_iter.next().map(|(e, c)| (*e, c)),
        }
    }
}

pub trait IntoCoeffsIterator: IntoIterator where Self: Sized, Self::IntoCoeffsIter: Iterator<Item=Self::Coeff>{
    type Coeff;
    type IntoCoeffsIter;

    fn coeffs(self) -> Self::IntoCoeffsIter;

    #[inline]
    fn nonzero_coeffs(self) -> Self::IntoIter { self.into_iter() }
}

impl<C> IntoCoeffsIterator for Polynomial<C> where C: Semiring {
    type Coeff = C;
    type IntoCoeffsIter = IntoCoeffsIter<C>;

    /// Return <code>impl Iterator<Item=C></code>.
    fn coeffs(self) -> Self::IntoCoeffsIter {
        match self {
            Polynomial::Zero() => IntoCoeffsIter::Zero(),
            Polynomial::Constant(cc) => IntoCoeffsIter::Constant(Some(cc.0)),
            Polynomial::Dense(dc) => IntoCoeffsIter::Dense(dc.0.into_iter()),
            Polynomial::Sparse(sc) => {
                let mut map_iter = sc.0.into_iter();
                IntoCoeffsIter::Sparse { index: 0, current: map_iter.next(), map_iter }
            },
        }
    }
}

pub enum IntoCoeffsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(vec::IntoIter<C>),
    Sparse{
        index: usize,
        current: Option<(usize, C)>,
        map_iter: btree_map::IntoIter<usize, C>
    },
}

impl<C> Iterator for IntoCoeffsIter<C> where C: Semiring {

    type Item = C;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoCoeffsIter::Zero() => None,
            IntoCoeffsIter::Constant(value) => value.take(),
            IntoCoeffsIter::Dense(vec_iter) =>  vec_iter.next(),
            IntoCoeffsIter::Sparse{ index, current, map_iter } => {
                let result = match current {
                    Some(c) => 
                        if c.0 == *index {
                            let r = match map_iter.next() {
                                Some(nxt) => current.replace(nxt),
                                None => current.take(),
                            };
                            Some(r.unwrap().1)
                        } else { 
                            Some(C::zero())
                        }
                    _ => None,
                };
        
                *index += 1;
                return result;
            },
        }
    }
}

impl<'a, C> IntoCoeffsIterator for &'a Polynomial<C> where C: Semiring {
    
    type Coeff = Option<&'a C>;

    type IntoCoeffsIter = CoeffsIter<'a, C>;

    /// Return <code>impl Iterator<Item=Option<Option<&C>>></code>.
    fn coeffs(self) -> Self::IntoCoeffsIter {
        match self {
            Polynomial::Zero() => CoeffsIter::Zero(),
            Polynomial::Constant(cc) => CoeffsIter::Constant(Some(Some(&cc.0))),
            Polynomial::Dense(dc) => CoeffsIter::Dense(dc.0.iter()),
            Polynomial::Sparse(sc) => {
                let mut map_iter = sc.0.iter();
                CoeffsIter::Sparse{ index: 0, current: map_iter.next(), map_iter}
            },
        }
    }
}