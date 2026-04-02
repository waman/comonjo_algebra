use std::{collections::btree_map, iter, slice, vec};

use crate::{algebra::Semiring, poly::Polynomial};

fn const_iter_size_hint<'a, C>(c: &'a Option<C>) -> (usize, Option<usize>) {
    if c.is_some() { (1, Some(1)) } else { (0, Some(0)) }
}

//********** Into Nonzero Coeffs Iterator (consumes self and returns value itself) **********/
pub(crate) fn into_nonzero_coeffs_iter<C>(p: Polynomial<C>) -> IntoNonzeroCoeffsIter<C> where C: Semiring {
    match p {
        Polynomial::Zero()       => IntoNonzeroCoeffsIter::Zero(),
        Polynomial::Constant(cc) => IntoNonzeroCoeffsIter::Constant(Some(cc.0)),
        Polynomial::Dense(dc)    => IntoNonzeroCoeffsIter::Dense(dc.0.into_iter().enumerate()),
        Polynomial::Sparse(sc)   => IntoNonzeroCoeffsIter::Sparse(sc.0.into_iter()),
    }
}

pub enum IntoNonzeroCoeffsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(iter::Enumerate<vec::IntoIter<C>>),
    Sparse(btree_map::IntoIter<usize, C>),
}

impl<C> Iterator for IntoNonzeroCoeffsIter<C> where C: Semiring {

    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoNonzeroCoeffsIter::Zero() => None,
            IntoNonzeroCoeffsIter::Constant(c) => 
                if let Some(c) = c.take() { Some((0, c)) } else { None },
            IntoNonzeroCoeffsIter::Dense(it) => 
                loop {
                    match it.next() {
                        Some((i, c)) =>  if !c.is_zero() { return Some((i, c)); },
                        None => return None,
                    }
            },
            IntoNonzeroCoeffsIter::Sparse(it) => it.next(),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoNonzeroCoeffsIter::Zero()      => (0, Some(0)),
            IntoNonzeroCoeffsIter::Constant(c) => const_iter_size_hint(c),
            IntoNonzeroCoeffsIter::Dense(it)   => it.size_hint(),
            IntoNonzeroCoeffsIter::Sparse(it)  => it.size_hint(),
        }
    }
}

//********** Nonzero Coeffs Iterator (NOT consume self) **********/
pub(crate) fn nonzero_coeffs_iter<'a, C>(p: &'a Polynomial<C>) -> NonzeroCoeffsIter<'a, C> where C: Semiring {
    match p {
        Polynomial::Zero()       => NonzeroCoeffsIter::Zero(),
        Polynomial::Constant(cc) => NonzeroCoeffsIter::Constant(Some(&cc.0)),
        Polynomial::Dense(dc)    => NonzeroCoeffsIter::Dense(dc.0.iter().enumerate()),
        Polynomial::Sparse(sc)   => NonzeroCoeffsIter::Sparse(sc.0.iter()),
    }
}

pub enum NonzeroCoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<&'a C>),
    Dense(iter::Enumerate<slice::Iter<'a, C>>),
    Sparse(btree_map::Iter<'a, usize, C>),
}

impl<'a, C> Iterator for NonzeroCoeffsIter<'a, C> where C: Semiring {

    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            NonzeroCoeffsIter::Zero() => None,
            NonzeroCoeffsIter::Constant(c) => 
                if let Some(c) = c.take() { Some((0, c)) } else { None },
            NonzeroCoeffsIter::Dense(it) => 
                loop {
                    match it.next() {
                        Some((i, c)) =>  if !c.is_zero() { return Some((i, c)); },
                        None => return None,
                    }
                },
            NonzeroCoeffsIter::Sparse(it) => it.next().map(|(i, c)| (*i, c)),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            NonzeroCoeffsIter::Zero()      => (0, Some(0)),
            NonzeroCoeffsIter::Constant(c) => const_iter_size_hint(c),
            NonzeroCoeffsIter::Dense(it)   => it.size_hint(),
            NonzeroCoeffsIter::Sparse(it)  => it.size_hint(),
        }
    }
}

//********** Into Coeffs Iterator (consumes self and returns value itself) **********/
pub(crate) fn into_coeffs_iter<C>(p: Polynomial<C>) -> IntoCoeffsIter<C> where C: Semiring {
    match p {
        Polynomial::Zero()       => IntoCoeffsIter::Zero(),
        Polynomial::Constant(cc) => IntoCoeffsIter::Constant(Some(cc.0)),
        Polynomial::Dense(dc)    => IntoCoeffsIter::Dense(dc.0.into_iter()),
        Polynomial::Sparse(sc)   => {
            let mut map_iter = sc.0.into_iter();
            IntoCoeffsIter::Sparse{ index: 0, current: map_iter.next(), map_iter }
        }, 
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
            IntoCoeffsIter::Zero()          => None,
            IntoCoeffsIter::Constant(value) => value.take(),
            IntoCoeffsIter::Dense(it)       => it.next(),
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
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoCoeffsIter::Zero()      => (0, Some(0)),
            IntoCoeffsIter::Constant(c) => const_iter_size_hint(c),
            IntoCoeffsIter::Dense(it)   => it.size_hint(),
            IntoCoeffsIter::Sparse{ index: _, current: _, map_iter} => map_iter.size_hint(),
        }
    }
}

//********** Coeffs Iterator (NOT consume self) **********/
pub(crate) fn coeffs_iter<'a, C>(p: &'a Polynomial<C>) -> CoeffsIter<C> where C: Semiring {
    match p {
        Polynomial::Zero()       => CoeffsIter::Zero(),
        Polynomial::Constant(cc) => CoeffsIter::Constant(Some(Some(&cc.0))),
        Polynomial::Dense(dc)    => CoeffsIter::Dense(dc.0.iter()),
        Polynomial::Sparse(sc)   => {
            let mut map_iter = sc.0.iter();
            CoeffsIter::Sparse{ index: 0, current: map_iter.next(), map_iter }
        }, 
    }
}

/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some(Some(0))</code> if the coefficient is zero.
pub enum CoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<Option<&'a C>>),
    Dense(slice::Iter<'a, C>),
    Sparse {
        index: usize,
        current: Option<(&'a usize, &'a C)>,
        map_iter: btree_map::Iter<'a, usize, C>,
    },
}

impl<'a, C> Iterator for CoeffsIter<'a, C> where C: Semiring {

    type Item = Option<&'a C>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            CoeffsIter::Zero() => None,
            CoeffsIter::Constant(value) => value.take(),
            CoeffsIter::Dense(it) =>
                if let Some(c) = it.next() { Some(Some(c)) } else { None },
            CoeffsIter::Sparse{ index, current, map_iter } => {
                if let Some(c) = current {
                    let result = if c.0 == index { 
                        let r = Some(c.1);
                        *current = map_iter.next();
                        r
                    } else {
                        None
                    };
                    *index += 1;
                    Some(result)
                } else {
                    return None;
                }
            },
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            CoeffsIter::Zero()      => (0, Some(0)),
            CoeffsIter::Constant(_) => (1, Some(1)),
            CoeffsIter::Dense(it)   => it.size_hint(),
            CoeffsIter::Sparse{ index: _, current: _, map_iter } => map_iter.size_hint(),
        }
    }
}