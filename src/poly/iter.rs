use std::{collections::btree_map, iter::Enumerate, slice::Iter, vec};

use crate::algebra::Semiring;

fn const_iter_size_hint<'a, C>(c: &'a Option<C>) -> (usize, Option<usize>) {
    if c.is_some() { (1, Some(1)) } else { (0, Some(0)) }
}

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
                        Some((i, c)) =>  if !c.is_zero() { return Some((i, c)); },
                        None => return None,
                    }
                },
            IntoNonzeroCoeffsIter::Sparse(map_iter) => map_iter.next(),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoNonzeroCoeffsIter::Zero() => (0, Some(0)),
            IntoNonzeroCoeffsIter::Constant(c) => const_iter_size_hint(c),
            IntoNonzeroCoeffsIter::Dense(ite) => ite.size_hint(),
            IntoNonzeroCoeffsIter::Sparse(ite) => ite.size_hint(),
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
                        Some((i, c)) =>  if !c.is_zero() { return Some((i, c)); },
                        None => return None,
                    }
                },
            NonzeroCoeffsIter::Sparse(map_iter) => map_iter.next().map(|(i, c)| (*i, c)),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            NonzeroCoeffsIter::Zero() => (0, Some(0)),
            NonzeroCoeffsIter::Constant(c) => const_iter_size_hint(c),
            NonzeroCoeffsIter::Dense(ite) => ite.size_hint(),
            NonzeroCoeffsIter::Sparse(ite) => ite.size_hint(),
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
        map_iter: btree_map::IntoIter<usize, C>,
    },
}

impl<C> Iterator for IntoCoeffsIter<C> where C: Semiring {

    type Item = C;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoCoeffsIter::Zero() => None,
            IntoCoeffsIter::Constant(value) => value.take(),
            IntoCoeffsIter::Dense(vec_iter) =>  vec_iter.next(),
            IntoCoeffsIter::Sparse{ index, current, map_iter, .. } => {
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
            IntoCoeffsIter::Zero() => (0, Some(0)),
            IntoCoeffsIter::Constant(c) => const_iter_size_hint(c),
            IntoCoeffsIter::Dense(ite) => ite.size_hint(),
            IntoCoeffsIter::Sparse { map_iter, .. } => map_iter.size_hint(),
        }
    }
}

/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some(Some(0))</code> if the coefficient is zero.
pub enum CoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<Option<&'a C>>),
    Dense(Iter<'a, C>),
    Sparse{
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
            CoeffsIter::Dense(vec_iter) => 
                match vec_iter.next() {
                    Some(c) => Some(Some(c)),
                    None => None,
                } ,
            CoeffsIter::Sparse{ index, current, map_iter, .. } => 
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
                },
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            CoeffsIter::Zero() => (0, Some(0)),
            CoeffsIter::Constant(_) => (1, Some(1)),
            CoeffsIter::Dense(ite) => ite.size_hint(),
            CoeffsIter::Sparse { map_iter, .. } => map_iter.size_hint(),
        }
    }
}