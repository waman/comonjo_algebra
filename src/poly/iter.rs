use crate::{algebra::Semiring, poly::{dense::{DCoeffsIter, DIntoCoeffsIter, DIntoNonzeroCoeffsIter, DNonzeroCoeffsIter}, sparse::{SCoeffsIter, SIntoCoeffsIter, SIntoNonzeroCoeffsIter, SNonzeroCoeffsIter}}};

fn const_iter_size_hint<'a, C>(c: &'a Option<C>) -> (usize, Option<usize>) {
    if c.is_some() { (1, Some(1)) } else { (0, Some(0)) }
}

pub enum IntoNonzeroCoeffsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(DIntoNonzeroCoeffsIter<C>),
    Sparse(SIntoNonzeroCoeffsIter<C>),
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
            IntoNonzeroCoeffsIter::Dense(di) => di.next(),
            IntoNonzeroCoeffsIter::Sparse(si) => si.next(),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoNonzeroCoeffsIter::Zero()      => (0, Some(0)),
            IntoNonzeroCoeffsIter::Constant(c) => const_iter_size_hint(c),
            IntoNonzeroCoeffsIter::Dense(di)   => di.size_hint(),
            IntoNonzeroCoeffsIter::Sparse(si)  => si.size_hint(),
        }
    }
}

pub enum NonzeroCoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<&'a C>),
    Dense(DNonzeroCoeffsIter<'a, C>),
    Sparse(SNonzeroCoeffsIter<'a, C>),
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
            NonzeroCoeffsIter::Dense(di) => di.next(),
            NonzeroCoeffsIter::Sparse(si) => si.next(),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            NonzeroCoeffsIter::Zero()      => (0, Some(0)),
            NonzeroCoeffsIter::Constant(c) => const_iter_size_hint(c),
            NonzeroCoeffsIter::Dense(di)   => di.size_hint(),
            NonzeroCoeffsIter::Sparse(si)  => si.size_hint(),
        }
    }
}

pub enum IntoCoeffsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(DIntoCoeffsIter<C>),
    Sparse(SIntoCoeffsIter<C>),
}

impl<C> Iterator for IntoCoeffsIter<C> where C: Semiring {

    type Item = C;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoCoeffsIter::Zero()          => None,
            IntoCoeffsIter::Constant(value) => value.take(),
            IntoCoeffsIter::Dense(di)       =>  di.next(),
            IntoCoeffsIter::Sparse(si)      => si.next(),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoCoeffsIter::Zero()      => (0, Some(0)),
            IntoCoeffsIter::Constant(c) => const_iter_size_hint(c),
            IntoCoeffsIter::Dense(di)   => di.size_hint(),
            IntoCoeffsIter::Sparse(si)  => si.size_hint(),
        }
    }
}

/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some(Some(0))</code> if the coefficient is zero.
pub enum CoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<Option<&'a C>>),
    Dense(DCoeffsIter<'a, C>),
    Sparse(SCoeffsIter<'a, C>),
}

impl<'a, C> Iterator for CoeffsIter<'a, C> where C: Semiring {

    type Item = Option<&'a C>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            CoeffsIter::Zero()          => None,
            CoeffsIter::Constant(value) => value.take(),
            CoeffsIter::Dense(di)       => di.next(),
            CoeffsIter::Sparse(si)      => si.next(),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            CoeffsIter::Zero()      => (0, Some(0)),
            CoeffsIter::Constant(_) => (1, Some(1)),
            CoeffsIter::Dense(di)   => di.size_hint(),
            CoeffsIter::Sparse(si)  => si.size_hint(),
        }
    }
}