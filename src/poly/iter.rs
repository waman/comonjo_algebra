use std::{collections::btree_map, iter, slice, vec};

use crate::{algebra::Semiring, poly::Polynomial};

fn const_iter_size_hint<'a, C>(c: &'a Option<C>) -> (usize, Option<usize>) {
    if c.is_some() { (1, Some(1)) } else { (0, Some(0)) }
}

fn sparse_next_including_zero<C, R, F, G>(
        index: &mut usize,
        current: &mut Option<(usize, C)>,
        map_iter: &mut btree_map::IntoIter<usize, C>,
        f_nonzero: F, f_zero: G) -> Option<R> 
        where F: Fn(usize, C) -> R, G: Fn(usize) -> R {
    match current.take() {
        Some(c) => {
            let result = if c.0 == *index {
                let r = f_nonzero(c.0, c.1);
                *current = map_iter.next();
                r
            } else {
                *current = Some(c);
                f_zero(*index)
            };

            *index += 1;
            Some(result)
        }
        _ => None,
    }
}

fn sparse_next_including_zero_with_ref_index<'a, C, R, F, G>(
        index: &mut usize,
        current: &mut Option<(&'a usize, &'a C)>,
        map_iter: &mut btree_map::Iter<'a, usize, C>,
        f_nonzero: F, f_zero: G) -> Option<R> 
        where F: Fn(usize, &'a C) -> R, G: Fn(usize) -> R {
    match current {
        Some(c) => {
            let result = if c.0 == index {
                let r = f_nonzero(*c.0, c.1);
                *current = map_iter.next();
                r
            } else { 
                f_zero(*index)
            };

            *index += 1;
            Some(result)
        }
        _ => None,
    }
}

//********** CoeffsIter *********/
//***** Into Coeffs Iterator (consumes self and returns value itself)
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
            IntoCoeffsIter::Sparse{ index, current, map_iter } =>
                sparse_next_including_zero(index, current, map_iter, |_, c| c, |_| C::zero()),
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

//***** Coeffs Iterator (NOT consume self)
pub(crate) fn coeffs_iter<'a, C>(p: &'a Polynomial<C>) -> CoeffsIter<C> where C: Semiring {
    match p {
        Polynomial::Zero()       => CoeffsIter::Zero(),
        Polynomial::Constant(cc) => CoeffsIter::Constant(Some(&cc.0)),
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
    Constant(Option<&'a C>),
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
            CoeffsIter::Constant(value) => value.take().map(|c| Some(c)),
            CoeffsIter::Dense(it) => it.next().map(|c| Some(c)),
            CoeffsIter::Sparse{ index, current, map_iter } =>
                sparse_next_including_zero_with_ref_index(index, current, map_iter, |_, c| Some(c), |_| None),
                // match current {
                //     Some(c) => {
                //         let result = if c.0 == index { 
                //             let r = Some(c.1);
                //             *current = map_iter.next();
                //             r
                //         } else {
                //             None
                //         };

                //         *index += 1;
                //         Some(result)
                //     },
                //     _ => None,
                // },
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

//********** NonzeroCoeffsIter **********/
//***** Into Nonzero Coeffs Iterator (consumes self and returns value itself)
pub(crate) fn into_nonzero_coeffs_iter<C>(p: Polynomial<C>) -> IntoNonzeroCoeffsIter<C>
        where C: Semiring {
    match p {
        Polynomial::Zero()       => IntoNonzeroCoeffsIter::Zero(),
        Polynomial::Constant(cc) => IntoNonzeroCoeffsIter::Constant(Some(cc.0)),
        Polynomial::Dense(dc)    => IntoNonzeroCoeffsIter::Dense(dc.0.into_iter()),
        Polynomial::Sparse(sc)   => IntoNonzeroCoeffsIter::Sparse(sc.0.into_iter()),
    }
}

pub enum IntoNonzeroCoeffsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(vec::IntoIter<C>),
    Sparse(btree_map::IntoIter<usize, C>),
}

impl<C> Iterator for IntoNonzeroCoeffsIter<C> where C: Semiring {

    type Item = C;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoNonzeroCoeffsIter::Zero()          => None,
            IntoNonzeroCoeffsIter::Constant(value) => value.take(),
            IntoNonzeroCoeffsIter::Dense(it)       => it.skip_while(|c|c.is_zero()).next(),
            IntoNonzeroCoeffsIter::Sparse(it)      => it.next().map(|(_, c)| c),
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

//***** Nonzero Coeffs Iterator (NOT consume self)
pub(crate) fn nonzero_coeffs_iter<'a, C>(p: &'a Polynomial<C>) -> NonzeroCoeffsIter<C> where C: Semiring {
    match p {
        Polynomial::Zero()       => NonzeroCoeffsIter::Zero(),
        Polynomial::Constant(cc) => NonzeroCoeffsIter::Constant(Some(&cc.0)),
        Polynomial::Dense(dc)    => NonzeroCoeffsIter::Dense(dc.0.iter()),
        Polynomial::Sparse(sc)   => NonzeroCoeffsIter::Sparse(sc.0.iter()),
    }
}

/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some(Some(0))</code> if the coefficient is zero.
pub enum NonzeroCoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<&'a C>),
    Dense(slice::Iter<'a, C>),
    Sparse(btree_map::Iter<'a, usize, C>),
}

impl<'a, C> Iterator for NonzeroCoeffsIter<'a, C> where C: Semiring {

    type Item = &'a C;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            NonzeroCoeffsIter::Zero() => None,
            NonzeroCoeffsIter::Constant(value) => value.take(),
            NonzeroCoeffsIter::Dense(it) => it.skip_while(|c| c.is_zero()).next(),
            NonzeroCoeffsIter::Sparse(it) => it.next().map(|(_, c)| c),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            NonzeroCoeffsIter::Zero()      => (0, Some(0)),
            NonzeroCoeffsIter::Constant(_) => (1, Some(1)),
            NonzeroCoeffsIter::Dense(it)   => it.size_hint(),
            NonzeroCoeffsIter::Sparse(it)  => it.size_hint(),
        }
    }
}

//********** TermsIter *********/
//***** Into Terms Iterator (consumes self and returns value itself)
pub(crate) fn into_terms_iter<C>(p: Polynomial<C>) -> IntoTermsIter<C> where C: Semiring {
    match p {
        Polynomial::Zero()       => IntoTermsIter::Zero(),
        Polynomial::Constant(cc) => IntoTermsIter::Constant(Some(cc.0)),
        Polynomial::Dense(dc)    => IntoTermsIter::Dense(dc.0.into_iter().enumerate()),
        Polynomial::Sparse(sc)   => {
            let mut map_iter = sc.0.into_iter();
            IntoTermsIter::Sparse{ index: 0, current: map_iter.next(), map_iter }
        }, 
    }
}

pub enum IntoTermsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(iter::Enumerate<vec::IntoIter<C>>),
    Sparse{
        index: usize,
        current: Option<(usize, C)>,
        map_iter: btree_map::IntoIter<usize, C>
    },
}

impl<C> Iterator for IntoTermsIter<C> where C: Semiring {

    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoTermsIter::Zero()          => None,
            IntoTermsIter::Constant(value) => value.take().map(|c| (0, c)),
            IntoTermsIter::Dense(it)       => it.next(),
            IntoTermsIter::Sparse{ index, current, map_iter } =>
                sparse_next_including_zero(index, current, map_iter, |i, c| (i, c), |i| (i, C::zero())),
                // match current {
                //     Some(c) => {
                //         let result = if c.0 == *index {
                //             let r = (c.0, c.1);
                //             *current = map_iter.next();
                //             Some(r)
                //         } else { 
                //             Some((*index, C::zero()))
                //         };

                //         *index += 1;
                //         result
                //     },
                //     _ => None,
                // },
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoTermsIter::Zero()      => (0, Some(0)),
            IntoTermsIter::Constant(c) => const_iter_size_hint(c),
            IntoTermsIter::Dense(it)   => it.size_hint(),
            IntoTermsIter::Sparse{ index: _, current: _, map_iter} => map_iter.size_hint(),
        }
    }
}

//***** Terms Iterator (NOT consume self)
pub(crate) fn terms_iter<'a, C>(p: &'a Polynomial<C>) -> TermsIter<C> where C: Semiring {
    match p {
        Polynomial::Zero()       => TermsIter::Zero(),
        Polynomial::Constant(cc) => TermsIter::Constant(Some(&cc.0)),
        Polynomial::Dense(dc)    => TermsIter::Dense(dc.0.iter().enumerate()),
        Polynomial::Sparse(sc)   => {
            let mut map_iter = sc.0.iter();
            TermsIter::Sparse{ index: 0, current: map_iter.next(), map_iter }
        }, 
    }
}

/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some((usize, Some(0)))</code> if the coefficient is zero.
pub enum TermsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<&'a C>),
    Dense(iter::Enumerate<slice::Iter<'a, C>>),
    Sparse {
        index: usize,
        current: Option<(&'a usize, &'a C)>,
        map_iter: btree_map::Iter<'a, usize, C>,
    },
}

impl<'a, C> Iterator for TermsIter<'a, C> where C: Semiring {

    type Item = (usize, Option<&'a C>);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            TermsIter::Zero() => None,
            TermsIter::Constant(value) => value.take().map(|c| (0, Some(c))),
            TermsIter::Dense(it) => it.next().map(|(i, c)| (i, Some(c))),
            TermsIter::Sparse{ index, current, map_iter } => 
                sparse_next_including_zero_with_ref_index(index, current, map_iter, 
                    |i, c| (i, Some(c)), |i| (i, None)),
                // match current {
                //     Some(c) => {
                //         let result = if c.0 == index { 
                //             let r = Some(c.1);
                //             *current = map_iter.next();
                //             r
                //         } else {
                //             None
                //         };
                //         *index += 1;
                //         Some(result)
                //     },
                //     _ => None
                // },
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            TermsIter::Zero()      => (0, Some(0)),
            TermsIter::Constant(_) => (1, Some(1)),
            TermsIter::Dense(it)   => it.size_hint(),
            TermsIter::Sparse{ index: _, current: _, map_iter } => map_iter.size_hint(),
        }
    }
}

//********** NonzeroTermsIter **********/
//***** Into Nonzero Terms Iterator (consumes self and returns value itself)
pub(crate) fn into_nonzero_terms_iter<C>(p: Polynomial<C>) -> IntoNonzeroTermsIter<C> where C: Semiring {
    match p {
        Polynomial::Zero()       => IntoNonzeroTermsIter::Zero(),
        Polynomial::Constant(cc) => IntoNonzeroTermsIter::Constant(Some(cc.0)),
        Polynomial::Dense(dc)    => IntoNonzeroTermsIter::Dense(dc.0.into_iter().enumerate()),
        Polynomial::Sparse(sc)   => IntoNonzeroTermsIter::Sparse(sc.0.into_iter()),
    }
}

pub enum IntoNonzeroTermsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(iter::Enumerate<vec::IntoIter<C>>),
    Sparse(btree_map::IntoIter<usize, C>),
}

impl<C> Iterator for IntoNonzeroTermsIter<C> where C: Semiring {

    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoNonzeroTermsIter::Zero() => None,
            IntoNonzeroTermsIter::Constant(c) => c.take().map(|c| (0, c)),
            IntoNonzeroTermsIter::Dense(it) => it.skip_while(|c| c.1.is_zero()).next(),
            IntoNonzeroTermsIter::Sparse(it) => it.next(),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoNonzeroTermsIter::Zero()      => (0, Some(0)),
            IntoNonzeroTermsIter::Constant(c) => const_iter_size_hint(c),
            IntoNonzeroTermsIter::Dense(it)   => it.size_hint(),
            IntoNonzeroTermsIter::Sparse(it)  => it.size_hint(),
        }
    }
}

//***** Nonzero Terms Iterator (NOT consume self)
pub(crate) fn nonzero_terms_iter<'a, C>(p: &'a Polynomial<C>) -> NonzeroTermsIter<'a, C> where C: Semiring {
    match p {
        Polynomial::Zero()       => NonzeroTermsIter::Zero(),
        Polynomial::Constant(cc) => NonzeroTermsIter::Constant(Some(&cc.0)),
        Polynomial::Dense(dc)    => NonzeroTermsIter::Dense(dc.0.iter().enumerate()),
        Polynomial::Sparse(sc)   => NonzeroTermsIter::Sparse(sc.0.iter()),
    }
}

pub enum NonzeroTermsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<&'a C>),
    Dense(iter::Enumerate<slice::Iter<'a, C>>),
    Sparse(btree_map::Iter<'a, usize, C>),
}

impl<'a, C> Iterator for NonzeroTermsIter<'a, C> where C: Semiring {

    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            NonzeroTermsIter::Zero() => None,
            NonzeroTermsIter::Constant(c) => c.take().map(|c| (0, c)),
            NonzeroTermsIter::Dense(it) => it.skip_while(|c| c.1.is_zero()).next(),
            NonzeroTermsIter::Sparse(it) => it.next().map(|(i, c)| (*i, c)),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            NonzeroTermsIter::Zero()      => (0, Some(0)),
            NonzeroTermsIter::Constant(c) => const_iter_size_hint(c),
            NonzeroTermsIter::Dense(it)   => it.size_hint(),
            NonzeroTermsIter::Sparse(it)  => it.size_hint(),
        }
    }
}