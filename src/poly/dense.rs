use std::{collections::BTreeMap, iter::Enumerate, slice::Iter, vec::IntoIter};
use crate::{algebra::{EuclideanRing, Field, Ring, Semiring}, poly::{CoeffsIterator, Differentiable, EuclideanRingPolyOps, Integrable, Polynomial, RingPolyOps, SemiringPolyOps, iter::{CoeffsIter, IntoCoeffsIter, IntoNonzeroCoeffsIter, NonzeroCoeffsIter}, mul_div_uint}};

#[derive(Clone)]
pub struct DenseCoeffs<C>(pub(crate) Vec<C>);

impl<C> DenseCoeffs<C> where C: Semiring {

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(n)
    }
    
    pub fn is_x(&self) -> bool {
        self.0.len() == 2 && self.0[0].is_zero() && self.0[1].is_one()
    }


// def apply(x: C)(implicit ring: Semiring[C]): C = {
//   if (isZero) return ring.zero

//   var even = coeffs.length - 1
//   var odd = coeffs.length - 2
//   if ((even & 1) == 1) { even = odd; odd = coeffs.length - 1 }

//   var c0 = coeffs(even)
//   val x2 = x.pow(2)
//   cfor(even - 2)(_ >= 0, _ - 2) { i =>
//     c0 = coeffs(i) + c0 * x2
//   }

//   if (odd >= 1) {
//     var c1 = coeffs(odd)
//     cfor(odd - 2)(_ >= 1, _ - 2) { i =>
//       c1 = coeffs(i) + c1 * x2
//     }
//     c0 + c1 * x
//   } else {
//     c0
//   }
// }
}

impl<C> SemiringPolyOps<C> for DenseCoeffs<C> where C: Semiring {
    
    fn to_vec(self) -> Vec<C> { self.0 }
    
    fn to_map(mut self) -> BTreeMap<usize, C> {
        let mut map = BTreeMap::new();
        while let Some(c) = self.0.pop() {
            if !c.is_zero() {
                map.insert(self.0.len(), c);
            }
        }
        debug_assert!(self.0.is_empty());
        map
    }

    fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let v: Vec<D> = self.0.into_iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::dense_from_vec(v)
    }

    fn reductum(mut self) -> Polynomial<C> {
        self.0.pop();
        Polynomial::dense_from_vec(self.0)
    }

    fn remove_zero_roots(self) -> Polynomial<C> {
        let vec: Vec<C> = self.0.into_iter().skip_while(|c|c.is_zero()).collect();
        Polynomial::dense_from_vec(vec)
    }

    fn reciprocal(self) -> Polynomial<C> {
        let mut vec = self.0;
        vec.reverse();
        Polynomial::dense_from_vec(vec)
    }
}

impl<'a, C> SemiringPolyOps<C> for &'a DenseCoeffs<C> where C: Semiring + Clone {
    
    fn to_vec(self) -> Vec<C> { self.0.clone() }
    
    fn to_map(self) -> BTreeMap<usize, C> {
        self.nonzero_coeffs_iter().map(|(i, c)|(i, c.clone())).collect()
    }

    fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let v: Vec<D> = self.0.iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c.clone()) }
        ).collect();
        Polynomial::dense_from_vec(v)
    }

    fn reductum(self) -> Polynomial<C> {
        let n = self.0.len();
        let vec: Vec<C> = self.0.iter().take(n-1).map(|c|c.clone()).collect();
        Polynomial::dense_from_vec(vec)
    }

    fn remove_zero_roots(self) -> Polynomial<C> {
        let vec: Vec<C> = self.0.iter().skip_while(|c|c.is_zero()).map(|c|c.clone()).collect();
        Polynomial::dense_from_vec(vec)
    }

    fn reciprocal(self) -> Polynomial<C> {
        let mut vec = Vec::with_capacity(self.0.len());
        self.0.iter().rev().for_each(|c| vec.push(c.clone()));
        Polynomial::dense_from_vec(vec)
    }
}

impl<C> RingPolyOps<C> for DenseCoeffs<C> where C: Ring {

    fn flip(mut self) -> Polynomial<C> {
        self.0.iter_mut().enumerate().filter(|(i, _)| i % 2 != 0).for_each(|(_, c)| *c = c.ref_neg());
        Polynomial::Dense(self)
    }
}

impl<'a, C> RingPolyOps<C> for &'a DenseCoeffs<C> where C: Ring + Clone {

    fn flip(self) -> Polynomial<C> {
        let vec: Vec<C> = 
            self.0.iter().enumerate().map(|(i, c)|{
                if i % 2 == 0 { c.clone() } else { c.ref_neg() }
            }).collect();
        Polynomial::Dense(DenseCoeffs(vec))
    }
}

impl<'a, C> EuclideanRingPolyOps<C> for &'a DenseCoeffs<C> where C: EuclideanRing + num::FromPrimitive + Clone {

    // ** From spire code *****
    // The trick here came from this answer:
    //   http://math.stackexchange.com/questions/694565/polynomial-shift
    // This is a heavily optimized version of the same idea. This is fairly
    // critical method to be fast, since it is the most expensive part of the
    // VAS root isolation algorithm.
    fn shift(self, h: C) -> Polynomial<C> {

        let mut coeffs: Vec<C> = self.0.clone();
        for (deg, c) in self.nonzero_coeffs_iter() {
            if deg == 0 { continue }
            let mut i: u128 = 1;
            let mut d: usize = deg;
            let mut m: u128 = 1;
            let mut k: C = c.clone();
            while d > 0 {
                m = mul_div_uint(m, d, i);  // (m * d) / i;
                k = k * &h;
                if let Some(coeff) = coeffs.get_mut(d-1) {
                    *coeff = coeff.ref_add(C::from_u128(m).unwrap() * &k);
                }
                d = d - 1;
                i = i + 1;
            }
        }
        Polynomial::dense_from_vec(coeffs)
    }
}

//********** Iterator **********/
impl<C> DenseCoeffs<C> where C: Semiring {

    pub fn into_nonzero_coeffs_iter(self) -> IntoNonzeroCoeffsIter<C> {
        IntoNonzeroCoeffsIter::Dense(DIntoNonzeroCoeffsIter(self.0.into_iter().enumerate()))
    }

    pub fn nonzero_coeffs_iter<'a>(&'a self) -> NonzeroCoeffsIter<'a, C> {
        NonzeroCoeffsIter::Dense(DNonzeroCoeffsIter(self.0.iter().enumerate()))
    }

    pub fn into_coeffs_iter(self) -> IntoCoeffsIter<C> {
        IntoCoeffsIter::Dense(DIntoCoeffsIter(self.0.into_iter()))
    }

    pub fn coeffs_iter<'a>(&'a self) -> CoeffsIter<'a, C> {
        CoeffsIter::Dense(DCoeffsIter(self.0.iter()))
    }
}

pub struct DIntoNonzeroCoeffsIter<C>(pub(crate) Enumerate<IntoIter<C>>) where C: Semiring;

impl<C> Iterator for DIntoNonzeroCoeffsIter<C> where C: Semiring {

    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.0.next() {
                Some((i, c)) =>  if !c.is_zero() { return Some((i, c)); },
                None => return None,
            }
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
}

pub struct DNonzeroCoeffsIter<'a, C>(pub(crate) Enumerate<Iter<'a, C>>) where C: Semiring;

impl<'a, C> Iterator for DNonzeroCoeffsIter<'a, C> where C: Semiring {

    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.0.next() {
                Some((i, c)) =>  if !c.is_zero() { return Some((i, c)); },
                None => return None,
            }
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
}

pub struct DIntoCoeffsIter<C>(pub(crate) IntoIter<C>) where C: Semiring;

impl<C> Iterator for DIntoCoeffsIter<C> where C: Semiring {

    type Item = C;

    fn next(&mut self) -> Option<Self::Item> { self.0.next() }

    fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
}

/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some(Some(0))</code> if the coefficient is zero.
pub struct DCoeffsIter<'a, C>(pub(crate) Iter<'a, C>) where C: Semiring;

impl<'a, C> Iterator for DCoeffsIter<'a, C> where C: Semiring {

    type Item = Option<&'a C>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.0.next() {
            Some(c) => Some(Some(c)),
            None => None,
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
}

//********** Add **********/
/// Return (max, min, left_is_higher)
fn max_min(x: usize, y: usize) -> (usize, usize, bool) {
    if x >= y { (x, y, true) } else { (y, x, false) }
}

pub(crate) fn add_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Semiring {

    let (d_max, d_min, lhs_is_higher) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(x), Some(y)) => v.push(x + y),
            _ => panic!(),
        }
    }

    let rest_iter = if lhs_is_higher { lhs_iter } else { rhs_iter };
    v.extend(rest_iter);

    Polynomial::dense_from_vec(v)
}

fn clone_coeff<'a, C>(op_c: Option<&'a C>) -> C where C: Semiring + Clone {
    match op_c {
        Some(c) => c.clone(),
        _ => C::zero(),
    }
}

pub(crate) fn add_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let (d_max, d_min, lhs_is_higher) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(x), Some(op_y)) => match op_y {
                Some(y) => v.push(x + y),
                _ => v.push(x),
            },
            _ => panic!(),
        }
    }

    if lhs_is_higher { 
        v.extend(lhs_iter); 
    } else {
        v.extend(rhs_iter.map(clone_coeff));
    }

    Polynomial::dense_from_vec(v)
}

pub(crate) fn add_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let (d_max, d_min, lhs_is_higher) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(op_x), Some(y)) => match op_x {
                Some(x) => v.push(x.ref_add(y)),
                _ => v.push(y),
            },
            _ => panic!(),
        }
    }

    if lhs_is_higher { 
        v.extend(lhs_iter.map(clone_coeff));
    } else {
        v.extend(rhs_iter); 
    }

    Polynomial::dense_from_vec(v)
}

pub(crate) fn add_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
        where C: Semiring + Clone {

    let (d_max, d_min, lhs_is_higher) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(op_x), Some(op_y)) => match (op_x, op_y) {
                (Some(x), Some(y)) => v.push(x.ref_add(y)),
                (Some(x), None) => v.push(x.clone()),
                (None, Some(y)) => v.push(y.clone()),
                (None, None) => v.push(C::zero()),
            },
            _ => panic!(),
        }
    }

    let rest_iter = if lhs_is_higher { lhs_iter } else { rhs_iter };
    v.extend(rest_iter.map(clone_coeff));

    Polynomial::dense_from_vec(v)
}

//********** Sub **********/
pub(crate) fn sub_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> 
        where C: Ring {

    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(x), Some(y)) => v.push(x - y),
            _ => panic!(),
        }
    }

    if lhs_is_longer {
        v.extend(lhs_iter);
    } else {
        v.extend(rhs_iter.map(|c|-c));
    }

    Polynomial::dense_from_vec(v)
}

fn neg_coeff<'a, C>(op_c: Option<&'a C>) -> C where C: Ring + Clone {
    match op_c {
        Some(c) => c.ref_neg(),
        _ => C::zero(),
    }
}

pub(crate) fn sub_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
        where C: Ring + Clone {

    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(x), Some(op_y)) => match op_y {
                Some(y) => v.push(x.ref_sub(y)),
                _ => v.push(x),
            },
            _ => panic!(),
        }
    }

    if lhs_is_longer {
        v.extend(lhs_iter);
    } else {
        v.extend(rhs_iter.map(neg_coeff));
    }

    Polynomial::dense_from_vec(v)
}

pub(crate) fn sub_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> 
        where C: Ring + Clone {

    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(op_x), Some(y)) => match op_x {
                Some(x) => v.push(x.ref_sub(y)),
                _ => v.push(-y),
            },
            _ => panic!(),
        }
    }

    if lhs_is_longer {
        v.extend(lhs_iter.map(|c|match c {
            Some(x) => x.clone(),
            None => C::zero(),
        }));
    } else {
        v.extend(rhs_iter.map(|c|-c));
    }

    Polynomial::dense_from_vec(v)
}

pub(crate) fn sub_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
        where C: Ring + Clone {

    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(some_x), Some(some_y)) => match (some_x, some_y) {
                (Some(x), Some(y)) => v.push(x.ref_sub(y)),
                (Some(x), None) => v.push(x.clone()),
                (None, Some(y)) => v.push(y.clone()),
                (None, None) => v.push(C::zero()),
            },
            _ => panic!(),
        }
    }

    if lhs_is_longer {
        v.extend(lhs_iter.map(|c| match c {
            Some(x) => x.clone(),
            None => C::zero(),
        }));
    } else {
        v.extend(rhs_iter.map(|c| match c {
            Some(x) => x.ref_neg(),
            None => C::zero(),
        }));
    }

    Polynomial::dense_from_vec(v)
}

//********** Mul **********/
pub(crate) fn mul<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let deg = lhs.degree() + rhs.degree();
    let mut v: Vec<C> = Vec::with_capacity(deg + 1);

    for (i, x) in lhs.nonzero_coeffs() {
        for (j, y) in rhs.nonzero_coeffs() {
            let k = i + j;
            let z = x.ref_mul(y);
            match v.get_mut(k) {
                Some(c) => *c = c.ref_add(z),
                None => {
                    let n = v.len();
                    for _ in 0..(k-n) { v.push(C::zero()) }
                    v.push(z);
                },
            }
        }
    }

    Polynomial::dense_from_vec(v)
}
//********** Div & Rem **********/
pub(crate) fn div_rem<C>(mut u: Vec<C>, rhs: &Polynomial<C>) -> (Polynomial<C>, Polynomial<C>) 
        where C: Field + Clone {

    let d_rhs = rhs.degree();
    let n = u.len() - d_rhs;  // = lhs.degree() + 1 - rhs.degree()
    let v0: &C = rhs.max_order_term_coeff().unwrap();
    let mut q: Vec<C> = Vec::with_capacity(n);

    for _ in 0..n {
        let u_last: C = u.pop().unwrap();
        if u_last.is_zero() {
            q.push(C::zero());
        } else {
            let q0: C = u_last.ref_div(v0);
            let offset = u.len() - d_rhs;  // the last of u is already popped
            for (i, c_rhs) in rhs.nonzero_coeffs() {
                if i == d_rhs { break; }
                let c_lhs = u.get_mut(offset + i).unwrap();
                *c_lhs = c_lhs.ref_sub(c_rhs.ref_mul(&q0));
            }
            q.push(q0);
        }
    }

    q.reverse();
    (Polynomial::dense_from_vec(q), Polynomial::dense_from_vec(u))
}

//********** Analysis **********/
impl<C> Differentiable<C> for DenseCoeffs<C> where C: Semiring + num::FromPrimitive {
    
    fn derivative(self) -> Polynomial<C> {
        let mut ite = self.0.into_iter().enumerate();
        ite.next();
        let vec: Vec<C> = ite.map(|(i, c)| C::from_usize(i).unwrap() * c).collect();
        Polynomial::dense_from_vec(vec)
    }
}

impl<'a, C> Differentiable<C> for &'a DenseCoeffs<C> where C: Semiring + num::FromPrimitive  + Clone {
    
    fn derivative(self) -> Polynomial<C> {
        let mut ite = self.0.iter().enumerate();
        ite.next();
        let vec: Vec<C> = ite.map(|(i, c)| C::from_usize(i).unwrap() * c).collect();
        Polynomial::dense_from_vec(vec)
    }
}
    
impl<C> Integrable<C> for DenseCoeffs<C> where C: Field + num::FromPrimitive {

    fn integral(self) -> Polynomial<C> {
        let n = self.0.len();
        let ite = self.0.into_iter().enumerate();
        let mut vec: Vec<C> = Vec::with_capacity(n+1);
        vec.push(C::zero());
        vec.extend(ite.map(|(i, c)| c / C::from_usize(i+1).unwrap()));
        Polynomial::dense_from_vec(vec)
    }
}

impl<'a, C> Integrable<C> for &'a DenseCoeffs<C> where C: Field + num::FromPrimitive  + Clone {

    fn integral(self) -> Polynomial<C> {
        let n = self.0.len();
        let ite = self.0.iter().enumerate();
        let mut vec: Vec<C> = Vec::with_capacity(n+1);
        vec.push(C::zero());
        vec.extend(ite.map(|(i, c)| c.clone() / C::from_usize(i+1).unwrap()));
        Polynomial::dense_from_vec(vec)
    }
}