use std::{collections::BTreeMap, fmt::Debug, iter::Enumerate, slice::Iter, vec::IntoIter};

use num::Integer;

use crate::{algebra::{EuclideanRing, Field, Ring, Semiring}, poly::{CoeffsIterator, Polynomial, PolynomialOps, factorial, iter::{CoeffsIter, IntoCoeffsIter, IntoNonzeroCoeffsIter, NonzeroCoeffsIter}, mul_div_uint, remove_tail_zeros}};

#[derive(Clone)]
pub struct DenseCoeffs<C>(pub(crate) Vec<C>);

impl<C> DenseCoeffs<C> where C: Semiring {

    pub(crate) fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub(crate) fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(n)
    }

    pub(crate) fn max_order_term(&self) -> Option<(usize, &C)> {
        self.0.last().map(|c|(self.degree(), c))
    }

    pub(crate) fn min_order_term(&self) -> Option<(usize, &C)> {
        self.0.iter().enumerate().skip_while(|e|e.1.is_zero()).next()
    }
    
    pub(crate) fn is_x(&self) -> bool {
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

impl<C> PolynomialOps<C> for DenseCoeffs<C> where C: Semiring {
    
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
}

impl<'a, C> PolynomialOps<C> for &'a DenseCoeffs<C> where C: Semiring + Clone {
    
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
}

//********** Polynomial Operations **********/
fn normalize_vec_or_new_poly<C>(vec: &mut Vec<C>) -> Option<Polynomial<C>> where C: Semiring {

    remove_tail_zeros(vec);

    match vec.len() {
        0 => Some(Polynomial::Zero()),
        1 if vec.len() == 1 =>
            Some(Polynomial::new_raw_const(vec.pop().unwrap())),
        _ => {
            vec.shrink_to_fit();
            None
        },
    }
}

impl<C> DenseCoeffs<C> where C: Semiring {

    pub(crate) fn reciprocal(&mut self) -> Option<Polynomial<C>> {
        self.0.reverse();
        normalize_vec_or_new_poly(&mut self.0)
    }

    pub(crate) fn remove_zero_roots(&mut self) -> Option<Polynomial<C>> {
        if let Some((i, _)) = self.min_order_term() {
            if i == 0 { return None }
            let vec = &mut self.0;
            vec.drain(0..i);
            normalize_vec_or_new_poly(vec)
        } else {
            None
        }
    }

    pub(crate) fn reductum(&mut self) -> Option<Polynomial<C>> {
        self.0.pop();
        normalize_vec_or_new_poly(&mut self.0)
    }
}

impl<C> DenseCoeffs<C> where C: Semiring + Clone {

    pub(crate) fn new_reciprocal(&self) -> Polynomial<C> {
        let mut vec = Vec::with_capacity(self.0.len());
        self.0.iter().rev().for_each(|c| vec.push(c.clone()));
        Polynomial::dense_from_vec(vec)
    }

    pub(crate) fn new_zero_roots_removed(&self) -> Polynomial<C> {
        let vec: Vec<C> = self.0.iter().skip_while(|c|c.is_zero()).map(|c|c.clone()).collect();
        Polynomial::dense_from_vec(vec)
    }

    pub(crate) fn new_reductum(&self) -> Polynomial<C> {
        let n = self.0.len();
        let vec: Vec<C> = self.0.iter().take(n-1).map(|c|c.clone()).collect();
        Polynomial::dense_from_vec(vec)
    }
}

impl<C> DenseCoeffs<C> where C: Semiring + num::FromPrimitive {
    
    pub(crate) fn differentiate(&mut self) -> Option<Polynomial<C>> {
        if self.0.len() == 2 {
            return Some(Polynomial::new_raw_const(self.0.pop().unwrap()));
        }

        let d = self.degree();
        for i in 1..=d {
            let v = C::from_usize(i).unwrap() * self.0.get(i).unwrap();
            self.0[i-1] = v;
        }
        self.0.pop();
        self.0.shrink_to_fit();
        None
    }
    
    pub(crate) fn new_derivative(&self) -> Polynomial<C> {
        let mut ite = self.0.iter().enumerate();
        ite.next();
        let vec: Vec<C> = ite.map(|(i, c)| C::from_usize(i).unwrap() * c).collect();
        Polynomial::dense_from_vec(vec)
    }

    /// k * self.0[[n]]
    fn deriv_coeff<'a>(&self, n: usize, k: u128) -> C
            where C: Semiring + num::FromPrimitive {
        C::from_u128(k).unwrap().ref_mul(self.0.get(n).unwrap())
    }
    
    pub(crate) fn n_differentiate(&mut self, n: usize) -> Option<Polynomial<C>> {
        debug_assert!(n > 1);

        match self.degree() {
            d if d < n => Some(Polynomial::Zero()),
            d if d == n => {
                let f: C = factorial(n);
                Some(Polynomial::new_raw_const(f * self.0.pop().unwrap()))
            },
            d => {
                let mut k: u128 = factorial(n);
                self.0[0] = self.deriv_coeff(n, k);

                for i in (n+1)..=d {
                    k = mul_div_u128(k, i, i-n);
                    self.0[i-n] = self.deriv_coeff(i, k);
                }
                self.0.truncate(d-n+1);
                self.0.shrink_to_fit();
                None
            }
        }
    }
    
    pub(crate) fn new_nth_derivative(&self, n: usize) -> Polynomial<C> {
        debug_assert!(n > 1);

        match self.degree() {
            d if d < n => Polynomial::Zero(),
            d if d == n => {
                let f: C = factorial(n);
                Polynomial::new_raw_const(f * self.0.last().unwrap())
            },
            d => {
                let mut vec: Vec<C> = Vec::with_capacity(d-n+1);

                let mut k: u128 = factorial(n);
                vec.push(self.deriv_coeff(n, k));

                for i in (n+1)..=d {
                    k = mul_div_u128(k, i, i-n);
                    vec.push(self.deriv_coeff(i, k));
                }
                Polynomial::dense_from_vec(vec)
            }
        }
    }
}

fn mul_div_u128(x: u128, y: usize, z: usize) -> u128 {
    let z128 = z as u128;
    let gcd = x.gcd(&z128);
    let x_red = x / gcd;
    let z_red = (z as u128) / gcd;
    x_red * ((y as u128) / z_red)
}

impl<C> DenseCoeffs<C> where C: Ring {

    pub(crate) fn flip(&mut self){
        self.0.iter_mut().enumerate()
            .filter(|(i, _)| i % 2 != 0)
            .for_each(|(_, c)| *c = c.ref_neg());
    }
}

impl<C> DenseCoeffs<C> where C: Ring + Clone {

    pub(crate) fn new_flipped(&self) -> Polynomial<C> {
        let vec: Vec<C> = 
            self.0.iter().enumerate().map(|(i, c)|{
                if i % 2 == 0 { c.clone() } else { c.ref_neg() }
            }).collect();
        Polynomial::new_raw_dense(vec)
    }
}

impl<C> DenseCoeffs<C> where C: EuclideanRing + num::FromPrimitive + num::Integer + Clone {

    // ** From spire code *****
    // The trick here came from this answer:
    //   http://math.stackexchange.com/questions/694565/polynomial-shift
    // This is a heavily optimized version of the same idea. This is fairly
    // critical method to be fast, since it is the most expensive part of the
    // VAS root isolation algorithm.
    pub(crate) fn shift(&self, h: C) -> Polynomial<C> {

        let mut coeffs: Vec<C> = self.0.clone();
        for (deg, c) in self.nonzero_coeffs_iter() {
            if deg == 0 { continue; }
            let mut i: C = C::one();
            let mut d: usize = deg;
            let mut m: C = C::one();
            let mut k: C = c.clone();
            while d > 0 {
                m = mul_div_uint(m, d, i.clone());  // (m * d) / i;
                k = k * &h;
                if let Some(coeff) = coeffs.get_mut(d-1) {
                    *coeff = coeff.ref_add(m.ref_mul(&k));
                }
                d = d - 1;
                i = i + C::one();
            }
        }
        Polynomial::new_raw_dense(coeffs)
    }
}

impl<C> DenseCoeffs<C> where C: Field + num::FromPrimitive + Clone {

    pub(crate) fn shift_f(&self, h: C) -> Polynomial<C> {

        let mut coeffs: Vec<C> = self.0.clone();
        for (deg, c) in self.nonzero_coeffs_iter() {
            if deg == 0 { continue; }
            let mut i: C = C::one();
            let mut d: usize = deg;
            let mut m: C = C::one();
            let mut k: C = c.clone();
            while d > 0 {
                m = m * C::from_usize(d).unwrap() / i.clone();
                k = k * &h;
                if let Some(coeff) = coeffs.get_mut(d-1) {
                    *coeff = coeff.ref_add(m.ref_mul(&k));
                }
                d = d - 1;
                i = i + C::one();
            }
        }
        Polynomial::new_raw_dense(coeffs)
    }
}

impl<C> DenseCoeffs<C> where C: Field {

    pub(crate) fn monic(&mut self) {
        let a = self.0.pop().unwrap();
        self.0.iter_mut().for_each(|c| *c = c.ref_div(&a));
        self.0.push(C::one());
    }
}

impl<C> DenseCoeffs<C> where C: Field + Clone {

    pub(crate) fn new_monic(&self) -> Polynomial<C> {
        let a = self.0.last().unwrap();
        let vec: Vec<C> = self.0.iter().map(|c| c.ref_div(a)).collect();
        Polynomial::new_raw_dense(vec)
    }
}
    
impl<C> DenseCoeffs<C> where C: Field + num::FromPrimitive + Debug {

    pub(crate) fn integrate(&mut self) {
        let d = self.degree();
        self.0.reserve(1);

        let v_last = self.0.get(d).unwrap().ref_div(C::from_usize(d+1).unwrap());
        self.0.push(v_last);

        for i in 1..=d {
            let j = d - i + 1;
            self.0[j] = self.0.get(j-1).unwrap().ref_div(C::from_usize(j).unwrap());
        }

        self.0[0] = C::zero();
    }

    pub(crate) fn new_integral(&self) -> Polynomial<C> {
        let mut vec: Vec<C> = Vec::with_capacity(self.0.len() + 1);
        vec.push(C::zero());
        vec.extend(self.0.iter().enumerate().map(|(i, c)| c.ref_div(C::from_usize(i+1).unwrap())));
        Polynomial::new_raw_dense(vec)
    }

    // /// self.0[[n]] / k;
    // fn integ_coeff(&self, n: usize, k: &C) -> C {
    //     self.0.get(n).unwrap().ref_div(k)
    // } 

    // pub(crate) fn n_integrate(&mut self, n: usize) {
    //     debug_assert!(n > 1);

    //     todo!()
    // }

    // pub(crate) fn new_nth_integral(&self, n: usize) -> Polynomial<C> {
        // debug_assert!(n > 1);

        // let d = self.degree();
        // let mut vec: Vec<C> = Vec::with_capacity(d+n+1);
        // vec.fill_with(|| C::one());
        // for _ in 0..n { vec.push(C::zero()); } 

        // let mut k: C = factorial(n);
        // vec.push(self.integ_coeff(0, &k));

        // for i in 1..d {
        //     k = k * C::from_usize(n+i).unwrap() / C::from_usize(i).unwrap();
        //     vec.push(self.integ_coeff(i, &k));
        // }

        // Polynomial::new_raw_dense(vec)
    // }
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

    Polynomial::new_raw_dense(v)
}
//********** Div & Rem **********/
pub(crate) fn div_rem<C>(mut u: Vec<C>, rhs: &Polynomial<C>) -> (Polynomial<C>, Polynomial<C>) 
        where C: Field + Clone {

    let d_rhs = rhs.degree();
    let n = u.len() - d_rhs;  // = lhs.degree() + 1 - rhs.degree()
    let v0: &C = rhs.max_order_term().unwrap().1;
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
                if let Some(c_lhs) = u.get_mut(offset + i) {
                    *c_lhs = c_lhs.ref_sub(c_rhs.ref_mul(&q0));
                }
            }
            q.push(q0);
        }
    }

    q.reverse();
    (Polynomial::dense_from_vec(q), Polynomial::dense_from_vec(u))
}