use std::{collections::BTreeMap};

use num::Integer;

use crate::{algebra::{EuclideanRing, Field, Ring, Semiring}, poly::{CoeffsIterator, Polynomial, factorial, mul_div_uint, remove_tail_zeros}};

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

    pub(crate) fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let v: Vec<D> = self.0.into_iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::from(v)
    }

    pub(crate) fn map_nonzero_ref<'a, D, F>(&'a self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, &'a C) -> D {
        let v: Vec<D> = self.0.iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::from(v)
    }

    pub(crate) fn to_map(self) -> BTreeMap<usize, C> {
        self.0.into_iter().enumerate().filter(|(_, c)| !c.is_zero()).collect()
    }

    pub(crate) fn nonzero_coeffs_iter(&self) -> impl Iterator<Item=(usize, &C)> {
        self.0.iter().enumerate().filter(|(_, c)| !c.is_zero())
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
        Polynomial::from(vec)
    }

    pub(crate) fn new_zero_roots_removed(&self) -> Polynomial<C> {
        let vec: Vec<C> = self.0.iter().skip_while(|c|c.is_zero()).map(|c|c.clone()).collect();
        Polynomial::from(vec)
    }

    pub(crate) fn new_reductum(&self) -> Polynomial<C> {
        let n = self.0.len();
        let vec: Vec<C> = self.0.iter().take(n-1).map(|c|c.clone()).collect();
        Polynomial::from(vec)
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

    /// k * self.0[[n]]
    fn deriv_coeff<'a>(&self, n: usize, k: u128) -> C
            where C: Semiring + num::FromPrimitive {
        C::from_u128(k).unwrap().ref_mul(self.0.get(n).unwrap())
    }
    
    pub(crate) fn differentiate_n(&mut self, n: usize) -> Option<Polynomial<C>> {
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
}

impl<C> DenseCoeffs<C> where C: Semiring + num::FromPrimitive + Clone {
    
    pub(crate) fn new_derivative(&self) -> Polynomial<C> {
        let mut ite = self.0.iter().enumerate();
        ite.next();
        let vec: Vec<C> = ite.map(|(i, c)| C::from_usize(i).unwrap() * c).collect();
        Polynomial::from(vec)
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
                Polynomial::from(vec)
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
    fn new_shifted_coeffs(&self, h: C) -> Vec<C> {
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
        
        coeffs
    }

    pub(crate) fn shift(&mut self, h: C) {
        self.0 = self.new_shifted_coeffs(h);
    }

    pub(crate) fn new_shifted(&self, h: C) -> Polynomial<C> {
        Polynomial::new_raw_dense(self.new_shifted_coeffs(h))
    }
}

impl<C> DenseCoeffs<C> where C: Field + num::FromPrimitive + Clone {

    pub fn new_shifted_coeffs_f(&self, h: C) -> Vec<C> {
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

        coeffs
    }

    pub(crate) fn shift_f(&mut self, h: C) {
        self.0 = self.new_shifted_coeffs_f(h);
    }

    pub(crate) fn new_shifted_f(&self, h: C) -> Polynomial<C> {
        Polynomial::new_raw_dense(self.new_shifted_coeffs_f(h))
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
    
impl<C> DenseCoeffs<C> where C: Field + num::FromPrimitive {

    fn new_integral_coeffs(&self) -> Vec<C> {
        let mut vec: Vec<C> = Vec::with_capacity(self.0.len() + 1);
        vec.push(C::zero());
        vec.extend(self.0.iter().enumerate().map(|(i, c)| c.ref_div(C::from_usize(i+1).unwrap())));
        vec
    }

    pub(crate) fn integrate(&mut self) {
        self.0 = self.new_integral_coeffs();
    }

    pub(crate) fn new_integral(&self) -> Polynomial<C> {
        Polynomial::new_raw_dense(self.new_integral_coeffs())
    }

    fn new_n_integral_coeffs(&self, n: usize) -> Vec<C> {
        let d = self.degree();
        let mut vec: Vec<C> = Vec::with_capacity(d+n+1);
        for _ in 0..n { vec.push(C::zero()); } 

        let mut k: C = factorial(n);
        vec.push(self.0.get(0).unwrap().ref_div(&k));

        for i in 1..=d {
            k = k * C::from_usize(n+i).unwrap() / C::from_usize(i).unwrap();
            vec.push(self.0.get(i).unwrap().ref_div(&k));
        }

        vec
    }

    pub(crate) fn integrate_n(&mut self, n: usize) {
        debug_assert!(n > 1);
        self.0 = self.new_n_integral_coeffs(n);
    }

    pub(crate) fn new_nth_integral(&self, n: usize) -> Polynomial<C> {
        debug_assert!(n > 1);
        Polynomial::new_raw_dense(self.new_n_integral_coeffs(n))
    }
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

    Polynomial::from(v)
}

fn clone_coeff<'a, C>(op_c: Option<&'a C>) -> C where C: Semiring + Clone {
    op_c.map_or_else(||C::zero(), |c|c.clone())
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

    Polynomial::from(v)
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

    Polynomial::from(v)
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

    Polynomial::from(v)
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

    Polynomial::from(v)
}

fn neg_coeff<'a, C>(op_c: Option<&'a C>) -> C where C: Ring + Clone {
    op_c.map_or_else(||C::zero(), |c|c.ref_neg())
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

    Polynomial::from(v)
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
        v.extend(lhs_iter.map(clone_coeff));
    } else {
        v.extend(rhs_iter.map(|c|-c));
    }

    Polynomial::from(v)
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
        v.extend(lhs_iter.map(clone_coeff));
    } else {
        v.extend(rhs_iter.map(neg_coeff));
    }

    Polynomial::from(v)
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
    (Polynomial::from(q), Polynomial::from(u))
}