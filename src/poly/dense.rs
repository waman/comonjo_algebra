use std::collections::BTreeMap;

use crate::{algebra::{Field, Ring, Semiring}, poly::{CoeffsIterator, Polynomial}};

#[derive(Clone)]
pub struct DenseContent<C>(pub(crate) Vec<C>);

impl<C> DenseContent<C> where C: Semiring {

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(n)
    }

    pub fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let v: Vec<D> = self.0.into_iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::dense_from_vec(v)
    }

    pub fn ref_map_nonzero<D, F>(&self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, &C) -> D {
        let v: Vec<D> = self.0.iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::dense_from_vec(v)
    }

    pub fn to_map(mut self) -> BTreeMap<usize, C> {
        let mut map = BTreeMap::new();
        while let Some(c) = self.0.pop() {
            if !c.is_zero() {
                map.insert(self.0.len(), c);
            }
        }
        debug_assert!(self.0.is_empty());
        map
    }
}

impl<C> DenseContent<C> where C: Semiring {

    pub fn reductum(mut self) -> Polynomial<C> {
        self.0.pop();
        Polynomial::dense_from_vec(self.0)
    }

    // pub fn shift(&mut self, h: C) {
    //     todo!()
    // }

    // pub fn differentiate(&mut self) {
    //     todo!()
    // }

    // pub fn integrate(&mut self) {
    //     todo!()
    // }

    // pub fn remove_zero_roots(&mut self) {
    //     todo!()
    // }

    // pub fn flip(&mut self) {
    //     todo!()
    // }

    // pub fn reciprocal(&mut self){
    //     // self.0.reverse();
    //     // self.0.d
    //     todo!()
    // }

    

// def *(rhs: Polynomial[C])(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] = {
//   if (rhs.isZero) return rhs
//   if (lhs.isZero) return lhs
//   val lcs = lhs.coeffsArray
//   val rcs = rhs.coeffsArray
//   val cs = new Array[C](lcs.length + rcs.length - 1)
//   cfor(0)(_ < cs.length, _ + 1) { i => cs(i) = ring.zero }
//   cfor(0)(_ < lcs.length, _ + 1) { i =>
//     val c = lcs(i)
//     var k = i
//     cfor(0)(_ < rcs.length, _ + 1) { j =>
//       cs(k) += c * rcs(j)
//       k += 1
//     }
//   }
//   Polynomial.dense(cs)
// }

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

// def derivative(implicit ring: Ring[C], eq: Eq[C]): Polynomial[C] = {
//   if (isZero) return this
//   val cs = new Array[C](degree)
//   var j = coeffs.length - 1
//   cfor(cs.length - 1)(_ >= 0, _ - 1) { i =>
//     cs(i) = ring.fromInt(j) * coeffs(j)
//     j -= 1
//   }
//   Polynomial.dense(cs)
// }

// def integral(implicit field: Field[C], eq: Eq[C]): Polynomial[C] = {
//   val cs = new Array[C](coeffs.length + 1)
//   cs(0) = field.zero
//   cfor(0)(_ < coeffs.length, _ + 1) { i => cs(i + 1) = coeffs(i) / field.fromInt(i + 1) }
//   Polynomial.dense(cs)
// }

// def *:(k: C)(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] =
//   if (k === ring.zero) {
//     Polynomial.dense(new Array[C](0))
//   } else {
//     val cs = new Array[C](coeffs.length)
//     cfor(0)(_ < cs.length, _ + 1) { i =>
//       cs(i) = k * coeffs(i)
//     }
//     Polynomial.dense(cs)
//   }
// }
}

impl<C> DenseContent<C> where C: Semiring + Clone {

    pub fn new_reductum(&self) -> Polynomial<C> {
        let n = self.0.len();
        let vec: Vec<C> = self.0.iter().take(n-1).map(|c|c.clone()).collect();
        Polynomial::dense_from_vec(vec)
    }

    // pub fn new_shifted(&self, h: C) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn new_derivative(&self) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn new_integral(&self) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn new_zero_roots_removed(&self) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn new_flipped(&self) -> Polynomial<C> {
    //     todo!()
    // }
    // pub fn new_reciprocal(&self) -> Polynomial<C> {
    //     todo!()
    // }
}

/// Return (max, min, left_is_higher)
fn max_min(x: usize, y: usize) -> (usize, usize, bool) {
    if x >= y { (x, y, true) } else { (y, x, false) }
}


//********** Add **********/
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