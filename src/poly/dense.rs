use std::collections::BTreeMap;

use crate::{algebra::{Field, Ring, Semiring}, poly::Polynomial};

use super::CoeffsAccessor;

#[derive(Clone)]
pub struct DenseContent<C>(pub(crate) Vec<C>);

impl<C> DenseContent<C> where C: Semiring {

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(n)
    }

    pub fn map<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let v: Vec<D> = self.0.into_iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::dense_from_vec(v)
    }

    pub fn map_ref<D, F>(&self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, &C) -> D {
        let v: Vec<D> = self.0.iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::dense_from_vec(v)
    }

    pub(crate) fn to_map(mut self) -> BTreeMap<usize, C> {
        let mut map = BTreeMap::new();
        while let Some(c) = self.0.pop() {
            if !c.is_zero() {
                map.insert(self.0.len(), c);
            }
        }
        debug_assert!(self.0.is_empty());
        map
    }

    pub fn reductum(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn monic(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn derivative(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn integral(&self) -> Polynomial<C> {
        todo!()
    }
}

impl<C> DenseContent<C> where C: Ring {

    pub fn neg(self) -> Polynomial<C> {
        let v: Vec<C> = self.0.into_iter().map(|c|-c).collect();
        Polynomial::Dense(DenseContent(v))
    }
}

impl<C> DenseContent<C> where C: Ring + Clone {

    pub fn neg_ref(&self) -> Polynomial<C> {
        let v: Vec<C> = self.0.iter().map(|c|c.neg_ref()).collect();
        Polynomial::Dense(DenseContent(v))
    }
}


/// Return (max, min, first_arg_is_higher)
fn max_min(x: usize, y: usize) -> (usize, usize, bool) {
    if x >= y { (x, y, true) } else { (y, x, false) }
}

pub(crate) fn add_val<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
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

pub(crate) fn add_ref<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
        where C: Semiring + Clone {

    let (d_max, d_min, lhs_is_higher) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(some_x), Some(some_y)) => match (some_x, some_y) {
                (Some(x), Some(y)) => v.push(x.add_ref(y)),
                (Some(x), None) => v.push(x.clone()),
                (None, Some(y)) => v.push(y.clone()),
                (None, None) => v.push(C::zero()),
            },
            _ => panic!(),
        }
    }

    let rest_iter = if lhs_is_higher { lhs_iter } else { rhs_iter };
    v.extend(rest_iter.map(|c| match c {
        Some(x) => x.clone(),
        None => C::zero(),
    }));

    Polynomial::dense_from_vec(v)
}

pub(crate) fn sub_val<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> 
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

pub(crate) fn sub_ref<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
        where C: Ring + Clone {

    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs();
    let mut rhs_iter = rhs.coeffs();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(some_x), Some(some_y)) => match (some_x, some_y) {
                (Some(x), Some(y)) => v.push(x.clone() - y.clone()),
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
            Some(x) => -x.clone(),
            None => C::zero(),
        }));
    }

    Polynomial::dense_from_vec(v)
}

pub(crate) fn mul_val<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> 
        where C: Semiring + Clone {
    mul_ref(&lhs, &rhs)
}

pub(crate) fn mul_ref<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let deg = lhs.degree() + rhs.degree();
    let mut v: Vec<C> = Vec::with_capacity(deg + 1);

    for (i, x) in lhs.nonzero_coeffs() {
        for (j, y) in rhs.nonzero_coeffs() {
            let k = i + j;
            let z = x.clone() * y.clone();
            match v.get_mut(k) {
                Some(c) => *c = c.clone() + z,
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
            let q0: C = u_last.div_ref(v0);
            let offset = u.len() - d_rhs;  // the last of u is already popped
            for (i, c_rhs) in rhs.nonzero_coeffs() {
                if i == d_rhs { break; }
                let c_lhs = u.get_mut(offset + i).unwrap();
                *c_lhs = c_lhs.sub_ref(&(c_rhs.mul_ref(&q0)));
            }
            q.push(q0);
        }
    }

    q.reverse();
    (Polynomial::dense_from_vec(q), Polynomial::dense_from_vec(u))
}

// pub(crate) fn div_rem_val<C>(lhs: Polynomial<C>, rhs: &Polynomial<C>) -> (Polynomial<C>, Polynomial<C>) 
//         where C: Field + Clone {

//     let n = lhs.degree() - rhs.degree();
//     let v0: &C = rhs.max_order_term_coeff().unwrap();
//     let mut q: Vec<C> = Vec::with_capacity(n);
//     let mut u: Vec<C> = lhs.to_vec();

//     for _ in 0..n {
//         let u_last: &C = &u.pop().unwrap();
//         let q0: C = u_last.div_ref(v0);
//         let u_offset = u.len() - rhs.degree();  // the last of u is already popped
//         for (i, c) in rhs.nonzero_coeffs().filter(|e|e.0 != rhs.degree()) {
//             let del = c.mul_ref(&q0);
//             match u.get_mut(u_offset + i) {
//                 Some(c) => {
//                     *c = c.sub_ref(&del);
//                 },
//                 None => panic!(),
//             }
//         }
//         q.push(q0);
//     }

//     q.reverse();
//     (Polynomial::dense_from_vec(q), Polynomial::dense_from_vec(u))
// }

// pub(crate) fn div_rem_ref<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> (Polynomial<C>, Polynomial<C>) 
//         where C: Field + Clone {

//     let n = lhs.degree() - rhs.degree();
//     let v0: &C = rhs.max_order_term_coeff().unwrap();
//     let mut q: Vec<C> = Vec::with_capacity(n);
//     let mut u: Vec<C> = lhs.to_vec();

//     for _ in 0..n {
//         let u_last: &C = &u.pop().unwrap();
//         let q0: C = u_last.div_ref(v0);
//         let u_offset = u.len() - rhs.degree();  // the last of u is already popped
//         for (i, c) in rhs.nonzero_coeffs().filter(|e|e.0 != rhs.degree()) {
//             let del = c.mul_ref(&q0);
//             match u.get_mut(u_offset + i) {
//                 Some(c) => {
//                     *c = c.sub_ref(&del);
//                 },
//                 None => panic!(),
//             }
//         }
//         q.push(q0);
//     }

//     q.reverse();
//     (Polynomial::dense_from_vec(q), Polynomial::dense_from_vec(u))
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

// def reductum(implicit e: Eq[C], ring: Semiring[C], ct: ClassTag[C]): Polynomial[C] = {
//   var i = coeffs.length - 2
//   while (i >= 0 && coeffs(i) === ring.zero) i -= 1
//   if (i < 0) {
//     new PolyDense(new Array[C](0))
//   } else {
//     val arr = new Array[C](i + 1)
//     System.arraycopy(coeffs, 0, arr, 0, i + 1)
//     new PolyDense(arr)
//   }
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