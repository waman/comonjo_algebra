use std::{collections::BTreeMap, ops::Mul};

use num::{complex::{Complex32, Complex64}, BigInt, BigRational, Rational32, Rational64};

use crate::{algebra::{RefMul, Semiring}, poly::{CoeffsIterator, Polynomial}};


impl<C> Mul<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: Self) -> Self::Output {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), rhs) => {
                if lhs.0.is_one() { return rhs }
                rhs.map_nonzero(|_, c| lhs.0.ref_mul(c))
            },
            (lhs, Polynomial::Constant(rhs)) => {
                if rhs.0.is_one() { return lhs }
                lhs.map_nonzero(|_, c| c * &rhs.0)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_mul(&lhs, &rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_mul(&lhs, &rhs),
        }
    }
}

impl<'b, C> Mul<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: &'b Self) -> Self::Output {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), rhs) => {
                if lhs.0.is_one() { return rhs.clone() }
                rhs.ref_map_nonzero(|_, c| lhs.0.ref_mul(c))
            },
            (lhs, Polynomial::Constant(rhs)) => {
                if rhs.0.is_one() { return lhs }
                lhs.map_nonzero(|_, c| c * &rhs.0)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_mul(&lhs, rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_mul(&lhs, rhs),
        }
    }
}

impl<'a, C> Mul<Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: Polynomial<C>) -> Self::Output {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), rhs) => {
                if lhs.0.is_one() { return rhs }
                rhs.map_nonzero(|_, c| lhs.0.ref_mul(c))
            },
            (lhs, Polynomial::Constant(rhs)) => {
                if rhs.0.is_one() { return lhs.clone() }
                lhs.ref_map_nonzero(|_, c| c.ref_mul(&rhs.0))
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_mul(lhs, &rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_mul(lhs, &rhs),
        }
    }
}

impl<'a, 'b, C> Mul<&'b Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;
    
    fn mul(self, other: &'b Polynomial<C>) -> Polynomial<C> {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), rhs) => {
                if lhs.0.is_one() { return rhs.clone() }
                rhs.ref_map_nonzero(|_, c| lhs.0.ref_mul(c))
            },
            (lhs, Polynomial::Constant(rhs)) => {
                if rhs.0.is_one() { return lhs.clone() }
                lhs.ref_map_nonzero(|_, c| c.ref_mul(&rhs.0))
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_mul(lhs, rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_mul(lhs, rhs),
        }
    }
}

impl<C> RefMul<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_mul(&self, other: Polynomial<C>) -> Polynomial<C> { self * other }
}

impl<'b, C> RefMul<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_mul(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self * other }
}

//********** Mul (Scalar) **********/
impl<C> Mul<C> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, k: C) -> Self::Output { self * &k }
}

impl<'b, C> Mul<&'b C> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, k: &'b C) -> Self::Output {
        if k.is_zero() { return Polynomial::Zero() }
        if k.is_one() { return self }
        self.map_nonzero(|_, c| c.ref_mul(k))
    }
}

impl<'a, C> Mul<C> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, k: C) -> Self::Output { self * &k }
}

impl<'a, 'b, C> Mul<&'b C> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;
    
    fn mul(self, k: &'b C) -> Polynomial<C> {
        if k.is_zero() { return Polynomial::Zero() }
        if k.is_one() { return self.clone() }
        self.ref_map_nonzero(|_, c| c.ref_mul(k))
    }
}

macro_rules! impl_scale_by_left {
    ( $( $t:ident ),* ) => {
        $(

            impl Mul<Polynomial<$t>> for $t {

                type Output = Polynomial<$t>;

                fn mul(self, poly: Polynomial<$t>) -> Self::Output { poly.scale_by_left(&self) }
            }

            impl<'b> Mul<&'b Polynomial<$t>> for $t {

                type Output = Polynomial<$t>;

                fn mul(self, poly: &'b Polynomial<$t>) -> Self::Output { poly.ref_scale_by_left(&self) }
            }

            impl<'a> Mul<Polynomial<$t>> for &'a $t {

                type Output = Polynomial<$t>;

                fn mul(self, poly: Polynomial<$t>) -> Self::Output { poly.scale_by_left(self) }
            }

            impl<'a, 'b> Mul<&'b Polynomial<$t>> for &'a $t {

                type Output = Polynomial<$t>;

                fn mul(self, poly: &'b Polynomial<$t>) -> Self::Output { poly.ref_scale_by_left(self) }
            }
        )*
    };
}

impl_scale_by_left!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64,
        BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

fn dense_mul<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
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

fn sparse_mul<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {
    let mut map: BTreeMap<usize, C> = BTreeMap::new();

    for (i, x) in lhs.nonzero_coeffs() {
        for (j, y) in rhs.nonzero_coeffs() {
            let k = i + j;
            let z = x.ref_mul(y);
            match map.get_mut(&k){
                Some(c) => *c = c.ref_add(z),
                None => { map.insert(k, z); },
            }
        }
    }

    Polynomial::sparse_from_map(map)
}