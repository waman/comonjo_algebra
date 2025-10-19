use std::{collections::BTreeMap, ops::{Div, Rem}};

use num::traits::Euclid;

use crate::{algebra::{Field, RefDiv, RefRem}, poly::{CoeffsIterator, Polynomial}};

impl<C> Polynomial<C> where C: Field + Clone {

    pub fn euclidean_fn(&self) -> usize { self.degree() }

    #[inline]
    fn div_rem_val<'a>(self, other: &'a Polynomial<C>) -> (Self, Self) {
        match (self, other) {
            (_, Polynomial::Zero()) => panic!("Can't divide by polynomial of zero!"),
            (Polynomial::Zero(), _) => (Polynomial::Zero(), Polynomial::Zero()),
            (lhs, Polynomial::Constant(rhs)) => (lhs.map_nonzero(|_, c| c / (&rhs.0)), Polynomial::Zero()),
            (lhs @ Polynomial::Constant(_), _) => (Polynomial::Zero(), lhs),
            (lhs @ Polynomial::Dense(_), rhs) => dense_div_rem(lhs.to_vec(), rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_div_rem(lhs.to_map(), rhs),
        }
    }

    #[inline]
    fn div_rem_ref(&self, other: &Self) -> (Self, Self) {
        match (self, other) {
            (_, Polynomial::Zero()) => panic!("Can't divide by polynomial of zero!"),
            (Polynomial::Zero(), _) => (Polynomial::Zero(), Polynomial::Zero()),
            (lhs, Polynomial::Constant(rhs)) => (lhs.ref_map_nonzero(|_, c| c.ref_div(&rhs.0)), Polynomial::Zero()),
            (lhs @ Polynomial::Constant(_), _) => (Polynomial::Zero(), lhs.clone()),
            (lhs @ Polynomial::Dense(_), rhs) => dense_div_rem(lhs.clone_to_vec(), rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_div_rem(lhs.clone_to_map(), rhs),
        }
    }
}

impl<C> Euclid for Polynomial<C> where C: Field + Clone {

    fn div_euclid(&self, other: &Self) -> Self { self / other }
    fn rem_euclid(&self, other: &Self) -> Self { self % other }
    fn div_rem_euclid(&self, other: &Self) -> (Self, Self) { self.div_rem_ref(other) }
}


impl<C> Div for Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn div(self, other: Self) -> Self::Output { self.div_rem_val(&other).0 }
}

impl<'b, C> Div<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn div(self, other: &'b Polynomial<C>) -> Self::Output { (&self).div_rem_ref(other).0 }
}

impl<'a, C> Div<Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn div(self, other: Polynomial<C>) -> Self::Output { self.div_rem_ref(&other).0 }
}

impl<'a, 'b, C> Div<&'b Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn div(self, other: &'b Polynomial<C>) -> Self::Output { self.div_rem_ref(other).0 }
}

impl<C> RefDiv<Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_div(&self, other: Polynomial<C>) -> Polynomial<C> { self / other }
}

impl<'b, C> RefDiv<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_div(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self / other }
}


impl<C> Rem for Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;
    
    fn rem(self, other: Self) -> Self::Output { self.div_rem_val(&other).1 }
}

impl<'b, C> Rem<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn rem(self, other: &'b Polynomial<C>) -> Self::Output { (&self).div_rem_ref(other).1 }
}

impl<'a, C> Rem<Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn rem(self, other: Polynomial<C>) -> Self::Output { self.div_rem_ref(&other).1 }
}

impl<'a, 'b, C> Rem<&'b Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn rem(self, other: &'b Polynomial<C>) -> Self::Output { self.div_rem_ref(other).1 }
}

impl<C> RefRem<Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_rem(&self, other: Polynomial<C>) -> Polynomial<C> { self % other }
}

impl<'b, C> RefRem<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_rem(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self % other }
}

fn dense_div_rem<C>(mut u: Vec<C>, rhs: &Polynomial<C>) -> (Polynomial<C>, Polynomial<C>) 
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

fn sparse_div_rem<C>(mut u: BTreeMap<usize, C>, rhs: &Polynomial<C>) -> (Polynomial<C>, Polynomial<C>) 
        where C: Field + Clone {

    let d_rhs = rhs.degree();
    let v0: &C = rhs.max_order_term_coeff().unwrap();
    let mut q: BTreeMap<usize, C> = BTreeMap::new();

    while !u.is_empty() {
        let i_last = *u.last_key_value().unwrap().0;
        if i_last < d_rhs { break; }
        let c_last = u.remove(&i_last).unwrap();

        let q0: C = c_last.ref_div(v0);
        let offset = i_last - d_rhs;  // the last of u is already removed
        for (i, c_rhs) in rhs.nonzero_coeffs() {
            if i == d_rhs { break; }
            let e = u.entry(offset + i).or_insert(C::zero());
            *e = e.ref_sub(c_rhs.ref_mul(&q0));
        }
        q.insert(i_last - d_rhs, q0);
    }

    (Polynomial::sparse_from_map(q), Polynomial::sparse_from_map(u))
}