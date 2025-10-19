use std::{collections::BTreeMap, ops::Sub};

use crate::{algebra::{AdditiveGroup, RefSub, Ring}, poly::{max_min, CoeffsIterator, Polynomial}};

impl<C> Sub for Polynomial<C> where C: Ring {

    type Output = Polynomial<C>;

    fn sub(self, other: Self) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => -rhs,

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 - rhs.0),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_sub_vv(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_sub_vv(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = (*c_lhs).ref_sub(c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_sub_vv(lhs, rhs),

            // Sparse
            (Polynomial::Sparse(mut sc), Polynomial::Constant(c_rhs)) => {
                if let Some(c_lhs) = sc.0.get_mut(&0) {
                    let dif = c_lhs.ref_sub(c_rhs.0);
                    if dif.is_zero() { sc.0.remove(&0); } else { *c_lhs = dif; }
                } else {
                    sc.0.insert(0, -c_rhs.0);
                }
                Polynomial::Sparse(sc)
            },
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_sub_vv(lhs, rhs),
        }
    }
}

impl<'b, C> Sub<&'b Polynomial<C>> for Polynomial<C> where C: Ring + Clone {

    type Output = Polynomial<C>;

    fn sub(self, other: &'b Self) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => rhs.ref_neg(),

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 - &rhs.0),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_sub_vr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_sub_vr(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = c_lhs.ref_sub(&c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_sub_vr(lhs, rhs),

            // Sparse
            (Polynomial::Sparse(mut sc), Polynomial::Constant(c_rhs)) => {
                if let Some(c_lhs) = sc.0.get_mut(&0) {
                    let dif = c_lhs.ref_sub(&c_rhs.0);
                    if dif.is_zero() { sc.0.remove(&0); } else { *c_lhs = dif; }
                } else {
                    sc.0.insert(0, c_rhs.0.ref_neg());
                }
                Polynomial::Sparse(sc)
            },
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_sub_vr(lhs, rhs),
        }
    }
}

impl<'a, C> Sub<Polynomial<C>> for &'a Polynomial<C> where C: Ring + Clone {

    type Output = Polynomial<C>;

    fn sub(self, other: Polynomial<C>) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => -rhs,

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_sub(rhs.0)),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_sub_rv(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_sub_rv(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense_sub_rv(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_sub_rv(lhs, rhs),
        }
    }
}

impl<'a, 'b, C> Sub<&'b Polynomial<C>> for &'a Polynomial<C> where C: Ring + Clone {

    type Output = Polynomial<C>;

    fn sub(self, other: &'b Polynomial<C>) -> Polynomial<C> {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => rhs.ref_neg(),

            // Const
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_sub(&rhs.0)),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_sub_rr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_sub_rr(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense_sub_rr(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_sub_rr(lhs, rhs),
        }
    }
}

impl<C> RefSub<Polynomial<C>> for Polynomial<C> where C: Ring + Clone {
    #[inline]
    fn ref_sub(&self, other: Polynomial<C>) -> Polynomial<C> { self + other }
}

impl<'b, C> RefSub<&'b Polynomial<C>> for Polynomial<C> where C: Ring + Clone {
    #[inline]
    fn ref_sub(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self + other }
}

pub(crate) fn dense_sub_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> 
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

pub(crate) fn dense_sub_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
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

pub(crate) fn dense_sub_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> 
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

fn dense_sub_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
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

//********** Sparse **********/
fn sparse_sub_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Ring {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1 - y.1;
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, -y.1);
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(lhs_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, -y.1);
        map.extend(rhs_iter.map(|(i, c)| (i, -c)));
    }

    Polynomial::sparse_from_map(map)
} 

fn sparse_sub_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1 - y.1;
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.ref_neg());
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(lhs_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.ref_neg());
        map.extend(rhs_iter.map(|(i, c)| (i, c.ref_neg())));
    }

    Polynomial::sparse_from_map(map)
} 

fn sparse_sub_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1.ref_sub(y.1);
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, -y.1);
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(lhs_iter.map(|(i, c)| (i, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, -y.1);
        map.extend(rhs_iter.map(|(i, c)| (i, -c)));
    }

    Polynomial::sparse_from_map(map)
} 

fn sparse_sub_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'a Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1.ref_sub(y.1);
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.ref_neg());
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(lhs_iter.map(|(e, c)|(e, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.ref_neg());
        map.extend(rhs_iter.map(|(e, c)|(e, c.ref_neg())));
    }

    Polynomial::sparse_from_map(map)
} 