use std::{collections::BTreeMap, ops::Add};

use crate::{algebra::{RefAdd, Semiring}, poly::{max_min, CoeffsIterator, Polynomial}};

impl<C> Add for Polynomial<C> where C: Semiring {

    type Output = Polynomial<C>;

    fn add(self, other: Self) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => rhs,

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 + rhs.0),
            (Polynomial::Constant(c_lhs), Polynomial::Dense(mut dc)) => {
                dc.0.get_mut(0).map(|c_rhs| *c_rhs = c_lhs.0 + &*c_rhs);
                Polynomial::Dense(dc)
            },
            (Polynomial::Constant(c_lhs), Polynomial::Sparse(mut sc)) => {
                if let Some(c_rhs) = sc.0.get_mut(&0) {
                    let sum = c_lhs.0 + &*c_rhs;
                    if sum.is_zero() { sc.0.remove(&0); } else { *c_rhs = sum; }
                } else {
                    sc.0.insert(0, c_lhs.0);
                }
                Polynomial::Sparse(sc)
            },

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = c_lhs.ref_add(c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_add_vv(lhs, rhs),

            // Sparse
            (Polynomial::Sparse(mut sc), Polynomial::Constant(c_rhs)) => {
                if let Some(c_lhs) = sc.0.get_mut(&0) {
                    let sum = c_lhs.ref_add(c_rhs.0);
                    if sum.is_zero() { sc.0.remove(&0); } else { *c_lhs = sum; }
                } else {
                    sc.0.insert(0, c_rhs.0);
                }
                Polynomial::Sparse(sc)
            },
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_add_vv(lhs, rhs),
        }
    }
}

impl<'b, C> Add<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn add(self, other: &'b Self) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => rhs.clone(),

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 + &rhs.0),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_add_vr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_add_vr(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = c_lhs.ref_add(&c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_add_vr(lhs, rhs),

            // Sparse
            (Polynomial::Sparse(mut sc), Polynomial::Constant(c_rhs)) => {
                if let Some(c_lhs) = sc.0.get_mut(&0) {
                    let sum = c_lhs.ref_add(&c_rhs.0);
                    if sum.is_zero() { sc.0.remove(&0); } else { *c_lhs = sum; }
                } else {
                    sc.0.insert(0, c_rhs.0.clone());
                }
                Polynomial::Sparse(sc)
            },
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_add_vr(lhs, rhs),
        }
    }
}

impl<'a, C> Add<Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn add(self, other: Polynomial<C>) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => rhs,

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_add(rhs.0)),
            (Polynomial::Constant(c_lhs), Polynomial::Dense(mut dc)) => {
                dc.0.get_mut(0).map(|c_rhs| *c_rhs = c_lhs.0.ref_add(&*c_rhs));
                Polynomial::Dense(dc)
            },
            (Polynomial::Constant(c_lhs), Polynomial::Sparse(mut sc)) => {
                if let Some(c_rhs) = sc.0.get_mut(&0) {
                    let sum = c_lhs.0.ref_add(&*c_rhs);
                    if sum.is_zero() { sc.0.remove(&0); } else { *c_rhs = sum; }
                } else {
                    sc.0.insert(0, c_lhs.0.clone());
                }
                Polynomial::Sparse(sc)
            },

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense_add_rv(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_add_rv(lhs, rhs),
        }
    }
}

impl<'a, 'b, C> Add<&'b Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn add(self, other: &'b Polynomial<C>) -> Polynomial<C> {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => rhs.clone(),

            // Const
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_add(&rhs.0)),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_add_rr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_add_rr(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense_add_rr(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_add_rr(lhs, rhs),
        }
    }
}

impl<C> RefAdd<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_add(&self, other: Polynomial<C>) -> Polynomial<C> { self + other }
}

impl<'b, C> RefAdd<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_add(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self + other }
}

//********** Dense **********/
fn dense_add_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
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

fn dense_add_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
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

fn dense_add_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
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

fn dense_add_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
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

//********** Sparse **********/
fn sparse_add_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Semiring {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1 + y.1;
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1);
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(lhs_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1);
        map.extend(rhs_iter);
    }

    Polynomial::sparse_from_map(map)
} 

fn sparse_add_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1 + y.1;
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.clone());
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(lhs_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.clone());
        map.extend(rhs_iter.map(|(i, c)| (i, c.clone())));
    }

    Polynomial::sparse_from_map(map)
} 

fn sparse_add_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1.ref_add(y.1);
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1);
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(lhs_iter.map(|(i, c)| (i, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1);
        map.extend(rhs_iter);
    }

    Polynomial::sparse_from_map(map)
} 

fn sparse_add_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'a Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1.ref_add(y.1);
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.clone());
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(lhs_iter.map(|(e, c)|(e, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.clone());
        map.extend(rhs_iter.map(|(e, c)|(e, c.clone())));
    }

    Polynomial::sparse_from_map(map)
} 