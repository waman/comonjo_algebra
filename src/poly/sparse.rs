use std::collections::{BTreeMap, btree_map::{IntoIter, Iter}};

use crate::{algebra::{EuclideanRing, Field, Ring, Semiring}, poly::{CoeffsIterator, Polynomial, SemiringPolyOps, iter::{CoeffsIter, IntoCoeffsIter, IntoNonzeroCoeffsIter, NonzeroCoeffsIter}, mul_div_uint}};

#[derive(Clone)]
pub struct SparseCoeffs<C>(pub(crate) BTreeMap<usize, C>);

impl<C> SparseCoeffs<C> where C: Semiring {

    pub(crate) fn degree(&self) -> usize {
        *self.0.last_key_value().unwrap().0
    }

    pub(crate) fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(&n)
    }

    pub(crate) fn max_order_term(&self) -> Option<(usize, &C)> {
        let last = self.0.last_key_value().unwrap();
        Some((*last.0, last.1))
    }

    pub(crate) fn min_order_term(&self) -> Option<(usize, &C)> {
        let first = self.0.first_key_value().unwrap();
        Some((*first.0, first.1))
    }

    pub(crate) fn is_x(&self) -> bool {
        self.0.len() == 1 && self.0.get(&1).is_some_and(|c|c.is_one())
    }
 
//    final private def expBits(x: C)(implicit ring: Semiring[C]): Array[C] = {
//      val bits = new Array[C](math.max(2, 32 - numberOfLeadingZeros(degree)))
//      bits(0) = x
//      // we use pow(2) here for the benefit of Interval[_], where
//      // x.pow(2) has better error bounds than than (x * x).
//      if (bits.length > 1) bits(1) = x.pow(2)
//      cfor(2)(_ < bits.length, _ + 1) { i =>
//        val prev = bits(i - 1)
//        bits(i) = prev * prev
//      }
//      bits
//    }
 
//    @tailrec
//    final private def fastExp(bits: Array[C], e: Int, i: Int, acc: C)(implicit ring: Semiring[C]): C = {
//      if (e == 0) acc
//      else {
//        val lb = numberOfTrailingZeros(e) + 1
//        val j = i + lb
//        fastExp(bits, e >>> lb, j, acc * bits(j - 1))
//      }
//    }
 
//    final private def fastExp(bits: Array[C], e: Int)(implicit ring: Semiring[C]): C = {
//      val lb = numberOfTrailingZeros(e) + 1
//      fastExp(bits, e >>> lb, lb, bits(lb - 1))
//    }

    pub(crate) fn reciprocal(&mut self) -> Option<Polynomial<C>>{
        let d = self.degree();
        let map0 = &mut self.0;
        let keys: Vec<usize> = map0.keys().map(|i|*i).collect();
        let mut map: BTreeMap<usize, C> = BTreeMap::new();

        for i in keys {
            map.insert(d - i, map0.remove(&i).unwrap());
        }
        Some(Polynomial::sparse_from_map(map))
    }

    pub(crate) fn remove_zero_roots(&mut self) -> Option<Polynomial<C>> {
        let i0 = self.min_order_term().unwrap().0;
        if i0 == 0 {
            None
        } else {
            let map0 = &mut self.0;
            let keys: Vec<usize> = map0.keys().map(|i|*i).collect();
            let mut map: BTreeMap<usize, C> = BTreeMap::new();

            for i in keys {
                map.insert(i - i0, map0.remove(&i).unwrap());
            }
            Some(Polynomial::sparse_from_map(map))
        }
    }

    pub(crate) fn reductum(&mut self) -> Option<Polynomial<C>> {
        let dim = self.degree();
        self.0.remove(&dim);
        match self.0.len() {
            0 => Some(Polynomial::Zero()),
            1 if self.0.contains_key(&0) => 
                Some(Polynomial::new_raw_const(self.0.remove(&0).unwrap())),
            _ => None,
        }
    }
}

impl<C> SparseCoeffs<C> where C: Semiring + Clone {
    
    pub(crate) fn new_reciprocal(&self) -> Polynomial<C> {
        let d = self.degree();
        let map: BTreeMap<usize, C> = self.0.iter().map(|(i, c)| (d-i, c.clone())).collect();
        Polynomial::sparse_from_map(map)
    }

    pub(crate) fn new_zero_roots_removed(&self) -> Polynomial<C> {
        let i0 = *self.0.first_key_value().unwrap().0;
        if i0 == 0 {
            Polynomial::Sparse(SparseCoeffs(self.0.clone()))
        } else {
            let map: BTreeMap<usize, C> = self.0.iter().map(|(i, c)|(i-i0, c.clone())).collect();
            Polynomial::sparse_from_map(map)
        }
    }

    pub(crate) fn new_reductum(&self) -> Polynomial<C> {
        let n = self.0.len();
        let map: BTreeMap<usize, C> = self.0.iter().take(n-1).map(|(i, c)| (*i, c.clone())).collect();
        Polynomial::sparse_from_map(map)
    }
}

//********** Iterator **********/
impl<C> SparseCoeffs<C> where C: Semiring {
    
    pub(crate) fn into_nonzero_coeffs_iter(self) -> IntoNonzeroCoeffsIter<C> {
        IntoNonzeroCoeffsIter::Sparse(SIntoNonzeroCoeffsIter(self.0.into_iter()))
    }

    pub(crate) fn nonzero_coeffs_iter<'a>(&'a self) -> NonzeroCoeffsIter<'a, C> {
        NonzeroCoeffsIter::Sparse(SNonzeroCoeffsIter(self.0.iter()))
    }

    pub(crate) fn into_coeffs_iter(self) -> IntoCoeffsIter<C> {
        let mut map_iter = self.0.into_iter();
        IntoCoeffsIter::Sparse(SIntoCoeffsIter{ index: 0, current: map_iter.next(), map_iter })
    }

    pub(crate) fn coeffs_iter<'a>(&'a self) -> CoeffsIter<'a, C> {
        let mut map_iter = self.0.iter();
        CoeffsIter::Sparse(SCoeffsIter{ index: 0, current: map_iter.next(), map_iter})
    }
}

pub struct SIntoNonzeroCoeffsIter<C>(pub(crate) IntoIter<usize, C>) where C: Semiring;

impl<C> Iterator for SIntoNonzeroCoeffsIter<C> where C: Semiring {
    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> { self.0.next() }

    fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
}

pub struct SNonzeroCoeffsIter<'a, C>(pub(crate) Iter<'a, usize, C>) where C: Semiring;

impl<'a, C> Iterator for SNonzeroCoeffsIter<'a, C> where C: Semiring {

    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> { self.0.next().map(|(i, c)| (*i, c)) }

    fn size_hint(&self) -> (usize, Option<usize>) { self.0.size_hint() }
}

pub struct SIntoCoeffsIter<C> where C: Semiring {
    pub(crate) index: usize,
    pub(crate) current: Option<(usize, C)>,
    pub(crate) map_iter: IntoIter<usize, C>,
}

impl<C> Iterator for SIntoCoeffsIter<C> where C: Semiring {

    type Item = C;

    fn next(&mut self) -> Option<Self::Item> {
        let result = match &mut self.current {
            Some(c) => 
                if c.0 == self.index {
                    let r = match self.map_iter.next() {
                        Some(nxt) => self.current.replace(nxt),
                        None => self.current.take(),
                    };
                    Some(r.unwrap().1)
                } else { 
                    Some(C::zero())
                }
            _ => None,
        };

        self.index += 1;
        return result;
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { self.map_iter.size_hint() }
}

/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some(Some(0))</code> if the coefficient is zero.
pub struct SCoeffsIter<'a, C> where C: Semiring {
    pub(crate) index: usize,
    pub(crate) current: Option<(&'a usize, &'a C)>,
    pub(crate) map_iter: Iter<'a, usize, C>,
}

impl<'a, C> Iterator for SCoeffsIter<'a, C> where C: Semiring {

    type Item = Option<&'a C>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(c) = self.current {
            let result = if *c.0 == self.index { 
                let r = Some(c.1);
                self.current = self.map_iter.next();
                r
            } else {
                None
            };
            self.index += 1;
            Some(result)
        } else {
            return None;
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { self.map_iter.size_hint() }
}

impl<C> SemiringPolyOps<C> for SparseCoeffs<C> where C: Semiring {

    fn to_vec(mut self) -> Vec<C> {
        let n = self.degree() + 1;
        let mut v: Vec<C> = Vec::with_capacity(n);
        for i in 0..n {
            match self.0.remove(&i) {
                Some(e) => v.push(e),
                None => v.push(C::zero()),
            } 
        }
        debug_assert!(self.0.is_empty());
        v
    }
    
    fn to_map(self) -> BTreeMap<usize, C> { self.0 }

    fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let m: BTreeMap<usize, D> = self.0.into_iter().map(|(i, c)| (i, f(i, c))).collect();
        Polynomial::sparse_from_map(m)
    }
}

impl<'a, C> SemiringPolyOps<C> for &'a SparseCoeffs<C> where C: Semiring + Clone {
    
    fn to_vec(self) -> Vec<C> {
        self.coeffs_iter().map(|c|{
            match c {
                Some(v) => v.clone(),
                _ => C::zero(),
            }
        }).collect()
    }
    
    fn to_map(self) -> BTreeMap<usize, C> { self.0.clone() }

    fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let m: BTreeMap<usize, D> = self.0.iter().map(|(i, c)| (*i, f(*i, c.clone()))).collect();
        Polynomial::sparse_from_map(m)
    }
}

impl<C> SparseCoeffs<C> where C: Semiring + num::FromPrimitive {
    
    pub(crate) fn new_derivative(&self) -> Polynomial<C> {
        let ite = self.0.iter();
        let map: BTreeMap<usize, C> = 
            ite.filter(|(i, _)| **i != 0)
            .map(|(i,c)| (i-1, C::from_usize(*i).unwrap() * c))
            .collect();
        Polynomial::sparse_from_map(map)
    }
}

impl<C> SparseCoeffs<C> where C: Ring {

    pub(crate) fn flip(&mut self){
        for (i, c) in self.0.iter_mut() {
            if i % 2 != 0 {
                *c = c.ref_neg();
            }
        }
    }
}

impl<C> SparseCoeffs<C> where C: Ring + Clone {

    pub(crate) fn new_flipped(&self) -> Polynomial<C> {
        let map: BTreeMap<usize, C> =
            self.0.iter().map(|(i, c)|{
                if i % 2 == 0 { (*i, c.clone()) } else { (*i, c.ref_neg()) }
            }).collect();
        Polynomial::Sparse(SparseCoeffs(map))
    }
}
    
impl<C> SparseCoeffs<C> where C: Field + num::FromPrimitive {

    pub(crate) fn integral(&self) -> Polynomial<C> {
        let map: BTreeMap<usize, C> = 
            self.0.iter().map(|(i, c)| {
                let j = i+1;
                let d = c.ref_div(C::from_usize(j).unwrap());
                (j, d)
            }).collect();
        Polynomial::sparse_from_map(map)
    }
}

impl<C> SparseCoeffs<C> where C: EuclideanRing + num::FromPrimitive + num::Integer + Clone {
    
    // ** From spire code *****
    // The trick here came from this answer:
    //   http://math.stackexchange.com/questions/694565/polynomial-shift
    // This is a heavily optimized version of the same idea. This is fairly
    // critical method to be fast, since it is the most expensive part of the
    // VAS root isolation algorithm.
    pub(crate) fn shift(&self, h: C) -> Polynomial<C> {
        let mut coeffs: BTreeMap<usize, C> = self.0.clone();
        for (deg, c) in self.nonzero_coeffs_iter() {
            if deg == 0 { continue; }
            let mut i: C = C::one();
            let mut d: usize = deg;
            let mut m: C = C::one();
            let mut k: C = c.clone();
            while d > 0 {
                m = mul_div_uint(m, d, i.clone());  // m * d / i
                k = k * &h;
                let dif = m.ref_mul(&k);
                coeffs.entry(d-1).and_modify(|v| *v = v.ref_add(&dif)).or_insert(dif);
                d = d - 1;
                i = i + C::one();
            }
        }
        Polynomial::sparse_from_map(coeffs)
    }
}

impl<C> SparseCoeffs<C> where C: Field + num::FromPrimitive + Clone {
    
    pub(crate) fn shift_f(&self, h: C) -> Polynomial<C> {
        let mut coeffs: BTreeMap<usize, C> = self.0.clone();
        for (deg, c) in self.nonzero_coeffs_iter() {
            if deg == 0 { continue; }
            let mut i: C = C::one();
            let mut d: usize = deg;
            let mut m: C = C::one();
            let mut k: C = c.clone();
            while d > 0 {
                m = m * C::from_usize(d).unwrap() / i.clone();  // m * d / i
                k = k * &h;
                let dif = m.ref_mul(&k);
                coeffs.entry(d-1).and_modify(|v| *v = v.ref_add(&dif)).or_insert(dif);
                d = d - 1;
                i = i + C::one();
            }
        }
        Polynomial::sparse_from_map(coeffs)
    }
}

//********** Add **********/
pub(crate) fn add_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Semiring {

    let mut map = BTreeMap::new();

    let mut x_iter = lhs.nonzero_coeffs();
    let mut y_iter = rhs.nonzero_coeffs();

    let mut x_next = x_iter.next();
    let mut y_next = y_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1 + y.1;
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = x_iter.next();
            y_next = y_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = x_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1);
            
            x_next = Some(x);
            y_next = y_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(x_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1);
        map.extend(y_iter);
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn add_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut x_iter = lhs.nonzero_coeffs();
    let mut y_iter = rhs.nonzero_coeffs();

    let mut x_next = x_iter.next();
    let mut y_next = y_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1 + y.1;
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = x_iter.next();
            y_next = y_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = x_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.clone());
            
            x_next = Some(x);
            y_next = y_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(x_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.clone());
        map.extend(y_iter.map(|(i, c)| (i, c.clone())));
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn add_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut x_iter = lhs.nonzero_coeffs();
    let mut y_iter = rhs.nonzero_coeffs();

    let mut x_next = x_iter.next();
    let mut y_next = y_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1.ref_add(y.1);
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = x_iter.next();
            y_next = y_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = x_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1);
            
            x_next = Some(x);
            y_next = y_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(x_iter.map(|(i, c)| (i, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1);
        map.extend(y_iter);
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn add_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'a Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut x_iter = lhs.nonzero_coeffs();
    let mut y_iter = rhs.nonzero_coeffs();

    let mut x_next = x_iter.next();
    let mut y_next = y_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1.ref_add(y.1);
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = x_iter.next();
            y_next = y_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = x_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.clone());
            
            x_next = Some(x);
            y_next = y_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(x_iter.map(|(e, c)|(e, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.clone());
        map.extend(y_iter.map(|(e, c)|(e, c.clone())));
    }

    Polynomial::sparse_from_map(map)
} 

//********** Sub **********/
pub(crate) fn sub_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Ring {

    let mut map = BTreeMap::new();

    let mut x_iter = lhs.nonzero_coeffs();
    let mut y_iter = rhs.nonzero_coeffs();

    let mut x_next = x_iter.next();
    let mut y_next = y_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1 - y.1;
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = x_iter.next();
            y_next = y_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = x_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, -y.1);
            
            x_next = Some(x);
            y_next = y_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(x_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, -y.1);
        map.extend(y_iter.map(|(i, c)| (i, -c)));
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn sub_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut x_iter = lhs.nonzero_coeffs();
    let mut y_iter = rhs.nonzero_coeffs();

    let mut x_next = x_iter.next();
    let mut y_next = y_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1 - y.1;
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = x_iter.next();
            y_next = y_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = x_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.ref_neg());
            
            x_next = Some(x);
            y_next = y_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(x_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.ref_neg());
        map.extend(y_iter.map(|(i, c)| (i, c.ref_neg())));
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn sub_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut x_iter = lhs.nonzero_coeffs();
    let mut y_iter = rhs.nonzero_coeffs();

    let mut x_next = x_iter.next();
    let mut y_next = y_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1.ref_sub(y.1);
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = x_iter.next();
            y_next = y_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = x_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, -y.1);
            
            x_next = Some(x);
            y_next = y_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(x_iter.map(|(i, c)| (i, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, -y.1);
        map.extend(y_iter.map(|(i, c)| (i, -c)));
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn sub_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'a Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut x_iter = lhs.nonzero_coeffs();
    let mut y_iter = rhs.nonzero_coeffs();

    let mut x_next = x_iter.next();
    let mut y_next = y_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1.ref_sub(y.1);
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = x_iter.next();
            y_next = y_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = x_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.ref_neg());
            
            x_next = Some(x);
            y_next = y_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(x_iter.map(|(e, c)|(e, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.ref_neg());
        map.extend(y_iter.map(|(e, c)|(e, c.ref_neg())));
    }

    Polynomial::sparse_from_map(map)
} 

//********** Mul **********/
pub(crate) fn mul<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
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

//********** Div & Rem **********/
pub(crate) fn div_rem<C>(mut u: BTreeMap<usize, C>, rhs: &Polynomial<C>) -> (Polynomial<C>, Polynomial<C>) 
        where C: Field + Clone {

    let d_rhs = rhs.degree();
    let v0: &C = rhs.max_order_term().unwrap().1;
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