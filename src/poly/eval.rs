use std::collections::BTreeMap;

use num::{BigInt, BigRational, BigUint, FromPrimitive, Rational32, Rational64, complex::{Complex32, Complex64}, pow::Pow};

use crate::{algebra::{Field, Semiring}, poly::{Polynomial, dense_coeffs::DenseCoeffs, sparse_coeffs::SparseCoeffs}};

pub struct Eval;

pub trait PolynomialEvaluator<C> where C: Semiring {
    fn eval(p: &Polynomial<C>, x: C) -> C;
}

fn eval_dense<C: Semiring>(dc: &DenseCoeffs<C>, x: C) -> C {
    let x2: &C = &(x.ref_mul(&x));
    let mut sum0: C = C::zero();
    let mut sum1: C = C::zero();

    let mut rchunks = dc.0.rchunks_exact(2);
    for cs in rchunks.by_ref() {
        sum0 = sum0 * x2 + &cs[0];
        sum1 = sum1 * x2 + &cs[1];
    }

    let rem = rchunks.remainder();
    match rem.len() {
        0 => sum0 + sum1 * x,
        1 => sum0 * x + sum1 * x2 + &rem[0],
        _ => panic!()
    }
}

fn eval_sparse<C>(coeffs: &SparseCoeffs<C>, x: C, pow: fn(C, usize) -> C) -> C 
        where C: Semiring + Clone {

    if coeffs.0.len() == 1 {
        match coeffs.max_order_term() {
            Some((i, c)) => return c.ref_mul(pow(x, i)),
            _ => panic!(),
        }
    }

    let mut ite = coeffs.0.iter().rev();

    let last = ite.next().unwrap();
    let mut prev_i = *last.0;
    let mut sum: C = last.1.clone();

    let calc = &mut PowerCalculator::new(x);
    let cache = &mut PowerCache::new();

    for (i, c) in ite {
        sum = sum * cache.get(prev_i - *i, calc) + c;
        prev_i = *i;
    }

    if prev_i == 0 {
        sum
    } else {
        sum * cache.get(prev_i, calc)
    }
}

struct PowerCache<C>(BTreeMap<usize, C>) where C: Semiring;

impl<C> PowerCache<C> where C: Semiring {

    fn new() -> PowerCache<C> { PowerCache(BTreeMap::new()) }

    fn get(&mut self, i: usize, calc: &mut PowerCalculator<C>) -> &C {
        self.0.entry(i).or_insert_with(|| calc.get(i))
    }
}

struct PowerCalculator<C> where C: Semiring {
    /// *x^{2^i}*
    exp_bits: Vec<C>
}

impl<C> PowerCalculator<C> where C: Semiring {

    fn new(x: C) -> PowerCalculator<C> {
        let x2 = x.ref_mul(&x);
        PowerCalculator { exp_bits: vec![x, x2] }
    }

    /// Returns *x^i*
    fn get(&mut self, i: usize) -> C {
        let mut exp: usize = i;
        let mut result: C = C::one();
        let mut j: usize = 0;

        while exp > 0 {
            let lb = exp.trailing_zeros() as usize + 1;
            j = j + lb;
            result = result * self.get_base_power(j-1);
            exp >>= lb;
        }
        
        result
    }

    /// Returns *x^{2^i}*.
    fn get_base_power(&mut self, i: usize) -> &C {
        let base_pows = &mut self.exp_bits;
        let len = base_pows.len();
        if i >= len {
            for _ in len..=i {
                let next_last: C;
                {  // scope of the ref of the last
                    let last: &C = base_pows.last().unwrap();
                    next_last = last.ref_mul(last);
                }
                base_pows.push(next_last);
            }
            base_pows.last().unwrap()
        } else {
            base_pows.get(i).unwrap()
        }
    }
}

//***** For Integer (usize-pow) *****/
macro_rules! eval_impl_for_pow_usize {
    ( $( $t:ident ),* ) => {
        $(
            impl PolynomialEvaluator<$t> for Eval {

                fn eval(p: &Polynomial<$t>, x: $t) -> $t {
                    eval_int(p, x)
                }
            }
        )*
    };
}

eval_impl_for_pow_usize!(usize, u8, u16, u32, u64, u128, BigUint,
                         isize, i8, i16, i32, i64, i128, BigInt,
                         Rational32, Rational64, BigRational);

fn eval_int<C>(p: &Polynomial<C>, x: C) -> C where C: Semiring + Pow<usize, Output=C> + Clone {
    match p {
        Polynomial::Zero() => C::zero(),
        Polynomial::Constant(cc) => cc.0.clone(),
        Polynomial::Dense(dc) => eval_dense(dc, x),
        Polynomial::Sparse(sc) => eval_sparse(sc, x, |x, i|{ x.pow(i) }),
    }
}

//***** For Float *****/
macro_rules! eval_impl_for_f {
    ( $( $t:ident ),* ) => {
        $(
            impl PolynomialEvaluator<$t> for Eval {

                fn eval(p: &Polynomial<$t>, x: $t) -> $t {
                    eval_f(p, x)
                }
            }
        )*
    };
}

eval_impl_for_f!(f32, f64, Complex32, Complex64);

fn eval_f<C>(p: &Polynomial<C>, x: C) -> C where C: Field + Pow<C, Output=C> + FromPrimitive + Clone {
    match p {
        Polynomial::Zero() => C::zero(),
        Polynomial::Constant(cc) => cc.0.clone(),
        Polynomial::Dense(dc) => eval_dense(dc, x),
        Polynomial::Sparse(sc) => eval_sparse(sc, x, |x, i|{ x.pow(C::from_usize(i).unwrap()) }),
    }
}