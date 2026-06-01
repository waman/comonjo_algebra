use num::{traits::Zero, BigInt, BigRational, BigUint, Rational32, Rational64, complex::{Complex32, Complex64}};

use crate::{algebra::*, poly::Polynomial};

pub struct Shift;

pub trait PolynomialShifter<C> where C: Semiring {

    fn shift(p: &mut Polynomial<C>, h: C);
    fn new_shifted(p: &Polynomial<C>, h: C) -> Polynomial<C>;
}

macro_rules! shift_impl {
    ( $mul_div:ident; $( $t:ident ),* ) => {
        $(
            impl PolynomialShifter<$t> for Shift {

                fn shift(p: &mut Polynomial<$t>, h: $t) {

                    if h.is_zero() { return; }
                    match p {
                        Polynomial::Dense(dc) => dc.shift(h, $mul_div),
                        Polynomial::Sparse(sc) => sc.shift(h, $mul_div),
                        _ => (),
                    }
                }

                fn new_shifted(p: &Polynomial<$t>, h: $t) -> Polynomial<$t> {

                    if h.is_zero() { return p.clone(); }
                    match p {
                        Polynomial::Dense(dc) => dc.new_shifted(h, $mul_div),
                        Polynomial::Sparse(sc) => sc.new_shifted(h, $mul_div),
                        _ => p.clone(),
                    }
                }
            }
        )*
    };
}

/// x * y / z
fn mul_div_int<C>(x: C, y: usize, z: C) -> C
        where C: Semiring + num::FromPrimitive + num::Integer + Clone {

    let gcd_xz = x.gcd(&z);
    let x_red = x / gcd_xz.clone();
    let z_red = z / gcd_xz;
    let y_red = C::from_usize(y).unwrap() / z_red;
    x_red * y_red
}

shift_impl!(mul_div_int; usize, u8, u16, u32, u64, u128, BigUint, isize, i8, i16, i32, i64, i128, BigInt);

/// x * y / z
fn mul_div_field<C>(x: C, y: usize, z: C) -> C where C: Field + num::FromPrimitive {
    x * C::from_usize(y).unwrap() / z
}

shift_impl!(mul_div_field; Rational32, Rational64, BigRational, f32, f64, Complex32, Complex64);


