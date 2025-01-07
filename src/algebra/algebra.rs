use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use num::{traits::{One, Zero}, BigInt, BigRational, BigUint, Complex, Rational32, Rational64};

macro_rules! impl_algebra {
    ( $al:tt ; $( $t:ident ),* ) => {
        $(
            impl $al for $t {}
        )*
    };
}

type Complex32 = Complex<f32>;
type Complex64 = Complex<f64>;

pub trait Semigroup: Mul<Self, Output=Self> where Self: Sized {}
impl_algebra!(Semigroup; usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64, 
    BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait Monoid: Semigroup + PartialEq + One {}
impl_algebra!(Monoid; usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64, 
    BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait Group: Monoid + Div<Self, Output=Self> {}
impl_algebra!(Group; f32, f64, Rational32, Rational64, BigRational, Complex32, Complex64);


pub trait AdditiveSemigroup: Add<Self, Output=Self> where Self: Sized {}
impl_algebra!(AdditiveSemigroup; usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64, 
    BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait AdditiveMonoid: AdditiveSemigroup + Zero {}
impl_algebra!(AdditiveMonoid; usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64, 
    BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait AdditiveGroup: AdditiveMonoid + Neg<Output=Self> + Sub<Self, Output=Self> {}
impl_algebra!(AdditiveGroup; isize, i8, i16, i32, i64, i128, f32, f64, 
    BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);


pub trait Semiring: Monoid + AdditiveMonoid {}
impl_algebra!(Semiring; usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64,
        BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait Ring: Semiring + AdditiveGroup {}
impl_algebra!(Ring; isize, i8, i16, i32, i64, i128, f32, f64,
    BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);


pub trait EuclideanRing: Ring + Eq + Div<Self, Output=Self> + Rem<Self, Output=Self>{}
impl_algebra!(EuclideanRing; isize, i8, i16, i32, i64, i128, BigInt);


pub trait Field: Ring + Group{}
impl_algebra!(Field; f32, f64, Rational32, Rational64, BigRational, Complex32, Complex64);