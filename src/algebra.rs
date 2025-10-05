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

//********** Traits for Arithmetic Operation on Reference  **********/
pub trait RefAdd<Rhs> {
    fn ref_add(&self, rhs: Rhs) -> Self;
} 

pub trait RefSub<Rhs> {
    fn ref_sub(&self, rhs: Rhs) -> Self;
} 

pub trait RefMul<Rhs> {
    fn ref_mul(&self, rhs: Rhs) -> Self;
} 

pub trait RefDiv<Rhs> {
    fn ref_div(&self, rhs: Rhs) -> Self;
} 

pub trait RefRem<Rhs> {
    fn ref_rem(&self, rhs: Rhs) -> Self;
} 

pub trait Semigroup:
        Mul<Self, Output=Self> +
        for<'b> Mul<&'b Self, Output=Self> +
        RefMul<Self> + 
        for<'b> RefMul<&'b Self> where Self: Sized {}

macro_rules! impl_ref_traits {
    ( $t:ident, $ref_type:tt, $ref_op:ident, $op:ident ) => {
        impl $ref_type<$t> for $t {
            #[inline]
            fn $ref_op(&self, other: Self) -> Self { self.$op(other) }
        }

        impl<'b> $ref_type<&'b $t> for $t {
            #[inline]
            fn $ref_op(&self, other: &Self) -> Self { self.$op(other) }
        }
    };
}

macro_rules! impl_semigroup {
    ( $( $t:ident ),* ) => {
        $(
            impl_ref_traits!($t, RefMul, ref_mul, mul);

            impl Semigroup for $t {}
        )*
    };
}

impl_semigroup!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64, 
    BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait Monoid: Semigroup + PartialEq + One {}
impl_algebra!(Monoid; usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64, 
    BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait Group: Monoid +
        Div<Self, Output=Self> +
        for<'b> Div<&'b Self, Output=Self> +
        RefDiv<Self> + 
        for<'b> RefDiv<&'b Self> where Self: Sized {}

macro_rules! impl_group {
    ( $( $t:ident ),* ) => {
        $(
            impl_ref_traits!($t, RefDiv, ref_div, div);

            impl Group for $t {}
        )*
    };
}

impl_group!(f32, f64, Rational32, Rational64, BigRational, Complex32, Complex64);


pub trait AdditiveSemigroup: 
        Add<Self, Output=Self> +
        for<'b> Add<&'b Self, Output=Self> +
        RefAdd<Self> + 
        for<'b> RefAdd<&'b Self> where Self: Sized {}

macro_rules! impl_additive_semigroup {
    ( $( $t:ident ),* ) => {
        $(
            impl_ref_traits!($t, RefAdd, ref_add, add);

            impl AdditiveSemigroup for $t {}
        )*
    };
}

impl_additive_semigroup!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64, 
    BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait AdditiveMonoid: AdditiveSemigroup + Zero {}
impl_algebra!(AdditiveMonoid; usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64, 
    BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait AdditiveGroup: AdditiveMonoid + 
        Neg<Output=Self> +
        Sub<Self, Output=Self> +
        for<'b> Sub<&'b Self, Output=Self> +
        RefSub<Self> + 
        for<'b> RefSub<&'b Self> where Self: Sized {
    fn ref_neg(&self) -> Self;

}

macro_rules! impl_additive_group {
    ( $( $t:ident ),* ) => {
        $(
            impl_ref_traits!($t, RefSub, ref_sub, sub);

            impl AdditiveGroup for $t {
                #[inline]
                fn ref_neg(&self) -> Self { -self }
            }
        )*
    };
}

impl_additive_group!(isize, i8, i16, i32, i64, i128, f32, f64, 
    BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);


pub trait Semiring: Monoid + AdditiveMonoid {}
impl_algebra!(Semiring; usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64,
        BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

pub trait Ring: Semiring + AdditiveGroup {}
impl_algebra!(Ring; isize, i8, i16, i32, i64, i128, f32, f64,
    BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);


pub trait EuclideanRing: Ring + 
        Div<Self, Output=Self> +
        for<'b> Div<&'b Self, Output=Self> +
        RefDiv<Self> + 
        for<'b> RefDiv<&'b Self> +

        Rem<Self, Output=Self> +
        for<'b> Rem<&'b Self, Output=Self> +
        RefRem<Self> + 
        for<'b> RefRem<&'b Self> where Self: Sized {

    fn div_rem(self, other: Self) -> (Self, Self);
    fn div_rem_ref(&self, other: &Self) -> (Self, Self);
}

macro_rules! impl_euclidean_ring {
    ( $( $t:ident ),* ) => {
        $(
            impl_ref_traits!($t, RefDiv, ref_div, div);
            impl_ref_traits!($t, RefRem, ref_rem, rem);

            impl EuclideanRing for $t {
                #[inline]
                fn div_rem(self, other: Self) -> (Self, Self) { 
                    (&self).div_rem_ref(&other)
                }

                #[inline]
                fn div_rem_ref(&self, other: &Self) -> (Self, Self) {
                    (self / other, self % other) 
                }
            }
        )*
    };
}

impl_euclidean_ring!(isize, i8, i16, i32, i64, i128, BigInt);


pub trait Field: Ring + Group{}
impl_algebra!(Field; f32, f64, Rational32, Rational64, BigRational, Complex32, Complex64);