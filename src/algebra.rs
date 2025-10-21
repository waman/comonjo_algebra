use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use num::{traits::{One, Zero, Euclid}, BigInt, BigRational, BigUint, Complex, Rational32, Rational64};

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
    ( $t:ident ) => {
        impl_ref_traits!($t, RefMul, ref_mul, mul);

        impl Semigroup for $t {}
    };
}

pub trait Monoid: Semigroup + PartialEq + One {}

macro_rules! impl_monoid {
    ( $t:ident ) => {
        impl_semigroup!($t);

        impl Monoid for $t {}
    };
}

pub trait Group: Monoid +
        Div<Self, Output=Self> +
        for<'b> Div<&'b Self, Output=Self> +
        RefDiv<Self> + 
        for<'b> RefDiv<&'b Self> where Self: Sized {}

// macro_rules! impl_group {
//     ( $t:ident ) => {
//         impl_monoid!($t);

//         impl_ref_traits!($t, RefDiv, ref_div, div);

//         impl Group for $t {}
//     };
// }

pub trait AdditiveSemigroup: 
        Add<Self, Output=Self> +
        for<'b> Add<&'b Self, Output=Self> +
        RefAdd<Self> + 
        for<'b> RefAdd<&'b Self> where Self: Sized {}

macro_rules! impl_additive_semigroup {
    ( $t:ident ) => {
        impl_ref_traits!($t, RefAdd, ref_add, add);

        impl AdditiveSemigroup for $t {}
    };
}

pub trait AdditiveMonoid: AdditiveSemigroup + Zero {}

macro_rules! impl_additive_monoid {
    ( $t:ident ) => {
        impl_additive_semigroup!($t);

        impl AdditiveMonoid for $t {}
    };
}

pub trait AdditiveGroup: AdditiveMonoid + 
        Neg<Output=Self> +
        Sub<Self, Output=Self> +
        for<'b> Sub<&'b Self, Output=Self> +
        RefSub<Self> + 
        for<'b> RefSub<&'b Self> where Self: Sized {
    fn ref_neg(&self) -> Self;
}

// macro_rules! impl_additive_group {
//     ( t:ident ) => {
//         impl_additive_monoid!($t);

//         impl_ref_traits!($t, RefSub, ref_sub, sub);

//         impl AdditiveGroup for $t {
//             #[inline]
//             fn ref_neg(&self) -> Self { -self }
//         }
//     };
// }

pub trait Semiring: Monoid + AdditiveMonoid {}

macro_rules! impl_semiring {
    ( $t:ident ) => {
        impl_monoid!($t);

        impl_additive_monoid!($t);

        impl Semiring for $t {}
    };
}

macro_rules! impl_semirings {
    ( $( $t:ident ),* ) => {
        $(
            impl_semiring!($t);
        )*
    };
}

impl_semirings!(usize, u8, u16, u32, u64, u128, BigUint);

pub trait Ring: Semiring + AdditiveGroup {}

macro_rules! impl_ring {
    ( $t:ident ) => {
        impl_semiring!($t);
        
        impl_ref_traits!($t, RefSub, ref_sub, sub);

        impl AdditiveGroup for $t {
            #[inline]
            fn ref_neg(&self) -> Self { -self }
        }

        impl Ring for $t {}
    };
}


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
    fn ref_div_rem(&self, other: &Self) -> (Self, Self);
}

macro_rules! impl_euclidean_ring {
    ( $t:ident ) => {
        impl_ring!($t);

        impl_ref_traits!($t, RefDiv, ref_div, div);
        impl_ref_traits!($t, RefRem, ref_rem, rem);

        impl EuclideanRing for $t {

            fn div_rem(self, other: Self) -> (Self, Self) { 
                self.ref_div_rem(&other)
            }

            fn ref_div_rem(&self, other: &Self) -> (Self, Self) {
                Euclid::div_rem_euclid(self, other) 
            }
        }
    };
}

macro_rules! impl_euclidean_rings {
    ( $( $t:ident ),* ) => {
        $(
            impl_euclidean_ring!($t);
        )*
    };
}

impl_euclidean_rings!(isize, i8, i16, i32, i64, i128, BigInt);


pub trait Field: Ring + Group{}

macro_rules! impl_field {
    ( $t:ident ) => {
        impl_ring!($t);

        impl_ref_traits!($t, RefDiv, ref_div, div);

        impl Group for $t {}

        impl Field for $t {}
    };
}

macro_rules! impl_fields {
    ( $( $t:ident ),* ) => {
        $(
            impl_field!($t);
        )*
    };
}

impl_fields!(f32, f64, Rational32, Rational64, BigRational, Complex32, Complex64);