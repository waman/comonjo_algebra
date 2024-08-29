use std::{i64, ops::{Add, Div, Mul, Neg, Rem, Shl, Shr, Sub}};
use num::{bigint::{ParseBigIntError, Sign}, integer::Roots, pow::Pow, traits::{ConstOne, ConstZero}, BigInt, Integer, Num, One, Zero, ToPrimitive};

use once_cell::sync::Lazy;

/// Refer to spire's <a href="https://github.com/typelevel/spire/blob/main/core/src/main/scala/spire/math/SafeLong.scala">SafeLong</a>.
#[derive(Clone, PartialEq, Eq)]
pub enum SafeI64 {
    I64(i64),
    BI(Box<BigInt>)
}

impl SafeI64 {

    fn new_raw(x: BigInt) -> SafeI64 {
        debug_assert_eq!(x.to_i64(), None, "An i64-range value is passed to 'new_raw': {:?}", x);
        SafeI64::BI(Box::new(x))
    }

    pub fn from(x: i64) -> SafeI64 {
        SafeI64::I64(x)  // raw creation
    }

    pub fn from_integer(x: BigInt) -> SafeI64 {
        match x.to_i64() {
            Some(i) => SafeI64::I64(i),  // raw creation
            _ => SafeI64::new_raw(x),  // raw creation
        }
    }

    pub fn is_primitive(&self) -> bool {
        match self {
            SafeI64::I64(_) => true,
            _ => false,
        }
    }

    pub fn sign(&self) -> Sign {
        match self {
            SafeI64::I64(x) => if *x > 0 {
                Sign::Plus
            } else if *x < 0 {
                Sign::Minus
            } else {
                Sign::NoSign
            },
            SafeI64::BI(x) => x.sign(),
        }
    }

    pub fn abs(&self) -> SafeI64 {
        match self {
            SafeI64::I64(x) => SafeI64::from(x.abs()),
            SafeI64::BI(x) => match x.sign() {
                Sign::Plus => self.clone(),
                Sign::NoSign => SafeI64::zero(),
                Sign::Minus => self.clone() * (-1),  // TODO
            },
        }
    }

    pub fn sqrt(&self) -> SafeI64 {
        match self {
            SafeI64::I64(x) => SafeI64::from(x.sqrt()),
            SafeI64::BI(x) => SafeI64::from_integer(x.sqrt()),
        }
    }

    pub fn to_bigint(&self) -> BigInt {
        match self {
            SafeI64::I64(x) => BigInt::from(*x),
            SafeI64::BI(x) => (**x).clone(),
        }
    }

    pub fn to_f64(&self) -> Option<f64> {
        match self {
            SafeI64::I64(x) => x.to_f64(),
            SafeI64::BI(x) => x.to_f64(),
        }
    }

    pub fn to_str_radix(&self, radix: u32) -> String {
        match self {
            SafeI64::I64(x) => BigInt::from(*x).to_str_radix(radix),
            SafeI64::BI(x) => x.to_str_radix(radix),
        }
    }
}

macro_rules! num_op_impl {
    ($typ:ident, $op:ident, $ch_op:ident, $gen_fn:ident) => {

        impl $typ for SafeI64 {

            type Output = SafeI64;
        
            fn $op(self, rhs: Self) -> Self::Output {
                match rhs {
                    SafeI64::I64(y) => self.$op(y),
                    SafeI64::BI(y) => self.$op(*y),
                }
            }
        }
        
        impl $typ<i64> for SafeI64 {
        
            type Output = SafeI64;
        
            fn $op(self, rhs: i64) -> Self::Output {
                match self {
                    SafeI64::I64(x) => {
                        match x.$ch_op(rhs) {
                            Some(z) => SafeI64::from(z),
                            _ => SafeI64::$gen_fn(BigInt::from(x).$op(rhs)),
                        }
                    },
                    SafeI64::BI(x) => SafeI64::from_integer((*x).$op(rhs)),
                }
        
            }
        }
        
        impl $typ<BigInt> for SafeI64 {
        
            type Output = SafeI64;
        
            fn $op(self, rhs: BigInt) -> Self::Output {
                match self {
                    SafeI64::I64(x) => SafeI64::from_integer(x.$op(rhs)),
                    SafeI64::BI(x) => SafeI64::from_integer((*x).$op(rhs)),
                }
            }
        }


        impl $typ for &SafeI64 {

            type Output = SafeI64;
        
            fn $op(self, rhs: Self) -> Self::Output {
                match rhs {
                    SafeI64::I64(y) => self.$op(*y),
                    SafeI64::BI(y) => self.$op(&**y),
                }
            }
        }
        
        impl $typ<i64> for &SafeI64 {
        
            type Output = SafeI64;
        
            fn $op(self, rhs: i64) -> Self::Output {
                match self {
                    SafeI64::I64(x) => {
                        match x.$ch_op(rhs) {
                            Some(z) => SafeI64::from(z),
                            _ => SafeI64::$gen_fn(BigInt::from(*x).$op(rhs)),
                        }
                    },
                    SafeI64::BI(x) => SafeI64::from_integer((&**x).$op(rhs)),
                }
        
            }
        }
        
        impl $typ<&BigInt> for &SafeI64 {
        
            type Output = SafeI64;
        
            fn $op(self, rhs: &BigInt) -> Self::Output {
                match self {
                    SafeI64::I64(x) => SafeI64::from_integer(x.$op(rhs)),
                    SafeI64::BI(x) => SafeI64::from_integer((&**x).$op(rhs)),
                }
            }
        }
    };
}

num_op_impl!(Add, add, checked_add, new_raw);
num_op_impl!(Sub, sub, checked_sub, new_raw);
num_op_impl!(Mul, mul, checked_mul, new_raw);
num_op_impl!(Div, div, checked_div, new_raw);
num_op_impl!(Rem, rem, checked_rem, from_integer);  // i64::MIN.checked_rem(-1) returns None while the result of BigInt is 0.

impl Pow<u32> for SafeI64{

    type Output = SafeI64;
    
    fn pow(self, p: u32) -> Self::Output {
        if p == 0 { return SafeI64::ONE; }
        match self {
            SafeI64::I64(x) => {
                match x.checked_pow(p) {
                    Some(z) => SafeI64::from(z),
                    _ => SafeI64::new_raw(BigInt::from(x).pow(p)),
                }
            },
            SafeI64::BI(x) => SafeI64::new_raw((*x).pow(p)),
        }
    }
}

impl Pow<u32> for &SafeI64{

    type Output = SafeI64;
    
    fn pow(self, p: u32) -> Self::Output {
        if p == 0 { return SafeI64::ONE; }
        match self {
            SafeI64::I64(x) => {
                match x.checked_pow(p) {
                    Some(z) => SafeI64::from(z),
                    _ => SafeI64::new_raw(BigInt::from(*x).pow(p)),
                }
            },
            SafeI64::BI(x) => SafeI64::new_raw(x.clone().pow(p)),
        }
    }
}

fn positive_shl(x: i64, p: u32) -> SafeI64 {
    debug_assert!(x > 0);
    debug_assert!(p != 0);
    if x.leading_zeros() > p {
        SafeI64::from(x << p)
    } else {
        SafeI64::new_raw(BigInt::from(x) << p)
    }
}

impl Shl<u32> for SafeI64{

    type Output = SafeI64;
    
    fn shl(self, p: u32) -> Self::Output {
        if p == 0 { return self; }
        match self {
            SafeI64::I64(x) => {
                if x > 0 {
                    positive_shl(x, p)
                } else if x == i64::MIN {
                    SafeI64::new_raw(BigInt::from(x) << p)  // note p != 0
                } else if x < 0 {
                    -positive_shl(-x, p)
                } else {
                    SafeI64::ZERO
                }
            },
            SafeI64::BI(x) => SafeI64::new_raw((*x) << p),
        }
    }
}

impl Shl<u32> for &SafeI64{

    type Output = SafeI64;
    
    fn shl(self, p: u32) -> Self::Output {
        if p == 0 { return self.clone(); }
        match self {
            SafeI64::I64(x) => {
                let y = *x;
                if y > 0 {
                    positive_shl(y, p)
                } else if y == i64::MIN {
                    SafeI64::new_raw(BigInt::from(y) << p)  // note p != 0
                } else if y < 0 {
                    -positive_shl(-y, p)
                } else {
                    SafeI64::ZERO
                }
            },
            SafeI64::BI(x) => SafeI64::new_raw(*x.clone() << p),
        }
    }
}

impl Shr<u32> for SafeI64{

    type Output = SafeI64;
    
    fn shr(self, rhs: u32) -> Self::Output {
        match self {
            SafeI64::I64(x) => {
                match x.checked_shr(rhs) {
                    Some(z) => SafeI64::from(z),
                    _ => SafeI64::from_integer(BigInt::from(x).shr(rhs)),
                }
            },
            SafeI64::BI(x) => SafeI64::from_integer((*x).shr(rhs)),
        }
    }
}

impl Shr<u32> for &SafeI64{

    type Output = SafeI64;
    
    fn shr(self, rhs: u32) -> Self::Output {
        match self {
            SafeI64::I64(x) => {
                match x.checked_shr(rhs) {
                    Some(z) => SafeI64::from(z),
                    _ => SafeI64::from_integer(BigInt::from(*x).shr(rhs)),
                }
            },
            SafeI64::BI(x) => SafeI64::from_integer(x.clone().shr(rhs)),
        }
    }
}

impl Zero for SafeI64{
    
    fn is_zero(&self) -> bool {
        match self {
            SafeI64::I64(x) => x.is_zero(),
            _ => false,
        }
    }
    
    fn zero() -> Self { SafeI64::ZERO }
}

impl ConstZero for SafeI64 {

    const ZERO: Self = SafeI64::I64(0_i64);  // raw creation
}

impl One for SafeI64{
    fn one() -> Self { Self::ONE }
}

impl ConstOne for SafeI64 {

    const ONE: Self = SafeI64::I64(1_i64);  // raw creation
}

/// BigInt value of i64::MAX + 1
pub(crate) static I64_MAX_P1: Lazy<BigInt> = Lazy::new(|| BigInt::from(i64::MAX) + 1 );

impl Neg for SafeI64{

    type Output = SafeI64;

    fn neg(self) -> Self::Output {
        match self {
            SafeI64::I64(x) => 
                if x == i64::MIN { 
                    SafeI64::new_raw(I64_MAX_P1.clone())  // raw creation
                }else{
                    SafeI64::I64(-x)  // raw creation
                },
            SafeI64::BI(x) => 
                if *x == *I64_MAX_P1 {
                    SafeI64::I64(i64::MIN)  // raw creation
                }else{
                    SafeI64::new_raw(-*x)  // raw creation
                },
        }
    }
}

impl Neg for &SafeI64{

    type Output = SafeI64;

    fn neg(self) -> Self::Output {
        match self {
            SafeI64::I64(x) => 
                if *x == i64::MIN { 
                    SafeI64::new_raw(I64_MAX_P1.clone())  // raw creation
                }else{
                    SafeI64::I64(-x)  // raw creation
                },
            SafeI64::BI(x) => 
                if **x == *I64_MAX_P1 {
                    SafeI64::I64(i64::MIN)  // raw creation
                }else{
                    SafeI64::new_raw(-*x.clone()) // raw creation
                },
        }
    }
}

impl Num for SafeI64{

    type FromStrRadixErr = ParseBigIntError;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        match i64::from_str_radix(str, radix){
            Ok(x) => Ok(SafeI64::I64(x)),  // raw creation
            Err(_) => {
                match BigInt::from_str_radix(str, radix) {
                    Ok(y) => Ok(SafeI64::new_raw(y)),  // raw creation
                    Err(err) => Err(err),
                }
            },
        }
    }
}

impl PartialOrd for SafeI64{

    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self {
            SafeI64::I64(x) => match other {
                SafeI64::I64(y) => x.partial_cmp(y),
                SafeI64::BI(y) => BigInt::from(*x).partial_cmp(y),
            },
            SafeI64::BI(x) => match other {
                SafeI64::I64(y) => (**x).partial_cmp(&BigInt::from(*y)),
                SafeI64::BI(y) => x.partial_cmp(y),
            },
        }
    }
}

impl Ord for SafeI64{

    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self {
            SafeI64::I64(x) => match other {
                SafeI64::I64(y) => x.cmp(y),
                SafeI64::BI(y) => BigInt::from(*x).cmp(y),
            },
            SafeI64::BI(x) => match other {
                SafeI64::I64(y) => (**x).cmp(&BigInt::from(*y)),
                SafeI64::BI(y) => x.cmp(y),
            },
        }
    }
}

macro_rules! unary_fn_impl {
    ($method:ident) => {
        fn $method(&self) -> bool {
            match self {
                SafeI64::I64(x) => x.$method(),
                SafeI64::BI(x) => x.$method(),
            }
        }
    };
}

macro_rules! binary_fn_impl {
    ($method:ident) => {
        fn $method(&self, other: &Self) -> Self {
            match self {
                SafeI64::I64(x) => match other {
                    SafeI64::I64(y) => {
                        if x == &i64::MIN || y == &i64::MIN {
                            SafeI64::from_integer(BigInt::from(*x).$method(&BigInt::from(*y)))
                        } else {
                            SafeI64::from(x.$method(y))
                        }
                    },
                    SafeI64::BI(y) => SafeI64::from_integer(BigInt::from(*x).$method(y)),
                },
                SafeI64::BI(x) => match other {
                    SafeI64::I64(y) => SafeI64::from_integer(x.$method(&BigInt::from(*y))),
                    SafeI64::BI(y) => SafeI64::from_integer(x.$method(y)),
                },
            }
        }
    };
}

impl Integer for SafeI64{

    unary_fn_impl!(is_even);
    unary_fn_impl!(is_odd);

    binary_fn_impl!(div_floor);
    binary_fn_impl!(mod_floor);
    binary_fn_impl!(gcd);

    fn lcm(&self, other: &Self) -> Self {
        let x = &self.to_bigint();
        let y = &other.to_bigint();
        SafeI64::from_integer(x.lcm(y))
    }

    fn is_multiple_of(&self, other: &Self) -> bool {
        match self {
            SafeI64::I64(x) => match other {
                SafeI64::I64(y) => {
                    if x == &i64::MIN || y == &i64::MIN {
                        BigInt::from(*x).is_multiple_of(&BigInt::from(*y))
                    } else {
                        x.is_multiple_of(y)
                    }
                },
                SafeI64::BI(y) => BigInt::from(*x).is_multiple_of(y),
            },
            SafeI64::BI(x) => match other {
                SafeI64::I64(y) => x.is_multiple_of(&BigInt::from(*y)),
                SafeI64::BI(y) => x.is_multiple_of(y),
            },
        }
    }

    fn div_rem(&self, other: &Self) -> (Self, Self) {
        fn to_safe_i64_tuple(r: (BigInt, BigInt)) -> (SafeI64, SafeI64) {
            (SafeI64::from_integer(r.0), SafeI64::from_integer(r.1))
        }

        match self {
            SafeI64::I64(x) => match other {
                SafeI64::I64(y) => {
                    if x == &i64::MIN || y == &i64::MIN {
                        let r = BigInt::from(*x).div_rem(&BigInt::from(*y));
                        (SafeI64::from_integer(r.0), SafeI64::from_integer(r.1))
                    } else {
                        let r = x.div_rem(y);
                        (SafeI64::from(r.0), SafeI64::from(r.1))
                    }
                },
                SafeI64::BI(y) => to_safe_i64_tuple(BigInt::from(*x).div_rem(y))
            },
            SafeI64::BI(x) => match other {
                SafeI64::I64(y) => to_safe_i64_tuple(x.div_rem(&BigInt::from(*y))),
                SafeI64::BI(y) => to_safe_i64_tuple(x.div_rem(y)),
            }
        }
    }
}

impl std::fmt::Display for SafeI64 {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SafeI64::I64(x) => x.fmt(f),
            SafeI64::BI(x) => x.fmt(f),
        }
    }
}

impl std::fmt::Debug for SafeI64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::I64(x) => f.debug_tuple("SafeI64::I64").field(x).finish(),
            Self::BI(x) => f.debug_tuple("SafeI64::BI").field(x).finish(),
        }
    }
}