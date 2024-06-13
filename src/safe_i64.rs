use std::{fmt::Display, ops::{Add, Div, Mul, Rem, Shl, Shr, Sub}};

use num::{bigint::{ParseBigIntError, Sign}, integer::Roots, pow::Pow, BigInt, Integer, Num, One, Zero};
use once_cell::sync::Lazy;

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum SafeI64 {
    I64(i64),
    BI(Box<BigInt>)
}

#[test]
fn test_size(){
    assert_eq!(std::mem::size_of::<i64>(), 8);  // <- use
    assert_eq!(std::mem::size_of::<i128>(), 16);
    assert_eq!(std::mem::size_of::<BigInt>(), 32);
    assert_eq!(std::mem::size_of::<&BigInt>(), 8);
    assert_eq!(std::mem::size_of::<Box<BigInt>>(), 8);  // <- use
}

static I64_MAX: Lazy<BigInt> = Lazy::new(|| BigInt::from(i64::MAX));
static I64_MIN: Lazy<BigInt> = Lazy::new(|| BigInt::from(i64::MIN));

impl SafeI64 {

    pub fn from(x: i64) -> SafeI64 {
        SafeI64::I64(x)
    }

    pub fn from_integer(x: BigInt) -> SafeI64 {
        match try_into_i64(&x) {
            Some(i) => SafeI64::I64(i),
            _ => SafeI64::BI(Box::new(x)),
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
                Sign::Minus => self.clone() * -1,
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

    pub fn to_str_radix(&self, radix: u32) -> String {
        self.to_bigint().to_str_radix(radix)
    }
}

fn try_into_i64(x: &BigInt) -> Option<i64>{
    if *x <= *I64_MAX && *I64_MIN <= *x {
        let (sign, d) = x.to_u64_digits();
        match sign {
            Sign::Plus => Some(d[0] as i64),
            Sign::NoSign => Some(0),
            Sign::Minus => Some((u64::MAX - d[0] + 1) as i64),
        }
    }else{
        None
    }
}

#[test]
fn test_try_into_i64(){
    assert_eq!(try_into_i64(&BigInt::from(0)), Some(0_i64));
    assert_eq!(try_into_i64(&BigInt::from(1)), Some(1_i64));
    assert_eq!(try_into_i64(&BigInt::from(234)), Some(234_i64));
    assert_eq!(try_into_i64(&BigInt::from(-1)), Some(-1_i64));
    assert_eq!(try_into_i64(&BigInt::from(-567)), Some(-567_i64));

    assert_eq!(try_into_i64(&(BigInt::from(i64::MAX) + BigInt::from(12))), None);
    assert_eq!(try_into_i64(&(BigInt::from(i64::MIN) - BigInt::from(34))), None);
}

macro_rules! num_op_impl {
    ($typ:ident, $op:ident, $ch_op:ident) => {

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
                            _ => SafeI64::BI(Box::new(BigInt::from(x).$op(rhs))),
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
                    SafeI64::I64(y) => self.$op(y),
                    SafeI64::BI(y) => self.$op(&**y),
                }
            }
        }
        
        impl $typ<&i64> for &SafeI64 {
        
            type Output = SafeI64;
        
            fn $op(self, rhs: &i64) -> Self::Output {
                match self {
                    SafeI64::I64(x) => {
                        match x.$ch_op(*rhs) {
                            Some(z) => SafeI64::from(z),
                            _ => SafeI64::BI(Box::new(BigInt::from(*x).$op(rhs))),
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

num_op_impl!(Add, add, checked_add);
num_op_impl!(Sub, sub, checked_sub);
num_op_impl!(Mul, mul, checked_mul);
num_op_impl!(Div, div, checked_div);
num_op_impl!(Rem, rem, checked_rem);

impl Pow<u32> for SafeI64 {

    type Output = Self;

    fn pow(self, p: u32) -> Self::Output {
        match self {
            SafeI64::I64(x) => match x.checked_pow(p) {
                Some(z) => SafeI64::from(z),
                _ => SafeI64::from_integer(BigInt::from(x).pow(p)),
            },
            SafeI64::BI(x) => SafeI64::from_integer(x.pow(p)),
        }
    }
}

macro_rules! bit_op_impl {
    ($typ:ident, $op:ident, $ch_op:ident) => {
        impl $typ<u32> for SafeI64{

            type Output = Self;
            
            fn $op(self, rhs: u32) -> Self::Output {
                match self {
                    SafeI64::I64(x) => {
                        match x.$ch_op(rhs) {
                            Some(z) => SafeI64::from(z),
                            _ => SafeI64::from_integer(BigInt::from(x).$op(rhs)),
                        }
                    },
                    SafeI64::BI(x) => SafeI64::from_integer((*x).$op(rhs)),
                }
            }
        }
    };
}

bit_op_impl!(Shl, shl, checked_shl);
bit_op_impl!(Shr, shr, checked_shr);


#[test]
fn test_add(){
    for _ in 1..1000 {
        let x = rand::random::<i64>();
        let y = rand::random::<i64>();

        let sx = SafeI64::from(x);
        let sy = SafeI64::from(y);
        let sz = sx + sy;

        match x.checked_add(y) {
            Some(z) => {
                assert!(sz.is_primitive());
                assert_eq!(sz, SafeI64::from(z));
            },
            _ => {
                assert!(!sz.is_primitive());
                let z = BigInt::from(x) + BigInt::from(y);
                assert_eq!(sz, SafeI64::from_integer(z));
            },
        }

        // BigInt random number
        let bx = BigInt::from(x) * 2 - 1;
        let by = BigInt::from(y) * 2 - 1;

        assert_eq!(sz * 2 - 2, SafeI64::from_integer(bx + by));
    }
}

impl Zero for SafeI64{
    
    fn is_zero(&self) -> bool {
        match self {
            SafeI64::I64(x) => x.is_zero(),
            _ => false,
        }
    }
    
    fn zero() -> Self { SafeI64::I64(0) }
}

impl One for SafeI64{
    fn one() -> Self { SafeI64::I64(1) }
}

impl Num for SafeI64{

    type FromStrRadixErr = ParseBigIntError;

    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        match i64::from_str_radix(str, radix){
            Ok(x) => Ok(SafeI64::I64(x)),
            Err(_) => {
                match BigInt::from_str_radix(str, radix) {
                    Ok(y) => Ok(SafeI64::BI(Box::new(y))),
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
                    SafeI64::I64(y) => SafeI64::from(x.$method(y)),
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
    binary_fn_impl!(lcm);

    fn is_multiple_of(&self, other: &Self) -> bool {
        match self {
            SafeI64::I64(x) => match other {
                SafeI64::I64(y) => x.is_multiple_of(y),
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
                    let r = x.div_rem(y);
                    (SafeI64::from(r.0), SafeI64::from(r.1))
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

impl Display for SafeI64 {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SafeI64::I64(x) => x.fmt(f),
            SafeI64::BI(x) => x.fmt(f),
        }
    }
}

#[test]
fn test_display(){
    let x = SafeI64::from(123);
    assert_eq!(format!("{}", x), "123");

    let y = SafeI64::from(i64::MAX) + SafeI64::from(123);
    println!("{}", i64::MAX);
    assert_eq!(format!("{}", y), "9223372036854775930");  //9223372036854775807 + 123
}

#[test]
fn test_debug(){
    let x = SafeI64::from(123);
    assert_eq!(format!("{:?}", x), "I64(123)");

    let y = SafeI64::from(i64::MAX) + SafeI64::from(123);
    assert_eq!(format!("{:?}", y), "BI(9223372036854775930)");  //9223372036854775807 + 123
}