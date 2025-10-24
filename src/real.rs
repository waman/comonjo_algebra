use std::{fmt::{Debug, Display, Formatter}, ops::{Add, Mul}};
use regex::Regex;
use num::{bigint::Sign, pow::Pow, rational::Ratio, traits::{ConstOne, ConstZero}, BigInt, One, Zero};

use crate::safe_i64::SafeI64;

pub(crate) type SafeRational = Ratio<SafeI64>;

/// Refer to spire's <a href="https://github.com/typelevel/spire/blob/main/core/src/main/scala/spire/math/Real.scala">Real</a>
pub enum Real<'a>{
    Zero(),
    One(),
    Exact(SafeRational),
    Inexact{ memo: Option<(u32, SafeI64)>, f: Box<dyn FnMut(u32) -> SafeI64 + 'a> }
    // PositiveInfinity(),
    // NegativeInfinity(),
}

impl<'a> Real<'a> {

    fn from<'b, N>(value: N, f: &dyn Fn(N) -> SafeRational) -> Real<'b> 
            where N: Eq + Zero + One{
        if value.is_zero() {
            Real::Zero()
        } else if value.is_one() {
            Real::One()
        } else{
            Real::Exact(f(value))
        }
    }

    pub fn from_i64<'b>(value: i64) -> Real<'b> {
        Real::from(value, &|value| SafeRational::from_integer(SafeI64::from(value)))
    }

    pub fn from_bigint<'b>(value: BigInt) -> Real<'b> {
        Real::from(value, &|value| SafeRational::from_integer(SafeI64::from_integer(value)))
    }

    pub fn from_safe_i64<'b>(value: SafeI64) -> Real<'b> {
        Real::from(value, &|value| SafeRational::from_integer(value))
    }

    pub fn from_rational<'b>(value: SafeRational) -> Real<'b> {
        Real::from(value, &|value| value)
    }

    fn from_values<'b, N>(num: N, deno: N, f: &dyn Fn(N, N) -> SafeRational) -> Real<'b> 
            where N: Eq + Zero{
        if num.is_zero() {
            Real::Zero()
        } else if num == deno {
            Real::One()
        } else{
            Real::Exact(f(num, deno))
        }
    }

    pub fn from_i64s<'b>(num: i64, deno: i64) -> Real<'b> {
        Real::from_values(num, deno, &|num, deno| SafeRational::new(SafeI64::I64(num), SafeI64::I64(deno)))
    }

    pub fn from_bigints<'b>(num: BigInt, deno: BigInt) -> Real<'b> {
        Real::from_values(num, deno, 
            &|num, deno| SafeRational::new(SafeI64::from_integer(num), SafeI64::from_integer(deno)))
    }

    pub fn from_safe_i64s<'b>(num: SafeI64, deno: SafeI64) -> Real<'b> {
        Real::from_values(num, deno,  &|num, deno| SafeRational::new(num, deno))
    }

    pub fn new_inexact<'b>(f: Box<dyn FnMut(u32) -> SafeI64 + 'b>) -> Real<'b> {
        Real::Inexact{ memo: None, f }
    }

    pub fn eval(&mut self, p: u32) -> SafeI64 {
        match self {
            Real::Zero() => SafeI64::zero(),
            Real::One() => SafeI64::one() << p,
            Real::Exact(ref r) => round_up(r * (SafeI64::one() << p)),

            Real::Inexact{ memo, f } => {
                match memo {
                    Some((bits, value)) if *bits >= p => 
                        round_up(SafeRational::new(value.clone(), SafeI64::one() << (*bits - p))),
                    _ => {
                        let result = f(p);
                        *memo = Some((p, result.clone()));
                        result
                    },
                }
            },
        }
    }

    pub fn eval_digits(&mut self, d: u32) -> SafeI64 {
        self.eval_digits_radix(d, 10)
    }

    pub fn eval_digits_radix(&mut self, d: u32, radix: u32) -> SafeI64 {
        self.eval(digits_to_bits(d, radix))
    }

    pub fn sqrt<'b>(&'b mut self) -> Real<'b> {
        Real::new_inexact(Box::new(|p| self.eval(p * 2).sqrt()))
    }

    pub fn try_into_rational(&self) -> Option<SafeRational> {
        match self {
            Real::Zero() => Some(SafeRational::zero()),
            Real::One() => Some(SafeRational::one()),
            Real::Exact(r) => Some(r.clone()),
            _ => None,
        }
    }
}

pub(crate) fn round_up(r: SafeRational) -> SafeI64 {
    let n = r.numer();
    let d = r.denom();
    if *n >= SafeI64::zero() {
        let m = n % d;
        if m * 2 >= *d { n / d + SafeI64::one() } else { n / d } 
    } else {
        let m = (n % d) * (-1);
        if m * 2 >= *d { n / d - SafeI64::one() } else { n / d }
    }
}

// impl<'a> From<SafeRational> for Real<'a> {
//     fn from(value: SafeRational) -> Self { Real::Exact(value) }
// }

// impl<'a> From<i64> for Real<'a> {
//     fn from(value: i64) -> Self { Real::Exact(SafeRational::from(SafeI64::from(value))) }
// }

// impl<'a> From<SafeI64> for Real<'a> {
//     fn from(value: SafeI64) -> Self { Real::Exact(SafeRational::from(value)) }
// }

// impl<'a> From<(i64, i64)> for Real<'a> {
//     fn from(value: (i64, i64)) -> Self { 
//         Real::Exact(SafeRational::new(SafeI64::from(value.0), SafeI64::from(value.1))) 
//     }
// }

impl<'a> Zero for Real<'a> {

    fn zero() -> Self { Real::Zero() }

    fn is_zero(&self) -> bool {
        match self {
            Real::Zero() => true,
            Real::Exact(r) => r.is_zero(),
            _ => false
        }
    }
}

impl<'a> ConstZero for Real<'a>{

    const ZERO: Self = Real::Zero();
}

impl<'a> One for Real<'a> {

    fn one() -> Self { Real::One() }
}

impl<'a> ConstOne for Real<'a> {
    const ONE: Self = Real::One();
}

impl<'a> Display for Real<'a> {

    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Real::Zero() => f.write_str("0"),
            Real::One() => f.write_str("1"),
            Real::Exact(n) => f.write_str(&format!("{}", n)),

            Real::Inexact { memo, f:_ } => {
                match memo {
                    Some((bit, value)) => {
                        let s = inexact_to_display_string(*bit, value.clone(), 10);
                        f.write_str(&s)
                    },
                    _ => f.write_str("Real::Inexact(Unevaluated)"),
                }
            },
        }
    }
}

/// Return min u32 value <i>x</i> that satisfies <i>radix</i>^<i>d</i> <= 2^<i>x</i>.
pub(crate) fn digits_to_bits(d: u32, radix: u32) -> u32 {
    ((d as f64) * (radix as f64).ln() / 2_f64.ln()).ceil() as u32
}

pub(crate) fn bits_to_digits(b: u32, radix: u32) -> u32 {
    ((b as f64) * 2_f64.ln() / (radix as f64).ln()).floor() as u32
}

fn inexact_to_display_string(bit: u32, value: SafeI64, radix: u32) -> String {
    let d = bits_to_digits(bit, radix);

    let rad = SafeI64::from(radix as i64);
    let r = SafeRational::new(value * rad.pow(d), SafeI64::one() << bit);
    let m = round_up(r);

    let (sign, abs) = match m.sign() {
        Sign::Minus  => ("-", m.abs().to_str_radix(radix)),
        Sign::NoSign => ("", "0".to_string()),
        Sign::Plus   => ("", m.to_str_radix(radix)),
    };

    let i = abs.len() - d as usize;
    let s = if i > 0 {
        format!("{}{}.{}", sign, &abs[0..i], &abs[i..])
    } else{
        let zeros = "0".repeat(d as usize - abs.len());
        format!("{}0.{}{}", sign, zeros, abs)
    };

    let s1 = Regex::new(r"0+$").unwrap().replace_all(&s, "");
    let s2 = Regex::new(r"\.$").unwrap().replace_all(&s1, "");
    s2.to_string()
}

impl<'a> Debug for Real<'a> {

    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Real::Zero() => f.write_str("Real::0"),
            Real::One() => f.write_str("Real::1"),
            Real::Exact(n) => 
                f.write_str(&format!("Real::Exact({})", n)),

            Real::Inexact { memo, f:_ } => {
                match memo {
                    Some((bit, value)) => {
                        let s = inexact_to_display_string(*bit, value.clone(), 10);
                        f.write_str(&format!("Real::Inexact({}; bit={}, value={})", s, bit, value))
                    },
                    _ => f.write_str(&format!("Real::Inexact(Unevaluated)")),
                }
            },
        }
    }
}

impl<'a> Add<Real<'a>> for Real<'a>{

    type Output = Real<'a>;

    fn add(self, rhs: Real<'a>) -> Self::Output {
        fn add_rationals<'b>(x: SafeRational, y: SafeRational) -> Real<'b> {
            Real::from_rational(x + y)
        }

        fn add_inexact<'b>(mut x: Real<'b>, mut y: Real<'b>) -> Real<'b> {
            Real::new_inexact(Box::new(move |p|
                round_up(SafeRational::new(x.eval(p+2) + y.eval(p+2), SafeI64::from(4_i64)))))
        }

        match rhs {
            Real::Zero() => self,
            Real::One() => match self {
                Real::Zero() => rhs,
                Real::One() => Real::from_i64(2_i64),
                Real::Exact(x) => add_rationals(x, SafeRational::one()),
                _ => add_inexact(self, rhs),
            },
            Real::Exact(y) => match self {
                Real::Zero() => Real::Exact(y),
                Real::One() => add_rationals(SafeRational::one(), y),
                Real::Exact(x) => add_rationals(x, y),
                _ => add_inexact(self, Real::Exact(y)),
            }
            _ => match self {
                Real::Zero() => rhs,
                _ => add_inexact(self, rhs),
            },
        }
    }
}

fn size_in_base(n: SafeI64, base: u32) -> u32 {

    fn div(n: SafeI64, acc: u32, base: i64) -> u32 {
        if n <= SafeI64::one() { acc + 1 } else { div(n / base, acc + 1, base) }
    }

    div(n.abs(), 0, base as i64)
}

impl<'a> Mul for Real<'a> {

    type Output = Real<'a>;

    fn mul(self, rhs: Self) -> Self::Output {
        fn mul_rationals<'b>(x: SafeRational, y: SafeRational) -> Real<'b> {
            Real::from_rational(x * y)
        }

        fn mul_inexact<'b>(mut x: Real<'b>, mut y: Real<'b>) -> Real<'b> {
            Real::new_inexact(Box::new(move |p|{
                let x0 = x.eval(0).abs() + 2;
                let y0 = y.eval(0).abs() + 2;
                let sx = size_in_base(x0, 2) + 3;
                let sy = size_in_base(y0, 2) + 3;
                round_up(SafeRational::new(x.eval(p + sy) * y.eval(p + sx), SafeI64::one() << (p + sx + sy)))
            }))
        }

        match rhs {
            Real::Zero() => Real::Zero(),
            Real::One() => self,
            Real::Exact(y) => match self {
                Real::Zero() => Real::Zero(),
                Real::One() => Real::Exact(y),
                Real::Exact(x) => mul_rationals(x, y),
                _ => mul_inexact(self, Real::Exact(y)),
            }
            _ => match self {
                Real::Zero() => Real::Zero(),
                Real::One() => rhs,
                _ => mul_inexact(self, rhs),
            },
        }
    }
}