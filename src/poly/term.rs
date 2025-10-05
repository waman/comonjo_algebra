use std::{borrow::Cow, collections::HashMap, fmt::Display};

use once_cell::sync::Lazy;

use crate::algebra::Semiring;

static SUPERSCRIPTS: &'static str = "⁰¹²³⁴⁵⁶⁷⁸⁹";

// [('0', '⁰'), ('1', '¹'), ...]
static MAPPING_TO_SUPERSCRIPTS: Lazy<HashMap<char, char>> = Lazy::new(||{
    SUPERSCRIPTS.chars().enumerate()
        .map(|e| (e.0.to_string().chars().next().unwrap(), e.1))
        .collect::<HashMap<char, char>>()
});

// 123_usize -> "¹²³"
fn to_superscript(p: usize) -> String {
    p.to_string().chars().map(|c|(*MAPPING_TO_SUPERSCRIPTS)[&c]).collect()
}

pub(crate) fn term_to_string<C>(exp: usize, coeff: &C) -> String where C: Semiring + Display {

    let coeff_str = format!("{}", coeff);

    if coeff.is_zero() || coeff_str == "0" {
        "".to_string()

    } else if coeff.is_one() || coeff_str == "1" {
        match exp {
            0 => " + 1".to_string(),
            1 => " + x".to_string(),
            _ => format!(" + x{}", to_superscript(exp)),
        }

    } else if coeff_str == "-1" {
        match exp {
            0 => " - 1".to_string(),
            1 => " - x".to_string(),
            _ => format!(" - x{}", to_superscript(exp)),
        }

    } else {
        if coeff_str.chars().all(|c| "0123456789-.Ee".contains(c)) {
            let exp_str = match exp {
                0 => "".to_string(),
                1 => "x".to_string(),
                _ => format!("x{}", to_superscript(exp))
            };

            if coeff_str.starts_with("-") {
                format!("{}{}", coeff_str.replace("-", " - "), exp_str)

            } else {
                format!(" + {}{}", coeff_str, exp_str)
            }

        } else {
            match exp {
                0 => if coeff_str.starts_with("-") {
                    format!(" + ({})", coeff)
                } else {
                    format!(" + {}", coeff)  // this sign may be removed
                },
                1 => format!(" + ({})x", coeff),
                _ => format!(" + ({})x{}", coeff, to_superscript(exp)),
            }       
        }
    }
}

/// Refer to spire's <a href="https://github.com/typelevel/spire/blob/main/core/src/main/scala/spire/math/poly/Term.scala">Term</a>
pub struct Term<'a, C> where C: Semiring + Clone {
    pub(crate) exp: usize,
    pub(crate) coeff: Cow<'a, C>,
}

impl<'a, C> Term<'a, C> where C: Semiring + Clone {

    pub fn from_ref(exp: usize, coeff: &'a C) -> Term<'a, C> {
        Term { exp, coeff: Cow::Borrowed(coeff) }
    }

    pub fn from_value(exp: usize, coeff: C) -> Term<'a, C> {
        Term { exp, coeff: Cow::Owned(coeff) }
    }

    pub fn is_index_zero(&self) -> bool {
        self.exp == 0
    }

    // 
    // pub fn eval(self, x: C) -> C {
    //     // if self.exp == 0 { return self.coeff.into_owned(); }
    //     // let mut x_n = x.clone();
    //     // for _ in 0..(self.exp-1) { x_n = x_n * x.clone() }
    //     // self.coeff.into_owned() * x_n

    //     fn calc_pow<C>(base: C, p: u32, extra: C) -> C where C: Semiring + Clone {
    //         if p == 1 {
    //             base * extra
    //         } else {
    //             let next_extra = if (p & 1) == 1 { base * extra } else { extra.clone() };
    //             calc_pow(base * base, p >> 1, next_extra)
    //         }
    //     }

    //     match self.exp {
    //         0 => C::one(),
    //         1 => self.coeff,
    //         _ => calc_pow(self.coeff, self.exp-1, self.coeff)
    //     }
    // }
}

impl<'a, C> Term<'a, C> where C: Semiring + Clone + Display {

    pub fn to_string(&self) -> String {
        term_to_string(self.exp, self.coeff.as_ref())
    }
}

// impl<'a, C> Term<'a, C> where C: Clone + FromPrimitive {

//     pub fn derivate(self) -> Term<'a, C> {
//         match FromPrimitive::from_usize(self.exp) {
//             Some(x) => Term::from_value(self.exp - 1, self.coeff.into_owned() * x),
//             None => panic!("can't derivate term of exp {}", self.exp),
//         }
//     }

//     pub fn integrate(self) -> Term<'a, C> {
//         match FromPrimitive::from_usize(self.exp + 1) {
//             Some(x) => Term::from_value(self.exp + 1, self.coeff.into_owned() / x),
//             None => panic!("can't integrate term of exp {}", self.exp),
//         }
//     }
// }

// impl<C> Term<C> where C: Ring + Pow<usize, Output=Term<C>>{

//     pub fn eval(&self, x: C) -> C {
//         if self.exp != 0 { self.coeff * (x.pow(self.exp)) } else { self.coeff }
//     }
// }

impl<'a, C> From<(usize, &'a C)> for Term<'a, C> where C: Semiring + Clone {

    fn from(value: (usize, &'a C)) -> Self {
        Term::from_ref(value.0, value.1)
    }
}

impl<'a, C> From<(usize, C)> for Term<'a, C> where C: Semiring + Clone {

    fn from(value: (usize, C)) -> Self {
        Term::from_value(value.0, value.1)
    }
}

// impl<'a, C> Zero for Term<'a, C> where C: Semiring + Clone {

//     fn zero() -> Self {
//         Term::from_value(0, C::zero())
//     }

//     fn is_zero(&self) -> bool {
//         self.coeff.is_zero()
//     }
// }

// impl<'a, C> One for Term<'a, C> where C: Semiring + Clone {

//     fn one() -> Self {
//         Term::from_value(0, C::one())
//     }
// }

// impl<'a, C> Neg for Term<'a, C> where C: Semiring + Clone {

//     type Output = Term<'a, C>;

//     fn neg(self) -> Self::Output {
//         Term::from_value(self.exp, -self.coeff.into_owned())
//     }
// }

// impl<'a, C> Add for Term<'a, C> where C: AdditiveSemigroup + Clone {

//     type Output = Term<'a, C>;

//     fn add(self, other: Self) -> Self::Output {
//         if self.exp == other.exp {
//             Term::from_value(self.exp, self.coeff.into_owned() + other.coeff.into_owned())
//         } else {
//             panic!("can't add terms of degree {} and ${}", self.exp, other.exp)
//         }
        
//     }
// }

// impl<'a, C> Mul for Term<'a, C> where C: Semigroup + Clone{

//     type Output = Term<'a, C>;

//     fn mul(self, other: Self) -> Self::Output {
//         Term::from_value(self.exp + other.exp, self.coeff.into_owned() * other.coeff.into_owned())
//     }
// }

// impl<'a, C> Div<C> for Term<'a, C> where C: Field + Clone {

//     type Output = Term<'a, C>;

//     fn div(self, rhs: C) -> Self::Output {
//         Term::from_value(self.exp, self.coeff.into_owned() / rhs)
//     }
// }