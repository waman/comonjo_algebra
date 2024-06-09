use std::fmt::{Debug, Display, Formatter};
use regex::Regex;
use num::{bigint::Sign, BigInt, BigRational, One, Zero};


pub enum Real<'a>{
    // Zero(),
    // One(),
    // MinusOne(),
    // PositiveInfinity(),
    // NegativeInfinity(),
    Exact(BigRational),
    Inexact{ memo: Option<(usize, BigInt)>, f: Box<dyn FnMut(usize) -> BigInt + 'a> }
}

impl<'a> Real<'a> {


    pub fn new_inexact<'b>(f: Box<dyn FnMut(usize) -> BigInt + 'b>) -> Real<'b> {
        // Real::Inexact(f)
        Real::Inexact{ memo: None, f }
    }

    pub fn eval(&mut self, p: usize) -> BigInt {
        match self {
            Real::Exact(ref r) => 
                round_up(r * (BigInt::one() << p)),

            Real::Inexact{ memo, f } => {
                match memo {
                    Some((bits, value)) if *bits >= p => 
                        round_up(BigRational::new(value.clone(), BigInt::one() << (*bits - p))),
                    _ => {
                        let result = f(p);
                        *memo = Some((p, result.clone()));
                        result
                    },
                }
            },
        }
    }

    // fn digits() -> usize { 40 }

    // fn bits() -> usize { Real::digits_to_bits(Real::digits()) }

    pub fn sqrt<'b>(&'b mut self) -> Real<'b> {
        Real::new_inexact(Box::new(|p| self.eval(p * 2).sqrt()))
    }
}

fn round_up(r: BigRational) -> BigInt {
    let n = r.numer();
    let d = r.denom();
    if n >= &BigInt::zero() {
        let m = n % d;
        if m * 2 >= *d { n / d + BigInt::one() } else { n / d } 
    } else {
        let m = -(n % d);
        if m * 2 >= *d { n / d - BigInt::one() } else { n / d }
    }
}
  
// fn digits_to_bits(d: usize) -> usize {
//     ((d as f64) * 10_f64.ln() / 2_f64.ln()).ceil() as usize + 4
// }

fn bits_to_digits(b: usize, radix: u32) -> usize {
    ((b as f64) * 2_f64.ln() / (radix as f64).ln() + 0.001).floor() as usize
}

fn inexact_to_display_string(bit: usize, value: &BigInt, radix: u32) -> String {
    let d = bits_to_digits(bit, radix);
    let rad = BigInt::from(radix);
    let r = BigRational::new(value * rad.pow(d as u32), BigInt::one() << bit);
    let m = round_up(r);

    let (sign, abs) = match m.sign() {
        Sign::Minus  => ("-", m.magnitude().to_str_radix(radix)),
        Sign::NoSign => ("", "0".to_string()),
        Sign::Plus   => ("", m.to_str_radix(radix))
    };

    let s = if abs.len() > d {
        let i = abs.len() - d;
        format!("{}{}.{}", sign, &abs[0..i], &abs[i..])
    } else{
        let zeros = "0".repeat(d - abs.len());
        format!("{}0.{}{}", sign, zeros, abs)
    };

    let s1 = Regex::new(r"0+$").unwrap().replace_all(&s, "");
    let s2 = Regex::new(r"\.$").unwrap().replace_all(&s1, "");
    s2.to_string()

}

// impl<'a> Add<Real<'a>> for Real<'a>{
//     type Output = Real<'a>;

//     fn add(self, rhs: Real<'a>) -> Self::Output {
//         todo!()
//     }
// }

impl<'a> From<BigRational> for Real<'a> {
    fn from(value: BigRational) -> Self { Real::Exact(value) }
}

impl<'a> From<i64> for Real<'a> {
    fn from(value: i64) -> Self { Real::Exact(BigRational::from(BigInt::from(value))) }
}

impl<'a> From<BigInt> for Real<'a> {
    fn from(value: BigInt) -> Self { Real::Exact(BigRational::from(value)) }
}

impl<'a> From<(i64, i64)> for Real<'a> {
    fn from(value: (i64, i64)) -> Self { 
        Real::Exact(BigRational::new(BigInt::from(value.0), BigInt::from(value.1))) 
    }
}

#[test]
fn test_sizes_of_some_types(){
    assert_eq!(std::mem::size_of::<BigRational>(), 64);
    assert_eq!(std::mem::size_of::<Box<BigRational>>(), 8);
    assert_eq!(std::mem::size_of::<Option<(usize, BigInt)>>(), 40);
    assert_eq!(std::mem::size_of::<Box<dyn Fn(i64) -> BigInt>>(), 16);
}

#[test] #[allow(non_snake_case)]
fn test_round_BigRational_to_BigInt(){
    fn test_fn(x: (i64, i64), expected: i64){
        let x = BigRational::new(BigInt::from(x.0), BigInt::from(x.1));
        let ex = BigInt::from(expected);
        assert_eq!(round_up(x), ex);
    }

    test_fn((0, 1), 0);
    test_fn((1, 1), 1); test_fn((-1, 1), -1);
    test_fn((2, 1), 2); test_fn((-2, 1), -2);

    test_fn((0, 2), 0);
    test_fn((1, 2), 1); test_fn((-1, 2), -1);
    test_fn((2, 2), 1); test_fn((-2, 2), -1);
    test_fn((3, 2), 2); test_fn((-3, 2), -2);
    test_fn((4, 2), 2); test_fn((-4, 2), -2);

    test_fn((0, 3), 0);
    test_fn((1, 3), 0); test_fn((-1, 3), 0);
    test_fn((2, 3), 1); test_fn((-2, 3), -1);
    test_fn((3, 3), 1); test_fn((-3, 3), -1);
    test_fn((4, 3), 1); test_fn((-4, 3), -1);
    test_fn((5, 3), 2); test_fn((-5, 3), -2);
    test_fn((6, 3), 2); test_fn((-6, 3), -2);
    test_fn((7, 3), 2); test_fn((-7, 3), -2);
    test_fn((8, 3), 3); test_fn((-8, 3), -3);
    test_fn((9, 3), 3); test_fn((-9, 3), -3);
}

#[test]
fn test_eval(){
    let mut x: Real = 2.into();
    assert_eq!(x.eval(3), BigInt::from(2 << 3));
    assert_eq!(x.eval(5), BigInt::from(2 << 5));
    assert_eq!(x.eval(10), BigInt::from(2 << 10));

    let mut y: Real = (1, 3).into();  // 1/3 == 0b0.010101010101...
    assert_eq!(y.eval(3), BigInt::from(0b11));
    assert_eq!(y.eval(5), BigInt::from(0b1011));
    assert_eq!(y.eval(10), BigInt::from(0b101010101));

    let mut z = x.sqrt();  // √2 == 0b1.0110101000001001
    assert_eq!(z.eval(3), BigInt::from(0b1011));
    assert_eq!(z.eval(5), BigInt::from(0b101101));
    assert_eq!(z.eval(12), BigInt::from(0b1011010100000));
}

impl<'a> Display for Real<'a> {

    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Real::Exact(n) => f.write_str(&format!("{}", n)),

            Real::Inexact { memo, f:_ } => {
                match memo {
                    Some((bit, value)) => {
                        let s = inexact_to_display_string(*bit, value, 10);
                        f.write_str(&s)
                    },
                    _ => f.write_str("Real::Inexact(Unevaluated)")
                }
            },
        }
    }
}

impl<'a> Debug for Real<'a> {

    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Real::Exact(n) => 
                f.write_str(&format!("Real::Exact({})", n)),

            Real::Inexact { memo, f:_ } => {
                match memo {
                    Some((bit, value)) => {
                        let s = inexact_to_display_string(*bit, value, 10);
                        f.write_str(&format!("Real::Inexact({}; bit={}, value={})", s, bit, value))
                    },
                    _ => f.write_str(&format!("Real::Inexact(Unevaluated)")),
                }
            },
        }
    }
}

#[test]
fn test_display(){
    let x: Real = 123.into();
    assert_eq!(format!("{}", x), "123");

    let y: Real = (123, 7).into();
    assert_eq!(format!("{}", y), "123/7");

    let mut z: Real = 2.into();
    let mut w = z.sqrt();
    assert_eq!(format!("{}", w), "Real::Inexact(Unevaluated)");

    w.eval(27);  // 8 digits ~ 27 bits
    let s = format!("{}", w);
    assert_eq!(s, "1.41421356");
}

#[test]
fn test_debug(){
    let x: Real = 123.into();
    assert_eq!(format!("{:?}", x), "Real::Exact(123)");

    let y: Real = (123, 7).into();
    assert_eq!(format!("{:?}", y), "Real::Exact(123/7)");

    let mut z: Real = 2.into();
    let mut w = z.sqrt();
    assert_eq!(format!("{:?}", w), "Real::Inexact(Unevaluated)");

    w.eval(27);  // 8 digits ~ 27 bits
    let s = format!("{:?}", w);
    assert_eq!(s, "Real::Inexact(1.41421356; bit=27, value=189812531)");
}