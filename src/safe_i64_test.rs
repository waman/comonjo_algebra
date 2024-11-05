use std::i64;

use num::pow::Pow;
use num::{BigInt, Integer, ToPrimitive};
use num::traits::Zero;
use num_bigint::RandomBits;
use rand::prelude::Distribution;

use crate::safe_i64::SafeI64;

const RAND_BITS: u64 = 66;

fn append_big_ints_around(vec: &mut Vec<BigInt>, center: BigInt, width: u64){
    let c = &center;
    let w = width as i64;
    for i in -w ..= w { 
        vec.push(c + i);
    }
}

fn new_boundary_values() -> Vec<BigInt> {
    let mut values = Vec::new();
    append_big_ints_around(&mut values, BigInt::ZERO, 2);
    append_big_ints_around(&mut values, BigInt::from(i64::MIN), 2);
    append_big_ints_around(&mut values, BigInt::from(i64::MAX), 2);
    values
}

fn assert_eq_safe_i64(x: SafeI64, ex: &BigInt){
    assert_eq!(&x.to_bigint(), ex, "{:?} != {:?}", x, ex);
    match ex.to_i64() { 
        Some(_) => assert!( x.is_primitive(), "{:?} must be primitive.", x),
        _       => assert!(!x.is_primitive(), "{:?} must not be primitive.", x),
    }
}

#[test]
fn test_size(){
    assert_eq!(std::mem::size_of::<i64>(), 8);  // <- use
    assert_eq!(std::mem::size_of::<i128>(), 16);
    assert_eq!(std::mem::size_of::<BigInt>(), 32);
    assert_eq!(std::mem::size_of::<&BigInt>(), 8);
    assert_eq!(std::mem::size_of::<Box<BigInt>>(), 8);  // <- use
}

#[test]
fn test_is_primitive(){
    fn test(x: BigInt){
        match x.to_i64() {
            Some(i) => {
                assert!(SafeI64::from_integer(x).is_primitive());
                assert!(SafeI64::from(i).is_primitive());
            }
            None => 
                assert!(!SafeI64::from_integer(x).is_primitive()), 
        }
    }
    
    // boundary values
    for x in new_boundary_values() { 
        test(x);
    }

    // random values
    let rng = &mut rand::thread_rng();
    let rand_bits = RandomBits::new(RAND_BITS);
    for _ in 0..100 {
        test(rand_bits.sample(rng));
    }
}

#[test]
fn test_num_binary_op(){

    fn test(x: &BigInt, y: &BigInt){
        let x_safe = &SafeI64::from_integer(x.clone());
        let y_safe = &SafeI64::from_integer(y.clone());

        let ex_add = x + y;
        assert_eq_safe_i64(x_safe + y_safe, &ex_add);  // Add for &SafeI64
        assert_eq_safe_i64(x_safe.clone() + y_safe.clone(), &ex_add);  // Add for SafeI64

        let ex_sub = x - y;
        assert_eq_safe_i64(x_safe - y_safe, &ex_sub);
        assert_eq_safe_i64(x_safe.clone() - y_safe.clone(), &ex_sub);

        let ex_mul = x * y;
        assert_eq_safe_i64(x_safe * y_safe, &ex_mul);
        assert_eq_safe_i64(x_safe.clone() * y_safe.clone(), &ex_mul);

        if !y.is_zero() {
            let ex_div = x / y;
            assert_eq_safe_i64(x_safe / y_safe, &ex_div);
            assert_eq_safe_i64(x_safe.clone() / y_safe.clone(), &ex_div);

            let ex_rem = x % y;
            assert_eq_safe_i64(x_safe % y_safe, &ex_rem);
            assert_eq_safe_i64(x_safe.clone() % y_safe.clone(), &ex_rem);
        }
    }

    for x in &new_boundary_values() {
        for y in &new_boundary_values() {
            test(x, y);
        }
    }

    let rng = &mut rand::thread_rng();
    let rand_bits = RandomBits::new(RAND_BITS);

    for _ in 0..100 {
        test(&rand_bits.sample(rng), &rand_bits.sample(rng));
    }
}

#[test]
fn test_pow_and_bit_shift(){

    fn test(x: &BigInt, p: u32){
        let x_safe = &SafeI64::from_integer(x.clone());

        let ex_pow = x.pow(p);
        assert_eq_safe_i64(x_safe.pow(p), &ex_pow);  // Pow for &SafeI64
        assert_eq_safe_i64(x_safe.clone().pow(p), &ex_pow);  // Pow for SafeI64

        let ex_shl = x << p;
        assert_eq_safe_i64(x_safe << p, &ex_shl);
        assert_eq_safe_i64(x_safe.clone() << p, &ex_shl);

        let ex_shr = x >> p;
        assert_eq_safe_i64(x_safe >> p, &ex_shr);
        assert_eq_safe_i64(x_safe.clone() >> p, &ex_shr);
    }

    let u32_boundary: &[u32] = &[0, 1, 2, 3, 10, 62, 63, 64, 65, 100];

    for x in &new_boundary_values() {
        for p in u32_boundary {
            test(x, *p);
        }
    }

    let rng = &mut rand::thread_rng();
    let rand_bits = RandomBits::new(RAND_BITS);

    for _ in 0..20 {
        for p in u32_boundary {
            test(&rand_bits.sample(rng), *p);
        }
    }
}

#[test]
fn test_neg(){

    fn test(x: &BigInt){
        let x_safe = &SafeI64::from_integer(x.clone());

        let ex_neg = -x;
        assert_eq_safe_i64(-x_safe, &ex_neg);  // Neg for &SafeI64
        assert_eq_safe_i64(-(x_safe.clone()), &ex_neg);  // Neg for SafeI64
    }

    // boundary values
    for x in &new_boundary_values() {
        test(x);
    }

    // random values
    let rng = &mut rand::thread_rng();
    let rand_bits = RandomBits::new(RAND_BITS);

    for _ in 0..100 {
        test(&rand_bits.sample(rng));
    }
}

#[test]
fn test_unary_integer_methods(){
    fn test(x: &BigInt){
        let x_safe = &SafeI64::from_integer(x.clone());

        let ex_is_odd = x.is_odd();
        assert_eq!(x_safe.is_odd() , ex_is_odd);  // is_odd() for &SafeI64
        assert_eq!(x_safe.clone().is_odd() , ex_is_odd);  // is_odd() for SafeI64

        let ex_is_even = x.is_even();
        assert_eq!(x_safe.is_even() , ex_is_even);
        assert_eq!(x_safe.clone().is_even() , ex_is_even);
    }

    // boundary values
    for x in &new_boundary_values() {
        test(x);
    }

    // random values
    let rng = &mut rand::thread_rng();
    let rand_bits = RandomBits::new(RAND_BITS);

    for _ in 0..20 {
        test(&rand_bits.sample(rng));
    }
}

#[test]
fn test_integer_binary_op(){

    fn test(x: &BigInt, y: &BigInt){
        let x_safe = &SafeI64::from_integer(x.clone());
        let y_safe = &SafeI64::from_integer(y.clone());

        assert_eq_safe_i64(x_safe.gcd(y_safe), &x.gcd(y));
        assert_eq_safe_i64(x_safe.lcm(y_safe), &x.lcm(y));

        // println!("{}.is_multiple_of({}) ?= {}", x_safe, y_safe, x.is_multiple_of(y));
        assert_eq!(x_safe.is_multiple_of(y_safe), x.is_multiple_of(y));

        if !y_safe.is_zero() {
            assert_eq_safe_i64(x_safe.div_floor(y_safe), &x.div_floor(y));
            assert_eq_safe_i64(x_safe.mod_floor(y_safe), &x.mod_floor(y));
            
            // div_rem()
            let dr = x_safe.div_rem(y_safe);
            let ex_dr = x.div_rem(y);
            assert_eq_safe_i64(dr.0, &ex_dr.0);
            assert_eq_safe_i64(dr.1, &ex_dr.1);
        }
    }

    for x in &new_boundary_values() {
        for y in &new_boundary_values() {
            test(x, y);
        }
    }

    let rng = &mut rand::thread_rng();
    let rand_bits = RandomBits::new(RAND_BITS);

    for _ in 0..100 {
        test(&rand_bits.sample(rng), &rand_bits.sample(rng));
    }
}


#[test]
fn test_display(){
    // setup
    let x = SafeI64::from(123);
    let exp_x = "123";
    // exercise
    let sut = format!("{}", x);
    // verify
    assert_eq!(sut, exp_x);

    // setup
    let y = SafeI64::from(i64::MAX) + SafeI64::from(123);
    let exp_y = "9223372036854775930";   //9223372036854775807 + 123
    // exercise
    let sut = format!("{}", y);
    // verify
    assert_eq!(sut, exp_y);
}

#[test]
fn test_debug(){
    // setup
    let x = SafeI64::from(123);
    let exp_x = "SafeI64::I64(123)";
    // exercise
    let sut = format!("{:?}", x);
    // verify
    assert_eq!(sut, exp_x);

    // setup
    let y = SafeI64::from(i64::MAX) + SafeI64::from(123);
    let exp_y = "SafeI64::BI(9223372036854775930)";  //9223372036854775807 + 123
    // exercise
    let sut = format!("{:?}", y);
    // verify
    assert_eq!(sut, exp_y);
}