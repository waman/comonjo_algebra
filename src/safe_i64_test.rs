use std::i64;

use num::{BigInt, ToPrimitive};
use num::traits::Zero;
use num_bigint::RandomBits;
use rand::prelude::Distribution;

use crate::safe_i64::SafeI64;
// use once_cell::sync::Lazy;

// static I64_MAX: Lazy<BigInt> = Lazy::new(|| BigInt::from(i64::MAX));
// static I64_MIN: Lazy<BigInt> = Lazy::new(|| BigInt::from(i64::MIN));

// fn append_i64s_around(vec: &mut Vec<i64>, center: i64, width: i64){
//     for i in -width ..= width { 
//         match center.checked_add(i) {
//             Some(x) => vec.push(x),
//             None => (),
//         }
//     }
// }

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

fn assert_eq_safe_i64(x: SafeI64, ex: BigInt){
    assert_eq!(x.to_bigint(), ex);
    match ex.to_i64() { 
        Some(_) => assert!(x.is_primitive()),
        _ => assert!(!x.is_primitive()),
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
fn test_neg(){
    fn test(x: BigInt){
        // setup
        let x_safe = SafeI64::from_integer(x.clone());
        let exp = SafeI64::from_integer(-x);
        // exercise
        let sut = -x_safe;
        // verify
        assert_eq!(sut, exp);
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
fn test_integer_impl(){

    fn test(x: &BigInt, y: &BigInt){
        let x_safe = &SafeI64::from_integer(x.clone());
        let y_safe = &SafeI64::from_integer(y.clone());

        assert_eq_safe_i64(x_safe + y_safe, x + y);  // Add for &SafeI64
        assert_eq_safe_i64(x_safe.clone() + y_safe.clone(), x + y);  // Add for SafeI64

        assert_eq_safe_i64(x_safe - y_safe, x - y);
        assert_eq_safe_i64(x_safe.clone() - y_safe.clone(), x - y);

        assert_eq_safe_i64(x_safe * y_safe, x * y);
        assert_eq_safe_i64(x_safe.clone() * y_safe.clone(), x * y);

        if !y.is_zero() {
            assert_eq_safe_i64(x_safe / y_safe, x / y);
            assert_eq_safe_i64(x_safe.clone() / y_safe.clone(), x / y);

            assert_eq_safe_i64(x_safe % y_safe, x % y);
            assert_eq_safe_i64(x_safe.clone() % y_safe.clone(), x % y);
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
    let exp_x = "I64(123)";
    // exercise
    let sut = format!("{:?}", x);
    // verify
    assert_eq!(sut, exp_x);

    // setup
    let y = SafeI64::from(i64::MAX) + SafeI64::from(123);
    let exp_y = "BI(9223372036854775930)";  //9223372036854775807 + 123
    // exercise
    let sut = format!("{:?}", y);
    // verify
    assert_eq!(sut, exp_y);
}