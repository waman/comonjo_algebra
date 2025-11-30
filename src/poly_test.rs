use std::{collections::HashMap, vec};

use num::{complex::c64, pow::Pow, One, Rational64, Zero};

use crate::{algebra::Semiring, dense, poly::{CoeffsIterator, Compose, Differentiable, EuclideanRingPolyOps, FieldPolyOps, Integrable, Polynomial, RingPolyOps, SemiringPolyOps}, sparse};

type PolyR64 = Polynomial<Rational64>;

fn zero<C: Semiring>() -> Polynomial<C> { Polynomial::zero() }
fn one<C: Semiring + Clone>() -> Polynomial<C> { Polynomial::one() }
fn cst<C: Semiring>(v: C) -> Polynomial<C> { Polynomial::constant(v) }

fn to_poly_r64(p: Polynomial<i64>) -> PolyR64 { p.map_nonzero(|_, c| ri(c)) }

fn r(n: i64, d: i64) -> Rational64 { Rational64::new(n, d) }
fn ri(n: i64) -> Rational64 { Rational64::new(n, 1) }

fn get_impls<'a, C>(x: &'a Polynomial<C>) -> Vec<Polynomial<C>> where C: Semiring + Clone {
    if x.degree() == 0 {
        vec![x.clone()]
    } else {
        vec![x.dense_clone(), x.sparse_clone()]
    }
}

// ⁰¹²³⁴⁵⁶⁷⁸⁹
#[test]
fn test_dense_macro(){
    let p0 = dense![1, 2, 3];
    assert_eq!(format!("{}", p0), "1 + 2x + 3x²");

    let p1 = dense![0, 4, 5, 6, 0, 7, 0, 0];
    assert_eq!(format!("{}", p1), "4x + 5x² + 6x³ + 7x⁵");

    let p2 = dense![0, 0, 8];
    assert_eq!(format!("{}", p2), "8x²");
}

#[test]
fn test_sparse_macro(){
    let p0 = sparse![(0, 1), (1, 2), (2, 3)];
    assert_eq!(format!("{}", p0), "1 + 2x + 3x²");

    let p1 = sparse![(0, 0), (1, 4), (2, 5), (5, 7), (6, 0), (3, 6)];
    assert_eq!(format!("{}", p1), "4x + 5x² + 6x³ + 7x⁵");


    let p2 = sparse![(2, 8)];
    assert_eq!(format!("{}", p2), "8x²");
}

#[test]
fn test_display(){
    fn test(p: Polynomial<i64>, exp: &str){
        assert_eq!(format!("{}", p), exp);
    }

    let table = [
        (zero(), "(0)"),
        (one(), "(1)"),

        (cst(0),  "(0)"),
        (cst(1),  "(1)"),
        (cst(-1), "(-1)"),
        (cst(2),  "(2)"),

        (dense![],      "(0)"),
        (dense![0, 0], "(0)"),
        (dense![1],    "(1)"),
        (dense![-1],   "(-1)"),
        (dense![2],    "(2)"),

        (dense![1, 2, 3],        "1 + 2x + 3x²"),
        (dense![-1, -2, -3],     "-1 - 2x - 3x²"),
        (dense![0, 2, 3],        "2x + 3x²"),
        (dense![1, 2, 3, 0, 0],  "1 + 2x + 3x²"),
        (dense![1, 1, 1, 1],     "1 + x + x² + x³"),
        (dense![-1, -1, -1, -1], "-1 - x - x² - x³"),


        (sparse![], "(0)"),
        (sparse![(0, 0), (1, 0)], "(0)"),
        (sparse![(0, 1)],  "(1)"),
        (sparse![(0, -1)], "(-1)"),
        (sparse![(0, 2)],   "(2)"),

        (sparse![(0, 1), (2, 3)],   "1 + 3x²"),
        (sparse![(0, -1), (2, -3)], "-1 - 3x²"),
        (sparse![(1, 2), (2, 3)],   "2x + 3x²"),
        (sparse![(0, 1), (2, 3), (4, 0)], "1 + 3x²"),
        (sparse![(0, 1), (1, 1), (2, 1), (3, 1)],     "1 + x + x² + x³"),
        (sparse![(0, -1), (1, -1), (2, -1), (3, -1)], "-1 - x - x² - x³"),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_display_for_rational_and_complex(){
    // rational
    let p_rat = dense![r(1, 2), r(1, 3), r(2, 3)];
    assert_eq!(format!("{}", p_rat), "1/2 + (1/3)x + (2/3)x²");

    // complex
    let p_cx = dense![c64(1., 2.), c64(3., 4.), c64(5., 6.)];
    assert_eq!(format!("{}", p_cx), "1+2i + (3+4i)x + (5+6i)x²");
}

#[test]
fn test_debug(){
    fn test(p: Polynomial<i64>, exp: &str){
        assert_eq!(format!("{:?}", p), exp);
    }

    let table = [
        (zero(), "Polynomial::<i64>::Zero[(0)]"),
        (one(), "Polynomial::<i64>::Constant[(1)]"),
        (cst(2),  "Polynomial::<i64>::Constant[(2)]"),
        (dense![1, 2, 3],  "Polynomial::<i64>::Dense[1 + 2x + 3x²]"),
        (dense![-1, -2, -3],  "Polynomial::<i64>::Dense[-1 - 2x - 3x²]"),
        (sparse![(0, 1), (1, 2), (2, 3)],   "Polynomial::<i64>::Sparse[1 + 2x + 3x²]"),
        (sparse![(0, -1), (1, -2), (2, -3)],   "Polynomial::<i64>::Sparse[-1 - 2x - 3x²]"),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_degree_and_is_xxx_methods(){

    fn test(p: Polynomial<i64>, exp_deg: usize, exp_zero: bool, exp_one: bool, exp_const: bool, exp_x: bool){
        assert_eq!(p.degree(), exp_deg);
        assert_eq!(p.is_zero(), exp_zero);
        assert_eq!(p.is_one(), exp_one);
        assert_eq!(p.is_constant(), exp_const);
        assert_eq!(p.is_x(), exp_x)
    }

    let table = [ // deg, zero, one, const, x
        (zero(), 0, true,  false, true, false),
        (one(), 0, false, true,  true, false),

        (cst(0),  0, true,  false, true, false),
        (cst(1),  0, false, true,  true, false),
        (cst(-1), 0, false, false, true, false),
        (cst(2),  0, false, false, true, false),


        (dense![],     0, true,  false, true, false),  // zero & const
        (dense![0, 0], 0, true,  false, true, false),  // zero & const
        (dense![1],    0, false, true,  true, false),  // const
        (dense![-1],   0, false, false, true, false),  // const
        (dense![2],    0, false, false, true, false),  // const

        (dense![0, 1],          1, false, false, false, true),  // x
        (dense![0, 1, 0],       1, false, false, false, true),  // x
        (dense![0, 2],          1, false, false, false, false),  // 2x
        (dense![1, 2, 3],       2, false, false, false, false),
        (dense![-1, -2, -3],    2, false, false, false, false),
        (dense![0, 2, 3],       2, false, false, false, false),
        (dense![1, 2, 3, 0, 0], 2, false, false, false, false),


        (sparse![],               0, true,  false, true, false),  // zero & const
        (sparse![(0, 0), (2, 0)], 0, true,  false, true, false),  // zero & const
        (sparse![(0, 1)],         0, false, true,  true, false),  // const
        (sparse![(0, -1)],        0, false, false, true, false),  // const
        (sparse![(0, 2)],         0, false, false, true, false),  // const

        (sparse![(1, 1)],                                  1, false, false, false, true),  // x
        (sparse![(1, 2)],                                  1, false, false, false, false),  // 2x
        (sparse![(2, 3)],                                  2, false, false, false, false),
        (sparse![(0, 1), (1, 2), (2, 3)],                  2, false, false, false, false),
        (sparse![(0, -1), (1, -2), (2, -3)],               2, false, false, false, false),
        (sparse![(0, 0), (1, 2), (2, 3)],                  2, false, false, false, false),
        (sparse![(0, 1), (1, 2), (2, 3), (5, 0), (10, 0)], 2, false, false, false, false),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2, entry.3, entry.4, entry.5);
    }
}

#[test]
fn test_is_dense_and_is_sparse(){
    
    fn test(p: Polynomial<i64>, exp_be_dense: bool, exp_be_sparse: bool){
        assert_eq!(p.is_dense(), exp_be_dense);
        assert_eq!(p.is_sparse(), exp_be_sparse);
    }

    let table = [ // deg, zero, one, const, x
        (zero(), false, false),
        (cst(2), false, false),
        (dense![1, 2, 3], true, false),
        (sparse![(0, 1), (1, 2), (2, 3)], false, true),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_coeffs_iter_methods(){
    fn test(p: Polynomial<i64>, exp: Vec<i64>){
        let nonzero_exp: HashMap<usize, i64> = exp.clone().into_iter().enumerate().filter(|(_, c)|*c != 0).collect();

        // into_iter() for Polynomial
        let cs0: HashMap<usize, i64> = p.clone().nonzero_coeffs().collect();
        assert_eq!(cs0, nonzero_exp);

        // nonzero_coeffs() for Polynomial (= into_iter())
        let cs1: HashMap<usize, i64> = p.clone().nonzero_coeffs().collect();
        assert_eq!(cs1, nonzero_exp);

        // coeffs() for Polynomial
        let cs2: Vec<i64> = p.clone().coeffs().collect();
        assert_eq!(cs2, exp);

        let p_ref: &Polynomial<i64> = &p;

        // into_iter() for &Polynomial
        let cs3: HashMap<usize, i64> = p_ref.into_iter().map(|(e, c)|(e, *c)).collect();
        assert_eq!(cs3, nonzero_exp);

        // nonzero_coeffs() for &Polynomial (= into_iter())
        let cs4: HashMap<usize, i64> = p_ref.nonzero_coeffs().map(|(e, c)|(e, *c)).collect();
        assert_eq!(cs4, nonzero_exp);

        // coeffs() for &Polynomial
        let cs5: Vec<i64> = p_ref.coeffs().map(|c| match c {
            Some(c) => *c,
            _ => 0,
        }).collect();
        assert_eq!(cs5, exp);
    }

    let v = vec![1, 2, 0, 3, 0, 4, 0, 0, 5];
    let table = [
        (Polynomial::zero(), vec![]),
        (Polynomial::one(), vec![1]),
        (cst(5), vec![5]),
        (Polynomial::dense_from_vec(v.clone()), v.clone()),
        (sparse![(0, 1), (1, 2), (3, 3), (5, 4), (8, 5)], v.clone())
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_eq(){
    let p_dense = dense![1, 0, 2, 4, 0, 0, 1];
    assert_eq!(p_dense, dense![1, 0, 2, 4, 0, 0, 1]);

    let p_sparse = sparse![(0, 1), (2, 2), (3, 4), (6, 1)];
    assert_eq!(p_sparse, sparse![(0, 1), (2, 2), (3, 4), (6, 1)]);

    assert_eq!(p_dense, p_sparse);
}

#[test]
fn test_clone_methods(){  // clone(), dense_clone(), sparse_clone(), to_dense() and to_sparse()
    fn test_with_clones<F>(p: Polynomial<i64>, test: F) where F: Fn(Polynomial<i64>) {
        for p_clone in [
            p.clone(),
            p.dense_clone(),
            p.sparse_clone(),
            p.clone().to_dense(),
            p.clone().to_sparse()
        ]{
            test(p_clone)
        }
    }

    test_with_clones(zero(), |p| assert!(p.is_zero()));
    test_with_clones(one(), |p| assert!(p.is_one()));
    test_with_clones(cst(5), |p| {
        assert!(p.is_constant());
        assert_eq!(p.nth(0), Some(&5));
    });
    
    // 1 + 2x² + 4x³ + x⁶
    let pd = dense![1, 0, 2, 4, 0, 0, 1];
    let ps = sparse![(0, 1), (2, 2), (3, 4), (6, 1)];

    test_with_clones(dense![1, 0, 2, 4, 0, 0, 1], |p|{
        assert_eq!(p, pd);
        assert_eq!(p, ps);
    });
    test_with_clones(sparse![(0, 1), (2, 2), (3, 4), (6, 1)], |p|{
        assert_eq!(p, pd);
        assert_eq!(p, ps);
    });
}

/// 1 + 2x + 3x² (= dense![[1, 2, 3]])
fn p0() -> Polynomial<i64> { dense![1, 2, 3] }

/// 4 + 5x + 6x³ + 7x⁴ (= dense![[4, 5, 0, 6, 7]])
fn p1() -> Polynomial<i64> { dense![4, 5, 0, 6, 7] }

/// 4 + 5x³ + 6x⁷ (= dense![[4, 0, 0, 5, 0, 0, 0, 6]])
fn p2() -> Polynomial<i64> { dense![4, 0, 0, 5, 0, 0, 0, 6] }

/// 1 + 2x⁴ (= dense![[1, 0, 0, 0, 2]])
fn p3() -> Polynomial<i64> { dense![1, 0, 0, 0, 2] }

/// 6x³ (= dense![[0, 0, 0, 6]])
fn p4() -> Polynomial<i64> { dense![0, 0, 0, 6] }

/// 5x + 7x⁴ (= dense![[0, 5, 0, 0, 7]])
fn p5() -> Polynomial<i64> { dense![0, 5, 0, 0, 7] }
    
/// 1 + 2x + 3x² (= dense![[1, 2, 3]])
fn pr0() -> PolyR64 { to_poly_r64(p0()) }

/// 4 + 5x + 6x³ + 7x⁴ (= dense![[4, 5, 0, 6, 7]])
fn pr1() -> PolyR64 { to_poly_r64(p1()) }

/// 4 + 5x³ + 6x⁷ (= dense![[4, 0, 0, 5, 0, 0, 0, 6]])
fn pr2() -> PolyR64 { to_poly_r64(p2()) }

/// 1 + 2x⁴ (= dense![[1, 0, 0, 0, 2]])
fn pr3() -> PolyR64 { to_poly_r64(p3()) }

/// 6x³ (= dense![[0, 0, 0, 6]])
fn pr4() -> PolyR64 { to_poly_r64(p4()) }

#[test]
fn test_neg(){

    fn test(x: Polynomial<i64>, exp: Polynomial<i64>){
    
        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, exp: &'b Polynomial<i64>){
            assert_eq!(-x, *exp);
            assert_eq!(-(x.clone()), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, &exp);
        }
    }

    let table = [
        (zero(), zero()),
        (one(),  cst(-1)),
        (cst(5), cst(-5)),
        (cst(-4), cst(4)),

        (p0(), dense![-1, -2, -3]),
        (p1(), dense![-4, -5, 0, -6, -7]),
        (p2(), dense![-4, 0, 0, -5, 0, 0, 0, -6]),
        (p3(), dense![-1, 0, 0, 0, -2]),
        (p4(), dense![0, 0, 0, -6]),
        (p5(), dense![0, -5, 0, 0, -7]),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_add(){

    fn test(x: Polynomial<i64>, y: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b, 'c>(x: &'a Polynomial<i64>, y: &'b Polynomial<i64>, exp: &'c Polynomial<i64>){
            assert_eq!(x + y, *exp);
            assert_eq!(x.clone() + y, *exp);
            assert_eq!(x + y.clone(), *exp);
            assert_eq!(x.clone() + y.clone(), *exp);
        }

        for x_ in get_impls(&x) {
            for y_ in get_impls(&y) {
                test_op(&x_, &y_, &exp);
                test_op(&y_, &x_, &exp);
            }
        }
    }

    let table = [
        (zero(), zero(), zero()),
        (zero(), cst(3), cst(3)),
        (zero(), p0(), dense![1, 2, 3]),
        (zero(), p1(), dense![4, 5, 0, 6, 7]),
        (zero(), p2(), dense![4, 0, 0, 5, 0, 0, 0, 6]),
        (zero(), p3(), dense![1, 0, 0, 0, 2]),
        (zero(), p4(), dense![0, 0, 0, 6]),
        
        (cst(5), cst(3), cst(8)),
        (cst(5), cst(-5), zero()),  // result be zero
        (cst(5), p0(), dense![6, 2, 3]),
        (cst(5), p1(), dense![9, 5, 0, 6, 7]),
        (cst(5), p2(), dense![9, 0, 0, 5, 0, 0, 0, 6]),
        (cst(5), p3(), dense![6, 0, 0, 0, 2]),
        (cst(5), p4(), dense![5, 0, 0, 6]),
        
        (cst(-4), cst(3), cst(-1)),
        (cst(-4), p0(), dense![-3, 2, 3]),
        (cst(-4), p1(), dense![0, 5, 0, 6, 7]),
        (cst(-4), p2(), dense![0, 0, 0, 5, 0, 0, 0, 6]),
        (cst(-4), p3(), dense![-3, 0, 0, 0, 2]),
        (cst(-4), p4(), dense![-4, 0, 0, 6]),
        
        (p1(), -p1(), zero()),                    // result be zero
        (p1(), dense![-1, -5, 0, -6, -7], cst(3)), // result be const
        (p0(), p1(), dense![5, 7, 3, 6, 7]),
        (p1(), p2(), dense![8, 5, 0, 11, 7, 0, 0, 6]),
        (p3(), p4(), dense![1, 0, 0, 6, 2]),
        (p1(), -p4(), dense![4, 5, 0, 0, 7]),
        (p1(), -p5(), dense![4, 0, 0, 6]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_add_c(){

    fn test(x: Polynomial<i64>, y: i64, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, y: i64, exp: &'b Polynomial<i64>){
            assert_eq!(x + &y, *exp);
            assert_eq!(x.clone() + &y, *exp);
            assert_eq!(x + y, *exp);
            assert_eq!(x.clone() + y, *exp);

            assert_eq!(&y + x, *exp);
            assert_eq!(&y + x.clone(), *exp);
            assert_eq!(y + x, *exp);
            assert_eq!(y + x.clone(), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, y, &exp);
        }
    }

    let table = [
        (zero(), 0,  zero()),
        (zero(), 3,  cst(3)),
        (zero(), -4, cst(-4)),
        
        (cst(5), 0,  cst(5)),
        (cst(5), 3,  cst(8)),
        (cst(5), -5, zero()),
        
        (cst(-4), 0,  cst(-4)),
        (cst(-4), 4,  zero()),
        (cst(-4), -5, cst(-9)),
        
        (p1(), 0, p1()),
        (p1(), 3, dense![7, 5, 0, 6, 7]),
        (p1(), -4, dense![0, 5, 0, 6, 7]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}


#[test]
fn test_sub(){

    fn test(x: Polynomial<i64>, y: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b, 'c>(x: &'a Polynomial<i64>, y: &'b Polynomial<i64>, exp: &'c Polynomial<i64>){
            assert_eq!(x - y, *exp);
            assert_eq!(x.clone() - y, *exp);
            assert_eq!(x - y.clone(), *exp);
            assert_eq!(x.clone() - y.clone(), *exp);
        }

        for x_ in get_impls(&x) {
            for y_ in get_impls(&y) {
                test_op(&x_, &y_, &exp);
                test_op(&y_, &x_, &-(exp.clone()));
            }
        }
    }

    let table = [
        (zero(), zero(), zero()),
        (zero(), cst(3), cst(-3)),
        (zero(), cst(-4), cst(4)),
        (zero(), p0(), dense![-1, -2, -3]),
        (zero(), p1(), dense![-4, -5, 0, -6, -7]),
        (zero(), p2(), dense![-4, 0, 0, -5, 0, 0, 0, -6]),
        (zero(), p3(), dense![-1, 0, 0, 0, -2]),
        (zero(), p4(), dense![0, 0, 0, -6]),
        
        (cst(5), cst(3), cst(2)),
        (cst(5), cst(5), zero()),  // result be zero
        (cst(5), p0(), dense![4, -2, -3]),
        (cst(5), p1(), dense![1, -5, 0, -6, -7]),
        (cst(5), p2(), dense![1, 0, 0, -5, 0, 0, 0, -6]),
        (cst(5), p3(), dense![4, 0, 0, 0, -2]),
        (cst(5), p4(), dense![5, 0, 0, -6]),
        
        (cst(4), p1(), dense![0, -5, 0, -6, -7]),
        (cst(4), p2(), dense![0, 0, 0, -5, 0, 0, 0, -6]),
        
        (p0(), p0(), zero()),             // result be zero
        (p0(), dense![-2, 2, 3], cst(3)), // result be const
        (p0(), p1(), dense![-3, -3, 3, -6, -7]),
        (p1(), p2(), dense![0, 5, 0, 1, 7, 0, 0, -6]),
        (p3(), p4(), dense![1, 0, 0, -6, 2]),
        (p1(), p4(), dense![4, 5, 0, 0, 7]),
        (p1(), p5(), dense![4, 0, 0, 6]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_sub_c(){

    fn test(x: Polynomial<i64>, y: i64, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, y: i64, exp: &'b Polynomial<i64>){
            assert_eq!(x - &y, *exp);
            assert_eq!(x.clone() - &y, *exp);
            assert_eq!(x - y, *exp);
            assert_eq!(x.clone() - y, *exp);

            let exp_neg = &(-exp);
            assert_eq!(&y - x, *exp_neg);
            assert_eq!(&y - x.clone(), *exp_neg);
            assert_eq!(y - x, *exp_neg);
            assert_eq!(y - x.clone(), *exp_neg);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, y, &exp);
        }
    }

    let table = [
        (zero(), 0,  zero()),
        (zero(), 3,  cst(-3)),
        (zero(), -4, cst(4)),
        
        (cst(5), 0,  cst(5)),
        (cst(5), 3,  cst(2)),
        (cst(5), -4, cst(9)),
        
        (cst(-4), 0,  cst(-4)),
        (cst(-4), 3,  cst(-7)),
        (cst(-4), -4, zero()),
        
        (p1(), 0, p1()),
        (p1(), 3, dense![1, 5, 0, 6, 7]),
        (p1(), 4, dense![0, 5, 0, 6, 7]),
        (p1(), -4, dense![8, 5, 0, 6, 7]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}


#[test]
fn test_mul(){

    fn test(x: Polynomial<i64>, y: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b, 'c>(x: &'a Polynomial<i64>, y: &'b Polynomial<i64>, exp: &'c Polynomial<i64>){
            assert_eq!(x * y, *exp);
            assert_eq!(x.clone() * y, *exp);
            assert_eq!(x * y.clone(), *exp);
            assert_eq!(x.clone() * y.clone(), *exp);
        }

        for x_ in get_impls(&x) {
            for y_ in get_impls(&y) {
                test_op(&x_, &y_, &exp);
                test_op(&y_, &x_, &exp);
            }
        }
    }

    fn exp_p0_mul_p1() -> Polynomial<i64> { dense![4, 13, 22, 21, 19, 32, 21] }
    fn exp_p2_mul_p3() -> Polynomial<i64> { sparse![(0, 4), (3, 5), (4, 8), (7, 16), (11, 12)] }

    let table = [
        (zero(), zero(),  zero()),
        (zero(), one(),   zero()),
        (zero(), cst(3),  zero()),
        (zero(), cst(-4), zero()),
        (zero(), p0(),    zero()),
        
        (one(), one(),   one()),
        (one(), cst(3),  cst(3)),
        (one(), cst(-4), cst(-4)),
        (one(), p0(),    p0()),
        
        (cst(5), cst(3), cst(15)),
        (cst(5), cst(-4), cst(-20)),
        (cst(5), p0(),   dense![5, 10, 15]),
        (cst(5), p2(),   sparse![(0, 20), (3, 25), (7, 30)]),
        
        (cst(-4), cst(3), cst(-12)),
        (cst(-4), cst(-5), cst(20)),
        (cst(-4), p0(),   dense![-4, -8, -12]),
        (cst(-4), p2(),   sparse![(0, -16), (3, -20), (7, -24)]),
        
        (p1(), zero(), zero()),
        (p1(), one(),  p1()),
        (p1(), cst(3), dense![12, 15, 0, 18, 21]),
        (p1(), cst(-4), dense![-16, -20, 0, -24, -28]),
        (p0(), p1(),   exp_p0_mul_p1()),
        (p2(), p3(),   exp_p2_mul_p3()),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_mul_by_c(){

    fn test(x: Polynomial<i64>, y: i64, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, y: i64, exp: &'b Polynomial<i64>){
            assert_eq!(x * &y, *exp);
            assert_eq!(x.clone() * &y, *exp);
            assert_eq!(x * y, *exp);
            assert_eq!(x.clone() * y, *exp);

            assert_eq!(&y * x, *exp);
            assert_eq!(&y * x.clone(), *exp);
            assert_eq!(y * x, *exp);
            assert_eq!(y * x.clone(), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, y, &exp);
        }
    }

    let table = [
        (zero(), 0,  zero()),
        (zero(), 3,  zero()),
        (zero(), -4, zero()),
        
        (cst(5), 0,  zero()),
        (cst(5), 3,  cst(15)),
        (cst(5), -4, cst(-20)),
        
        (cst(-4), 0,  zero()),
        (cst(-4), 3,  cst(-12)),
        (cst(-4), -4, cst(16)),
        
        (p1(), 0, zero()),
        (p1(), 3, dense![12, 15, 0, 18, 21]),
        (p1(), -4, dense![-16, -20, 0, -24, -28]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_div_rem(){

    fn pr1_div_3() -> PolyR64 { dense![r(4, 3), r(5, 3), ri(0), r(6, 3), r(7, 3)] }
    fn pr1_div_m4() -> PolyR64 { dense![ri(-1), r(-5, 4), ri(0), r(-3, 2), r(-7, 4)] }

    fn pr1_div_pr0() -> PolyR64 { dense![r(-29, 27), r(4, 9), r(7, 3)] }
    fn pr1_rem_pr0() -> PolyR64 { dense![r(137, 27), r(181, 27)] }

    fn pr2_div_pr3() -> PolyR64 { to_poly_r64(dense![0, 0, 0, 3]) }
    fn pr2_rem_pr3() -> PolyR64 { to_poly_r64(dense![4, 0, 0, 2]) }
    
    fn pr2_div_pr4() -> PolyR64 { dense![r(5, 6), ri(0), ri(0), ri(0), ri(1)] }
    fn pr2_rem_pr4() -> PolyR64 { dense![ri(4)] }


    fn test(x: PolyR64, y: PolyR64, exp_div: PolyR64, exp_rem: PolyR64){

        fn test_op<'a, 'b, 'c, 'd>(x: &'a PolyR64, y: &'b PolyR64, exp_div: &'c PolyR64, exp_rem: &'d PolyR64){
            assert_eq!(x / y, *exp_div);
            assert_eq!(x % y, *exp_rem);

            assert_eq!(x / y.clone(), *exp_div);
            assert_eq!(x % y.clone(), *exp_rem);

            assert_eq!(x.clone() / y, *exp_div);
            assert_eq!(x.clone() % y, *exp_rem);

            assert_eq!(x.clone() / y.clone(), *exp_div);
            assert_eq!(x.clone() % y.clone(), *exp_rem);

            assert_eq!(x.clone(), y * exp_div + exp_rem);
        }

        for x_ in get_impls(&x) {
            for y_ in get_impls(&y) {
                test_op(&x_, &y_, &exp_div, &exp_rem);
            }
        }
    }

    let table = [
        (zero(), one(),       zero(), zero()),
        (zero(), cst(ri(3)),  zero(), zero()),
        (zero(), cst(ri(-4)), zero(), zero()),
        (zero(), pr0(),       zero(), zero()),
        
        (cst(ri(5)), one(),      cst(ri(5)),   zero()),  
        (cst(ri(5)), cst(ri(3)), cst(r(5, 3)), zero()),
        (cst(ri(5)), pr0(),      zero(),       cst(ri(5))),
        
        (cst(ri(-4)), one(),      cst(ri(-4)),   zero()),  
        (cst(ri(-4)), cst(ri(3)), cst(r(-4, 3)), zero()),
        (cst(ri(-4)), pr0(),      zero(),        cst(ri(-4))),
        
        (pr1(),         one(),       pr1(),         zero()),
        (pr1(),         cst(ri(3)),  pr1_div_3(),   zero()),
        (pr1(),         cst(ri(-4)), pr1_div_m4(),  zero()),
        (pr1(),         pr0(),       pr1_div_pr0(), pr1_rem_pr0()),
        (pr2() * pr3(), pr3(),       pr2(),         zero()),
        (pr2() * pr3(), pr2(),       pr3(),         zero()),
        (pr2(),         pr3(),       pr2_div_pr3(), pr2_rem_pr3()),
        (pr2(),         pr4(),       pr2_div_pr4(), pr2_rem_pr4()),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2, entry.3);
    }
}

#[test]
fn test_div_by_c(){

    fn test(x: PolyR64, y: Rational64, exp: PolyR64){

        fn test_op<'a, 'b>(x: &'a PolyR64, y: Rational64, exp: &'b PolyR64){
            assert_eq!(x / &y, *exp);
            assert_eq!(x.clone() / &y, *exp);
            assert_eq!(x / y, *exp);
            assert_eq!(x.clone() / y, *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, y, &exp);
        }
    }

    let table = [
        (zero(), ri(3),  zero()),
        (zero(), ri(-4), zero()),
        
        (cst(ri(5)), ri(3),  cst(r(5, 3))),
        (cst(ri(5)), ri(-4), cst(r(-5, 4))),
        
        (cst(ri(-4)), ri(3),  cst(r(-4, 3))),
        (cst(ri(-4)), ri(-4), cst(ri(1))),
        
        (pr1(), ri(3), dense![r(4, 3), r(5, 3), ri(0), ri(2), r(7, 3)]),
        (pr1(), ri(-4), dense![ri(-1), r(-5, 4), ri(0), r(-3, 2), r(-7, 4)]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

fn suppress_stack_traces() {
    std::panic::set_hook(Box::new(|_|{}));
}

macro_rules! test_div_by_zero {
    ($test_name0:ident, $test_name1:ident, $arg:expr) => {
        #[test]
        #[should_panic(expected="Can't divide by zero!")]
        fn $test_name0(){
            suppress_stack_traces();
            let _ = $arg / Polynomial::Zero::<Rational64>();
        }

        #[test]
        #[should_panic(expected="Can't divide by zero!")]
        fn $test_name1(){
            suppress_stack_traces();
            let _ = $arg / Rational64::zero();
        }
    };
}

test_div_by_zero!(test_div_zero_by_poly_zero,   test_div_zero_by_zero,   Polynomial::<Rational64>::zero());
test_div_by_zero!(test_div_one_by_poly_zero,    test_div_one_by_zero,    Polynomial::<Rational64>::one());
test_div_by_zero!(test_div_const_by_poly_zero,  test_div_const_by_zero,  Polynomial::constant(ri(5)));
test_div_by_zero!(test_div_dense_by_poly_zero,  test_div_dense_by_zero,  dense![ri(1), ri(2), ri(3)]);
test_div_by_zero!(test_div_sparse_by_poly_zero, test_div_sparse_by_zero, sparse![(0, ri(1)), (1, ri(2)), (2, ri(3))]);

#[test]
fn test_pow(){

    fn test(x: Polynomial<i64>, p: u32, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, p: u32, exp: &'b Polynomial<i64>){
            assert_eq!(x.pow(p), *exp);
            assert_eq!(x.clone().pow(p), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, p, &exp);
        }
    }

    let table = [
        (zero(), 0, one()),  // ok? In rust, 0_u32.pow(0) == 1
        (zero(), 1, zero()),
        (zero(), 2, zero()),
        (zero(), 5, zero()),

        (one(), 0, one()),
        (one(), 1, one()),
        (one(), 2, one()),
        (one(), 5, one()),

        (cst(3), 0, one()),
        (cst(3), 1, cst(3)),
        (cst(3), 2, cst(9)),
        (cst(3), 5, cst(243)),

        (cst(-4), 0, one()),
        (cst(-4), 1, cst(-4)),
        (cst(-4), 2, cst(16)),
        (cst(-4), 5, cst(-1024)),

        (p2(), 0, one()),
        (p2(), 1, p2()),
        (p2(), 2, p2() * p2()),
        (p2(), 5, p2() * p2() * p2() * p2() * p2()),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_compose(){

    fn test(x: Polynomial<i64>, y: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b, 'c>(x: &'a Polynomial<i64>, y: &'b Polynomial<i64>, exp: &'c Polynomial<i64>){
            assert_eq!(x.compose(y), *exp);
            assert_eq!(x.compose(y.clone()), *exp);
            assert_eq!(x.clone().compose(y), *exp);
            assert_eq!(x.clone().compose(y.clone()), *exp);
        }

        for x_ in get_impls(&x) {
            for y_ in get_impls(&y) {
                test_op(&x_, &y_, &exp);
            }
        }
    }

    let table = [
        (zero(), zero(),           zero()),
        (zero(), one(),            zero()),
        (zero(), cst(3),           zero()),
        (zero(), Polynomial::x(),  zero()),
        (zero(), Polynomial::x2(), zero()),
        (zero(), p0(),             zero()),
        
        (one(), zero(),           one()),
        (one(), one(),            one()),
        (one(), cst(3),           one()),
        (one(), Polynomial::x(),  one()),
        (one(), Polynomial::x2(), one()),
        (one(), p0(),             one()),
        
        (cst(5), zero(),           cst(5)),
        (cst(5), one(),            cst(5)),
        (cst(5), cst(3),           cst(5)),
        (cst(5), Polynomial::x(),  cst(5)),
        (cst(5), Polynomial::x2(), cst(5)),
        (cst(5), p0(),             cst(5)),

        (Polynomial::x(), zero(),          zero()),
        (Polynomial::x(), one(),           one()),
        (Polynomial::x(), cst(3),          cst(3)),
        (Polynomial::x(), Polynomial::x(), Polynomial::x()),
        (Polynomial::x(), Polynomial::x2(), Polynomial::x2()),
        (Polynomial::x(), p0(),            p0()),
        (Polynomial::x(), p4(),            p4()),

        (Polynomial::x2(), zero(),           zero()),
        (Polynomial::x2(), one(),            one()),
        (Polynomial::x2(), cst(3),           cst(9)),
        (Polynomial::x2(), Polynomial::x(),  Polynomial::x2()),
        (Polynomial::x2(), Polynomial::x2(), Polynomial::x4()),
        (Polynomial::x2(), p0(),             p0().pow(2)),
        (Polynomial::x2(), p4(),             p4().pow(2)),
        
        (p1(), zero(),           cst(4)),
        (p1(), one(),            cst(4 + 5 + 6 + 7)),
        (p1(), cst(3),           cst(4 + 5 * 3 + 6 * 3*3*3 + 7 * 3*3*3*3)),
        (p1(), Polynomial::x(),  p1()),
        (p1(), Polynomial::x2(), dense![4, 0, 5, 0, 0, 0, 6, 0, 7]),
        (p1(), p0(),             4 + 5 * p0() + 6 * p0().pow(3) + 7 * p0().pow(4)),
        (p1(), p4(),             4 + 5 * p4() + 6 * p4().pow(3) + 7 * p4().pow(4)),

        (p4(), zero(),           zero()),
        (p4(), one(),            cst(6)),
        (p4(), cst(3),           cst(162)),
        (p4(), Polynomial::x(),  p4()),
        (p4(), Polynomial::x2(), 6_i64 * Polynomial::x().pow(6)),
        (p4(), p0(),             6_i64 * p0().pow(3)),
        (p4(), p4(),             6_i64 * p4().pow(3)),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

//********** Utility Methods **********/
#[test]
fn test_reductum(){

    fn test(x: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, exp: &'b Polynomial<i64>){
            assert_eq!(x.reductum(), *exp);
            assert_eq!(x.clone().reductum(), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, &exp);
        }
    }

    let table = [
        (zero(),           zero()),
        (cst(5),           zero()),
        (cst(-4),          zero()),
        (Polynomial::x(),  zero()),
        (Polynomial::x2(), zero()),
        (p0(),             dense![1, 2]),
        (p1(),             dense![4, 5, 0, 6]),
        (p2(),             dense![4, 0, 0, 5]),
        (p3(),             dense![1]),
        (p4(),             zero()),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_remove_zero_roots(){

    fn test(x: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, exp: &'b Polynomial<i64>){
            assert_eq!(x.remove_zero_roots(), *exp);
            assert_eq!(x.clone().remove_zero_roots(), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, &exp);
        }
    }

    let table = [
        (zero(),           zero()),
        (cst(5),           cst(5)),
        (cst(-4),          cst(-4)),
        (Polynomial::x(),  one()),
        (Polynomial::x2(), one()),
        (p0(),             dense![1, 2, 3]),
        (p1(),             dense![4, 5, 0, 6, 7]),
        (p2(),             dense![4, 0, 0, 5, 0, 0, 0, 6]),
        (p3(),             dense![1, 0, 0, 0, 2]),
        (p4(),             cst(6)),
        (dense![0, 0, 0, 1, 2, 3], dense![1, 2, 3])
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_derivative(){

    fn test(x: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, exp: &'b Polynomial<i64>){
            assert_eq!(x.derivative(), *exp);
            assert_eq!(x.clone().derivative(), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, &exp);
        }
    }

    let table = [
        (zero(),           zero()),
        (cst(5),           zero()),
        (cst(-4),          zero()),
        (Polynomial::x(),  cst(1)),
        (Polynomial::x2(), dense![0, 2]),
        (p0(),             dense![2, 6]),
        (p1(),             dense![5, 0, 18, 28]),
        (p2(),             dense![0, 0, 15, 0, 0, 0, 42]),
        (p3(),             dense![0, 0, 0, 8]),
        (p4(),             dense![0, 0, 18]),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_integral(){

    fn test(x: PolyR64, exp: PolyR64){

        fn test_op<'a, 'b>(x: &'a PolyR64, exp: &'b PolyR64){
            assert_eq!(x.integral(), *exp);
            assert_eq!(x.clone().integral(), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, &exp);
        }
    }

    let table = [
        (zero(),           zero()),
        (cst(ri(5)),       dense![ri(0), ri(5)]),
        (cst(ri(-4)),      dense![ri(0), ri(-4)]),
        (Polynomial::x(),  dense![ri(0), ri(0), r(1, 2)]),
        (Polynomial::x2(), dense![ri(0), ri(0), ri(0), r(1, 3)]),
        (pr0(),            dense![ri(0), ri(1), ri(1), ri(1)]),
        (pr1(),            dense![ri(0), ri(4), r(5, 2), ri(0), r(3, 2), r(7, 5)]),
        (pr2(),            dense![ri(0), ri(4), ri(0), ri(0), r(5, 4), ri(0), ri(0), ri(0), r(3, 4)]),
        (pr3(),            dense![ri(0), ri(1), ri(0), ri(0), ri(0), r(2, 5)]),
        (pr4(),            dense![ri(0), ri(0), ri(0), ri(0), r(3, 2)]),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_reciprocal(){

    fn test(x: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, exp: &'b Polynomial<i64>){
            assert_eq!(x.reciprocal(), *exp);
            assert_eq!(x.clone().reciprocal(), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, &exp);
        }
    }

    let table = [
        (zero(),           zero()),
        (cst(5),           cst(5)),
        (cst(-4),          cst(-4)),
        (Polynomial::x(),  cst(1)),
        (Polynomial::x2(), cst(1)),
        (p0(),             dense![3, 2, 1]),
        (p1(),             dense![7, 6, 0, 5, 4]),
        (p2(),             dense![6, 0, 0, 0, 5, 0, 0, 4]),
        (p3(),             dense![2, 0, 0, 0, 1]),
        (p4(),             cst(6)),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_flip(){

    fn test(x: Polynomial<i64>, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, exp: &'b Polynomial<i64>){
            assert_eq!(x.flip(), *exp);
            assert_eq!(x.clone().flip(), *exp);

            // a flipped polynomial of p(x) equals p(-x) 
            let exp2 = x.compose(&-Polynomial::x());
            assert_eq!(x.flip(), exp2.clone());
            assert_eq!(x.clone().flip(), exp2);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, &exp);
        }
    }

    let table = [
        (zero(),           zero()),
        (cst(5),           cst(5)),
        (cst(-4),          cst(-4)),
        (Polynomial::x(),  -Polynomial::x()),
        (Polynomial::x2(), Polynomial::x2()),
        (p0(),             dense![1, -2, 3]),
        (p1(),             dense![4, -5, 0, -6, 7]),
        (p2(),             dense![4, 0, 0, -5, 0, 0, 0, -6]),
        (p3(),             dense![1, 0, 0, 0, 2]),
        (p4(),             dense![0, 0, 0, -6]),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_shift(){

    fn test(x: Polynomial<i64>, h: i64, exp: Polynomial<i64>){

        fn test_op<'a, 'b>(x: &'a Polynomial<i64>, h: i64, exp: &'b Polynomial<i64>){
            assert_eq!(x.shift(h), *exp);
            assert_eq!(x.clone().shift(h), *exp);

            // a shifted polynomial of p(x) by h equals p(x + h) 
            let exp2 = x.compose(Polynomial::x() + h);
            assert_eq!(x.shift(h), exp2.clone());
            assert_eq!(x.clone().shift(h), exp2);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, h, &exp);
        }
    }

    let table = [
        (zero(), -2, zero()),
        (zero(),  0, zero()),
        (zero(),  2, zero()),

        (cst(5), -2, cst(5)),
        (cst(5),  0, cst(5)),
        (cst(5),  2, cst(5)),

        (Polynomial::x(), -2, dense![-2, 1]),
        (Polynomial::x(),  0, dense![0, 1]),
        (Polynomial::x(),  2, dense![2, 1]),

        (Polynomial::x2(), -2, dense![4, -4, 1]),
        (Polynomial::x2(), 0, dense![0, 0, 1]),
        (Polynomial::x2(), 2, dense![4, 4, 1]),

        (p0(), -2, dense![9, -10, 3]),
        (p0(),  0, dense![1, 2, 3]),
        (p0(),  2, dense![17, 14, 3]),

        (p1(), -2, dense![58, -147, 132, -50, 7]),
        (p1(),  0, dense![4, 5, 0, 6, 7]),
        (p1(),  2, dense![174, 301, 204, 62, 7]),
        (p2(), -2, dense![-804, 2748, -4062, 3365, -1680, 504, -84, 6]),
        (p2(),  0, dense![4, 0, 0, 5, 0, 0, 0, 6]),
        (p2(),  2, dense![812, 2748, 4062, 3365, 1680, 504, 84, 6]),

        (p3(), -2, dense![33, -64, 48, -16, 2]),
        (p3(),  0, dense![1, 0, 0, 0, 2]),
        (p3(),  2, dense![33, 64, 48, 16, 2]),

        (p4(), -2, dense![-48, 72, -36, 6]),
        (p4(),  0, dense![0, 0, 0, 6]),
        (p4(),  2, dense![48, 72, 36, 6]),

        (dense![1], -1, dense![1]),
        (dense![0, 1], -1, dense![-1, 1]),
        (dense![0, 0, 1], -1, dense![1, -2, 1]),
        (dense![0, 0, 0, 1], -1, dense![-1, 3, -3, 1]),
        (dense![0, 0, 0, 0, 1], -1, dense![1, -4, 6, -4, 1]),
        (dense![0, 0, 0, 0, 0, 1], -1, dense![-1, 5, -10, 10, -5, 1]),
        (dense![0, 0, 0, 0, 0, 0, 1], -1, dense![1, -6, 15, -20, 15, -6, 1]),
        (dense![0, 0, 0, 0, 0, 0, 0, 1], -1, dense![-1, 7, -21, 35, -35, 21, -7, 1]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_monic(){

    fn test(x: PolyR64, exp: PolyR64){

        fn test_op<'a, 'b>(x: &'a PolyR64, exp: &'b PolyR64){
            assert_eq!(x.monic(), *exp);
            assert_eq!(x.clone().monic(), *exp);
        }

        for x_ in get_impls(&x) {
            test_op(&x_, &exp);
        }
    }

    let table = [
        (zero(),           zero()),
        (cst(ri(5)),       one()),
        (cst(r(5, 3)),     one()),
        (Polynomial::x(),  Polynomial::x()),
        (Polynomial::x2(), Polynomial::x2()),
        (pr0(),            dense![r(1, 3), r(2, 3), ri(1)]),
        (pr1(),            dense![r(4, 7), r(5, 7), ri(0), r(6, 7), ri(1)]),
        (pr2(),            dense![r(2, 3), ri(0), ri(0), r(5, 6), ri(0), ri(0), ri(0), ri(1)]),
        (pr3(),            dense![r(1, 2), ri(0), ri(0), ri(0), ri(1)]),
        (pr4(),            dense![ri(0), ri(0), ri(0), ri(1)]),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}