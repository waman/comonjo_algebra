use std::{collections::HashMap, vec};

use num::{complex::c64, pow::Pow, One, Rational64, Zero};

use crate::{algebra::Semiring, dense, poly::{CoeffsIterator, Polynomial}, sparse};

fn zero<C: Semiring>() -> Polynomial<C> { Polynomial::zero() }
fn one<C: Semiring + Clone>() -> Polynomial<C> { Polynomial::one() }
fn cst<C: Semiring>(v: C) -> Polynomial<C> { Polynomial::constant(v) }

fn r(n: i64, d: i64) -> Rational64 { Rational64::new(n, d) }
fn ri(n: i64) -> Rational64 { Rational64::new(n, 1) }

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
        (sparse![(0, 1), (1, 2), (2, 3)],   "Polynomial::<i64>::Sparse[1 + 2x + 3x²]"),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_degree_and_is_xxx_methods(){
    fn test(p: Polynomial<i64>, exp_deg: usize, exp_zero: bool, exp_one: bool, exp_const: bool){
        assert_eq!(p.degree(), exp_deg);
        assert_eq!(p.is_zero(), exp_zero);
        assert_eq!(p.is_one(), exp_one);
        assert_eq!(p.is_constant(), exp_const);
    }

    let table = [ // deg, zero, one, const
        (zero(), 0, true,  false, true),
        (one(), 0, false, true,  true),

        (cst(0),  0, true,  false, true),
        (cst(1),  0, false, true,  true),
        (cst(-1), 0, false, false, true),
        (cst(2),  0, false, false, true),


        (dense![],     0, true,  false, true),  // zero & const
        (dense![0, 0], 0, true,  false, true),  // zero & const
        (dense![1],    0, false, true,  true),  // const
        (dense![-1],   0, false, false, true),  // const
        (dense![2],    0, false, false, true),  // const

        (dense![1, 2, 3],       2, false, false, false),
        (dense![-1, -2, -3],    2, false, false, false),
        (dense![0, 2, 3],       2, false, false, false),
        (dense![1, 2, 3, 0, 0], 2, false, false, false),


        (sparse![],               0, true,  false, true),  // zero & const
        (sparse![(0, 0), (2, 0)], 0, true,  false, true),  // zero & const
        (sparse![(0, 1)],         0, false, true,  true),  // const
        (sparse![(0, -1)],        0, false, false, true),  // const
        (sparse![(0, 2)],         0, false, false, true),  // const

        (sparse![(2, 3)],                                  2, false, false, false),
        (sparse![(0, 1), (1, 2), (2, 3)],                  2, false, false, false),
        (sparse![(0, -1), (1, -2), (2, -3)],               2, false, false, false),
        (sparse![(0, 0), (1, 2), (2, 3)],                  2, false, false, false),
        (sparse![(0, 1), (1, 2), (2, 3), (5, 0), (10, 0)], 2, false, false, false),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2, entry.3, entry.4);
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

/// 4x³ (= dense![[0, 0, 0, 4]])
fn p4() -> Polynomial<i64> { dense![0, 0, 0, 4] }

fn get_impls<'a, C>(x: &'a Polynomial<C>) -> Vec<Polynomial<C>> where C: Semiring + Clone {
    if x.degree() == 0 {
        vec![x.clone()]
    } else {
        vec![x.dense_clone(), x.sparse_clone()]
    }
}

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

        (p0(), dense![-1, -2, -3]),
        (p1(), dense![-4, -5, 0, -6, -7]),
        (p2(), dense![-4, 0, 0, -5, 0, 0, 0, -6]),
        (p3(), dense![-1, 0, 0, 0, -2]),
        (p4(), dense![0, 0, 0, -4]),
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
            // assert_eq!(x + y.clone(), *exp);
            // assert_eq!(x.clone() + y.clone(), *exp);
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
        (zero(), one(), one()),
        (zero(), cst(3), cst(3)),
        (zero(), p0(), dense![1, 2, 3]),
        (zero(), p1(), dense![4, 5, 0, 6, 7]),
        (zero(), p2(), dense![4, 0, 0, 5, 0, 0, 0, 6]),
        (zero(), p3(), dense![1, 0, 0, 0, 2]),
        (zero(), p4(), dense![0, 0, 0, 4]),
        
        (one(), one(), cst(2)),
        (one(), cst(3), cst(4)),
        (one(), p0(), dense![2, 2, 3]),
        (one(), p1(), dense![5, 5, 0, 6, 7]),
        (one(), p2(), dense![5, 0, 0, 5, 0, 0, 0, 6]),
        (one(), p3(), dense![2, 0, 0, 0, 2]),
        (one(), p4(), dense![1, 0, 0, 4]),
        
        (cst(5), cst(3), cst(8)),
        (cst(5), cst(-5), zero()),  // result be zero
        (cst(5), p0(), dense![6, 2, 3]),
        (cst(5), p1(), dense![9, 5, 0, 6, 7]),
        (cst(5), p2(), dense![9, 0, 0, 5, 0, 0, 0, 6]),
        (cst(5), p3(), dense![6, 0, 0, 0, 2]),
        (cst(5), p4(), dense![5, 0, 0, 4]),
        
        (cst(-4), cst(3), cst(-1)),
        (cst(-4), p0(), dense![-3, 2, 3]),
        (cst(-4), p1(), dense![0, 5, 0, 6, 7]),
        (cst(-4), p2(), dense![0, 0, 0, 5, 0, 0, 0, 6]),
        (cst(-4), p3(), dense![-3, 0, 0, 0, 2]),
        (cst(-4), p4(), dense![-4, 0, 0, 4]),
        
        (p1(), -p1(), zero()),                    // result be zero
        (p1(), dense![-1, -5, 0, -6, -7], cst(3)), // result be const
        (p0(), p1(), dense![5, 7, 3, 6, 7]),
        (p1(), p2(), dense![8, 5, 0, 11, 7, 0, 0, 6]),
        (p3(), p4(), dense![1, 0, 0, 4, 2]),
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
        (zero(), one(), cst(-1)),
        (zero(), cst(3), cst(-3)),
        (zero(), p0(), dense![-1, -2, -3]),
        (zero(), p1(), dense![-4, -5, 0, -6, -7]),
        (zero(), p2(), dense![-4, 0, 0, -5, 0, 0, 0, -6]),
        (zero(), p3(), dense![-1, 0, 0, 0, -2]),
        (zero(), p4(), dense![0, 0, 0, -4]),
        
        (one(), one(), zero()),
        (one(), cst(3), cst(-2)),
        (one(), p0(), dense![0, -2, -3]),
        (one(), p1(), dense![-3, -5, 0, -6, -7]),
        (one(), p2(), dense![-3, 0, 0, -5, 0, 0, 0, -6]),
        (one(), p3(), dense![0, 0, 0, 0, -2]),
        (one(), p4(), dense![1, 0, 0, -4]),
        
        (cst(5), cst(3), cst(2)),
        (cst(5), cst(5), zero()),  // result be zero
        (cst(5), p0(), dense![4, -2, -3]),
        (cst(5), p1(), dense![1, -5, 0, -6, -7]),
        (cst(5), p2(), dense![1, 0, 0, -5, 0, 0, 0, -6]),
        (cst(5), p3(), dense![4, 0, 0, 0, -2]),
        (cst(5), p4(), dense![5, 0, 0, -4]),
        
        (cst(4), p1(), dense![0, -5, 0, -6, -7]),
        (cst(4), p2(), dense![0, 0, 0, -5, 0, 0, 0, -6]),
        
        (p0(), p0(), zero()),             // result be zero
        (p0(), dense![-2, 2, 3], cst(3)), // result be const
        (p0(), p1(), dense![-3, -3, 3, -6, -7]),
        (p1(), p2(), dense![0, 5, 0, 1, 7, 0, 0, -6]),
        (p3(), p4(), dense![1, 0, 0, -4, 2]),
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

    fn exp_0_mul_1() -> Polynomial<i64> { dense![4, 13, 22, 21, 19, 32, 21] }
    fn exp_2_mul_3() -> Polynomial<i64> { sparse![(0, 4), (3, 5), (4, 8), (7, 16), (11, 12)] }

    let table = [
        (zero(), zero(), zero()),
        (zero(), one(),  zero()),
        (zero(), cst(3), zero()),
        (zero(), p0(),   zero()),
        
        (one(), one(),  one()),
        (one(), cst(3), cst(3)),
        (one(), p0(),   p0()),
        
        (cst(5), cst(3),      cst(15)),
        (cst(5), p0(),  dense![5, 10, 15]),
        (cst(5), p2(), sparse![(0, 20), (3, 25), (7, 30)]),
        
        (p1(), zero(),      zero()),
        (p1(), one(),       p1()),
        (p1(), cst(3),      dense![12, 15, 0, 18, 21]),
        (p0(), p1(),  exp_0_mul_1()),
        (p2(), p3(),  exp_2_mul_3()),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_div_rem(){
    type PolyR64 = Polynomial<Rational64>;

    fn cst_r(n: i64, d: i64) -> PolyR64 { Polynomial::constant(r(n, d)) }
    fn cst_i(i: i64) -> PolyR64 { Polynomial::constant(ri(i)) }

    fn to_poly_r64(p: Polynomial<i64>) -> PolyR64 { p.map(|_, c| ri(c)) }
    
    /// 1 + 2x + 3x² (= dense![[1, 2, 3]])
    fn pr0() -> PolyR64 { to_poly_r64(p0()) }

    /// 4 + 5x + 6x³ + 7x⁴ (= dense![[4, 5, 0, 6, 7]])
    fn pr1() -> PolyR64 { to_poly_r64(p1()) }

    /// 4 + 5x³ + 6x⁷ (= dense![[4, 0, 0, 5, 0, 0, 0, 6]])
    fn pr2() -> PolyR64 { to_poly_r64(p2()) }

    /// 1 + 2x⁴ (= dense![[1, 0, 0, 0, 2]])
    fn pr3() -> PolyR64 { to_poly_r64(p3()) }

    /// 4x³ (= dense![[0, 0, 0, 4]])
    fn pr4() -> PolyR64 { to_poly_r64(p4()) }

    fn pr1_div_3() -> PolyR64 { dense![r(4, 3), r(5, 3), ri(0), r(6, 3), r(7, 3)] }

    fn pr1_div_pr0() -> PolyR64 { dense![r(-29, 27), r(4, 9), r(7, 3)] }
    fn pr1_rem_pr0() -> PolyR64 { dense![r(137, 27), r(181, 27)] }

    fn pr2_div_pr3() -> PolyR64 { to_poly_r64(dense![0, 0, 0, 3]) }
    fn pr2_rem_pr3() -> PolyR64 { to_poly_r64(dense![4, 0, 0, 2]) }
    
    fn pr2_div_pr4() -> PolyR64 { dense![r(5, 4), ri(0), ri(0), ri(0), r(3, 2)] }
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
        (zero(), one(),    zero(), zero()),
        (zero(), cst_i(3), zero(), zero()),
        (zero(), pr0(),    zero(), zero()),
           
        (one(), one(),    one(),       zero()),
        (one(), cst_i(3), cst_r(1, 3), zero()),
        (one(), pr0(),    zero(),      one()),
        
        (cst_i(5), one(),    cst_i(5),    zero()),  
        (cst_i(5), cst_i(3), cst_r(5, 3), zero()),
        (cst_i(5), pr0(),    zero(),      cst_i(5)),
        
        (pr1(),           one(),    pr1(),         zero()),
        (pr1(),           cst_i(3), pr1_div_3(),   zero()),
        (pr1(),           pr0(),    pr1_div_pr0(), pr1_rem_pr0()),
        (pr2() * pr3(),   pr3(),    pr2(),         zero()),
        (pr2() * pr3(),   pr2(),    pr3(),         zero()),
        (pr2(),           pr3(),    pr2_div_pr3(), pr2_rem_pr3()),
        (pr2()        ,   pr4(),    pr2_div_pr4(), pr2_rem_pr4()),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2, entry.3);
    }
}

macro_rules! test_divide_by_zero {
    ($test_name:ident, $arg:expr) => {
        #[test]
        #[allow(unused_must_use)]
        #[should_panic(expected="Can't divide by polynomial of zero!")]
        fn $test_name(){
            $arg / Polynomial::Zero::<Rational64>();
        }
    };
}

test_divide_by_zero!(test_divide_zero_by_zero, Polynomial::<Rational64>::zero());
test_divide_by_zero!(test_divide_one_by_zero, Polynomial::<Rational64>::one());
test_divide_by_zero!(test_divide_const_by_zero, Polynomial::constant(ri(5)));
test_divide_by_zero!(test_divide_dense_by_zero, dense![ri(1), ri(2), ri(3)]);
test_divide_by_zero!(test_divide_sparse_by_zero, sparse![(0, ri(1)), (1, ri(2)), (2, ri(3))]);

#[test]
fn test_pow(){
    fn test(x: Polynomial<i64>, p: u32, exp: Polynomial<i64>){
        let ref_x = &x;
        assert_eq!(ref_x.pow(p), exp);
        assert_eq!(x.pow(p), exp);
    }

    let p = dense![4, 5, 6];
    let p2 = p.clone() * p.clone();
    let p5 = p2.clone() * p2.clone() * p.clone();

    let table = [
        (zero(),                          0, one()),  // ok?
        (one(),                           0, one()),
        (cst(3),                          0, one()),
        (dense![4, 5, 6],                 0, one()),
        (sparse![(0, 4), (1, 5), (2, 6)], 0, one()),

        (zero(),                          1, zero()),
        (one(),                           1, one()),
        (cst(3),                          1, cst(3)),
        (dense![4, 5, 6],                 1, dense![4, 5, 6]),
        (sparse![(0, 4), (1, 5), (2, 6)], 1, dense![4, 5, 6]),

        (zero(),                          2, zero()),
        (one(),                           2, one()),
        (cst(3),                          2, cst(9)),
        (dense![4, 5, 6],                 2, dense![16, 40, 73, 60, 36]),
        (sparse![(0, 4), (1, 5), (2, 6)], 2, dense![16, 40, 73, 60, 36]),

        (zero(),                          5, zero()),
        (one(),                           5, one()),
        (cst(3),                          5, cst(243)),
        (dense![4, 5, 6],                 5, p5.clone()),
        (sparse![(0, 4), (1, 5), (2, 6)], 5, p5),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

// #[test]
// fn test_neg_performance(){

//     fn gen_vec(order: usize, rng: &mut ThreadRng) -> Vec<i64> {
//         let mut v: Vec<i64> = Vec::with_capacity(order + 1);
//         for _ in 0..=order {
//             v.push(rng.gen());
//         }
//         v
//     }

//     let mut rng = rand::thread_rng();
//     const N: usize = 1000;

//     let vv0: Vec<Vec<i64>> = (0..N).map(move|_| gen_vec(20, &mut rng)).collect();
//     let ps0: Vec<Polynomial<i64>> = vv0.iter().map(|v|Polynomial::dense_from_vec(v.clone())).collect();
//     let ps1: Vec<Polynomial<i64>> = vv0.into_iter().map(|v| Polynomial::dense_from_vec(v.into_iter().map(|e|-e).collect())).collect();
//     let mut ps2 = Vec::with_capacity(N);

//     let now = Instant::now();
//     for p in ps0 {
//         ps2.push(-p)
//     }
//     let elapsed_time = now.elapsed();

//     ps1.iter().zip(ps2.iter()).for_each(|(p1, p2)| assert_eq!(p1, p2));
//     println!("Running slow_function() took {} microseconds.", elapsed_time.as_micros());
// }