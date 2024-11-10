use num::{Complex, One, Rational64, Zero};

use crate::{dense, spears, polynomial::Polynomial};

// ⁰¹²³⁴⁵⁶⁷⁸⁹
#[test]
fn test_dense_macro(){
    let p0 = dense![i64; 1, 2, 3];
    assert_eq!(format!("{}", p0), "1 + 2x + 3x²");

    let p1 = dense![i64; 0, 4, 5, 6, 0, 7, 0, 0];
    assert_eq!(format!("{}", p1), "4x + 5x² + 6x³ + 7x⁵");
}

#[test]
fn test_spears_macro(){
    let p0 = spears![i64; (2, 6), (5, 1)];
    assert_eq!(format!("{}", p0), "6x² + x⁵");
}

#[test]
fn test_display(){
    fn test(p: Polynomial<i64>, exp: &str){
        assert_eq!(format!("{}", p), exp);
    }

    let table = [
        (Polynomial::<i64>::zero(), "(0)"),
        (Polynomial::<i64>::one(), "(1)"),

        (Polynomial::<i64>::constant(0),  "(0)"),
        (Polynomial::<i64>::constant(1),  "(1)"),
        (Polynomial::<i64>::constant(-1), "(-1)"),
        (Polynomial::<i64>::constant(2),  "(2)"),

        // dense!
        (dense![i64;],      "(0)"),
        (dense![i64; 0, 0], "(0)"),
        (dense![i64; 1],    "(1)"),
        (dense![i64; -1],   "(-1)"),
        (dense![i64; 2],    "(2)"),

        (dense![i64; 1, 2, 3],        "1 + 2x + 3x²"),
        (dense![i64; -1, -2, -3],     "-1 - 2x - 3x²"),
        (dense![i64; 0, 2, 3],        "2x + 3x²"),
        (dense![i64; 1, 2, 3, 0, 0],  "1 + 2x + 3x²"),
        (dense![i64; 1, 1, 1, 1],     "1 + x + x² + x³"),
        (dense![i64; -1, -1, -1, -1], "-1 - x - x² - x³"),

        // spears!
        (spears![i64;], "(0)"),
        (spears![i64; (0, 0), (1, 0)], "(0)"),
        (spears![i64; (0, 1)],  "(1)"),
        (spears![i64; (0, -1)], "(-1)"),
        (spears![i64; (0, 2)],   "(2)"),

        (spears![i64; (0, 1), (2, 3)],   "1 + 3x²"),
        (spears![i64; (0, -1), (2, -3)], "-1 - 3x²"),
        (spears![i64; (1, 2), (2, 3)],   "2x + 3x²"),
        (spears![i64; (0, 1), (2, 3), (4, 0)], "1 + 3x²"),
        (spears![i64; (0, 1), (1, 1), (2, 1), (3, 1)],     "1 + x + x² + x³"),
        (spears![i64; (0, -1), (1, -1), (2, -1), (3, -1)], "-1 - x - x² - x³"),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_display_for_rational_and_complex(){
    // rational
    fn r(n: i64, d: i64) -> Rational64 { Rational64::new(n, d) }
    let p_rat = dense![Rational64; r(1, 2), r(1, 3), r(2, 3)];
    assert_eq!(format!("{}", p_rat), "1/2 + (1/3)x + (2/3)x²");

    // complex
    fn c(re: f64, im: f64) -> Complex<f64> { Complex { re, im } }
    let p_complex = dense!(Complex<f64>; c(1., 2.), c(3., 4.), c(5., 6.));
    assert_eq!(format!("{}", p_complex), "1+2i + (3+4i)x + (5+6i)x²");
}

#[test]
fn test_debug(){
    fn test(p: Polynomial<i64>, exp: &str){
        assert_eq!(format!("{:?}", p), exp);
    }

    let table = [
        (Polynomial::<i64>::zero(), "Polynomial::Zero[(0)]"),

        (Polynomial::<i64>::one(), "Polynomial::Constant[(1)]"),
        (Polynomial::<i64>::constant(2),  "Polynomial::Constant[(2)]"),

        (dense![i64; 1, 2, 3],  "Polynomial::Dense[1 + 2x + 3x²]"),
        (spears![i64; (0, 1), (1, 2), (2, 3)],   "Polynomial::Spears[1 + 2x + 3x²]"),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_degree_zero_one_constant(){
    fn test(p: Polynomial<i64>, exp_deg: usize, exp_zero: bool, exp_one: bool, exp_const: bool){
        assert_eq!(p.degree(), exp_deg);
        assert_eq!(p.is_zero(), exp_zero);
        assert_eq!(p.is_one(), exp_one);
        assert_eq!(p.is_constant(), exp_const);
    }

    let table = [                     // deg zero   one    const
        (Polynomial::<i64>::zero(),       0, true,  false, true),
        (Polynomial::<i64>::one(),        0, false, true,  true),

        (Polynomial::<i64>::constant(0),  0, true,  false, true),
        (Polynomial::<i64>::constant(1),  0, false, true,  true),
        (Polynomial::<i64>::constant(-1), 0, false, false, true),
        (Polynomial::<i64>::constant(2),  0, false, false, true),

        // dense!
        (dense![i64; ],                   0, true,  false, true),  // zero & const
        (dense![i64; 0, 0],               0, true,  false, true),  // zero & const
        (dense![i64; 1],                  0, false, true,  true),  // const
        (dense![i64; -1],                 0, false, false, true),  // const
        (dense![i64; 2],                  0, false, false, true),  // const

        (dense![i64; 1, 2, 3],            2, false, false, false),
        (dense![i64; -1, -2, -3],         2, false, false, false),
        (dense![i64; 0, 2, 3],            2, false, false, false),
        (dense![i64; 1, 2, 3, 0, 0],      2, false, false, false),

        // spears!
        (spears![i64; ],                  0, true,  false, true),  // zero & const
        (spears![i64; (0, 0), (2, 0)],    0, true,  false, true),  // zero & const
        (spears![i64; (0, 1)],            0, false, true,  true),  // const
        (spears![i64; (0, -1)],           0, false, false, true),  // const
        (spears![i64; (0, 2)],            0, false, false, true),  // const

        (spears![i64; (0, 1), (1, 2), (2, 3)],                  2, false, false, false),
        (spears![i64; (0, -1), (1, -2), (2, -3)],               2, false, false, false),
        (spears![i64; (0, 0), (1, 2), (2, 3)],                  2, false, false, false),
        (spears![i64; (0, 1), (1, 2), (2, 3), (5, 0), (10, 0)], 2, false, false, false),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2, entry.3, entry.4);
    }
}

#[test]
fn test_eq(){
    let p0 = dense![i64; 1, 0, 2, 4, 0, 0, 1];
    let p1 = spears![i64; (0, 1), (2, 2), (3, 4), (6, 1)];
    assert_eq!(p0, p1);
}

#[test]
fn test_clone_category_methods(){  // clone(), dense_clone(), spears_clone(), to_dense() and to_spears()
    // 1 + 2x² + 4x³ + x⁶
    let pd = dense![i64; 1, 0, 2, 4, 0, 0, 1];
    let ps = spears![i64; (0, 1), (2, 2), (3, 4), (6, 1)];
    assert_eq!(pd, ps);

    // dense_clone()
    let p1 = pd.dense_clone();
    assert_eq!(p1, pd);
    assert_eq!(p1, ps);

    let p2 = ps.dense_clone();
    assert_eq!(p2, pd);
    assert_eq!(p2, ps);

    // spears_clone()
    let p3 = pd.spears_clone();
    assert_eq!(p3, pd);
    assert_eq!(p3, ps);

    let p4 = ps.spears_clone();
    assert_eq!(p4, pd);
    assert_eq!(p4, ps);

    // to_dense()
    let p5 = pd.clone().to_dense();
    assert_eq!(p5, pd);
    assert_eq!(p5, ps);

    let p6 = ps.clone().to_dense();
    assert_eq!(p6, pd);
    assert_eq!(p6, ps);

    // to_spears()
    let p7 = pd.clone().to_spears();
    assert_eq!(p7, pd);
    assert_eq!(p7, ps);

    let p8 = ps.clone().to_spears();
    assert_eq!(p8, pd);
    assert_eq!(p8, ps);
}