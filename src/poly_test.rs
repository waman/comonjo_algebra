use std::collections::HashMap;

use num::{complex::c64, pow::Pow, One, Rational64, Zero};

use crate::{dense, poly::{iter::IntoCoeffsIterator, Polynomial}, sparse};

fn zero() -> Polynomial<i64> { Polynomial::<i64>::zero() }
fn one() -> Polynomial<i64> { Polynomial::<i64>::one() }
fn cst(v: i64) -> Polynomial<i64> { Polynomial::<i64>::constant(v) }

// ⁰¹²³⁴⁵⁶⁷⁸⁹
#[test]
fn test_dense_macro(){
    let p0 = dense![1, 2, 3];
    assert_eq!(format!("{}", p0), "1 + 2x + 3x²");

    let p1 = dense![0, 4, 5, 6, 0, 7, 0, 0];
    assert_eq!(format!("{}", p1), "4x + 5x² + 6x³ + 7x⁵");
}

#[test]
fn test_sparse_macro(){
    let p0 = sparse![(2, 6), (5, 1)];
    assert_eq!(format!("{}", p0), "6x² + x⁵");
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

        // dense!
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

        // sparse!
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
    fn r(n: i64, d: i64) -> Rational64 { Rational64::new(n, d) }
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
fn test_degree_zero_one_constant(){
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

        // dense!
        (dense![],     0, true,  false, true),  // zero & const
        (dense![0, 0], 0, true,  false, true),  // zero & const
        (dense![1],    0, false, true,  true),  // const
        (dense![-1],   0, false, false, true),  // const
        (dense![2],    0, false, false, true),  // const

        (dense![1, 2, 3],       2, false, false, false),
        (dense![-1, -2, -3],    2, false, false, false),
        (dense![0, 2, 3],       2, false, false, false),
        (dense![1, 2, 3, 0, 0], 2, false, false, false),

        // sparse!
        (sparse![],               0, true,  false, true),  // zero & const
        (sparse![(0, 0), (2, 0)], 0, true,  false, true),  // zero & const
        (sparse![(0, 1)],         0, false, true,  true),  // const
        (sparse![(0, -1)],        0, false, false, true),  // const
        (sparse![(0, 2)],         0, false, false, true),  // const

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

        // into_iter() for &Polynomial
        let cs3: HashMap<usize, i64> = (&p).into_iter().map(|(e, c)|(e, *c)).collect();
        assert_eq!(cs3, nonzero_exp);

        // nonzero_coeffs() for &Polynomial (= into_iter())
        let cs4: HashMap<usize, i64> = (&p).nonzero_coeffs().map(|(e, c)|(e, *c)).collect();
        assert_eq!(cs4, nonzero_exp);

        // coeffs() for &Polynomial
        let cs5: Vec<i64> = (&p).coeffs().map(|c| match c {
            Some(c) => *c,
            None => 0,
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
    let p0 = dense![1, 0, 2, 4, 0, 0, 1];
    assert_eq!(p0, dense![1, 0, 2, 4, 0, 0, 1]);

    let p1 = sparse![(0, 1), (2, 2), (3, 4), (6, 1)];
    assert_eq!(p1, sparse![(0, 1), (2, 2), (3, 4), (6, 1)]);

    assert_eq!(p0, p1);
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

#[test]
fn test_neg(){
    fn test(x: Polynomial<i64>, exp: Polynomial<i64>){
        assert_eq!(-&x, exp.clone());
        assert_eq!(-x, exp);
    }

    let table = [
        (zero(), zero()),
        (one(), cst(-1)),
        (cst(5), cst(-5)),
        (dense![4, 5, 0, 6, 7], dense![-4, -5, 0, -6, -7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], dense![-4, -5, 0, -6, -7]),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_add(){
    fn test(x: Polynomial<i64>, y: Polynomial<i64>, exp: Polynomial<i64>){
        assert_eq!(&x + &y, exp.clone());
        assert_eq!(x + y, exp);
    }

    let table = [
        (zero(), zero(),                          zero()),
        (zero(), one(),                           one()),
        (zero(), cst(3),                          cst(3)),
        (zero(), dense![1, 2, 3],                 dense![1, 2, 3]),
        (zero(), sparse![(0, 1), (2, 3), (4, 5)], dense![1, 0, 3, 0, 5]),
        
        (one(), zero(),                          one()),
        (one(), one(),                           cst(2)),
        (one(), cst(3),                          cst(4)),
        (one(), dense![1, 2, 3],                 dense![2, 2, 3]),
        (one(), sparse![(0, 1), (2, 3), (4, 5)], dense![2, 0, 3, 0, 5]),
        
        (cst(5), zero(),                          cst(5)),
        (cst(5), one(),                           cst(6)),
        (cst(5), cst(3),                          cst(8)),
        (cst(5), dense![1, 2, 3],                 dense![6, 2, 3]),
        (cst(5), sparse![(0, 1), (2, 3), (4, 5)], dense![6, 0, 3, 0, 5]),
        
        (dense![4, 5, 0, 6, 7], zero(),                          dense![4, 5, 0, 6, 7]),
        (dense![4, 5, 0, 6, 7], one(),                           dense![5, 5, 0, 6, 7]),
        (dense![4, 5, 0, 6, 7], cst(3),                          dense![7, 5, 0, 6, 7]),
        (dense![4, 5, 0, 6, 7], dense![1, 2, 3],                 dense![5, 7, 3, 6, 7]),
        (dense![4, 5, 0, 6, 7], sparse![(0, 1), (2, 3), (4, 5)], dense![5, 5, 3, 6, 12]),
        
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], zero(),                          dense![4, 5, 0, 6, 7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], one(),                           dense![5, 5, 0, 6, 7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], cst(3),                          dense![7, 5, 0, 6, 7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], dense![1, 2, 3],                 dense![5, 7, 3, 6, 7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], sparse![(0, 1), (2, 3), (4, 5)], dense![5, 5, 3, 6, 12]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_sub(){
    fn test(x: Polynomial<i64>, y: Polynomial<i64>, exp: Polynomial<i64>){
        assert_eq!(&x - &y, exp.clone());
        assert_eq!(x - y, exp);
    }

    let table = [
        (zero(), zero(),                          zero()),
        (zero(), one(),                           cst(-1)),
        (zero(), cst(3),                          cst(-3)),
        (zero(), dense![1, 2, 3],                 dense![-1, -2, -3]),
        (zero(), sparse![(0, 1), (2, 3), (4, 5)], dense![-1, 0, -3, 0, -5]),
        
        (one(), zero(),                          one()),
        (one(), one(),                           zero()),
        (one(), cst(3),                          cst(-2)),
        (one(), dense![1, 2, 3],                 dense![0, -2, -3]),
        (one(), sparse![(0, 1), (2, 3), (4, 5)], dense![0, 0, -3, 0, -5]),
        
        (cst(5), zero(),                          cst(5)),
        (cst(5), one(),                           cst(4)),
        (cst(5), cst(3),                          cst(2)),
        (cst(5), dense![1, 2, 3],                 dense![4, -2, -3]),
        (cst(5), sparse![(0, 1), (2, 3), (4, 5)], dense![4, 0, -3, 0, -5]),
        
        (dense![4, 5, 0, 6, 7], zero(),                          dense![4, 5, 0, 6, 7]),
        (dense![4, 5, 0, 6, 7], one(),                           dense![3, 5, 0, 6, 7]),
        (dense![4, 5, 0, 6, 7], cst(3),                          dense![1, 5, 0, 6, 7]),
        (dense![4, 5, 0, 6, 7], dense![1, 2, 3],                 dense![3, 3, -3, 6, 7]),
        (dense![4, 5, 0, 6, 7], sparse![(0, 1), (2, 3), (4, 5)], dense![3, 5, -3, 6, 2]),
        
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], zero(),                          dense![4, 5, 0, 6, 7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], one(),                           dense![3, 5, 0, 6, 7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], cst(3),                          dense![1, 5, 0, 6, 7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], dense![1, 2, 3],                 dense![3, 3, 3, 6, 7]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], sparse![(0, 1), (2, 3), (4, 5)], dense![3, 5, 3, 6, 2]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_mul(){
    fn test(x: Polynomial<i64>, y: Polynomial<i64>, exp: Polynomial<i64>){
        assert_eq!(&x * &y, exp.clone());
        assert_eq!(x * y, exp);
    }

    let table = [
        (zero(), zero(),                          zero()),
        (zero(), one(),                           zero()),
        (zero(), cst(3),                          zero()),
        (zero(), dense![1, 2, 3],                 zero()),
        (zero(), sparse![(0, 1), (2, 3), (4, 5)], zero()),
        
        (one(), zero(),                          zero()),
        (one(), one(),                           one()),
        (one(), cst(3),                          cst(3)),
        (one(), dense![1, 2, 3],                 dense![1, 2, 3]),
        (one(), sparse![(0, 1), (2, 3), (4, 5)], sparse![(0, 1), (2, 3), (4, 5)]),
        
        (cst(5), zero(),                          zero()),
        (cst(5), one(),                           cst(5)),
        (cst(5), cst(3),                          cst(15)),
        (cst(5), dense![1, 2, 3],                 dense![5, 10, 15]),
        (cst(5), sparse![(0, 1), (2, 3), (4, 5)], sparse![(0, 5), (2, 15), (4, 25)]),
        
        (dense![4, 5, 0, 6, 7], zero(),          zero()),
        (dense![4, 5, 0, 6, 7], one(),           dense![4, 5, 0, 6, 7]),
        (dense![4, 5, 0, 6, 7], cst(3),          dense![12, 15, 0, 18, 21]),
        (dense![4, 5, 0, 6, 7], dense![1, 2, 3], dense![4, 13, 22, 21, 19, 32, 21]),
        (sparse![(0, 4), (3, 5), (7, 6)].to_dense(), sparse![(0, 1), (4, 2)].to_dense(),
                sparse![(0, 4), (3, 5), (4, 8), (7, 16), (11, 12)]),
        (dense![4, 5, 0, 6, 7], dense![1, 2, 3].to_sparse(), dense![4, 13, 22, 21, 19, 32, 21]),
        (sparse![(0, 4), (3, 5), (7, 6)].to_dense(), sparse![(0, 1), (4, 2)], 
                sparse![(0, 4), (3, 5), (4, 8), (7, 16), (11, 12)]),
        
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], zero(), zero()),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], one(),  sparse![(0, 4), (1, 5), (3, 6), (4, 7)]),
        (sparse![(0, 4), (1, 5), (3, 6), (4, 7)], cst(3), sparse![(0, 12), (1, 15), (3, 18), (4, 21)]),
        (dense![4, 5, 0, 6, 7].to_sparse(), dense![1, 2, 3],             dense![4, 13, 22, 21, 19, 32, 21]),
        (sparse![(0, 4), (3, 5), (7, 6)], sparse![(0, 1), (4, 2)].to_dense(),
                sparse![(0, 4), (3, 5), (4, 8), (7, 16), (11, 12)]),
        (dense![4, 5, 0, 6, 7].to_sparse(), dense![1, 2, 3].to_sparse(), dense![4, 13, 22, 21, 19, 32, 21]),
        (sparse![(0, 4), (3, 5), (7, 6)], sparse![(0, 1), (4, 2)], 
                sparse![(0, 4), (3, 5), (4, 8), (7, 16), (11, 12)]),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2);
    }
}

#[test]
fn test_pow(){
    fn test(x: Polynomial<i64>, p: u32, exp: Polynomial<i64>){
        let ref_x = &x;
        assert_eq!(ref_x.pow(p), exp.clone());
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