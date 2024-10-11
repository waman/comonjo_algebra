use crate::{dense, spears, polynomial::Polynomial};

#[test]
fn test_dense_macro(){
    let p0 = dense![f64; 1, 2, 3];
    assert_eq!(format!("{}", p0), "1+2x+3x^2");

    let p1 = dense![f64; 0, 4, 5, 6, 0, 7, 0, 0];
    assert_eq!(format!("{}", p1), "4x+5x^2+6x^3+7x^5");
}

#[test]
fn test_spears_macro(){
    let p0 = spears![f64; (2, 6), (5, 1)];
    assert_eq!(format!("{}", p0), "6x^2+x^5");
}

#[test]
fn test_display(){
    fn test(p: Polynomial<i64>, exp: &str){
        assert_eq!(format!("{}", p), exp);
    }

    let table = [
        (Polynomial::<f64>::zero(), "0"),

        (Polynomial::<f64>::constant(0),  "0"),
        (Polynomial::<f64>::constant(1),  "1"),
        (Polynomial::<f64>::constant(-1), "-1"),
        (Polynomial::<f64>::constant(2),  "2"),

        // dense!
        (dense![f64;],      "0"),
        (dense![f64; 0, 0], "0"),
        (dense![f64; 1],    "1"),
        (dense![f64; -1],   "-1"),
        (dense![f64; 2],    "2"),

        (dense![f64; 1, 2, 3],        "1+2x+3x^2"),
        (dense![f64; -1, -2, -3],     "-1-2x-3x^2"),
        (dense![f64; 0, 2, 3],        "2x+3x^2"),
        (dense![f64; 1, 2, 3, 0, 0],  "1+2x+3x^2"),
        (dense![f64; 1, 1, 1, 1],     "1+x+x^2+x^3"),
        (dense![f64; -1, -1, -1, -1], "-1-x-x^2-x^3"),

        // spears!
        (spears![f64;], "0"),
        (spears![f64; (0, 0), (1, 0)], "0"),
        (spears![f64; (0, 1)],  "1"),
        (spears![f64; (0, -1)], "-1"),
        (spears![f64; (0, 2)],   "2"),

        (spears![f64; (0, 1), (2, 3)],   "1+3x^2"),
        (spears![f64; (0, -1), (2, -3)], "-1-3x^2"),
        (spears![f64; (1, 2), (2, 3)],   "2x+3x^2"),
        (spears![f64; (0, 1), (2, 3), (4, 0)], "1+3x^2"),
        (spears![f64; (0, 1), (1, 1), (2, 1), (3, 1)],     "1+x+x^2+x^3"),
        (spears![f64; (0, -1), (1, -1), (2, -1), (3, -1)], "-1-x-x^2-x^3"),
    ];

    for entry in table {
        test(entry.0, entry.1);
    }
}

#[test]
fn test_degree_zero_and_constant(){
    fn test(p: Polynomial<i64>, exp_deg: usize, exp_zero: bool, exp_const: bool){
        assert_eq!(p.degree(), exp_deg);
        assert_eq!(p.is_zero(), exp_zero);
        assert_eq!(p.is_constant(), exp_const);
    }

    let table = [
        (Polynomial::<f64>::zero(), 0, true, true),

        (Polynomial::<f64>::constant(0),  0, true, true),
        (Polynomial::<f64>::constant(1),  0, false, true),
        (Polynomial::<f64>::constant(-1), 0, false, true),
        (Polynomial::<f64>::constant(2),  0, false, true),

        // dense!
        (dense![f64; ],     0, true, true),  // zero & const
        (dense![f64; 0, 0], 0, true, true),  // zero & const
        (dense![f64; 1],    0, false, true),  // const
        (dense![f64; -1],   0, false, true),  // const
        (dense![f64; 2],    0, false, true),  // const

        (dense![f64; 1, 2, 3],       2, false, false),
        (dense![f64; -1, -2, -3],    2, false, false),
        (dense![f64; 0, 2, 3],       2, false, false),
        (dense![f64; 1, 2, 3, 0, 0], 2, false, false),

        // spears!
        (spears![f64; ],               0, true, true),  // zero & const
        (spears![f64; (0, 0), (2, 0)], 0, true, true),  // zero & const
        (spears![f64; (0, 1)],         0, false, true),  // const
        (spears![f64; (0, -1)],        0, false, true),  // const
        (spears![f64; (0, 2)],         0, false, true),  // const

        (spears![f64; (0, 1), (1, 2), (2, 3)],                  2, false, false),
        (spears![f64; (0, -1), (1, -2), (2, -3)],               2, false, false),
        (spears![f64; (0, 0), (1, 2), (2, 3)],                  2, false, false),
        (spears![f64; (0, 1), (1, 2), (2, 3), (5, 0), (10, 0)], 2, false, false),
    ];

    for entry in table {
        test(entry.0, entry.1, entry.2, entry.3);
    }
}