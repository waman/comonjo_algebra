pub(crate) mod dense_coeffs;
pub(crate) mod sparse_coeffs;
pub(crate) mod iter;

use std::{collections::{BTreeMap, HashMap}, fmt::{Debug, Display}, ops::*};
use num::{self, BigInt, BigRational, BigUint, One, Rational32, Rational64, Zero, complex::{Complex32, Complex64}, pow::Pow, traits::{ConstOne, ConstZero, Euclid}};
use once_cell::sync::Lazy;

use crate::{algebra::*, poly::{dense_coeffs::DenseCoeffs, iter::{CoeffsIter, IntoCoeffsIter, IntoNonzeroCoeffsIter, NonzeroCoeffsIter}, sparse_coeffs::SparseCoeffs}};

/// Polynomial type.
/// Refer to spire's [Polynomial](https://github.com/typelevel/spire/blob/main/core/src/main/scala/spire/math/Polynomial.scala).
/// 
/// ## Method Name
/// `Polynomial` has several similar-named methods like `reciprocal()` and `new_reciprocal()`.
/// In general, `new_xxx()` methods create a new `Polynomial` instance and `self` remains unchanged,
/// while the others mutate `self` (and return nothing).
/// 
///     use comonjo_algebra::poly::Polynomial;
///     use comonjo_algebra::dense;
/// 
///     // new_reciprocal()
///     let p = dense![1, 2, 0, 3];  // 1 + 2x + 3x³
///     assert_eq!(p.new_reciprocal(), dense![3, 0, 2, 1]);
/// 
///     assert_eq!(p, dense![1, 2, 0, 3]);  // p is not mutated
/// 
///     // reciprocal()
///     let mut q = dense![1, 2, 0, 3];  // 1 + 2x + 3x³
///     q.reciprocal();
///     assert_eq!(q, dense![3, 0, 2, 1]);  // r is mutated
/// 
/// In a bit addition, `new_xxx()` methods demand `Clone` trait to coefficient type,
/// while the others do not in almost all cases.
/// 
/// ## Iterator
/// `Polynomial`'s `into_iter()` method returns an `Iterator`
/// that iterates *only nonzero* coefficients (with its term's degree).
/// `Polynomial` and `&Polynomial` implement the `IntoIterator` trait respectively,
/// so those methods behave differently:
/// 
///     use comonjo_algebra::poly::Polynomial;
///     use comonjo_algebra::dense;
/// 
///     let p = dense![1, 2, 0, 3];  // 1 + 2x + 3x³
/// 
///     // &Polynomial
///     let mut ite0 = (&p).into_iter();
///     assert_eq!(ite0.next(), Some((0, &1)));
///     assert_eq!(ite0.next(), Some((1, &2)));
///     assert_eq!(ite0.next(), Some((3, &3)));
///     assert_eq!(ite0.next(), None);
/// 
///     // Polynomial
///     let mut ite1 = p.into_iter();
///     assert_eq!(ite1.next(), Some((0, 1)));
///     assert_eq!(ite1.next(), Some((1, 2)));
///     assert_eq!(ite1.next(), Some((3, 3)));
///     assert_eq!(ite1.next(), None);
/// 
/// To iterate coefficients including zero, import `CoeffsIterator` and
/// use `coeffs()` methods:
/// 
///     use comonjo_algebra::poly::Polynomial;
///     use comonjo_algebra::dense;
///     use comonjo_algebra::poly::CoeffsIterator;
/// 
///     let p = dense![1, 2, 0, 3];  // 1 + 2x + 3x³
///     let mut ite = p.coeffs();
///     assert_eq!(ite.next(), Some(1));
///     assert_eq!(ite.next(), Some(2));
///     assert_eq!(ite.next(), Some(0));
///     assert_eq!(ite.next(), Some(3));
///     assert_eq!(ite.next(), None);
/// 
/// `CoeffsIterator::nonzero_coeffs()` methods are equivalent to `into_iter()`.
pub enum Polynomial<C> where C: Semiring {

    Zero(),

    /// Use `Polynomial::constant()` function to instantiate.
    Constant(ConstCoeff<C>),

    /// Use `dense!()` macro to instantiate.
    Dense(DenseCoeffs<C>),
    
    /// Use `sparse!()` macro to instantiate.
    Sparse(SparseCoeffs<C>)
}

pub struct ConstCoeff<C: Semiring>(pub(crate) C);

#[macro_export]
macro_rules! dense {
    // [ $t:ty; $( $x:expr ),* ] => {
    //     Polynomial::<$t>::dense_from_vec(vec![ $( $x ),* ])
    // };
    [ $( $x:expr ),* ] => {
        Polynomial::dense_from_vec(vec![ $( $x ),* ])
    };
    [ $( $x:expr ),+ , ] => {
        dense![ $( $x ),* ]
    };
}

#[macro_export]
macro_rules! sparse {
    [ $( ($order:expr, $coeff:expr) ),* ] => {
        Polynomial::sparse_from_map(std::collections::BTreeMap::from([ $( ($order, $coeff) ),* ]))
    };
    [ $( ($order:expr, $coeff:expr) ),+ , ] => {
        sparse![ $( ($order, $coeff) ),* ]
    };
}

//********** Factory Methods ******
impl<C> Polynomial<C> where C: Semiring {

    // Not to check the arg not to be zero.
    #[inline]
    pub(crate) fn new_raw_const(c: C) -> Polynomial<C> {
        debug_assert!(!c.is_zero());

        Polynomial::Constant(ConstCoeff(c))
    }

    // Not to remove the tailing zeros, to check to be constant, and not to execute `shrink_to_fit()`
    #[inline]
    pub(crate) fn new_raw_dense(vec: Vec<C>) -> Polynomial<C> {
        debug_assert!(vec.len() > 1);
        debug_assert!(!vec.last().unwrap().is_zero());
        debug_assert_eq!(vec.len(), vec.capacity());

        Polynomial::Dense(DenseCoeffs(vec))
    }

    // Not to remove zero-value entries, and to check to be const
    #[inline]
    pub(crate) fn new_raw_sparse(map: BTreeMap<usize, C>) -> Polynomial<C> {
        debug_assert!(if map.contains_key(&0) { map.len() > 1 } else { !map.is_empty() } );
        debug_assert!(map.values().all(|v|!v.is_zero()));

        Polynomial::Sparse(SparseCoeffs(map))
    }

    /// Returns a new `Constant` polynomial.
    /// If the argument is zero, return the `Zero` polynomial.
    pub fn constant(c: C) -> Polynomial<C> {
        if c.is_zero() {
            Polynomial::Zero()
        } else {
            Polynomial::new_raw_const(c)
        }
    }

    /// Returns a new `Dense` polynomial whose coefficients are specified by the argument.
    /// If the length of the argument `Vec` is 0 or 1, return `Zero` or `Constant` polynomial respectively.
    pub fn dense_from_vec(mut coeffs: Vec<C>) -> Polynomial<C>  {

        remove_tail_zeros(&mut coeffs);

        match coeffs.len() {
            0 => Polynomial::Zero(),
            1 => Polynomial::new_raw_const(coeffs.pop().unwrap()),
            _ => {
                coeffs.shrink_to_fit();
                Polynomial::new_raw_dense(coeffs)
            }
        }
    }

    /// Returns a new `Sparse` polynomial whose coefficients are specified by the argument.
    /// If the length of the argument `BTreeMap` is 0 or 1, return `Zero` or `Constant` polynomial respectively.
    pub fn sparse_from_map(mut coeffs: BTreeMap<usize, C>) -> Polynomial<C> {

        coeffs.retain(|_, v| !v.is_zero());

        match coeffs.len() {
            0 => Polynomial::Zero(),
            1 => {
                if coeffs.contains_key(&0) {
                    Polynomial::new_raw_const(coeffs.remove(&0).unwrap())
                } else {
                    Polynomial::new_raw_sparse(coeffs)
                }
            },
            _ => Polynomial::new_raw_sparse(coeffs)
        }
    }

//     pub fn parse(s: &str) -> Polynomial<C> {
    
// //    private[this] val termRe = "([0-9]+\\.[0-9]+|[0-9]+/[0-9]+|[0-9]+)?(?:([a-z])(?:\\^([0-9]+))?)?".r
 
// //    private[this] val operRe = " *([+-]) *".r
 
// //    private[spire] def parse(s: String): Polynomial[Rational] = {
 
// //      // represents a term, plus a named variable v
// //      case class T(c: Rational, v: String, e: Int)
 
// //      // parse all the terms and operators out of the string
// //      @tailrec def parse(s: String, ts: List[T]): List[T] =
// //        if (s.isEmpty) {
// //          ts
// //        } else {
// //          val (op, s2) = operRe.findPrefixMatchOf(s) match {
// //            case Some(m) => (m.group(1), s.substring(m.end))
// //            case None    => if (ts.isEmpty) ("+", s) else throw new IllegalArgumentException(s)
// //          }
 
// //          val m2 = termRe.findPrefixMatchOf(s2).getOrElse(throw new IllegalArgumentException(s2))
// //          val c0 = Option(m2.group(1)).getOrElse("1")
// //          val c = if (op == "-") "-" + c0 else c0
// //          val v = Option(m2.group(2)).getOrElse("")
// //          val e0 = Option(m2.group(3)).getOrElse("")
// //          val e = if (e0 != "") e0 else if (v == "") "0" else "1"
 
// //          val t =
// //            try {
// //              T(Rational(c), v, e.toInt)
// //            } catch {
// //              case _: Exception => throw new IllegalArgumentException(s"illegal term: $c*x^$e")
// //            }
// //          parse(s2.substring(m2.end), if (t.c == 0) ts else t :: ts)
// //        }
 
// //      // do some pre-processing to remove whitespace/outer parens
// //      val t = s.trim
// //      val u = if (t.startsWith("(") && t.endsWith(")")) t.substring(1, t.length - 1) else t
// //      val v = Term.removeSuperscript(u)
 
// //      // parse out the terms
// //      val ts = parse(v, Nil)
 
// //      // make sure we have at most one variable
// //      val vs = ts.view.map(_.v).toSet.filter(_ != "")
// //      if (vs.size > 1) throw new IllegalArgumentException("only univariate polynomials supported")
 
// //      // we're done!
// //      ts.foldLeft(Polynomial.zero[Rational])((a, t) => a + Polynomial(t.c, t.e))
// //    }
//         todo!()
//     }

//     pub fn parse_to_dense(s: &str) -> Polynomial<C> {
//         todo!()
//     }

//     pub fn parse_to_sparse(s: &str) -> Polynomial<C> {
//         todo!()
//     }
    
//     pub fn interpolate(points: &[(C, C)]) -> Polynomial<C> {
//     //   def interpolate[C: Field: Eq: ClassTag](points: (C, C)*): Polynomial[C] = {
//     //     def loop(p: Polynomial[C], xs: List[C], pts: List[(C, C)]): Polynomial[C] =
//     //       pts match {
//     //         case Nil =>
//     //           p
//     //         case (x, y) :: tail =>
//     //           val c = Polynomial.constant((y - p(x)) / xs.map(x - _).qproduct)
//     //           val prod = xs.foldLeft(Polynomial.one[C]) { (prod, xn) =>
//     //             prod * (Polynomial.x[C] - constant(xn))
//     //           }
//     //           loop(p + c * prod, x :: xs, tail)
//     //       }
//     //     loop(Polynomial.zero[C], Nil, points.toList)
//     //   }
//     // }
//         todo!()
//     }

    /// `x`
    pub fn x() -> Polynomial<C> { sparse![(1, C::one())] }

    /// `x²`
    pub fn x2() -> Polynomial<C> { sparse![(2, C::one())] }

    /// `x³`
    pub fn x3() -> Polynomial<C> { sparse![(3, C::one())] }

    /// `x⁴`
    pub fn x4() -> Polynomial<C> { sparse![(4, C::one())] }

    /// `x⁵`
    pub fn x5() -> Polynomial<C> { sparse![(5, C::one())] }

    /// `2x`
    pub fn two_x() -> Polynomial<C> {
        sparse![(1, C::one() + C::one())]
    }

    /// Creates a linear polynomial *ax*.
    pub fn linear_monomial(a: C) -> Polynomial<C> { sparse![(1, a)] }

    /// Creates a linear polynomial *ax + b*. Note the coefficients order.
    pub fn linear(a: C, b: C) -> Polynomial<C> { sparse![(1, a), (0, b)] }

    /// Creates a quadratic polynomial *ax²*.
    pub fn quadratic_monomial(a: C) -> Polynomial<C> { sparse![(2, a)] }

    /// Creates a quadratic polynomial *ax² + bx + c*. Note the coefficients order.
    pub fn quadratic(a: C, b: C, c: C) -> Polynomial<C> { sparse![(2, a), (1, b), (0, c)] }

    /// Creates a cubic polynomial *ax³ + bx² + cx + d*. Note the coefficients order.
    pub fn cubic_monomial(a: C) -> Polynomial<C> { sparse![(3, a)] }

    /// Creates a cubic polynomial *ax³ + bx² + cx + d*. Note the coefficients order.
    pub fn cubic(a: C, b: C, c: C, d: C) -> Polynomial<C> { sparse![(3, a), (2, b), (1, c), (0, d)] }
}

pub(crate) fn remove_tail_zeros<C>(vec: &mut Vec<C>) where C: Semiring {
    while let Some(c) = vec.last() {
        if c.is_zero() {
            vec.pop();
        } else {
            break;
        }
    }
}

pub(crate) fn factorial<C>(n: usize) -> C where C: Semiring + num::FromPrimitive {
    if n == 0 || n == 1 { return C::one(); }

    let mut result = C::one();
    for i in 2..=n {
        result = result * C::from_usize(i).unwrap();
    }
    result
}

//********** METHODS with Semiring coefficients ******
impl<C> Polynomial<C> where C: Semiring {
    
    //   /**
    //    * Evaluate the polynomial at `x`.
    //    */
    //   def apply(x: C)(implicit r: Semiring[C]): C
    
    //   def evalWith[A: Semiring: Eq: ClassTag](x: A)(f: C => A): A =
    //     this.map(f).apply(x)

    /// Returns the degree of the `self` polynomial.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.degree(), 3);
    /// 
    ///     // if zero polynomial
    ///     let zero: Polynomial<i64> = Polynomial::Zero();
    ///     assert_eq!(zero.degree(), 0);
    /// 
    ///     // if constant polynomial
    ///     let cst: Polynomial<i64> = Polynomial::constant(2);
    ///     assert_eq!(cst.degree(), 0)
    /// 
    pub fn degree(&self) -> usize {
        match self {
            Polynomial::Dense(dc) => dc.degree(),
            Polynomial::Sparse(sc) => sc.degree(),
            _ => 0,
        }
    }

    /// Returns the reference of n-th coefficient if exists.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.nth(0), Some(&1));
    ///     assert_eq!(p.nth(2), Some(&3));
    ///     assert_eq!(p.nth(10), None);
    /// 
    pub fn nth(&self, n: usize) -> Option<&C> {
        match self {
            Polynomial::Zero() => None,
            Polynomial::Constant(cc) => if n == 0 { Some(&cc.0) } else { None },
            Polynomial::Dense(dc) => dc.nth(n),
            Polynomial::Sparse(sc) => sc.nth(n),
        }
    }

    /// Returns the max order term's degree and coefficient as a tuple if exists.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.max_order_term(), Some((3, &4)));
    /// 
    ///     // if zero polynomial
    ///     let zero: Polynomial<i64> = Polynomial::Zero();
    ///     assert_eq!(zero.max_order_term(), None);
    /// 
    ///     // if constant polynomial
    ///     let cst: Polynomial<i64> = Polynomial::constant(2);
    ///     assert_eq!(cst.max_order_term(), Some((0, &2)))
    /// 
    pub fn max_order_term(&self) -> Option<(usize, &C)> {
        match self {
            Polynomial::Zero() => None,
            Polynomial::Constant(cc) => Some((0, &cc.0)),
            Polynomial::Dense(dc) => dc.max_order_term(),
            Polynomial::Sparse(sc) => sc.max_order_term(),
        }
    }

    /// Returns the min order term's degree and coefficient as a tuple if exists.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![0, 0, 1, 2, 3];  // x² + 2x³ + 3x⁴
    ///     assert_eq!(p.min_order_term(), Some((2, &1)));
    /// 
    ///     // if zero polynomial
    ///     let zero: Polynomial<i64> = Polynomial::Zero();
    ///     assert_eq!(zero.min_order_term(), None);
    /// 
    ///     // if constant polynomial
    ///     let cst: Polynomial<i64> = Polynomial::constant(2);
    ///     assert_eq!(cst.min_order_term(), Some((0, &2)))
    /// 
    pub fn min_order_term(&self) -> Option<(usize, &C)> {
        match self {
            Polynomial::Zero() => None,
            Polynomial::Constant(cc) => Some((0, &cc.0)),
            Polynomial::Dense(dc) => dc.min_order_term(),
            Polynomial::Sparse(sc) => sc.min_order_term(),
        }
    }

    /// Returns true if `self` is constant, or false for otherwise.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![0, 0, 1, 2, 3];  // x² + 2x³ + 3x⁴
    ///     assert_eq!(p.min_order_term(), Some((2, &1)));
    /// 
    ///     // if zero polynomial
    ///     let zero: Polynomial<i64> = Polynomial::Zero();
    ///     assert_eq!(zero.min_order_term(), None);
    /// 
    ///     // if constant polynomial
    ///     let cst: Polynomial<i64> = Polynomial::constant(2);
    ///     assert_eq!(cst.min_order_term(), Some((0, &2)))
    /// 
    pub fn is_constant(&self) -> bool {
        match self {
            Polynomial::Zero() | Polynomial::Constant(_) => true,
            _ => false,
        }
    }

    /// Returns true if `self` is x, that is, a polynomial its 1st-order coefficient is 1 and others are 0.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p0: Polynomial<i64> = dense![0, 1];  // x
    ///     assert!(p0.is_x());
    /// 
    ///     let p1: Polynomial<i64> = dense![0, 2];  // 2x
    ///     assert!(!p1.is_x());
    /// 
    ///     let p2: Polynomial<i64> = dense![0, 0, 1];  // x²
    ///     assert!(!p2.is_x());
    /// 
    pub fn is_x(&self) -> bool {
        match self {
            Polynomial::Zero() | Polynomial::Constant(_) => false,
            Polynomial::Dense(dc) => dc.is_x(),
            Polynomial::Sparse(sc) => sc.is_x(),
        }
    }

    /// Returns true if `self` is dense.
    /// Note that this returns false if self is Zero or Constant.
    pub fn is_dense(&self) -> bool {
        match self {
            Polynomial::Dense(_) => true,
            _ => false,
        }
    }

    /// Returns true if `self` is sparse.
    /// Note that this returns false if `self` is `Zero` or `Constant`.
    pub fn is_sparse(&self) -> bool {
        match self {
            Polynomial::Sparse(_) => true,
            _ => false,
        }
    }
     
    /// Returns dense expression of this polynomial if this is sparse, otherwise (zero or constant) self.
    pub fn to_dense(self) -> Polynomial<C> {
        match self {
            Polynomial::Sparse(sc) => Polynomial::new_raw_dense(sc.to_vec()),
            _ => self,
        }
    }

    /// Returns sparse expression of this polynomial if this is dense, otherwise (zero or constant) self.
    pub fn to_sparse(self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => Polynomial::new_raw_sparse(dc.to_map()),
            _ => self
        }
    }

    pub fn sign_variations(&self) -> usize {
//   /**
//    * Returns the number of sign variations in the coefficients of this polynomial. Given 2 consecutive terms (ignoring 0
//    * terms), a sign variation is indicated when the terms have differing signs.
//    */
//   def signVariations(implicit ring: Semiring[C], order: Order[C], signed: Signed[C]): Int = {
//     var prevSign: Sign = Signed.Zero
//     var variations = 0
//     foreachNonzero { (_, c) =>
//       val sign = signed.sign(c)
//       if (Signed.Zero != prevSign && sign != prevSign) {
//         variations += 1
//       }
//       prevSign = sign
//     }
//     variations
//   }
        todo!()
    }

    /// Returns the reciprocal polynomial 
    /// (Reference: [Reciprocal polynomial](http://en.wikipedia.org/wiki/Reciprocal_polynomial)).
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     p.reciprocal();
    ///     assert_eq!(p, dense![4, 3, 2, 1]);  // 4 + 3x + 2x² + x³
    /// 
    pub fn reciprocal(&mut self) {
        match self {
            Polynomial::Dense(dc) => 
                if let Some(p) = dc.reciprocal() { *self = p; },
            Polynomial::Sparse(sc) => 
                if let Some(p) = sc.reciprocal() { *self = p; },
            _ => (),
        }
    }

    /// Removes all zero roots from the `self` polynomial.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p = dense![0, 0, 1, 2, 3];  // x² + 2x³ +3x⁵
    ///     p.remove_zero_roots();
    ///     assert_eq!(p, dense![1, 2, 3]);  // 1 + 2x + 3x²
    /// 
    pub fn remove_zero_roots(&mut self) {
        match self {
            Polynomial::Dense(dc) => 
                if let Some(p) = dc.remove_zero_roots() { *self = p; },
            Polynomial::Sparse(sc) =>
                if let Some(p) = sc.remove_zero_roots() { *self = p; },
            _ => (),
        }
    }

    /// Removes the max order term from the `self` polynomial.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     p.reductum();
    ///     assert_eq!(p, dense![1, 2, 3]);  // 1 + 2x + 3x²
    /// 
    pub fn reductum(&mut self) {
        match self {
            Polynomial::Zero() => (),
            Polynomial::Constant(_) => *self = Polynomial::Zero(),
            Polynomial::Dense(dc) =>
                if let Some(p) = dc.reductum() { *self = p; },
            Polynomial::Sparse(sc) =>
                if let Some(p) = sc.reductum() { *self = p; },
        }
    }
}

//********** Methods with Semiring + Clone Coefficients **********/
impl<C> Clone for Polynomial<C> where C: Semiring + Clone {

    fn clone(&self) -> Self {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::new_raw_const(cc.0.clone()),
            Polynomial::Dense(dc) => Polynomial::Dense(dc.clone()),
            Polynomial::Sparse(sc) => Polynomial::Sparse(sc.clone()),
        }
    }
}

impl<C> Polynomial<C> where C: Semiring + Clone {

    /// Creates a dense clone of this polynomial if the `self` is sparse, otherwise (zero, constant or dense) returns `self`.
    pub fn dense_clone(&self) -> Polynomial<C> {
        match self {
            Polynomial::Sparse(sc) => {
                let n: usize = sc.degree() + 1;
                let mut v: Vec<C> = Vec::with_capacity(n);
                for i in 0..n {
                    match sc.nth(i) {
                        Some(c) => v.push(c.clone()),
                        None => v.push(C::zero()),
                    } 
                }
                Polynomial::new_raw_dense(v)
            },
            _ => self.clone()
        }
    }

    /// Creates a sparse clone of this polynomial if `self` is dense, otherwise (zero, constant or sparse) returns `self`.
    pub fn sparse_clone(&self) -> Polynomial<C> {
        match self {
            d @ Polynomial::Dense(_) => {
                let mut map = BTreeMap::new();
                for (i, c) in d.nonzero_coeffs() {
                    map.insert(i, c.clone());  // note c is not zero
                }
                Polynomial::new_raw_sparse(map)
            },
            _ => self.clone()
        }
    }

    /// Returns a reciprocal polynomial of the `self` as a new instance.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.new_reciprocal(), dense![4, 3, 2, 1]);  // 4 + 3x + 2x² + x
    /// 
    pub fn new_reciprocal(&self) -> Polynomial<C> {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            c @ Polynomial::Constant(_)=> c.clone(),
            Polynomial::Dense(dc) => dc.new_reciprocal(),
            Polynomial::Sparse(sc) => sc.new_reciprocal(),
        }
    }

    /// Returns a polynomial that zero roots are removed from the `self` polynomial, as a new instance.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![0, 0, 1, 2, 3];  // x² + 2x³ + 3x⁴
    ///     assert_eq!(p.new_zero_roots_removed(), dense![1, 2, 3]);  // 1 + 2x + 3x²
    /// 
    pub fn new_zero_roots_removed(&self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => dc.new_zero_roots_removed(),
            Polynomial::Sparse(sc) => sc.new_zero_roots_removed(),
            _ => self.clone(),
        }
    }

    /// Returns a polynomial that the max-order term is removed from the `self` polynomial, as a new instance.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.new_reductum(), dense![1, 2, 3]);  // 1 + 2x + 3x²
    /// 
    pub fn new_reductum(&self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => dc.new_reductum(),
            Polynomial::Sparse(sc) => sc.new_reductum(),
            _ => Polynomial::Zero(),
        }
    }

    fn scale_by_left(self, k: &C) -> Polynomial<C> {
        if k.is_zero() { return Polynomial::Zero() }
        if k.is_one() { return self }
        self.map_nonzero(|_, c| k.ref_mul(c))
    }

    fn ref_scale_by_left(&self, k: &C) -> Polynomial<C> {
        if k.is_zero() { return Polynomial::Zero() }
        if k.is_one() { return self.clone() }
        self.map_nonzero(|_, c| k.ref_mul(c))
    }
}
 
//********** Display and Debug **********
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

        } else {  // Rational, Complex, etc.
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

impl<C> Display for Polynomial<C> where C: Semiring + Display {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        match self {
            Polynomial::Zero() => f.write_str("(0)"),
            Polynomial::Constant(c) => f.write_fmt(format_args!("({})", c.0)),
            _ => {
                let s: String = self.nonzero_coeffs().map(|(i, c)|term_to_string(i, c)).collect();
                let first_sign = if s.starts_with(" - ") { "-" } else { "" };
                f.write_fmt(format_args!("{}{}", first_sign, &s[3..]))
            },
        }
    }
}
 
impl<C> Debug for Polynomial<C> where C: Semiring + Display {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Polynomial::Zero() => "Zero",
            Polynomial::Constant(_) => "Constant",
            Polynomial::Dense(_) => "Dense",
            Polynomial::Sparse(_) => "Sparse",
        };
        f.write_fmt(format_args!("Polynomial::<{}>::{}[{}]", std::any::type_name::<C>(), s, self))
    }
}

//********** Eq, PartialEq, Zero and One **********
impl<C> PartialEq for Polynomial<C> where C: Semiring {

    fn eq(&self, other: &Polynomial<C>) -> bool {

        if self.degree() != other.degree() { return false; }

        fn has_the_same_content<C>(vec: &Vec<C>, map: &BTreeMap<usize, C>) -> bool where C: num::Zero + PartialEq {
            let vec_iter = vec.iter().enumerate().filter(|(_, c)| !c.is_zero());
            let map_iter = map.iter();
            vec_iter.zip(map_iter).all(|(v, m)| v.0 == *m.0 && v.1 == m.1)
        }

        match self {
            Polynomial::Zero() => other.is_zero(),
            
            Polynomial::Constant(lhs) => match other {
                Polynomial::Constant(rhs) => lhs.0 == rhs.0,
                _ => false,
            },

            Polynomial::Dense(lhs) => match other {
                Polynomial::Dense(rhs) => lhs.0 == rhs.0,
                Polynomial::Sparse(rhs) => has_the_same_content(&lhs.0, &rhs.0),
                _ => false,
            },

            Polynomial::Sparse(lhs) => match other {
                Polynomial::Dense(rhs) => has_the_same_content(&rhs.0, &lhs.0),
                Polynomial::Sparse(rhs) => lhs.0 == rhs.0,
                _ => false,
            },
        }
    }
}

impl<C> Eq for Polynomial<C> where C: Semiring + Eq {}

impl<C> Zero for Polynomial<C> where C: Semiring {

    fn zero() -> Self { Polynomial::Zero() }

    fn is_zero(&self) -> bool {
        match self {
            Polynomial::Zero() => true,
            _ => false,
        }
    }
}

impl<C> ConstZero for Polynomial<C> where C: Semiring + ConstZero + Clone {
    const ZERO: Self = Polynomial::Zero();
}

impl<C> One for Polynomial<C> where C: Semiring + Clone {

    fn one() -> Self { 
        Polynomial::new_raw_const(C::one())
    }
    
    fn is_one(&self) -> bool {
        match self {
            Polynomial::Constant(cc) => cc.0.is_one(),
            _ => false,
        }
    }
}

impl<C> ConstOne for Polynomial<C> where C: Semiring + Clone + ConstOne {
    const ONE: Self = Polynomial::Constant(ConstCoeff(C::ONE));
}

//********** Iterator **********
/// Iterate nonzero coefficients.
impl<C> IntoIterator for Polynomial<C> where C: Semiring {

    type Item = (usize, C);
    type IntoIter = IntoNonzeroCoeffsIter<C>;

    /// Return `impl Iterator<Item=(usize, C)>`.
    fn into_iter(self) -> Self::IntoIter {
        match self {
            Polynomial::Zero() => IntoNonzeroCoeffsIter::Zero(),
            Polynomial::Constant(cc) => IntoNonzeroCoeffsIter::Constant(Some(cc.0)),
            Polynomial::Dense(dc) => dc.into_nonzero_coeffs_iter(),
            Polynomial::Sparse(sc) => sc.into_nonzero_coeffs_iter(),
        }
    }
}

/// Iterate reference of nonzero coefficient.
impl<'a, C> IntoIterator for &'a Polynomial<C> where C: Semiring {

    type Item = (usize, &'a C);
    type IntoIter = NonzeroCoeffsIter<'a, C>;

    /// Iterate nonzero coefficients.
    /// Return `impl Iterator<Item=(usize, &C)>`.
    fn into_iter(self) -> Self::IntoIter {
        match self {
            Polynomial::Zero() => NonzeroCoeffsIter::Zero(),
            Polynomial::Constant(cc) => NonzeroCoeffsIter::Constant(Some(&cc.0)),
            Polynomial::Dense(dc) => dc.nonzero_coeffs_iter(),
            Polynomial::Sparse(sc) => sc.nonzero_coeffs_iter(),
        }
    }
}

pub trait CoeffsIterator<C>: IntoIterator where Self: Sized, C: Semiring, Self::IntoCoeffsIter: Iterator<Item=Self::Coeff>{
    type Coeff;
    type IntoCoeffsIter;

    fn coeffs(self) -> Self::IntoCoeffsIter;

    #[inline]
    fn nonzero_coeffs(self) -> Self::IntoIter { self.into_iter() }
}

impl<C> CoeffsIterator<C> for Polynomial<C> where C: Semiring {
    
    type Coeff = C;
    type IntoCoeffsIter = IntoCoeffsIter<C>;

    /// Return `impl Iterator<Item=C>`.
    fn coeffs(self) -> Self::IntoCoeffsIter {
        match self {
            Polynomial::Zero() => IntoCoeffsIter::Zero(),
            Polynomial::Constant(cc) => IntoCoeffsIter::Constant(Some(cc.0)),
            Polynomial::Dense(dc) => dc.into_coeffs_iter(),
            Polynomial::Sparse(sc) => sc.into_coeffs_iter(), 
        }
    }
}

impl<'a, C> CoeffsIterator<C> for &'a Polynomial<C> where C: Semiring {
    
    type Coeff = Option<&'a C>;
    type IntoCoeffsIter = CoeffsIter<'a, C>;

    /// Return `impl Iterator<Item=Option<Option<&C>>>`.
    fn coeffs(self) -> Self::IntoCoeffsIter {
        match self {
            Polynomial::Zero() => CoeffsIter::Zero(),
            Polynomial::Constant(cc) => CoeffsIter::Constant(Some(Some(&cc.0))),
            Polynomial::Dense(dc) => dc.coeffs_iter(),
            Polynomial::Sparse(sc) => sc.coeffs_iter(), 
        }
    }
}

//********** Polynomial Basic Operations **********/
pub trait PolynomialOps<C> where C: Semiring {
    
    fn to_vec(self) -> Vec<C>;
    fn to_map(self) -> BTreeMap<usize, C>;
    fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D;
}

impl<C> PolynomialOps<C> for Polynomial<C> where C: Semiring {
    
    fn to_vec(self) -> Vec<C> {
        match self {
            Polynomial::Zero() => Vec::new(),
            Polynomial::Constant(cc) => vec![cc.0],
            Polynomial::Dense(dc) => dc.0,
            Polynomial::Sparse(sc) => sc.to_vec(),
        }
    }
    
    fn to_map(self) -> BTreeMap<usize, C> {
        match self {
            Polynomial::Zero() => BTreeMap::new(),
            Polynomial::Constant(cc) => BTreeMap::from([(0, cc.0)]),
            Polynomial::Dense(dc) => dc.to_map(),
            Polynomial::Sparse(sc) => sc.0,
        }
    }
    
    /// Note that `f` is applied only to nonzero coefficients
    fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::constant(f(0, cc.0)),
            Polynomial::Dense(dc) => dc.map_nonzero(f),
            Polynomial::Sparse(sc) => sc.map_nonzero(f),
        }
    }
}

impl<'a, C> PolynomialOps<C> for &'a Polynomial<C> where C: Semiring + Clone {
    
    fn to_vec(self) -> Vec<C> {
        self.coeffs().map(|option_c| match option_c {
            Some(c) => c.clone(),
            None => C::zero(),
        }).collect()
    }
    
    fn to_map(self) -> BTreeMap<usize, C> {
        self.nonzero_coeffs().map(|e| (e.0, e.1.clone())).collect()
    }
    
    /// Note that `f` is applied only to nonzero coefficients
    fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::constant(f(0, cc.0.clone())),
            Polynomial::Dense(dc) => dc.map_nonzero(f),
            Polynomial::Sparse(sc) => sc.map_nonzero(f),
        }
    }
}

//********** Compose **********/
fn compose_polynomials<'a, 'b, C> (lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut acc = Polynomial::Zero();

    for (i, c) in lhs.nonzero_coeffs() {
        let z = rhs.pow(i as u32).scale_by_left(c);
        acc = acc + z;
    }

    acc
}

pub trait Compose<C, RHS> where C: Semiring + Clone {

    /// Composes this polynomial with another.
    fn compose(self, other: RHS) -> Polynomial<C>;
}

impl<'a, C> Compose<C, Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    fn compose(self, y: Polynomial<C>) -> Polynomial<C> { self.compose(&y) }
}

impl<'a, 'b, C> Compose<C, &'b Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    fn compose(self, y: &'b Polynomial<C>) -> Polynomial<C> {
        match (self, y) {
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (lhs @ Polynomial::Constant(_), _) => lhs.clone(),
            (lhs, Polynomial::Zero()) => match lhs.nth(0) {
                Some(c) => Polynomial::constant(c.clone()),
                _ => Polynomial::Zero(),
            },
            (lhs, rhs) if rhs.is_x() => lhs.clone(),
            (lhs, rhs) => compose_polynomials(lhs, rhs),
        }
    }
}
    
//********** Methods with Semiring Coefficients **********/
impl<C> Polynomial<C> where C: Semiring + num::FromPrimitive {

    /// Differentiates the `self` polynomial.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     p.differentiate();
    ///     assert_eq!(p, dense![2, 6, 12]);  // 2 + 6x + 12x²
    /// 
    pub fn differentiate(&mut self) {
        match self {
            Polynomial::Dense(dc) => 
                if let Some(p) = dc.differentiate() { *self = p },
            Polynomial::Sparse(sc) => 
                if let Some(p) = sc.differentiate() { *self = p },
            _ => *self = Polynomial::Zero(),
        }
    }

    /// differentiates the `self` polynomial `n` times.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     p.differentiate_n(2);
    ///     assert_eq!(p, dense![6, 24]);  // 6 + 24x
    /// 
    pub fn differentiate_n(&mut self, n: usize) {
        if n == 0 { return; }
        if n == 1 { 
            self.differentiate();
            return;
        }

        match self {
            Polynomial::Dense(dc) => 
                if let Some(p) = dc.differentiate_n(n) { *self = p; },
            Polynomial::Sparse(sc) => 
                if let Some(p) = sc.differentiate_n(n) { *self = p },
            _ => *self = Polynomial::Zero(),
        }
    }
}    

impl<C> Polynomial<C> where C: Semiring + num::FromPrimitive + Clone {
    
    /// Returns the derivative of the `self` polynomial.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.new_derivative(), dense![2, 6, 12]);  // 2 + 6x + 12x²
    /// 
    pub fn new_derivative(&self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => dc.new_derivative(),
            Polynomial::Sparse(sc) => sc.new_derivative(),
            _ => Polynomial::Zero(),
        }
    }
    
    /// Returns the n-th derivative of the `self` polynomial.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.new_nth_derivative(2), dense![6, 24]);  // 6 + 24x
    /// 
    pub fn new_nth_derivative(&self, n: usize) -> Polynomial<C> {
        if n == 0 { return self.clone(); }
        if n == 1 { return self.new_derivative(); }

        match self {
            Polynomial::Dense(dc) => dc.new_nth_derivative(n),
            Polynomial::Sparse(sc) => sc.new_nth_derivative(n),
            _ => Polynomial::Zero(),
        }
    }
}

//********** Methods with Ring Coefficients **********/
impl<C> Polynomial<C> where C: Ring {

    /// Flips the `self` polynomial, that is, changes the sign of the odd-order term's coefficients.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     p.flip();
    ///     assert_eq!(p, dense![1, -2, 3, -4]);  // 1 - 2x + 3x² - 4x³
    /// 
    /// The `flip()` method returns the same result as `compose()` with `-x`.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     use comonjo_algebra::poly::Compose;
    /// 
    ///     let mut p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     let q = p.clone();
    /// 
    ///     p.flip();
    ///     assert_eq!(p, q.compose(-Polynomial::x()));
    /// 
    pub fn flip(&mut self){
        match self {
            Polynomial::Dense(dc) => dc.flip(),
            Polynomial::Sparse(sc) => sc.flip(),
            _ => (),
        }
    }
}

impl<'a, C> Polynomial<C> where C: Ring + Clone {

    /// Returns the flipped polynomial of the `self` polynomial,
    /// that is, the polynomial whose odd-order terms' coefficients are changed.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.new_flipped(), dense![1, -2, 3, -4]);  // 1 - 2x + 3x² - 4x³
    /// 
    /// The `new_flipped()` method returns the same result as `compose()` with `-x`.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     use comonjo_algebra::poly::Compose;
    /// 
    ///     let p = dense![1, 2, 3, 4];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.new_flipped(), p.compose(-Polynomial::x()));
    /// 
    pub fn new_flipped(&self) -> Polynomial<C> {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(_) => self.clone(),
            Polynomial::Dense(dc) => dc.new_flipped(),
            Polynomial::Sparse(sc) => sc.new_flipped(),
        }
    }
}

//********** Methods with Euclidean Ring Coefficients **********/
/// x * y / z
pub(crate) fn mul_div_uint<C>(x: C, y: usize, z: C) -> C
        where C: EuclideanRing + num::FromPrimitive + num::Integer + Clone {

    let gcd_xz = x.gcd(&z);
    let x_red = x / gcd_xz.clone();
    let z_red = z / gcd_xz;
    let y_red = C::from_usize(y).unwrap() / z_red;
    x_red * y_red
}

impl<C> Polynomial<C> where C: EuclideanRing + num::FromPrimitive + num::Integer + Clone {

    /// Shifts the `self` polynomial by `h`, that is, substitutes `x + h` into `x`, and expands.
    /// The coefficient type must be Euclidean ring (and num::Integer).
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p = dense![0, 0, 0, 1];  // x³
    ///     p.shift(1);
    ///     assert_eq!(p, dense![1, 3, 3, 1]);  // (x + 1)³ = 1 + 3x + 3x² + x³
    /// 
    /// The `shift()` method returns the same result as `compose()` with `x + h`.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     use comonjo_algebra::poly::Compose;
    /// 
    ///     let mut p = dense![0, 0, 0, 3];  // x³
    ///     let q = p.clone();
    /// 
    ///     p.shift(1);
    ///     assert_eq!(p, q.compose(Polynomial::x() + 1));
    /// 
    pub fn shift(&mut self, h: C) {
        if h.is_zero() { return; }
        match self {
            Polynomial::Dense(dc) => dc.shift(h),
            Polynomial::Sparse(sc) => sc.shift(h),
            _ => (),
        }
    }

    /// Returns a polynomial that the `self` is shifted by `h`, that is, `x + h` is substituted into `x`, and expanded.
    /// The coefficient type must be Euclidean ring (and num::Integer).
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p = dense![0, 0, 0, 1];  // x³
    ///     assert_eq!(p.new_shifted(1), dense![1, 3, 3, 1]);  // (x + 1)³ = 1 + 3x + 3x² + x³
    /// 
    /// The `new_shifted()` method returns the same result as `compose()` with `x + h`.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     use comonjo_algebra::poly::Compose;
    /// 
    ///     let p = dense![0, 0, 0, 3];  // x³
    ///     assert_eq!(p.new_shifted(1), p.compose(Polynomial::x() + 1));
    /// 
    pub fn new_shifted(&self, h: C) -> Polynomial<C> {
        if h.is_zero() { return self.clone(); }
        match self {
            Polynomial::Dense(dc) => dc.new_shifted(h),
            Polynomial::Sparse(sc) => sc.new_shifted(h),
            _ => self.clone(),
        }
    }
}

impl<C> Polynomial<C> where C: Field + num::FromPrimitive + Clone {

    /// Shifts the `self` polynomial by `h`, that is, substitutes `x + h` into `x`, and expands.
    /// The coefficient type must implement `comonjo_algebra::algebra::Field`.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p: Polynomial<f64> = dense![0., 0., 0., 1.];  // x³
    ///     p.shift_f(1.);
    ///     assert_eq!(p, dense![1., 3., 3., 1.]);  // (x + 1)³ = 1 + 3x + 3x² + x³
    /// 
    /// The `shift_f()` method returns the same result as `compose()` with `x + h`.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     use comonjo_algebra::poly::Compose;
    /// 
    ///     let mut p: Polynomial<f64> = dense![0., 0., 0., 3.];  // x³
    ///     let q = p.clone();
    /// 
    ///     p.shift_f(1.);
    ///     assert_eq!(p, q.compose(Polynomial::x() + 1.));
    /// 
    pub fn shift_f(&mut self, h: C) {
        if h.is_zero() { return; }
        match self {
            Polynomial::Dense(dc) => dc.shift_f(h),
            Polynomial::Sparse(sc) => sc.shift_f(h),
            _ => (),
        }
    }

    /// Returns a polynomial that the `self` is shifted by `h`, that is, `x + h` is substituted into `x`, and expanded.
    /// The coefficient type must implement `comonjo_algebra::algebra::Field`.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p: Polynomial<f64> = dense![0., 0., 0., 1.];  // x³
    ///     assert_eq!(p.new_shifted_f(1.), dense![1., 3., 3., 1.]);  // (x + 1)³ = 1 + 3x + 3x² + x³
    /// 
    /// The `new_shifted_f()` method returns the same result as `compose()` with `x + h`.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     use comonjo_algebra::poly::Compose;
    /// 
    ///     let p: Polynomial<f64> = dense![0., 0., 0., 3.];  // x³
    ///     assert_eq!(p.new_shifted_f(1.), p.compose(Polynomial::x() + 1.));
    /// 
    pub fn new_shifted_f(&self, h: C) -> Polynomial<C> {
        if h.is_zero() { return self.clone(); }
        match self {
            Polynomial::Dense(dc) => dc.new_shifted_f(h),
            Polynomial::Sparse(sc) => sc.new_shifted_f(h),
            _ => self.clone(),
        }
    }
}

//********** Methods with Field Coefficients **********/


    //   /**
    //    * Returns the real roots of this polynomial.
    //    *
    //    * Depending on `C`, the `finder` argument may need to be passed "explicitly" via an implicit conversion. This is
    //    * because some types (eg `BigDecimal`, `Rational`, etc) require an error bound, and so provide implicit conversions
    //    * to `RootFinder`s from the error type. For instance, `BigDecimal` requires either a scale or MathContext. So, we'd
    //    * call this method with `poly.roots(MathContext.DECIMAL128)`, which would return a `Roots[BigDecimal` whose roots are
    //    * approximated to the precision specified in `DECIMAL128` and rounded appropriately.
    //    *
    //    * On the other hand, a type like `Double` doesn't require an error bound and so can be called without specifying the
    //    * `RootFinder`.
    //    *
    //    * @param finder
    //    *   a root finder to extract roots with
    //    * @return
    //    *   the real roots of this polynomial
    //    */
    //   def roots(implicit finder: RootFinder[C]): Roots[C] =
    //     finder.findRoots(this)
    // pub fn roots(self) -> Roots<C> {
    //     todo!()
    // }


impl<C> Polynomial<C> where C: Field {

    /// Makes the `self` polynomial be a monic polynomial, that is,
    /// a polynomial which is scaled for the leading term's coefficient to be equal to 1.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p: Polynomial<f64> = dense![1., 2., 3., 4.];  // 1 + 2x + 3x² + 4x³
    ///     p.monic();
    ///     assert_eq!(p, dense![1./4., 2./4., 3./4., 1.]);  // (1/4) + (1/2)x + (3/4)x² + x³
    /// 
    pub fn monic(&mut self) {
        match self {
            Polynomial::Zero() => (),
            Polynomial::Constant(_) => *self = Polynomial::new_raw_const(C::one()),
            Polynomial::Dense(dc) => dc.monic(),
            Polynomial::Sparse(sc) => sc.monic(),
        }
    }
}

impl<C> Polynomial<C> where C: Field + Clone {

    /// Returns a monic polynomial of the `self` polynomial, that is,
    /// a polynomial which is scaled for the leading term's coefficient to be equal to 1.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p: Polynomial<f64> = dense![1., 2., 3., 4.];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.new_integral(), dense![0., 1., 1., 1., 1.]);  // x + x² + x³ + x⁴
    /// 
    pub fn new_monic(&self) -> Polynomial<C> {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(_) => Polynomial::one(),
            Polynomial::Dense(dc) => dc.new_monic(),
            Polynomial::Sparse(sc) => sc.new_monic(),
        }
    }
}

impl<C> Polynomial<C> where C: Field + num::FromPrimitive + Clone + Debug{

    /// Integrates the `self` polynomial.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p: Polynomial<f64> = dense![1., 2., 3., 4.];  // 1 + 2x + 3x² + 4x³
    ///     p.integrate();
    ///     assert_eq!(p, dense![0., 1., 1., 1., 1.]);  // x + x² + x³ + x⁴
    /// 
    pub fn integrate(&mut self) {
        match self {
            Polynomial::Zero() => (),
            Polynomial::Constant(cc) => *self = Polynomial::linear_monomial(cc.0.clone()),
            Polynomial::Dense(dc) => dc.integrate(),
            Polynomial::Sparse(sc) => sc.integrate(),
        }
    }

    /// Returns the integral of the `self` polynomial.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p: Polynomial<f64> = dense![1., 2., 3., 4.];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(p.new_integral(), dense![0., 1., 1., 1., 1.]);  // x + x² + x³ + x⁴
    /// 
    pub fn new_integral(&self) -> Polynomial<C> {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::linear_monomial(cc.0.clone()),
            Polynomial::Dense(dc) => dc.new_integral(),
            Polynomial::Sparse(sc) => sc.new_integral(),
        }
    }

    /// integrates the `self` polynomial n times.
    ///
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let mut p: Polynomial<f64> = dense![1., 2., 3., 4.];  // 1 + 2x + 3x² + 4x³
    ///     p.integrate_n(2);
    ///     assert_eq!(
    ///         p,
    ///         dense![0., 0., 1./2., 1./3., 1./4., 1./5.]);  // (1/2)x² + (1/3)x³ + (1/4)x⁴ + (1/5)x⁵
    /// 
    pub fn integrate_n(&mut self, n: usize) {
        if n == 0 { return; }
        if n == 1 { 
            self.integrate();
            return;
        }

        match self {
            Polynomial::Zero() => (),
            Polynomial::Constant(cc) => {
                let f: C = factorial(n);
                *self = sparse![(n, cc.0.ref_div(f))]
            },
            Polynomial::Dense(dc) => dc.integrate_n(n),
            Polynomial::Sparse(sc) => sc.integrate_n(n),
        }
    }

    /// Returns the n-th integral of the `self` polynomial.
    /// 
    ///     # use comonjo_algebra::poly::Polynomial;
    ///     # use comonjo_algebra::dense;
    ///     let p: Polynomial<f64> = dense![1., 2., 3., 4.];  // 1 + 2x + 3x² + 4x³
    ///     assert_eq!(
    ///         p.new_nth_integral(2),
    ///         dense![0., 0., 1./2., 1./3., 1./4., 1./5.]);  // (1/2)x² + (1/3)x³ + (1/4)x⁴ + (1/5)x⁵
    /// 
    pub fn new_nth_integral(&self, n: usize) -> Polynomial<C> {
        if n == 0 { return self.clone(); }
        if n == 1 { return self.new_integral(); }

        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => {
                let f: C = factorial(n);
                sparse![(n, cc.0.clone() / f)]
            },
            Polynomial::Dense(dc) => dc.new_nth_integral(n),
            Polynomial::Sparse(sc) => sc.new_nth_integral(n),
        }
    }
}

//********** Operator Overloads **********
//********** Neg **********/
impl<C> Neg for Polynomial<C> where C: Ring {

    type Output = Polynomial<C>;

    fn neg(self) -> Self::Output {
        self.map_nonzero(|_, c| -c)
    }
}

impl<'a, C> Neg for &'a Polynomial<C> where C: Ring + Clone {

    type Output = Polynomial<C>;

    fn neg(self) -> Self::Output {
        self.map_nonzero(|_, c| c.ref_neg())
    }
}

//********** Add **********/
impl<C> Add for Polynomial<C> where C: Semiring {

    type Output = Polynomial<C>;

    fn add(self, other: Self) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => rhs,

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 + rhs.0),
            (Polynomial::Constant(c_lhs), Polynomial::Dense(mut dc)) => {
                dc.0.get_mut(0).map(|c_rhs| *c_rhs = c_lhs.0 + &*c_rhs);
                Polynomial::Dense(dc)
            },
            (Polynomial::Constant(c_lhs), Polynomial::Sparse(mut sc)) => {
                if let Some(c_rhs) = sc.0.get_mut(&0) {
                    let sum = c_lhs.0 + &*c_rhs;
                    if sum.is_zero() { sc.0.remove(&0); } else { *c_rhs = sum; }
                } else {
                    sc.0.insert(0, c_lhs.0);
                }
                Polynomial::Sparse(sc)
            },

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = c_lhs.ref_add(c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::add_vv(lhs, rhs),

            // Sparse
            (Polynomial::Sparse(mut sc), Polynomial::Constant(c_rhs)) => {
                if let Some(c_lhs) = sc.0.get_mut(&0) {
                    let sum = c_lhs.ref_add(c_rhs.0);
                    if sum.is_zero() { sc.0.remove(&0); } else { *c_lhs = sum; }
                } else {
                    sc.0.insert(0, c_rhs.0);
                }
                Polynomial::Sparse(sc)
            },
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::add_vv(lhs, rhs),
        }
    }
}

impl<'b, C> Add<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn add(self, other: &'b Self) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => rhs.clone(),

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 + &rhs.0),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_coeffs::add_vr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_coeffs::add_vr(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = c_lhs.ref_add(&c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::add_vr(lhs, rhs),

            // Sparse
            (Polynomial::Sparse(mut sc), Polynomial::Constant(c_rhs)) => {
                if let Some(c_lhs) = sc.0.get_mut(&0) {
                    let sum = c_lhs.ref_add(&c_rhs.0);
                    if sum.is_zero() { sc.0.remove(&0); } else { *c_lhs = sum; }
                } else {
                    sc.0.insert(0, c_rhs.0.clone());
                }
                Polynomial::Sparse(sc)
            },
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::add_vr(lhs, rhs),
        }
    }
}

impl<'a, C> Add<Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn add(self, other: Polynomial<C>) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => rhs,

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_add(rhs.0)),
            (Polynomial::Constant(c_lhs), Polynomial::Dense(mut dc)) => {
                dc.0.get_mut(0).map(|c_rhs| *c_rhs = c_lhs.0.ref_add(&*c_rhs));
                Polynomial::Dense(dc)
            },
            (Polynomial::Constant(c_lhs), Polynomial::Sparse(mut sc)) => {
                if let Some(c_rhs) = sc.0.get_mut(&0) {
                    let sum = c_lhs.0.ref_add(&*c_rhs);
                    if sum.is_zero() { sc.0.remove(&0); } else { *c_rhs = sum; }
                } else {
                    sc.0.insert(0, c_lhs.0.clone());
                }
                Polynomial::Sparse(sc)
            },

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::add_rv(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::add_rv(lhs, rhs),
        }
    }
}

impl<'a, 'b, C> Add<&'b Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn add(self, other: &'b Polynomial<C>) -> Polynomial<C> {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => rhs.clone(),

            // Const
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_add(&rhs.0)),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_coeffs::add_rr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_coeffs::add_rr(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::add_rr(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::add_rr(lhs, rhs),
        }
    }
}

//***** Add constant
impl<C> Add<C> for Polynomial<C> where C: Semiring {
    type Output = Polynomial<C>;
    fn add(self, rhs: C) -> Self::Output { self + Polynomial::constant(rhs) }
}

impl<'b, C> Add<&'b C> for Polynomial<C> where C: Semiring + Clone {
    type Output = Polynomial<C>;
    fn add(self, rhs: &'b C) -> Self::Output { self + Polynomial::constant(rhs.clone()) }
}

impl<'a, C> Add<C> for &'a Polynomial<C> where C: Semiring + Clone {
    type Output = Polynomial<C>;
    fn add(self, rhs: C) -> Self::Output { self + Polynomial::constant(rhs) }
}

impl<'a, 'b, C> Add<&'b C> for &'a Polynomial<C> where C: Semiring + Clone {
    type Output = Polynomial<C>;
    fn add(self, rhs: &'b C) -> Polynomial<C> { self + Polynomial::constant(rhs.clone())  }
}

macro_rules! impl_add_to_const {
    ( $( $t:ident ),* ) => {
        $(
            impl Add<Polynomial<$t>> for $t {
                type Output = Polynomial<$t>;
                fn add(self, rhs: Polynomial<$t>) -> Self::Output { Polynomial::constant(self) + rhs }
            }

            impl<'b> Add<&'b Polynomial<$t>> for $t {
                type Output = Polynomial<$t>;
                fn add(self, rhs: &'b Polynomial<$t>) -> Self::Output { Polynomial::constant(self) + rhs }
            }

            impl<'a> Add<Polynomial<$t>> for &'a $t {
                type Output = Polynomial<$t>;
                fn add(self, rhs: Polynomial<$t>) -> Self::Output { Polynomial::constant(self.clone()) + rhs }
            }

            impl<'a, 'b> Add<&'b Polynomial<$t>> for &'a $t {
                type Output = Polynomial<$t>;
                fn add(self, rhs: &'b Polynomial<$t>) -> Self::Output { Polynomial::constant(self.clone()) + rhs }
            }
        )*
    };
}

impl_add_to_const!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64,
        BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);


//********** Sub **********/
impl<C> Sub for Polynomial<C> where C: Ring {

    type Output = Polynomial<C>;

    fn sub(self, other: Self) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => -rhs,

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 - rhs.0),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_coeffs::sub_vv(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_coeffs::sub_vv(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = (*c_lhs).ref_sub(c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::sub_vv(lhs, rhs),

            // Sparse
            (Polynomial::Sparse(mut sc), Polynomial::Constant(c_rhs)) => {
                if let Some(c_lhs) = sc.0.get_mut(&0) {
                    let dif = c_lhs.ref_sub(c_rhs.0);
                    if dif.is_zero() { sc.0.remove(&0); } else { *c_lhs = dif; }
                } else {
                    sc.0.insert(0, -c_rhs.0);
                }
                Polynomial::Sparse(sc)
            },
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::sub_vv(lhs, rhs),
        }
    }
}

impl<'b, C> Sub<&'b Polynomial<C>> for Polynomial<C> where C: Ring + Clone {

    type Output = Polynomial<C>;

    fn sub(self, other: &'b Self) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => rhs.ref_neg(),

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 - &rhs.0),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_coeffs::sub_vr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_coeffs::sub_vr(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = c_lhs.ref_sub(&c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::sub_vr(lhs, rhs),

            // Sparse
            (Polynomial::Sparse(mut sc), Polynomial::Constant(c_rhs)) => {
                if let Some(c_lhs) = sc.0.get_mut(&0) {
                    let dif = c_lhs.ref_sub(&c_rhs.0);
                    if dif.is_zero() { sc.0.remove(&0); } else { *c_lhs = dif; }
                } else {
                    sc.0.insert(0, c_rhs.0.ref_neg());
                }
                Polynomial::Sparse(sc)
            },
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::sub_vr(lhs, rhs),
        }
    }
}

impl<'a, C> Sub<Polynomial<C>> for &'a Polynomial<C> where C: Ring + Clone {

    type Output = Polynomial<C>;

    fn sub(self, other: Polynomial<C>) -> Self::Output {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => -rhs,

            // Constant
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_sub(rhs.0)),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_coeffs::sub_rv(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_coeffs::sub_rv(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::sub_rv(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::sub_rv(lhs, rhs),
        }
    }
}

impl<'a, 'b, C> Sub<&'b Polynomial<C>> for &'a Polynomial<C> where C: Ring + Clone {

    type Output = Polynomial<C>;

    fn sub(self, other: &'b Polynomial<C>) -> Polynomial<C> {
        match (self, other) {
            // Zero
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => rhs.ref_neg(),

            // Const
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_sub(&rhs.0)),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense_coeffs::sub_rr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse_coeffs::sub_rr(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::sub_rr(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::sub_rr(lhs, rhs),
        }
    }
}


//***** Sub constant
impl<C> Sub<C> for Polynomial<C> where C: Ring {
    type Output = Polynomial<C>;
    fn sub(self, rhs: C) -> Self::Output { self - Polynomial::constant(rhs) }
}

impl<'b, C> Sub<&'b C> for Polynomial<C> where C: Ring + Clone {
    type Output = Polynomial<C>;
    fn sub(self, rhs: &'b C) -> Self::Output { self - Polynomial::constant(rhs.clone()) }
}

impl<'a, C> Sub<C> for &'a Polynomial<C> where C: Ring + Clone {
    type Output = Polynomial<C>;
    fn sub(self, rhs: C) -> Self::Output { self - Polynomial::constant(rhs) }
}

impl<'a, 'b, C> Sub<&'b C> for &'a Polynomial<C> where C: Ring + Clone {
    type Output = Polynomial<C>;
    fn sub(self, rhs: &'b C) -> Polynomial<C> { self - Polynomial::constant(rhs.clone())  }
}

macro_rules! impl_sub_from_const {
    ( $( $t:ident ),* ) => {
        $(
            impl Sub<Polynomial<$t>> for $t {
                type Output = Polynomial<$t>;
                fn sub(self, rhs: Polynomial<$t>) -> Self::Output { Polynomial::constant(self) - rhs }
            }

            impl<'b> Sub<&'b Polynomial<$t>> for $t {
                type Output = Polynomial<$t>;
                fn sub(self, rhs: &'b Polynomial<$t>) -> Self::Output { Polynomial::constant(self) - rhs }
            }

            impl<'a> Sub<Polynomial<$t>> for &'a $t {
                type Output = Polynomial<$t>;
                fn sub(self, rhs: Polynomial<$t>) -> Self::Output { Polynomial::constant(self.clone()) - rhs }
            }

            impl<'a, 'b> Sub<&'b Polynomial<$t>> for &'a $t {
                type Output = Polynomial<$t>;
                fn sub(self, rhs: &'b Polynomial<$t>) -> Self::Output { Polynomial::constant(self.clone()) - rhs }
            }
        )*
    };
}

impl_sub_from_const!(isize, i8, i16, i32, i64, i128, f32, f64, 
        BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

//********** Mul **********/
impl<C> Mul<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: Self) -> Self::Output {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), rhs) => {
                if lhs.0.is_one() { return rhs }
                rhs.map_nonzero(|_, c| lhs.0.ref_mul(c))
            },
            (lhs, Polynomial::Constant(rhs)) => {
                if rhs.0.is_one() { return lhs }
                lhs.map_nonzero(|_, c| c * &rhs.0)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::mul(&lhs, &rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::mul(&lhs, &rhs),
        }
    }
}

impl<'b, C> Mul<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: &'b Self) -> Self::Output {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), rhs) => {
                if lhs.0.is_one() { return rhs.clone() }
                rhs.map_nonzero(|_, c| lhs.0.ref_mul(c))
            },
            (lhs, Polynomial::Constant(rhs)) => {
                if rhs.0.is_one() { return lhs }
                lhs.map_nonzero(|_, c| c * &rhs.0)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::mul(&lhs, rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::mul(&lhs, rhs),
        }
    }
}

impl<'a, C> Mul<Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: Polynomial<C>) -> Self::Output {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), rhs) => {
                if lhs.0.is_one() { return rhs }
                rhs.map_nonzero(|_, c| lhs.0.ref_mul(c))
            },
            (lhs, Polynomial::Constant(rhs)) => {
                if rhs.0.is_one() { return lhs.clone() }
                lhs.map_nonzero(|_, c| c.ref_mul(&rhs.0))
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::mul(lhs, &rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::mul(lhs, &rhs),
        }
    }
}

impl<'a, 'b, C> Mul<&'b Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;
    
    fn mul(self, other: &'b Polynomial<C>) -> Polynomial<C> {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), rhs) => {
                if lhs.0.is_one() { return rhs.clone() }
                rhs.map_nonzero(|_, c| lhs.0.ref_mul(c))
            },
            (lhs, Polynomial::Constant(rhs)) => {
                if rhs.0.is_one() { return lhs.clone() }
                lhs.map_nonzero(|_, c| c.ref_mul(&rhs.0))
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::mul(lhs, rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::mul(lhs, rhs),
        }
    }
}

//***** Multiply by constant
impl<C> Mul<C> for Polynomial<C> where C: Semiring + Clone {
    type Output = Polynomial<C>;
    fn mul(self, k: C) -> Self::Output { self * &k }
}

impl<'b, C> Mul<&'b C> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, k: &'b C) -> Self::Output {
        if k.is_zero() { return Polynomial::Zero() }
        if k.is_one() { return self }
        self.map_nonzero(|_, c| c.ref_mul(k))
    }
}

impl<'a, C> Mul<C> for &'a Polynomial<C> where C: Semiring + Clone {
    type Output = Polynomial<C>;
    fn mul(self, k: C) -> Self::Output { self * &k }
}

impl<'a, 'b, C> Mul<&'b C> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;
    
    fn mul(self, k: &'b C) -> Polynomial<C> {
        if k.is_zero() { return Polynomial::Zero() }
        if k.is_one() { return self.clone() }
        self.map_nonzero(|_, c| c.ref_mul(k))
    }
}

macro_rules! impl_mul_to_const {
    ( $( $t:ident ),* ) => {
        $(
            impl Mul<Polynomial<$t>> for $t {
                type Output = Polynomial<$t>;
                fn mul(self, rhs: Polynomial<$t>) -> Self::Output { rhs.scale_by_left(&self) }
            }

            impl<'b> Mul<&'b Polynomial<$t>> for $t {
                type Output = Polynomial<$t>;
                fn mul(self, rhs: &'b Polynomial<$t>) -> Self::Output { rhs.ref_scale_by_left(&self) }
            }

            impl<'a> Mul<Polynomial<$t>> for &'a $t {
                type Output = Polynomial<$t>;
                fn mul(self, rhs: Polynomial<$t>) -> Self::Output { rhs.scale_by_left(self) }
            }

            impl<'a, 'b> Mul<&'b Polynomial<$t>> for &'a $t {
                type Output = Polynomial<$t>;
                fn mul(self, rhs: &'b Polynomial<$t>) -> Self::Output { rhs.ref_scale_by_left(self) }
            }
        )*
    };
}

impl_mul_to_const!(usize, u8, u16, u32, u64, u128, isize, i8, i16, i32, i64, i128, f32, f64,
        BigUint, BigInt, Rational32, Rational64, BigRational, Complex32, Complex64);

//********** Div & Rem **********/
fn panic_to_divide_by_zero<E>() -> E { panic!("Can't divide by zero!") }

impl<C> Polynomial<C> where C: Field + Clone {

    pub fn euclidean_fn(&self) -> usize { self.degree() }

    #[inline]
    fn div_rem_val<'a>(self, other: &'a Polynomial<C>) -> (Self, Self) {
        match (self, other) {
            (_, Polynomial::Zero()) => panic_to_divide_by_zero(),
            (Polynomial::Zero(), _) => (Polynomial::Zero(), Polynomial::Zero()),
            (lhs, Polynomial::Constant(rhs)) => (lhs.map_nonzero(|_, c| c / (&rhs.0)), Polynomial::Zero()),
            (lhs @ Polynomial::Constant(_), _) => (Polynomial::Zero(), lhs),
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::div_rem(lhs.to_vec(), rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::div_rem(lhs.to_map(), rhs),
        }
    }

    #[inline]
    fn div_rem_ref(&self, other: &Self) -> (Self, Self) {
        match (self, other) {
            (_, Polynomial::Zero()) => panic_to_divide_by_zero(),
            (Polynomial::Zero(), _) => (Polynomial::Zero(), Polynomial::Zero()),
            (lhs, Polynomial::Constant(rhs)) => (lhs.map_nonzero(|_, c| c.ref_div(&rhs.0)), Polynomial::Zero()),
            (lhs @ Polynomial::Constant(_), _) => (Polynomial::Zero(), lhs.clone()),
            (lhs @ Polynomial::Dense(_), rhs) => dense_coeffs::div_rem(lhs.to_vec(), rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse_coeffs::div_rem(lhs.to_map(), rhs),
        }
    }
}

impl<C> Euclid for Polynomial<C> where C: Field + Clone {

    fn div_euclid(&self, other: &Self) -> Self { self / other }
    fn rem_euclid(&self, other: &Self) -> Self { self % other }
    fn div_rem_euclid(&self, other: &Self) -> (Self, Self) { self.div_rem_ref(other) }
}


impl<C> Div for Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn div(self, other: Self) -> Self::Output { self.div_rem_val(&other).0 }
}

impl<'b, C> Div<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn div(self, other: &'b Polynomial<C>) -> Self::Output { self.div_rem_ref(other).0 }
}

impl<'a, C> Div<Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn div(self, other: Polynomial<C>) -> Self::Output { self.div_rem_ref(&other).0 }
}

impl<'a, 'b, C> Div<&'b Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn div(self, other: &'b Polynomial<C>) -> Self::Output { self.div_rem_ref(other).0 }
}


impl<C> Rem for Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn rem(self, other: Self) -> Self::Output { self.div_rem_val(&other).1 }
}

impl<'b, C> Rem<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn rem(self, other: &'b Polynomial<C>) -> Self::Output { self.div_rem_ref(other).1 }
}

impl<'a, C> Rem<Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn rem(self, other: Polynomial<C>) -> Self::Output { self.div_rem_ref(&other).1 }
}

impl<'a, 'b, C> Rem<&'b Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn rem(self, other: &'b Polynomial<C>) -> Self::Output { self.div_rem_ref(other).1 }
}

//***** Divide by constant
impl<C> Div<C> for Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn div(self, k: C) -> Self::Output { self / &k }
}

impl<'b, C> Div<&'b C> for Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn div(self, k: &'b C) -> Self::Output {
        if k.is_zero() { panic_to_divide_by_zero() }
        if k.is_one() { return self }
        self.map_nonzero(|_, c| c.ref_div(k))
    }
}

impl<'a, C> Div<C> for &'a Polynomial<C> where C: Field + Clone {
    type Output = Polynomial<C>;
    fn div(self, k: C) -> Self::Output { self / &k }
}

impl<'a, 'b, C> Div<&'b C> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;
    
    fn div(self, k: &'b C) -> Polynomial<C> {
        if k.is_zero() { panic_to_divide_by_zero() }
        if k.is_one() { return self.clone() }
        self.map_nonzero(|_, c| c.ref_div(k))
    }
}

//********* Pow **********/
fn calc_pow<C>(base: &Polynomial<C>, p: u32, extra: &Polynomial<C>) -> Polynomial<C> where C: Semiring + Clone {
    if p == 1 {
        base * extra
    } else {
        let next_extra = if (p & 1) == 1 { base * extra } else { extra.clone() };
        calc_pow(&(base * base), p >> 1, &next_extra)
    }
}

/// Note that 0^0 returns 1 for simplicity.
impl<C> Pow<u32> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn pow(self, power: u32) -> Self::Output {
        match power {
            0 => Polynomial::one(),
            1 => self,
            _ => calc_pow(&self, power-1, &self)
        }
    }
}

/// Note that 0^0 returns 1 for simplicity.
impl<'a, C> Pow<u32> for &'a Polynomial<C> where C: Semiring + Clone{

    type Output = Polynomial<C>;

    fn pow(self, power: u32) -> Polynomial<C> {
        match power {
            0 => Polynomial::one(),
            1 => self.clone(),
            _ => calc_pow(self, power-1, self)
        }
    }
}

//********** Implementation of Algebra *********/
impl<C> RefAdd<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_add(&self, other: Polynomial<C>) -> Polynomial<C> { self + other }
}

impl<'b, C> RefAdd<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_add(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self + other }
}

impl<C> AdditiveSemigroup for Polynomial<C> where C: Semiring + Clone {}
impl<C> AdditiveMonoid for Polynomial<C> where C: Semiring + Clone {}


impl<C> RefSub<Polynomial<C>> for Polynomial<C> where C: Ring + Clone {
    #[inline]
    fn ref_sub(&self, other: Polynomial<C>) -> Polynomial<C> { self + other }
}

impl<'b, C> RefSub<&'b Polynomial<C>> for Polynomial<C> where C: Ring + Clone {
    #[inline]
    fn ref_sub(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self + other }
}


impl<C> RefMul<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_mul(&self, other: Polynomial<C>) -> Polynomial<C> { self * other }
}

impl<C> AdditiveGroup for Polynomial<C> where C: Ring + Clone {
    fn ref_neg(&self) -> Self { -self }
}


impl<'b, C> RefMul<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_mul(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self * other }
}

impl<C> Semigroup for Polynomial<C> where C: Semiring + Clone {}
impl<C> Monoid for Polynomial<C> where C: Semiring + Clone {}

impl<C> Semiring for Polynomial<C> where C: Ring + Clone {}
impl<C> Ring for Polynomial<C> where C: Ring + Clone{}


impl<C> RefDiv<Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_div(&self, other: Polynomial<C>) -> Polynomial<C> { self / other }
}

impl<'b, C> RefDiv<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_div(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self / other }
}


impl<C> RefRem<Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_rem(&self, other: Polynomial<C>) -> Polynomial<C> { self % other }
}

impl<'b, C> RefRem<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_rem(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self % other }
}

impl<C> EuclideanRing for Polynomial<C> where C: Field + Clone {

    fn div_rem(self, other: Self) -> (Self, Self) { 
        self.div_rem_euclid(&other)
    }

    fn ref_div_rem(&self, other: &Self) -> (Self, Self) {
        self.div_rem_euclid(other)
    }
}