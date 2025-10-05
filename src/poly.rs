pub(crate) mod dense;
pub(crate) mod sparse;
pub(crate) mod iter;
pub mod term;

use std::{collections::BTreeMap, fmt::{Debug, Display}, ops::*};
use iter::{CoeffsIter, IntoCoeffsIter, IntoNonzeroCoeffsIter, NonzeroCoeffsIter};
use num::{pow::Pow, traits::{ConstOne, ConstZero, Euclid}, One, Zero};

use dense::DenseContent;
use sparse::SparseContent;
use term::Term;

use crate::algebra::*;

/// Polynomial type.
/// 
/// Refer to spire's 
/// <a href="https://github.com/typelevel/spire/blob/main/core/src/main/scala/spire/math/Polynomial.scala">Polynomial</a>.
/// 
/// <code>Polynomial</code>'s <code>into_iter()</code> method returns <code>Iterator</code> that iterates nonzero coefficients
/// (with its term's degree). 
/// <code>Polynomial</code> and <code>&Polynomial</code> implement the <code>IntoIterator</code> trait respectively,
/// so those methods behave differently:
/// 
///     # use comonjo_algebra::poly::Polynomial;
///     # use comonjo_algebra::dense;
/// 
///     let p = dense![4, 5, 0, 6];  // 4 + 5x + 6x³
/// 
///     // &Polynomial
///     let mut ite0 = (&p).into_iter();
///     assert_eq!(ite0.next(), Some((0, &4)));
///     assert_eq!(ite0.next(), Some((1, &5)));
///     assert_eq!(ite0.next(), Some((3, &6)));
///     assert_eq!(ite0.next(), None);
/// 
///     // Polynomial
///     let mut ite1 = p.into_iter();
///     assert_eq!(ite1.next(), Some((0, 4)));
///     assert_eq!(ite1.next(), Some((1, 5)));
///     assert_eq!(ite1.next(), Some((3, 6)));
///     assert_eq!(ite1.next(), None);
/// 
/// To iterate coefficients including zero, import <code>CoeffsIterator</code> and
/// use <code>coeffs()</code> methods:
/// 
///     # use comonjo_algebra::poly::Polynomial;
///     # use comonjo_algebra::dense;
///     use comonjo_algebra::poly::CoeffsIterator;
/// 
///     let p = dense![4, 5, 0, 6];  // 4 + 5x + 6x³
///     let mut ite = p.coeffs();
///     assert_eq!(ite.next(), Some(4));
///     assert_eq!(ite.next(), Some(5));
///     assert_eq!(ite.next(), Some(0));
///     assert_eq!(ite.next(), Some(6));
///     assert_eq!(ite.next(), None);
/// 
/// <code>nonzero_coeffs()</code> methods are equivalent to <code>into_iter()</code>.
pub enum Polynomial<C> where C: Semiring {

    Zero(),

    /// Use <code>Polynomial::constant()</code> function to instantiate.
    Constant(ConstContent<C>),

    /// Use <code>dense!()</code> macro to instantiate.
    Dense(DenseContent<C>),
    
    /// Use <code>sparse!()</code> macro to instantiate.
    Sparse(SparseContent<C>)
}

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

impl<C> Polynomial<C> where C: Semiring {

    pub fn constant(c: C) -> Polynomial<C> {
        if c.is_zero() {
            Polynomial::Zero()
        } else {
            Polynomial::Constant(ConstContent(c))
        }
    }

    pub fn dense_from_vec(mut coeffs: Vec<C>) -> Polynomial<C>  {

        while let Some(c) = coeffs.last() {
            if c.is_zero() {
                coeffs.pop();
            } else {
                break;
            }
        }

        match coeffs.len() {
            0 => Polynomial::Zero(),
            1 => Polynomial::constant(coeffs.pop().unwrap()),
            _ => {
                coeffs.shrink_to_fit();
                Polynomial::Dense(DenseContent(coeffs))
            }
        }
    }

    pub fn sparse_from_map(mut coeffs: BTreeMap<usize, C>) -> Polynomial<C> {

        coeffs.retain(|_, v| !v.is_zero());

        match coeffs.len() {
            0 => Polynomial::Zero(),
            1 => {
                if coeffs.contains_key(&0) {
                    Polynomial::constant(coeffs.remove(&0).unwrap())
                } else {
                    Polynomial::Sparse(SparseContent(coeffs))
                }
            },
            _ => Polynomial::Sparse(SparseContent(coeffs))
        }
    }

    pub fn parse(s: &str) -> Polynomial<C> {
    
//    private[this] val termRe = "([0-9]+\\.[0-9]+|[0-9]+/[0-9]+|[0-9]+)?(?:([a-z])(?:\\^([0-9]+))?)?".r
 
//    private[this] val operRe = " *([+-]) *".r
 
//    private[spire] def parse(s: String): Polynomial[Rational] = {
 
//      // represents a term, plus a named variable v
//      case class T(c: Rational, v: String, e: Int)
 
//      // parse all the terms and operators out of the string
//      @tailrec def parse(s: String, ts: List[T]): List[T] =
//        if (s.isEmpty) {
//          ts
//        } else {
//          val (op, s2) = operRe.findPrefixMatchOf(s) match {
//            case Some(m) => (m.group(1), s.substring(m.end))
//            case None    => if (ts.isEmpty) ("+", s) else throw new IllegalArgumentException(s)
//          }
 
//          val m2 = termRe.findPrefixMatchOf(s2).getOrElse(throw new IllegalArgumentException(s2))
//          val c0 = Option(m2.group(1)).getOrElse("1")
//          val c = if (op == "-") "-" + c0 else c0
//          val v = Option(m2.group(2)).getOrElse("")
//          val e0 = Option(m2.group(3)).getOrElse("")
//          val e = if (e0 != "") e0 else if (v == "") "0" else "1"
 
//          val t =
//            try {
//              T(Rational(c), v, e.toInt)
//            } catch {
//              case _: Exception => throw new IllegalArgumentException(s"illegal term: $c*x^$e")
//            }
//          parse(s2.substring(m2.end), if (t.c == 0) ts else t :: ts)
//        }
 
//      // do some pre-processing to remove whitespace/outer parens
//      val t = s.trim
//      val u = if (t.startsWith("(") && t.endsWith(")")) t.substring(1, t.length - 1) else t
//      val v = Term.removeSuperscript(u)
 
//      // parse out the terms
//      val ts = parse(v, Nil)
 
//      // make sure we have at most one variable
//      val vs = ts.view.map(_.v).toSet.filter(_ != "")
//      if (vs.size > 1) throw new IllegalArgumentException("only univariate polynomials supported")
 
//      // we're done!
//      ts.foldLeft(Polynomial.zero[Rational])((a, t) => a + Polynomial(t.c, t.e))
//    }
        todo!()
    }

    //********** METHODS ******
    pub fn is_constant(&self) -> bool {
        match self {
            Polynomial::Zero() => true,
            Polynomial::Constant(_) => true,
            Polynomial::Dense(_) => false,
            Polynomial::Sparse(_) => false,
        }
    }

    pub fn degree(&self) -> usize {
        match self {
            Polynomial::Dense(dc) => dc.degree(),
            Polynomial::Sparse(sc) => sc.degree(),
            _ => 0,
        }
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        match self {
            Polynomial::Zero() => None,
            Polynomial::Constant(cc) => if n == 0 { Some(&cc.0) } else { None },
            Polynomial::Dense(dc) => dc.nth(n),
            Polynomial::Sparse(sc) => sc.nth(n),
        }
    }
    
    /// Note that <code>f</code> is applied only to nonzero coefficients
    pub fn map<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::constant(f(0, cc.0)),
            Polynomial::Dense(dc) => dc.map(f),
            Polynomial::Sparse(sc) => sc.map(f),
        }
    }
    
    /// Note that <code>f</code> is applied only to nonzero coefficients
    pub fn map_ref<D, F>(&self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, &C) -> D {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::constant(f(0, &cc.0)),
            Polynomial::Dense(dc) => dc.map_ref(f),
            Polynomial::Sparse(sc) => sc.map_ref(f),
        }
    }
     
    /// Return dense expression of this polynomial if this is sparse, otherwise (zero or constant) self.
    pub fn to_dense(self) -> Polynomial<C> {
        match self {
            Polynomial::Sparse(sc) => Polynomial::dense_from_vec(sc.to_vec()),
            _ => self,
        }
    }

    /// Return sparse expression of this polynomial if this is dense, otherwise (zero or constant) self.
    pub fn to_sparse(self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => Polynomial::sparse_from_map(dc.to_map()),
            _ => self
        }
    }
    
    pub fn to_vec(self) -> Vec<C> {
        match self {
            Polynomial::Zero() => Vec::new(),
            Polynomial::Constant(cc) => vec![cc.0],
            Polynomial::Dense(dc) => dc.0,
            Polynomial::Sparse(sc) => sc.to_vec(),
        }
    }
    
    pub fn to_map(self) -> BTreeMap<usize, C> {
        match self {
            Polynomial::Zero() => BTreeMap::new(),
            Polynomial::Constant(cc) => BTreeMap::from([(0, cc.0)]),
            Polynomial::Dense(dc) => dc.to_map(),
            Polynomial::Sparse(sc) => sc.0,
        }
    }


    // the following methods are used?
    pub fn max_order_term_coeff(&self) -> Option<&C> {
        match self {
            Polynomial::Zero() => None,
            _ => self.nth(self.degree())
        }
    }

    // pub fn split(&self) -> (Vec<usize>, Vec<&C>) {
    //     let n = self.degree();
    //     let mut es = Vec::with_capacity(n);
    //     let mut cs = Vec::with_capacity(n);

    //     self.nonzero_coeffs().for_each(|(e, c)| {
    //         es.push(e);
    //         cs.push(c);
    //     });

    //     (es, cs)
    // }

    /// Returns a polynomial with the max term removed.
    pub fn reductum(&self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => dc.reductum(),
            Polynomial::Sparse(sc) => sc.reductum(),
            _ => Polynomial::Zero(),
        }
    }
    
    pub fn derivative(&self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => dc.derivative(),
            Polynomial::Sparse(sc) => sc.derivative(),
            _ => Polynomial::Zero(),
        }
    }

    //********** Some factory methods **********/
    pub fn x() -> Polynomial<C> {
        sparse![(1, C::one())]
    }

    pub fn two_x() -> Polynomial<C> {
        sparse![(1, C::one() + C::one())]
    }

    /// Create a linear polynomial <i>ax</i>.
    pub fn linear_monomial(a: C) -> Polynomial<C> {
        sparse![(1, a)]
    }

    /// Create a linear polynomial <i>ax + b</i>
    pub fn linear(c1: C, c0: C) -> Polynomial<C> {
        sparse![(1, c1), (0, c0)]
    }

    /// Create a quadratic polynomial <i>ax<sup>2</sup><i>
    pub fn quadratic_monomial(a: C) -> Polynomial<C> {
        sparse![(2, a)]
    }

    /// Create a quadratic polynomial <i>ax<sup>2</sup> + bx + c<i>
    pub fn quadratic(a: C, b: C, c: C) -> Polynomial<C> {
        sparse![(2, a), (1, b), (0, c)]
    }

    /// Create a cubic polynomial <i>ax<sup>3</sup> + bx<sup>2</sup> + cx + d<i>
    pub fn cubic_monomial(a: C) -> Polynomial<C> {
        sparse![(3, a)]
    }

    /// Create a cubic polynomial <i>ax<sup>3</sup> + bx<sup>2</sup> + cx + d<i>
    pub fn cubic(a: C, b: C, c: C, d: C) -> Polynomial<C> {
        sparse![(3, a), (2, b), (1, c), (0, d)]
    }
}

pub struct ConstContent<C: Semiring>(pub(crate) C);

impl<C> Clone for Polynomial<C> where C: Semiring + Clone {

    fn clone(&self) -> Self {
        match self {
            Self::Zero() => Self::Zero(),
            Self::Constant(cc) => Self::Constant(ConstContent(cc.0.clone())),  // raw creation
            Self::Dense(dc) => Self::Dense(dc.clone()),
            Self::Sparse(sc) => Self::Sparse(sc.clone()),
        }
    }
}

impl<C> Polynomial<C> where C: Semiring + Clone {

    /// Create dense clone of this polynomial if the self is sparse, otherwise (zero, constant or dense) return self.
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
                Polynomial::dense_from_vec(v)
            },
            _ => self.clone()
        }
    }

    /// Create sparse clone of this polynomial if self is dense, otherwise (zero, constant or sparse) return self.
    pub fn sparse_clone(&self) -> Polynomial<C> {
        match self {
            d @ Polynomial::Dense(_) => {
                let mut map = BTreeMap::new();
                for (i, c) in d.nonzero_coeffs() {
                    map.insert(i, c.clone());  // note c is not zero
                }
                Polynomial::sparse_from_map(map)
            },
            _ => self.clone()
        }
    }
    
    pub fn clone_to_vec(&self) -> Vec<C> {
        self.coeffs().map(|option_c| match option_c {
            Some(c) => c.clone(),
            None => C::zero(),
        }).collect()
    }
    
    pub fn clone_to_map(&self) -> BTreeMap<usize, C> {
        self.nonzero_coeffs().map(|e| (e.0, e.1.clone())).collect()
    }

    pub fn terms<'a>(&'a self) -> impl Iterator<Item=Term<'a, C>> {
        self.nonzero_coeffs().map(|(i, c)| Term::from_ref(i, c))
    }

    pub fn max_term<'a>(&'a self) -> Term<'a, C> {
        match self.max_order_term_coeff() {
            Some(c) => Term::from_ref(self.degree(), c),
            _ => Term::from_value(0, C::zero()),
        }
    }

    pub fn min_term<'a>(&'a self) -> Term<'a, C> {
        match self.terms().next() {
            Some(t) => t,
            _ => Term::from_value(0, C::zero()),
        }
    }

    /// Returns this polynomial as a monic polynomial, where the leading coefficient (ie. `maxOrderTermCoeff`) is 1.
    pub fn monic(&self) -> Polynomial<C> {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(_) => Polynomial::one(),
            Polynomial::Dense(dc) => dc.monic(),
            Polynomial::Sparse(sc) => sc.monic(),
        }
    }

    pub fn integral(&self) -> Polynomial<C> {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::linear_monomial(cc.0.clone()),
            Polynomial::Dense(dc) => dc.integral(),
            Polynomial::Sparse(sc) => sc.integral(),
        }
    }
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
    
    //   /**
    //    * Evaluate the polynomial at `x`.
    //    */
    //   def apply(x: C)(implicit r: Semiring[C]): C
    
    //   def evalWith[A: Semiring: Eq: ClassTag](x: A)(f: C => A): A =
    //     this.map(f).apply(x)
    
    //   /**
    //    * Compose this polynomial with another.
    //    */
    //   def compose(y: Polynomial[C])(implicit ring: Rig[C], eq: Eq[C]): Polynomial[C] = {
    //     var polynomial: Polynomial[C] = Polynomial.zero[C]
    //     foreachNonzero { (e, c) =>
    //       val z: Polynomial[C] = y.pow(e) :* c
    //       polynomial = polynomial + z
    //     }
    //     polynomial
    //   }


    
    //   def interpolate[C: Field: Eq: ClassTag](points: (C, C)*): Polynomial[C] = {
    //     def loop(p: Polynomial[C], xs: List[C], pts: List[(C, C)]): Polynomial[C] =
    //       pts match {
    //         case Nil =>
    //           p
    //         case (x, y) :: tail =>
    //           val c = Polynomial.constant((y - p(x)) / xs.map(x - _).qproduct)
    //           val prod = xs.foldLeft(Polynomial.one[C]) { (prod, xn) =>
    //             prod * (Polynomial.x[C] - constant(xn))
    //           }
    //           loop(p + c * prod, x :: xs, tail)
    //       }
    //     loop(Polynomial.zero[C], Nil, points.toList)
    //   }
    // }

}
 
//********** Display and Debug **********
impl<C> Display for Polynomial<C> where C: Semiring + Display {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        match self {
            Polynomial::Zero() => f.write_str("(0)"),
            Polynomial::Constant(c) => f.write_fmt(format_args!("({})", c.0)),
            _ => {
                let s: String = self.nonzero_coeffs().map(|(i, c)|term::term_to_string(i, c)).collect();
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

//********** Iterator **********
/// Iterate nonzero coefficients
impl<C> IntoIterator for Polynomial<C> where C: Semiring {

    type Item = (usize, C);
    type IntoIter = IntoNonzeroCoeffsIter<C>;

    /// Return <code>impl Iterator<Item=(usize, C)></code>.
    fn into_iter(self) -> Self::IntoIter {
        match self {
            Polynomial::Zero() => IntoNonzeroCoeffsIter::Zero(),
            Polynomial::Constant(cc) => IntoNonzeroCoeffsIter::Constant(Some(cc.0)),
            Polynomial::Dense(dc) => IntoNonzeroCoeffsIter::Dense(dc.0.into_iter().enumerate()),
            Polynomial::Sparse(sc) => IntoNonzeroCoeffsIter::Sparse(sc.0.into_iter()),
        }
    }
}

/// Iterate reference of nonzero coefficient
impl<'a, C> IntoIterator for &'a Polynomial<C> where C: Semiring {

    type Item = (usize, &'a C);
    type IntoIter = NonzeroCoeffsIter<'a, C>;

    /// Iterate nonzero coefficients.
    /// Return <code>impl Iterator<Item=(usize, &C)></code>.
    fn into_iter(self) -> Self::IntoIter {
        match self {
            Polynomial::Zero() => NonzeroCoeffsIter::Zero(),
            Polynomial::Constant(cc) => NonzeroCoeffsIter::Constant(Some(&cc.0)),
            Polynomial::Dense(dc) => NonzeroCoeffsIter::Dense(dc.0.iter().enumerate()),
            Polynomial::Sparse(sc) => NonzeroCoeffsIter::Sparse(sc.0.iter()),
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

    /// Return <code>impl Iterator<Item=C></code>.
    fn coeffs(self) -> Self::IntoCoeffsIter {
        match self {
            Polynomial::Zero() => IntoCoeffsIter::Zero(),
            Polynomial::Constant(cc) => IntoCoeffsIter::Constant(Some(cc.0)),
            Polynomial::Dense(dc) => IntoCoeffsIter::Dense(dc.0.into_iter()),
            Polynomial::Sparse(sc) => {
                let mut map_iter = sc.0.into_iter();
                IntoCoeffsIter::Sparse { index: 0, current: map_iter.next(), map_iter }
            },
        }
    }
}

impl<'a, C> CoeffsIterator<C> for &'a Polynomial<C> where C: Semiring {
    
    type Coeff = Option<&'a C>;
    type IntoCoeffsIter = CoeffsIter<'a, C>;

    /// Return <code>impl Iterator<Item=Option<Option<&C>>></code>.
    fn coeffs(self) -> Self::IntoCoeffsIter {
        match self {
            Polynomial::Zero() => CoeffsIter::Zero(),
            Polynomial::Constant(cc) => CoeffsIter::Constant(Some(Some(&cc.0))),
            Polynomial::Dense(dc) => CoeffsIter::Dense(dc.0.iter()),
            Polynomial::Sparse(sc) => {
                let mut map_iter = sc.0.iter();
                CoeffsIter::Sparse{ index: 0, current: map_iter.next(), map_iter }
            },
        }
    }
}

//********** Eq, PartialEq, Zero and One **********
impl<C> PartialEq for Polynomial<C> where C: Semiring {

    fn eq(&self, other: &Polynomial<C>) -> bool {

        if self.degree() != other.degree() { return false; }

        fn has_the_same_content<C>(vec: &Vec<C>, map: &BTreeMap<usize, C>) -> bool where C: Zero + PartialEq {
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
        Polynomial::Constant(ConstContent(C::one()))  // raw creation
    }
    
    fn is_one(&self) -> bool {
        match self {
            Polynomial::Constant(cc) => cc.0.is_one(),
            _ => false,
        }
    }
}

impl<C> ConstOne for Polynomial<C> where C: Semiring + Clone + ConstOne {
    const ONE: Self = Polynomial::Constant(ConstContent(C::ONE));
}

//********** Operator Overloads **********
//********** Neg **********/
impl<C> Neg for Polynomial<C> where C: Ring {

    type Output = Polynomial<C>;

    fn neg(self) -> Self::Output {
        self.map(|_, c| -c)
    }
}

impl<'a, C> Neg for &'a Polynomial<C> where C: Ring + Clone {

    type Output = Polynomial<C>;

    fn neg(self) -> Self::Output {
        self.map_ref(|_, c| c.ref_neg())
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
            (lhs @ Polynomial::Dense(_), rhs) => dense::add_vv(lhs, rhs),

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
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::add_vv(lhs, rhs),
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
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense::add_vr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse::add_vr(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = c_lhs.ref_add(&c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense::add_vr(lhs, rhs),

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
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::add_vr(lhs, rhs),
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
            (lhs @ Polynomial::Dense(_), rhs) => dense::add_rv(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::add_rv(lhs, rhs),
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
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense::add_rr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse::add_rr(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense::add_rr(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::add_rr(lhs, rhs),
        }
    }
}

impl<C> RefAdd<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_add(&self, other: Polynomial<C>) -> Polynomial<C> { self + other }
}

impl<'b, C> RefAdd<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_add(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self + other }
}

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
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense::sub_vv(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse::sub_vv(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = (*c_lhs).ref_sub(c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense::sub_vv(lhs, rhs),

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
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::sub_vv(lhs, rhs),
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
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense::sub_vr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse::sub_vr(lhs, rhs),

            // Dense
            (Polynomial::Dense(mut dc), Polynomial::Constant(c_rhs)) => {
                dc.0.get_mut(0).map(|c_lhs| *c_lhs = c_lhs.ref_sub(&c_rhs.0));
                Polynomial::Dense(dc)
            },
            (lhs @ Polynomial::Dense(_), rhs) => dense::sub_vr(lhs, rhs),

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
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::sub_vr(lhs, rhs),
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
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense::sub_rv(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse::sub_rv(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense::sub_rv(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::sub_rv(lhs, rhs),
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
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => dense::sub_rr(lhs, rhs),
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => sparse::sub_rr(lhs, rhs),

            // Dense
            (lhs @ Polynomial::Dense(_), rhs) => dense::sub_rr(lhs, rhs),

            // Sparse
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::sub_rr(lhs, rhs),
        }
    }
}

impl<C> RefSub<Polynomial<C>> for Polynomial<C> where C: Ring + Clone {
    #[inline]
    fn ref_sub(&self, other: Polynomial<C>) -> Polynomial<C> { self + other }
}

impl<'b, C> RefSub<&'b Polynomial<C>> for Polynomial<C> where C: Ring + Clone {
    #[inline]
    fn ref_sub(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self + other }
}


//********** Mul **********/
impl<C> Polynomial<C> where C: Semiring + Clone {

    #[inline]
    fn mul_ref<'b>(&self, other: &'b Polynomial<C>) -> Polynomial<C> {
        match (self, other) {
            (_, Polynomial::Zero()) |
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0.ref_mul(&rhs.0)),
            (Polynomial::Constant(lhs), rhs) => rhs.map_ref(|_, c| lhs.0.ref_mul(c)),
            (lhs, Polynomial::Constant(rhs)) => lhs.map_ref(|_, c| c.ref_mul(&rhs.0)),
            (lhs @ Polynomial::Dense(_), rhs) => dense::mul(lhs, rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::mul(lhs, rhs),
        }
    }
}

impl<C> Mul<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: Self) -> Self::Output {
        (&self).mul_ref(&other)
    }
}

impl<'b, C> Mul<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: &'b Self) -> Self::Output {
        (&self).mul_ref(other)
    }
}

impl<'a, C> Mul<Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;

    fn mul(self, other: Polynomial<C>) -> Self::Output {
        self.mul_ref(&other)
    }
}

impl<'a, 'b, C> Mul<&'b Polynomial<C>> for &'a Polynomial<C> where C: Semiring + Clone {

    type Output = Polynomial<C>;
    
    fn mul(self, other: &'b Polynomial<C>) -> Polynomial<C> {
        self.mul_ref(other)
    }
}

impl<C> RefMul<Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_mul(&self, other: Polynomial<C>) -> Polynomial<C> { self * other }
}

impl<'b, C> RefMul<&'b Polynomial<C>> for Polynomial<C> where C: Semiring + Clone {
    #[inline]
    fn ref_mul(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self * other }
}

//********** Div & Rem **********/
impl<C> Polynomial<C> where C: Field + Clone {

    pub fn euclidean_fn(&self) -> usize { self.degree() }

    #[inline]
    fn div_rem_val<'a>(self, other: &'a Polynomial<C>) -> (Self, Self) {
        match (self, other) {
            (_, Polynomial::Zero()) => panic!("Can't divide by polynomial of zero!"),
            (Polynomial::Zero(), _) => (Polynomial::Zero(), Polynomial::Zero()),
            (lhs, Polynomial::Constant(rhs)) => (lhs.map(|_, c| c / (&rhs.0)), Polynomial::Zero()),
            (lhs @ Polynomial::Constant(_), _) => (Polynomial::Zero(), lhs),
            (lhs @ Polynomial::Dense(_), rhs) => dense::div_rem(lhs.to_vec(), rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::div_rem(lhs.to_map(), rhs),
        }
    }

    #[inline]
    fn div_rem_ref(&self, other: &Self) -> (Self, Self) {
        match (self, other) {
            (_, Polynomial::Zero()) => panic!("Can't divide by polynomial of zero!"),
            (Polynomial::Zero(), _) => (Polynomial::Zero(), Polynomial::Zero()),
            (lhs, Polynomial::Constant(rhs)) => (lhs.map_ref(|_, c| c.ref_div(&rhs.0)), Polynomial::Zero()),
            (lhs @ Polynomial::Constant(_), _) => (Polynomial::Zero(), lhs.clone()),
            (lhs @ Polynomial::Dense(_), rhs) => dense::div_rem(lhs.clone_to_vec(), rhs),
            (lhs @ Polynomial::Sparse(_), rhs) => sparse::div_rem(lhs.clone_to_map(), rhs),
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

    fn div(self, other: &'b Polynomial<C>) -> Self::Output { (&self).div_rem_ref(other).0 }
}

impl<'a, C> Div<Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn div(self, other: Polynomial<C>) -> Self::Output { self.div_rem_ref(&other).0 }
}

impl<'a, 'b, C> Div<&'b Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn div(self, other: &'b Polynomial<C>) -> Self::Output { self.div_rem_ref(other).0 }
}

impl<C> RefDiv<Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_div(&self, other: Polynomial<C>) -> Polynomial<C> { self / other }
}

impl<'b, C> RefDiv<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_div(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self / other }
}


impl<C> Rem for Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;
    
    fn rem(self, other: Self) -> Self::Output { self.div_rem_val(&other).1 }
}

impl<'b, C> Rem<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn rem(self, other: &'b Polynomial<C>) -> Self::Output { (&self).div_rem_ref(other).1 }
}

impl<'a, C> Rem<Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn rem(self, other: Polynomial<C>) -> Self::Output { self.div_rem_ref(&other).1 }
}

impl<'a, 'b, C> Rem<&'b Polynomial<C>> for &'a Polynomial<C> where C: Field + Clone {

    type Output = Polynomial<C>;

    fn rem(self, other: &'b Polynomial<C>) -> Self::Output { self.div_rem_ref(other).1 }
}

impl<C> RefRem<Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_rem(&self, other: Polynomial<C>) -> Polynomial<C> { self % other }
}

impl<'b, C> RefRem<&'b Polynomial<C>> for Polynomial<C> where C: Field + Clone {
    #[inline]
    fn ref_rem(&self, other: &'b Polynomial<C>) -> Polynomial<C> { self % other }
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

// impl<C: Ring + Clone> BitXor<u32> for Polynomial<C> {

//     type Output = Polynomial<C>;

//     fn bitxor(self, power: u32) -> Self::Output {
//         self.pow(power)
//     }
// }

// impl<'a, C: Ring + Clone> BitXor<u32> for &'a Polynomial<C> {

//     type Output = Polynomial<C>;

//     fn bitxor(self, power: u32) -> Self::Output {
//         self.pow(power)
//     }
// }

//   /**
//    * Shift this polynomial along the x-axis by `h`, so that `this(x + h) == this.shift(h).apply(x)`. This is equivalent
//    * to calling `this.compose(Polynomial.x + h)`, but is likely to compute the shifted polynomial much faster.
//    */
//   def shift(h: C)(implicit ring: Ring[C], eq: Eq[C]): Polynomial[C] = {
//     // The trick here came from this answer:
//     //   http://math.stackexchange.com/questions/694565/polynomial-shift
//     // This is a heavily optimized version of the same idea. This is fairly
//     // critical method to be fast, since it is the most expensive part of the
//     // VAS root isolation algorithm.

//     def fromSafeLong(x: SafeLong): C =
//       if (x.isValidInt) {
//         ring.fromInt(x.toInt)
//       } else {
//         val d = ring.fromInt(1 << 30)
//         val mask = (1L << 30) - 1

//         @tailrec def loop(k: C, y: SafeLong, acc: C): C =
//           if (y.isValidInt) {
//             k * ring.fromInt(y.toInt) + acc
//           } else {
//             val z = y >> 30
//             val r = ring.fromInt((y & mask).toInt)
//             loop(d * k, z, k * r + acc)
//           }

//         loop(ring.one, x, ring.zero)
//       }

//     // The basic idea here is that instead of working with all the derivatives
//     // of the whole polynomial, we can just break the polynomial up and work
//     // with the derivatives of the individual terms. This let's us save a whole
//     // bunch of allocations in a clean way.
//     val coeffs: Array[C] = this.coeffsArray.clone()
//     this.foreachNonzero { (deg, c) =>
//       var i = 1 // Leading factor in factorial in denominator of Taylor series.
//       var d = deg - 1 // The degree of the current derivative of this term.
//       var m = SafeLong(1L) // The multiplier of our derivative
//       var k = c // The current delta (to some power) of the Taylor series.
//       while (d >= 0) {
//         // Note that we do division, but only on SafeLongs. This is not just
//         // for performance, but also required for us to only ask for a Ring,
//         // rather than a EuclideanRing. We always know that m * (d + 1) is
//         // divisible by i, so this is exact.
//         m = (m * (d + 1)) / i
//         k *= h
//         coeffs(d) = coeffs(d) + fromSafeLong(m) * k
//         d -= 1
//         i += 1
//       }
//     }
//     Polynomial.dense(coeffs)
//   }

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

//   /**
//    * Removes all zero roots from this polynomial.
//    */
//   def removeZeroRoots(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] = {
//     val Term(_, k) = minTerm
//     mapTerms { case Term(c, n) => Term(c, n - k) }
//   }

//   def mapTerms[D: Semiring: Eq: ClassTag](f: Term[C] => Term[D]): Polynomial[D] =
//     Polynomial(termsIterator.map(f))

//   /**
//    * This will flip/mirror the polynomial about the y-axis. It is equivalent to `poly.compose(-Polynomial.x)`, but will
//    * likely be faster to calculate.
//    */
//   def flip(implicit ring: Rng[C], eq: Eq[C]): Polynomial[C] =
//     mapTerms { case term @ Term(coeff, exp) =>
//       if (exp % 2 == 0) term
//       else Term(-coeff, exp)
//     }

//   /**
//    * Returns the reciprocal of this polynomial. Essentially, if this polynomial is `p` with degree `n`, then returns a
//    * polynomial `q(x) = x^n*p(1/x)`.
//    *
//    * @see
//    *   http://en.wikipedia.org/wiki/Reciprocal_polynomial
//    */
//   def reciprocal(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] = {
//     val d = degree
//     mapTerms { case term @ Term(coeff, exp) =>
//       Term(coeff, d - exp)
//     }
//   }

//   // VectorSpace ops.

//   def *:(k: C)(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C]
//   def :*(k: C)(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] = k *: lhs
//   def :/(k: C)(implicit field: Field[C], eq: Eq[C]): Polynomial[C] = this :* k.reciprocal

//   override def hashCode: Int = {
//     val it = lhs.termsIterator
//     @tailrec def loop(n: Int): Int =
//       if (it.hasNext) {
//         val term = it.next()
//         loop(n ^ (0xfeed1257 * term.exp ^ term.coeff.##))
//       } else n
//     loop(0)
//   }






impl<C> Semigroup for Polynomial<C> where C: Semiring + Clone {}
impl<C> Monoid for Polynomial<C> where C: Semiring + Clone {}

impl<C> AdditiveSemigroup for Polynomial<C> where C: Semiring + Clone {}

impl<C> AdditiveMonoid for Polynomial<C> where C: Semiring + Clone {}
impl<C> AdditiveGroup for Polynomial<C> where C: Ring + Clone {
    fn ref_neg(&self) -> Self { -self }
}

impl<C> Semiring for Polynomial<C> where C: Ring + Clone {}

impl<C> Ring for Polynomial<C> where C: Ring + Clone{}

impl<C> EuclideanRing for Polynomial<C> where C: Field + Clone {

    #[inline]
    fn div_rem(self, other: Self) -> (Self, Self) { 
        (&self).div_rem_ref(&other)
    }

    #[inline]
    fn div_rem_ref(&self, other: &Self) -> (Self, Self) {
        self.div_rem_euclid(other)
    }
} 