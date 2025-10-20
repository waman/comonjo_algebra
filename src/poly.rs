pub(crate) mod add;
pub(crate) mod sub;
pub(crate) mod mul;
pub(crate) mod div_rem;
pub mod term;

use std::{collections::{btree_map, BTreeMap, HashMap}, fmt::{Debug, Display}, iter::Enumerate, ops::*, slice::Iter, vec::{self, IntoIter}};
use num::{pow::Pow, traits::{ConstOne, ConstZero, Euclid}, One, Zero};

use once_cell::sync::Lazy;
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
    pub fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::constant(f(0, cc.0)),
            Polynomial::Dense(dc) => dc.map_nonzero(f),
            Polynomial::Sparse(sc) => sc.map_nonzero(f),
        }
    }
    
    /// Note that <code>f</code> is applied only to nonzero coefficients
    pub fn ref_map_nonzero<D, F>(&self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, &C) -> D {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => Polynomial::constant(f(0, &cc.0)),
            Polynomial::Dense(dc) => dc.ref_map_nonzero(f),
            Polynomial::Sparse(sc) => sc.ref_map_nonzero(f),
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
    /// <code>x</code>
    pub fn x() -> Polynomial<C> { sparse![(1, C::one())] }

    /// <code>x²</code>
    pub fn x2() -> Polynomial<C> { sparse![(2, C::one())] }

    /// <code>x³</code>
    pub fn x3() -> Polynomial<C> { sparse![(3, C::one())] }

    /// <code>x⁴</code>
    pub fn x4() -> Polynomial<C> { sparse![(4, C::one())] }

    /// <code>x⁵</code>
    pub fn x5() -> Polynomial<C> { sparse![(5, C::one())] }

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

#[derive(Clone)]
pub struct DenseContent<C>(pub(crate) Vec<C>);

impl<C> DenseContent<C> where C: Semiring {

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(n)
    }

    pub fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let v: Vec<D> = self.0.into_iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::dense_from_vec(v)
    }

    pub fn ref_map_nonzero<D, F>(&self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, &C) -> D {
        let v: Vec<D> = self.0.iter().enumerate().map(|(i, c)|
            if c.is_zero() { D::zero() } else { f(i, c) }
        ).collect();
        Polynomial::dense_from_vec(v)
    }

    pub(crate) fn to_map(mut self) -> BTreeMap<usize, C> {
        let mut map = BTreeMap::new();
        while let Some(c) = self.0.pop() {
            if !c.is_zero() {
                map.insert(self.0.len(), c);
            }
        }
        debug_assert!(self.0.is_empty());
        map
    }

    pub fn reductum(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn monic(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn derivative(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn integral(&self) -> Polynomial<C> {
        todo!()
    }

    // pub fn add_dense(mut self, mut other: DenseContent<C>) -> Polynomial<C> {
    //     // if self.degree() >= other.degree() {
    //     //     let mut iteref_rhs = other.0.iter().enumerate();
    //     //     for (i, c_lhs) in self.0.iteref_mut().enumerate() {

    //     //     }
    //     //     let mut lhs_iter = lhs.coeffs();
    //     //     let mut rhs_iter = rhs.coeffs();
        
    //     //     for _ in 0..=d_min {
    //     //         match (lhs_iter.next(), rhs_iter.next()) {
    //     //             (Some(x), Some(y)) => v.push(x + y),
    //     //             _ => panic!(),
    //     //         }
    //     //     }
        
    //     //     let rest_iter = if lhs_is_higher { lhs_iter } else { rhs_iter };
    //     //     v.extend(rest_iter);
        
    //     //     Polynomial::dense_from_vec(v)

    //     // }
    //     todo!()
    // }

    // pub fn add_val(mut self, mut other: Polynomial<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn add_ref<'b>(mut self, mut other: &'b Polynomial<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn ref_add_dense(&self, mut other: DenseContent<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn ref_add_sparse(&self, mut other: Polynomial<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn ref_add_ref<'b>(&self, mut other: &'b Polynomial<C>) -> Polynomial<C> {
    //     todo!()
    // }


    // pub fn sub_dense(mut self, mut other: DenseContent<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn sub_sparse(mut self, mut other: SparseContent<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn sub_ref<'b>(mut self, mut other: &'b Polynomial<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn ref_sub_dense(&self, mut other: DenseContent<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn ref_sub_sparse(&self, mut other: Polynomial<C>) -> Polynomial<C> {
    //     todo!()
    // }

    // pub fn ref_sub_ref<'b>(&self, mut other: &'b Polynomial<C>) -> Polynomial<C> {
    //     todo!()
    // }
    

// def *(rhs: Polynomial[C])(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] = {
//   if (rhs.isZero) return rhs
//   if (lhs.isZero) return lhs
//   val lcs = lhs.coeffsArray
//   val rcs = rhs.coeffsArray
//   val cs = new Array[C](lcs.length + rcs.length - 1)
//   cfor(0)(_ < cs.length, _ + 1) { i => cs(i) = ring.zero }
//   cfor(0)(_ < lcs.length, _ + 1) { i =>
//     val c = lcs(i)
//     var k = i
//     cfor(0)(_ < rcs.length, _ + 1) { j =>
//       cs(k) += c * rcs(j)
//       k += 1
//     }
//   }
//   Polynomial.dense(cs)
// }

// def reductum(implicit e: Eq[C], ring: Semiring[C], ct: ClassTag[C]): Polynomial[C] = {
//   var i = coeffs.length - 2
//   while (i >= 0 && coeffs(i) === ring.zero) i -= 1
//   if (i < 0) {
//     new PolyDense(new Array[C](0))
//   } else {
//     val arr = new Array[C](i + 1)
//     System.arraycopy(coeffs, 0, arr, 0, i + 1)
//     new PolyDense(arr)
//   }
// }

// def apply(x: C)(implicit ring: Semiring[C]): C = {
//   if (isZero) return ring.zero

//   var even = coeffs.length - 1
//   var odd = coeffs.length - 2
//   if ((even & 1) == 1) { even = odd; odd = coeffs.length - 1 }

//   var c0 = coeffs(even)
//   val x2 = x.pow(2)
//   cfor(even - 2)(_ >= 0, _ - 2) { i =>
//     c0 = coeffs(i) + c0 * x2
//   }

//   if (odd >= 1) {
//     var c1 = coeffs(odd)
//     cfor(odd - 2)(_ >= 1, _ - 2) { i =>
//       c1 = coeffs(i) + c1 * x2
//     }
//     c0 + c1 * x
//   } else {
//     c0
//   }
// }

// def derivative(implicit ring: Ring[C], eq: Eq[C]): Polynomial[C] = {
//   if (isZero) return this
//   val cs = new Array[C](degree)
//   var j = coeffs.length - 1
//   cfor(cs.length - 1)(_ >= 0, _ - 1) { i =>
//     cs(i) = ring.fromInt(j) * coeffs(j)
//     j -= 1
//   }
//   Polynomial.dense(cs)
// }

// def integral(implicit field: Field[C], eq: Eq[C]): Polynomial[C] = {
//   val cs = new Array[C](coeffs.length + 1)
//   cs(0) = field.zero
//   cfor(0)(_ < coeffs.length, _ + 1) { i => cs(i + 1) = coeffs(i) / field.fromInt(i + 1) }
//   Polynomial.dense(cs)
// }

// def *:(k: C)(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] =
//   if (k === ring.zero) {
//     Polynomial.dense(new Array[C](0))
//   } else {
//     val cs = new Array[C](coeffs.length)
//     cfor(0)(_ < cs.length, _ + 1) { i =>
//       cs(i) = k * coeffs(i)
//     }
//     Polynomial.dense(cs)
//   }
// }
}

#[derive(Clone)]
pub struct SparseContent<C>(pub(crate) BTreeMap<usize, C>);

impl<C> SparseContent<C> where C: Semiring {

    pub fn degree(&self) -> usize {
        *self.0.last_key_value().unwrap().0
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(&n)
    }

    pub fn map_nonzero<D, F>(self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, C) -> D {
        let m: BTreeMap<usize, D> = self.0.into_iter().map(|(i, c)| (i, f(i, c))).collect();
        Polynomial::sparse_from_map(m)
    }

    pub fn ref_map_nonzero<D, F>(&self, f: F) -> Polynomial<D> where D: Semiring, F: Fn(usize, &C) -> D {
        let m: BTreeMap<usize, D> = self.0.iter().map(|(i, c)| (*i, f(*i, c))).collect();
        Polynomial::sparse_from_map(m)
    }

    pub(crate) fn to_vec(mut self) -> Vec<C> {
        let n = self.degree() + 1;
        let mut v: Vec<C> = Vec::with_capacity(n);
        for i in 0..n {
            match self.0.remove(&i) {
                Some(e) => v.push(e),
                None => v.push(C::zero()),
            } 
        }
        debug_assert!(self.0.is_empty());
        v
    }

    pub fn reductum(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn monic(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn derivative(&self) -> Polynomial<C> {
        todo!()
    }

    pub fn integral(&self) -> Polynomial<C> {
        todo!()
    }
//    def reductum(implicit e: Eq[C], ring: Semiring[C], ct: ClassTag[C]): Polynomial[C] = {
//      var i = coeff.length - 2
//      while (i >= 0 && coeff(i) === ring.zero) i -= 1
//      if (i < 0) {
//        new PolySparse(new Array[Int](0), new Array[C](0))
//      } else {
//        val len = i + 1
//        val es = new Array[Int](len)
//        val cs = new Array[C](len)
//        System.arraycopy(coeff, 0, cs, 0, len)
//        System.arraycopy(exp, 0, es, 0, len)
//        new PolySparse(es, cs)
//      }
//    }
 
//    final private def expBits(x: C)(implicit ring: Semiring[C]): Array[C] = {
//      val bits = new Array[C](math.max(2, 32 - numberOfLeadingZeros(degree)))
//      bits(0) = x
//      // we use pow(2) here for the benefit of Interval[_], where
//      // x.pow(2) has better error bounds than than (x * x).
//      if (bits.length > 1) bits(1) = x.pow(2)
//      cfor(2)(_ < bits.length, _ + 1) { i =>
//        val prev = bits(i - 1)
//        bits(i) = prev * prev
//      }
//      bits
//    }
 
//    @tailrec
//    final private def fastExp(bits: Array[C], e: Int, i: Int, acc: C)(implicit ring: Semiring[C]): C = {
//      if (e == 0) acc
//      else {
//        val lb = numberOfTrailingZeros(e) + 1
//        val j = i + lb
//        fastExp(bits, e >>> lb, j, acc * bits(j - 1))
//      }
//    }
 
//    final private def fastExp(bits: Array[C], e: Int)(implicit ring: Semiring[C]): C = {
//      val lb = numberOfTrailingZeros(e) + 1
//      fastExp(bits, e >>> lb, lb, bits(lb - 1))
//    }
 
//    def apply(x: C)(implicit ring: Semiring[C]): C = if (isZero) {
//      ring.zero
//    } else if (exp.length == 1) {
//      if (exp(0) != 0) coeff(0) * (x.pow(exp(0))) else coeff(0)
//    } else {
//      // TODO: Rewrite this to be more like PolyDense.
//      val bits = expBits(x)
//      val e0 = exp(0)
//      val c0 = coeff(0)
//      var sum = if (e0 == 0) c0 else c0 * fastExp(bits, e0)
//      cfor(1)(_ < exp.length, _ + 1) { i =>
//        sum += coeff(i) * fastExp(bits, exp(i))
//      }
//      sum
//    }
 
//    def derivative(implicit ring: Ring[C], eq: Eq[C]): Polynomial[C] =
//      if (exp.length == 0) this
//      else {
//        val i0 = if (exp(0) == 0) 1 else 0
//        val es = new Array[Int](exp.length - i0)
//        val cs = new Array[C](es.length)
 
//        @tailrec
//        def loop(i: Int, j: Int): Unit = if (j < es.length) {
//          val e = exp(i)
//          es(j) = e - 1
//          cs(j) = e * coeff(i)
//          loop(i + 1, j + 1)
//        }
 
//        loop(i0, 0)
//        PolySparse.safe(es, cs)
//      }
 
//    def integral(implicit field: Field[C], eq: Eq[C]): Polynomial[C] = {
//      val es = new Array[Int](exp.length)
//      val cs = new Array[C](es.length)
 
//      cfor(0)(_ < es.length, _ + 1) { i =>
//        val e = exp(i) + 1
//        es(i) = e
//        cs(i) = coeff(i) / field.fromInt(e)
//      }
 
//      PolySparse.safe(es, cs)
//    }
 
//    def unary_-(implicit ring: Rng[C]): Polynomial[C] = {
//      val cs = new Array[C](coeff.length)
//      cfor(0)(_ < cs.length, _ + 1) { i => cs(i) = -coeff(i) }
//      new PolySparse(exp, cs)
//    }
 
//    def *(rhs0: Polynomial[C])(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] = {
//      val rhs: PolySparse[C] = PolySparse(rhs0)
//      PolySparse.multiplySparse(lhs, rhs)
//    }
 
//    def *:(k: C)(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C] = {
//      if (k === ring.zero) {
//        PolySparse.zero[C]
//      } else {
//        val cs = new Array[C](coeff.length)
//        cfor(0)(_ < cs.length, _ + 1) { i =>
//          cs(i) = k * coeff(i)
//        }
//        new PolySparse(exp, cs)
//      }
//    }
//  }
 
//  object PolySparse {
//    final private[math] def dense2sparse[@sp(Double) C: Semiring: Eq: ClassTag](poly: PolyDense[C]): PolySparse[C] = {
//      val cs = poly.coeffs
//      val es = new Array[Int](cs.length)
//      cfor(0)(_ < es.length, _ + 1) { i => es(i) = i }
//      PolySparse.safe(es, cs)
//    }
 
//    final private[math] def safe[@sp(Double) C: Semiring: Eq: ClassTag](exp: Array[Int],
//                                                                        coeff: Array[C]
//    ): PolySparse[C] = {
//      var len = 0
//      cfor(0)(_ < coeff.length, _ + 1) { i =>
//        if (coeff(i) =!= Semiring[C].zero)
//          len += 1
//      }
 
//      if (len == coeff.length) {
//        new PolySparse(exp, coeff)
//      } else {
//        val es = new Array[Int](len)
//        val cs = new Array[C](len)
//        @tailrec def loop(i: Int, j: Int): PolySparse[C] =
//          if (i < coeff.length) {
//            val c = coeff(i)
//            if (c =!= Semiring[C].zero) {
//              es(j) = exp(i)
//              cs(j) = c
//              loop(i + 1, j + 1)
//            } else {
//              loop(i + 1, j)
//            }
//          } else new PolySparse(es, cs)
//        loop(0, 0)
//      }
//    }
 
//    final def apply[@sp(Double) C: Semiring: Eq: ClassTag](data: IterableOnce[Term[C]]): PolySparse[C] = {
//      import spire.scalacompat.arrayBuilderMake
 
//      var expBldr = arrayBuilderMake[Int]
//      var coeffBldr = arrayBuilderMake[C]
//      val zero = Semiring[C].zero
//      var inReverseOrder = true
//      var inOrder = true
//      var lastDeg = -1
 
//      data.iterator.foreach { case Term(c, i) =>
//        if (c =!= zero) {
//          expBldr += i
//          coeffBldr += c
//          inOrder &&= (lastDeg < i)
//          inReverseOrder &&= (lastDeg < 0 || lastDeg > i)
//          lastDeg = i
//        }
//      }
 
//      val exp = expBldr.result()
//      val coeff = coeffBldr.result()
//      if (inOrder) {
//        PolySparse(exp, coeff)
//      } else if (inReverseOrder) {
//        reverse(exp); reverse(coeff)
//        PolySparse(exp, coeff)
//      } else {
//        val indices = Array.range(0, exp.length)
//        indices.qsortBy(exp(_))
//        expBldr = arrayBuilderMake[Int]
//        coeffBldr = arrayBuilderMake[C]
//        var i = 1
//        var j = indices(0)
//        var e = exp(j)
//        var c = coeff(j)
//        while (i < indices.length) {
//          val j0 = indices(i)
//          val e0 = exp(j0)
//          val c0 = coeff(j0)
//          if (e != e0) {
//            if (!c.isZero) {
//              expBldr += e
//              coeffBldr += c
//            }
//            c = c0
//          } else {
//            c += c0
//          }
//          e = e0
//          j = j0
//          i += 1
//        }
//        if (!c.isZero) {
//          expBldr += e
//          coeffBldr += c
//        }
//        val poly = PolySparse(expBldr.result(), coeffBldr.result())
//        poly
//      }
//    }
 
//    private def reverse[@sp(Double) A](arr: Array[A]): Unit = {
//      var i = 0
//      var j = arr.length - 1
//      while (i < j) {
//        val tmp = arr(i)
//        arr(i) = arr(j)
//        arr(j) = tmp
//        i += 1
//        j -= 1
//      }
//    }
 
//    final def apply[@sp(Double) C: Semiring: Eq: ClassTag](data: Map[Int, C]): PolySparse[C] = {
//      val data0 = data.toArray
//      data0.qsortBy(_._1)
//      val es = new Array[Int](data0.length)
//      val cs = new Array[C](data0.length)
//      cfor(0)(_ < data0.length, _ + 1) { i =>
//        val (e, c) = data0(i)
//        es(i) = e
//        cs(i) = c
//      }
//      safe(es, cs)
//    }
 
//    final def apply[@sp(Double) C: Semiring: Eq: ClassTag](poly: Polynomial[C]): PolySparse[C] = {
//      poly match {
//        case (poly: PolySparse[_]) =>
//          poly
 
//        case (_: PolyDense[_]) =>
//          dense2sparse(poly.asInstanceOf[PolyDense[C]]) // Yay...
 
//        case _ =>
//          var len = 0
//          poly.foreachNonzero { (_, _) => len += 1 }
//          val es = new Array[Int](len)
//          val cs = new Array[C](len)
//          var i = 0
//          poly.foreachNonzero { (e, c) =>
//            es(i) = e
//            cs(i) = c
//            i += 1
//          }
//          PolySparse.safe(es, cs)
//      }
//    }
 
//    final private def multiplyTerm[@sp(Double) C: Semiring: Eq: ClassTag](poly: PolySparse[C],
//                                                                          c: C,
//                                                                          e: Int
//    ): PolySparse[C] = {
//      val exp = poly.exp
//      val coeff = poly.coeff
//      val cs = new Array[C](coeff.length)
//      val es = new Array[Int](exp.length)
//      cfor(0)(_ < coeff.length, _ + 1) { i =>
//        cs(i) = c * coeff(i)
//        es(i) = exp(i) + e
//      }
//      new PolySparse(es, cs)
//    }
 
//    final private def multiplySparse[@sp(Double) C: Semiring: Eq: ClassTag](lhs: PolySparse[C],
//                                                                            rhs: PolySparse[C]
//    ): PolySparse[C] = {
//      val lexp = lhs.exp
//      val lcoeff = lhs.coeff
//      var sum = new PolySparse(new Array[Int](0), new Array[C](0))
//      cfor(0)(_ < lexp.length, _ + 1) { i =>
//        sum = addSparse(sum, multiplyTerm(rhs, lcoeff(i), lexp(i)))
//      }
//      sum
//    }
 
//    final private def countSumTerms[@sp(Double) C](lhs: PolySparse[C],
//                                                   rhs: PolySparse[C],
//                                                   lOffset: Int = 0,
//                                                   rOffset: Int = 0
//    ): Int = {
//      val PolySparse(lexp, lcoeff) = lhs
//      val PolySparse(rexp, rcoeff) = rhs
 
//      @tailrec
//      def loop(i: Int, j: Int, count: Int): Int =
//        if (i < lexp.length && j < rexp.length) {
//          val cmp = lexp(i) + lOffset - rexp(j) - rOffset
//          if (cmp == 0) loop(i + 1, j + 1, count + 1)
//          else if (cmp < 0) loop(i + 1, j, count + 1)
//          else loop(i, j + 1, count + 1)
//        } else {
//          count + (lexp.length - i) + (rexp.length - j)
//        }
 
//      loop(0, 0, 0)
//    }
}

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

    fn scale_by_left(self, k: &C) -> Polynomial<C> {
        if k.is_zero() { return Polynomial::Zero() }
        if k.is_one() { return self }
        self.map_nonzero(|_, c| k.ref_mul(c))
    }

    fn ref_scale_by_left(&self, k: &C) -> Polynomial<C> {
        if k.is_zero() { return Polynomial::Zero() }
        if k.is_one() { return self.clone() }
        self.ref_map_nonzero(|_, c| k.ref_mul(c))
    }

    /// Compose this polynomial with another.
    pub fn compose<'a>(&self, y: &'a Polynomial<C>) -> Polynomial<C> {
        match (self, y) {
            (Polynomial::Zero(), _) => Polynomial::Zero(),
            (lhs @ Polynomial::Constant(_), _) => lhs.clone(),
            (lhs, Polynomial::Zero()) => match lhs.nth(0) {
                Some(c) => Polynomial::constant(c.clone()),
                _ => Polynomial::Zero(),
            }
            (lhs, rhs) => {
                let mut acc = Polynomial::Zero();

                for (i, c) in lhs.nonzero_coeffs() {
                    let z = rhs.pow(i as u32).scale_by_left(c);
                    acc = acc + z;
                }

                acc
            }
        }
    }
    

    
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

fn term_to_string<C>(exp: usize, coeff: &C) -> String where C: Semiring + Display {

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

        } else {
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
}fn const_iter_size_hint<'a, C>(c: &'a Option<C>) -> (usize, Option<usize>) {
    if c.is_some() { (1, Some(1)) } else { (0, Some(0)) }
}

pub enum IntoNonzeroCoeffsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(Enumerate<IntoIter<C>>),
    Sparse(btree_map::IntoIter<usize, C>),
}

impl<C> Iterator for IntoNonzeroCoeffsIter<C> where C: Semiring {

    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoNonzeroCoeffsIter::Zero() => None,
            IntoNonzeroCoeffsIter::Constant(value) => 
                match value.take() {
                    Some(c) => Some((0, c)),
                    None => None,
                },
            IntoNonzeroCoeffsIter::Dense(vec_iter) => 
                loop {
                    match vec_iter.next() {
                        Some((i, c)) =>  if !c.is_zero() { return Some((i, c)); },
                        None => return None,
                    }
                },
            IntoNonzeroCoeffsIter::Sparse(map_iter) => map_iter.next(),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoNonzeroCoeffsIter::Zero() => (0, Some(0)),
            IntoNonzeroCoeffsIter::Constant(c) => const_iter_size_hint(c),
            IntoNonzeroCoeffsIter::Dense(ite) => ite.size_hint(),
            IntoNonzeroCoeffsIter::Sparse(ite) => ite.size_hint(),
        }
    }
}

pub enum NonzeroCoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<&'a C>),
    Dense(Enumerate<Iter<'a, C>>),
    Sparse(btree_map::Iter<'a, usize, C>),
}

impl<'a, C> Iterator for NonzeroCoeffsIter<'a, C> where C: Semiring {

    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            NonzeroCoeffsIter::Zero() => None,
            NonzeroCoeffsIter::Constant(value) => 
                match value.take() {
                    Some(c) => Some((0, c)),
                    None => None,
                },
            NonzeroCoeffsIter::Dense(vec_iter) => 
                loop {
                    match vec_iter.next() {
                        Some((i, c)) =>  if !c.is_zero() { return Some((i, c)); },
                        None => return None,
                    }
                },
            NonzeroCoeffsIter::Sparse(map_iter) => map_iter.next().map(|(i, c)| (*i, c)),
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            NonzeroCoeffsIter::Zero() => (0, Some(0)),
            NonzeroCoeffsIter::Constant(c) => const_iter_size_hint(c),
            NonzeroCoeffsIter::Dense(ite) => ite.size_hint(),
            NonzeroCoeffsIter::Sparse(ite) => ite.size_hint(),
        }
    }
}

pub enum IntoCoeffsIter<C> where C: Semiring {
    Zero(),
    Constant(Option<C>),
    Dense(vec::IntoIter<C>),
    Sparse{
        index: usize,
        current: Option<(usize, C)>,
        map_iter: btree_map::IntoIter<usize, C>,
    },
}

impl<C> Iterator for IntoCoeffsIter<C> where C: Semiring {

    type Item = C;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoCoeffsIter::Zero() => None,
            IntoCoeffsIter::Constant(value) => value.take(),
            IntoCoeffsIter::Dense(vec_iter) =>  vec_iter.next(),
            IntoCoeffsIter::Sparse{ index, current, map_iter, .. } => {
                let result = match current {
                    Some(c) => 
                        if c.0 == *index {
                            let r = match map_iter.next() {
                                Some(nxt) => current.replace(nxt),
                                None => current.take(),
                            };
                            Some(r.unwrap().1)
                        } else { 
                            Some(C::zero())
                        }
                    _ => None,
                };
        
                *index += 1;
                return result;
            },
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            IntoCoeffsIter::Zero() => (0, Some(0)),
            IntoCoeffsIter::Constant(c) => const_iter_size_hint(c),
            IntoCoeffsIter::Dense(ite) => ite.size_hint(),
            IntoCoeffsIter::Sparse { map_iter, .. } => map_iter.size_hint(),
        }
    }
}

/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some(Some(0))</code> if the coefficient is zero.
pub enum CoeffsIter<'a, C> where C: Semiring {
    Zero(),
    Constant(Option<Option<&'a C>>),
    Dense(Iter<'a, C>),
    Sparse{
        index: usize,
        current: Option<(&'a usize, &'a C)>,
        map_iter: btree_map::Iter<'a, usize, C>,
    },
}

impl<'a, C> Iterator for CoeffsIter<'a, C> where C: Semiring {

    type Item = Option<&'a C>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            CoeffsIter::Zero() => None,
            CoeffsIter::Constant(value) => value.take(),
            CoeffsIter::Dense(vec_iter) => 
                match vec_iter.next() {
                    Some(c) => Some(Some(c)),
                    None => None,
                } ,
            CoeffsIter::Sparse{ index, current, map_iter, .. } => 
                if let Some(c) = current {
                    let result = if c.0 == index { 
                        let r = Some(c.1);
                        *current = map_iter.next();
                        r
                    } else {
                        None
                    };
                    *index += 1;
                    Some(result)
                } else {
                    return None;
                },
        }
    }
    
    fn size_hint(&self) -> (usize, Option<usize>) { 
        match self {
            CoeffsIter::Zero() => (0, Some(0)),
            CoeffsIter::Constant(_) => (1, Some(1)),
            CoeffsIter::Dense(ite) => ite.size_hint(),
            CoeffsIter::Sparse { map_iter, .. } => map_iter.size_hint(),
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
        self.ref_map_nonzero(|_, c| c.ref_neg())
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

// operator priority is problem...
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



//********** Utility Functions *********/
/// Return (max, min, left_is_higher)
fn max_min(x: usize, y: usize) -> (usize, usize, bool) {
    if x >= y { (x, y, true) } else { (y, x, false) }
}


//********** Implementation of Algebra *********/
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
        self.div_rem_euclid(&other)
    }

    #[inline]
    fn div_rem_ref(&self, other: &Self) -> (Self, Self) {
        self.div_rem_euclid(other)
    }
} 