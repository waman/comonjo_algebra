use std::{collections::BTreeMap, fmt::{Debug, Display}, ops::{Add, Mul}};

use num::{traits::{ConstOne, ConstZero}, Num, One, Zero};

use crate::poly::{dense_content::{DenseContent, DenseIter}, spears_content::{SpearsContent, SpearsIter}, term::Term};

/// Refer to spire's <a href="https://github.com/typelevel/spire/blob/main/core/src/main/scala/spire/math/Polynomial.scala">Polynomial</a>

pub enum Polynomial<C: Num> {

    Zero(),

    /// Use <code>Polynomial::constant()</code> function to instantiate.
    Constant(ConstContent<C>),

    /// Use <code>dense!()</code> macro to instantiate.
    Dense(DenseContent<C>),
    
    /// Use <code>spears!()</code> macro to instantiate.
    Spares(SpearsContent<C>)
}

impl<C: Num> Polynomial<C> {

    pub fn constant<D: Num>(c: D) -> Polynomial<D> {
        if c.is_zero() {
            Polynomial::Zero()
        } else {
            Polynomial::Constant(ConstContent(c))
        }
    }

    pub fn dense_from_vec<D: Num>(mut coeffs: Vec<D>) -> Polynomial<D>  {

        while !coeffs.is_empty() && coeffs.last().unwrap().is_zero() {
            coeffs.pop();
        }

        match coeffs.len() {
            0 => Polynomial::<D>::constant(D::zero()),
            1 => Polynomial::<D>::constant(coeffs.pop().unwrap()),
            _ => {
                coeffs.shrink_to_fit();
                Polynomial::Dense(DenseContent(coeffs))
            }
        }
    }

    pub fn spears_from_map<D: Num>(mut coeffs: BTreeMap<usize, D>) -> Polynomial<D> {

        coeffs.retain(|_, v| !v.is_zero());

        match coeffs.len() {
            0 => Polynomial::<D>::constant(D::zero()),
            1 => {
                let mut key = 0_usize;
                for k in coeffs.keys() { key = *k; }
                let value = coeffs.remove(&key).unwrap();
                Polynomial::<D>::constant(value)
            },
            _ => Polynomial::Spares(SpearsContent(coeffs))
        }
    }

    pub fn parse<D: Num>(s: &str) -> Polynomial<D> {
    
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

    pub fn x<D: Num>() -> Polynomial<D> {
        Polynomial::<D>::spears_from_map(BTreeMap::from([(1, D::one())]))
    }

    pub fn two_x<D: Num>() -> Polynomial<D> {
        Polynomial::<D>::spears_from_map(BTreeMap::from([(1, D::one() + D::one())]))
    }

    /// Create a linear polynomial <i>ax</i>.
    pub fn linear_monomial<D: Num>(a: D) -> Polynomial<D> {
        Polynomial::<D>::spears_from_map(BTreeMap::from([(1, a)]))
    }

    /// Create a linear polynomial <i>ax + b</i>
    pub fn linear<D: Num>(c1: D, c0: D) -> Polynomial<D> {
        Polynomial::<D>::spears_from_map(BTreeMap::from([(1, c1), (0, c0)]))
    }

    /// Create a quadratic polynomial <i>ax<sup>2</sup><i>
    pub fn quadratic_monomial<D: Num>(a: D) -> Polynomial<D> {
        Polynomial::<D>::spears_from_map(BTreeMap::from([(2, a)]))
    }

    /// Create a quadratic polynomial <i>ax<sup>2</sup> + bx + c<i>
    pub fn quadratic<D: Num>(a: D, b: D, c: D) -> Polynomial<D> {
        Polynomial::<D>::spears_from_map(BTreeMap::from([(2, a), (1, b), (0, c)]))
    }

    /// Create a cubic polynomial <i>ax<sup>3</sup> + bx<sup>2</sup> + cx + d<i>
    pub fn cubic_monomial<D: Num>(a: D) -> Polynomial<D> {
        Polynomial::<D>::spears_from_map(BTreeMap::from([(3, a)]))
    }

    /// Create a cubic polynomial <i>ax<sup>3</sup> + bx<sup>2</sup> + cx + d<i>
    pub fn cubic<D: Num>(a: D, b: D, c: D, d: D) -> Polynomial<D> {
        Polynomial::<D>::spears_from_map(BTreeMap::from([(3, a), (2, b), (1, c), (0, d)]))
    }

    //********** METHODS ******
    pub fn is_constant(&self) -> bool {
        match self {
            Polynomial::Zero() => true,
            Polynomial::Constant(_) => true,
            Polynomial::Dense(_) => false,
            Polynomial::Spares(_) => false,
        }
    }

    pub fn degree(&self) -> usize {
        match self {
            Polynomial::Zero() => 0,
            Polynomial::Constant(_) => 0,
            Polynomial::Dense(dc) => dc.degree(),
            Polynomial::Spares(sc) => sc.degree(),
        }
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        match self {
            Polynomial::Zero() => None,
            Polynomial::Constant(const_content) => 
                if n == 0 && !const_content.0.is_zero() {
                    Some(&const_content.0)
                } else {
                    None
                },
            Polynomial::Dense(dense_content) => dense_content.nth(n),
            Polynomial::Spares(spears_content) => spears_content.nth(n),
        }
    }
     
    /// Return dense expression of this polynomial if this is spears, otherwise (zero or constant) self.
    pub fn to_dense(self) -> Polynomial<C> {
        match self {
            Polynomial::Spares(mut sc) => {
                let n = sc.degree() + 1;
                let mut v: Vec<C> = Vec::with_capacity(n);
                for i in 0..n {
                    match sc.0.remove(&i) {
                        Some(e) => v.push(e),
                        None => v.push(C::zero()),
                    } 
                }
                debug_assert!(sc.0.is_empty());
                Polynomial::<C>::dense_from_vec(v)
            },
            _ => self,
        }
    }

    /// Return spears expression of this polynomial if this is dense, otherwise (zero or constant) self.
    pub fn to_spears(self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(mut dc) => {
                let mut map = BTreeMap::new();
                while let Some(e) = dc.0.pop() {
                    map.insert(dc.0.len(), e);
                }
                Polynomial::<C>::spears_from_map(map)
            },
            _ => self
        }
    }


    // the following methods are used?

    // iterate non zero terms
    pub fn iter<'a>(&'a self) -> PolyIter<'a, C> {
        match self {
            Polynomial::Zero() => PolyIter::Zero(),
            Polynomial::Constant(cc) => cc.iter(),
            Polynomial::Dense(dc) => dc.iter(),
            Polynomial::Spares(sc) => sc.iter(),
        }
    }

    // pub fn coeffs(&self) -> Vec<&C> {

    // } 

    // pub fn data(&self) -> Map<usize, &C> {

    // }

    pub fn split(&self) -> (Vec<usize>, Vec<&C>) {
        let n = self.degree();
        let mut es = Vec::with_capacity(n);
        let mut cs = Vec::with_capacity(n);

        self.iter().for_each(|(e, c)| {
            es.push(e);
            cs.push(c);
        });

        (es, cs)
    }
}

impl<C: Num + Clone> Polynomial<C> {

    /// Create dense clone of this polynomial if the self is spears, otherwise (zero, constant or dense) return self.
    pub fn dense_clone(&self) -> Polynomial<C> {
        match self {
            Polynomial::Spares(sc) => {
                let n: usize = sc.degree() + 1;
                let mut v: Vec<C> = Vec::with_capacity(n);
                for i in 0..n {
                    match sc.nth(i) {
                        Some(e) => v.push(e.clone()),
                        None => v.push(C::zero()),
                    } 
                }
                Polynomial::<C>::dense_from_vec(v)
            },
            _ => self.clone()
        }
    }

    /// Create spears clone of this polynomial if self is dense, otherwise (zero, constant or spears) return self.
    pub fn spears_clone(&self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => {
                let mut map = BTreeMap::new();
                for (i, e) in dc.iter() {
                    map.insert(i, e.clone());  // note e is not zero
                }
                Polynomial::<C>::spears_from_map(map)
            },
            _ => self.clone()
        }
    }

    pub fn terms<'a>(&'a self) -> impl Iterator<Item=Term<'a, C>> {
        self.iter().map(|(e, c)| Term::<C>::from_ref(e, c))
    }
}

impl<C: Num + Clone> Clone for Polynomial<C> {

    fn clone(&self) -> Self {
        match self {
            Self::Zero() => Self::Zero(),
            Self::Constant(cc) => Self::Constant(cc.clone()),
            Self::Dense(dc) => Self::Dense(dc.clone()),
            Self::Spares(sc) => Self::Spares(sc.clone()),
        }
    }
}

//   /**
//    * Returns the coefficients in little-endian order. So, the i-th element is coeffsArray(i) * (x ** i).
//    */
//   def coeffsArray(implicit ring: Semiring[C]): Array[C]

//   /**
//    * Returns a map from exponent to coefficient of this polynomial.
//    */
//   def data(implicit ring: Semiring[C], eq: Eq[C]): Map[Int, C] = {
//     val bldr = scala.collection.mutable.Map.newBuilder[Int, C]
//     foreachNonZero { (e, c) =>
//       bldr += ((e, c))
//     }
//     bldr.result().toMap
//   }

#[derive(Clone)]
pub struct ConstContent<C: Num>(C);

impl<C: Num> ConstContent<C> {

    fn iter<'a>(&'a self) -> PolyIter<'a, C> {
        PolyIter::Constant(ConstIter { is_visited: false, ref_to_const: &self.0 })

    }
}

pub enum PolyIter<'a, C: Num> {
    Zero(),
    Constant(ConstIter<'a, C>),
    Dense(DenseIter<'a, C>),
    Spears(SpearsIter<'a, C>),
}

impl<'a, C: Num> Iterator for PolyIter<'a, C> {
    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            PolyIter::Zero() => None,
            PolyIter::Constant(cc) => cc.next(),
            PolyIter::Dense(dc) => dc.next(),
            PolyIter::Spears(sc) => sc.next(),
        }
    }
}

pub struct ConstIter<'a, C: Num> {
    is_visited: bool,
    ref_to_const: &'a C
}

impl<'a, C: Num> Iterator for ConstIter<'a, C> {
    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_visited {
            self.is_visited = true;
            Some((0, self.ref_to_const))
        } else {
            None
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

#[macro_export]
macro_rules! dense {
    [ $t:ty; $( $x:expr ),* ] => {
        Polynomial::<$t>::dense_from_vec(vec![ $( $x ),* ])
    }
}

#[macro_export]
macro_rules! spears {
    [ $t:ty; $( $x:expr ),* ] => {
        Polynomial::<$t>::spears_from_map(std::collections::BTreeMap::from([ $( $x ),* ]))
    }
}

impl<C: Num> Zero for Polynomial<C> {

    fn zero() -> Self {
        Polynomial::<C>::ZERO
    }

    fn is_zero(&self) -> bool {
        match self {
            Polynomial::Zero() => true,
            _ => false,
        }
    }
}

impl<C: Num> ConstZero for Polynomial<C> {
    const ZERO: Self = Polynomial::Zero();
}

impl<C: Num> PartialEq for Polynomial<C> {

    fn eq(&self, other: &Self) -> bool {

        if self.degree() != other.degree() {
            return false;
        }

        fn has_the_same_content<C: Num>(vec: &Vec<C>, map: &BTreeMap<usize, C>) -> bool {
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
                Polynomial::Spares(rhs) => has_the_same_content(&lhs.0, &rhs.0),
                _ => false,
            },

            Polynomial::Spares(lhs) => match other {
                Polynomial::Dense(rhs) => has_the_same_content(&rhs.0, &lhs.0),
                Polynomial::Spares(rhs) => lhs.0 == rhs.0,
                _ => false,
            },
        }
    }
}

impl<C: Num> One for Polynomial<C> {

    fn one() -> Self { 
        Polynomial::<C>::constant(C::one())
    }
    
    fn is_one(&self) -> bool {
        match self {
            Polynomial::Constant(c) => c.0.is_one(),
            _ => false,
        }
    }
}

impl<C: Num + ConstOne> ConstOne for Polynomial<C> {

    const ONE: Self = Polynomial::Constant(ConstContent(C::ONE));
}

impl<C: Num> Add for Polynomial<C> {

    type Output = Polynomial<C>;

    fn add(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

impl<C: Num> Mul for Polynomial<C> {

    type Output = Polynomial<C>;

    fn mul(self, rhs: Self) -> Self::Output {
        todo!()
    }
}
 
impl<C: Num + Clone + Display> Display for Polynomial<C> {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        match self {
            Polynomial::Zero() => f.write_str("(0)"),

            Polynomial::Constant(c) => f.write_fmt(format_args!("({})", c.0)),
            
            _ => {
                let s: String = self.terms().map(|t|t.to_string()).collect();
                let first_sign = if s.starts_with(" - ") { "-" } else { "" }; 
                f.write_fmt(format_args!("{}{}", first_sign, &s[3..]))
            },
        }
    }
}
 
impl<C> Debug for Polynomial<C> where C: Num + Clone + Display{

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Polynomial::Zero() => "Zero",
            Polynomial::Constant(_) => "Constant",
            Polynomial::Dense(_) => "Dense",
            Polynomial::Spares(_) => "Spears",
        };
        f.write_fmt(format_args!("Polynomial::{}[{}]", s, self))
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
//    * Returns the coefficient of the n-th degree term.
//    */
//   def nth(n: Int)(implicit ring: Semiring[C]): C

//   /**
//    * Returns the term of the highest degree in this polynomial.
//    */
//   def maxTerm(implicit ring: Semiring[C]): Term[C] = Term(maxOrderTermCoeff, degree)

//   /**
//    * Returns the non-zero term of the minimum degree in this polynomial, unless it is zero. If this polynomial is zero,
//    * then this returns a zero term.
//    */
//   @nowarn
//   def minTerm(implicit ring: Semiring[C], eq: Eq[C]): Term[C] = {
//     foreachNonZero { (n, c) =>
//       return Term(c, n)
//     }
//     Term(ring.zero, 0)
//   }

//   /**
//    * Returns `true` iff this polynomial is constant.
//    */
//   def isConstant: Boolean =
//     degree == 0

//   /**
//    * Returns the degree of this polynomial.
//    */
//   def degree: Int

//   /**
//    * Returns the coefficient of max term of this polynomial.
//    */
//   def maxOrderTermCoeff(implicit ring: Semiring[C]): C

//   /**
//    * Returns a polynomial with the max term removed.
//    */
//   def reductum(implicit e: Eq[C], ring: Semiring[C], ct: ClassTag[C]): Polynomial[C]

//   /**
//    * Returns `true` if this polynomial is `ring.zero`.
//    */
//   def isZero: Boolean

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
//     foreachNonZero { (e, c) =>
//       val z: Polynomial[C] = y.pow(e) :* c
//       polynomial = polynomial + z
//     }
//     polynomial
//   }

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
//     this.foreachNonZero { (deg, c) =>
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
//    * Returns this polynomial as a monic polynomial, where the leading coefficient (ie. `maxOrderTermCoeff`) is 1.
//    */
//   def monic(implicit f: Field[C], eq: Eq[C]): Polynomial[C] = this :/ maxOrderTermCoeff

//   def derivative(implicit ring: Ring[C], eq: Eq[C]): Polynomial[C]
//   def integral(implicit field: Field[C], eq: Eq[C]): Polynomial[C]

//   /**
//    * Returns the number of sign variations in the coefficients of this polynomial. Given 2 consecutive terms (ignoring 0
//    * terms), a sign variation is indicated when the terms have differing signs.
//    */
//   def signVariations(implicit ring: Semiring[C], order: Order[C], signed: Signed[C]): Int = {
//     var prevSign: Sign = Signed.Zero
//     var variations = 0
//     foreachNonZero { (_, c) =>
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

//   def map[D: Semiring: Eq: ClassTag](f: C => D): Polynomial[D] =
//     mapTerms { case Term(c, n) => Term(f(c), n) }

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

//   def unary_-(implicit ring: Rng[C]): Polynomial[C]
//   def +(rhs: Polynomial[C])(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C]
//   def -(rhs: Polynomial[C])(implicit ring: Rng[C], eq: Eq[C]): Polynomial[C] = lhs + (-rhs)
//   def *(rhs: Polynomial[C])(implicit ring: Semiring[C], eq: Eq[C]): Polynomial[C]

//   def **(k: Int)(implicit ring: Rig[C], eq: Eq[C]): Polynomial[C] = pow(k)

//   def pow(k: Int)(implicit ring: Rig[C], eq: Eq[C]): Polynomial[C] = {
//     def loop(b: Polynomial[C], k: Int, extra: Polynomial[C]): Polynomial[C] =
//       if (k == 1)
//         b * extra
//       else
//         loop(b * b, k >>> 1, if ((k & 1) == 1) b * extra else extra)

//     if (k < 0) {
//       throw new IllegalArgumentException("negative exponent")
//     } else if (k == 0) {
//       Polynomial.one[C]
//     } else if (k == 1) {
//       this
//     } else {
//       loop(this, k - 1, this)
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

//   override def equals(that: Any): Boolean = that match {
//     case rhs: Polynomial[_] if lhs.degree == rhs.degree =>
//       val it1 = lhs.termsIterator
//       val it2 = rhs.termsIterator
//       @tailrec def loop(): Boolean = {
//         val has1 = it1.hasNext
//         val has2 = it2.hasNext
//         if (has1 && has2) {
//           if (it1.next() == it2.next()) loop() else false
//         } else has1 == has2
//       }
//       loop()

//     case rhs: Polynomial[_] =>
//       false

//     case n if lhs.isZero =>
//       n == 0

//     case n if lhs.degree == 0 =>
//       val (_, lcs) = Polynomial.split(lhs)
//       lcs(0) == n

//     case _ =>
//       false
//   }

// trait PolynomialOverRig[@sp(Double) C] extends PolynomialOverSemiring[C] with Rig[Polynomial[C]] {
//   implicit override val scalar: Rig[C]

//   def one: Polynomial[C] = Polynomial.one[C]
// }

// trait PolynomialOverRng[@sp(Double) C] extends PolynomialOverSemiring[C] with Rng[Polynomial[C]] {
//   implicit override val scalar: Rng[C]

//   def timesl(r: C, v: Polynomial[C]): Polynomial[C] = r *: v
//   def negate(x: Polynomial[C]): Polynomial[C] = -x
// }

// // we skip the CSemiring, CRng, CRig instances

// trait PolynomialOverCRing[@sp(Double) C]
//     extends CRing[Polynomial[C]]
//     with PolynomialOverRing[C]
//     with RingAssociativeAlgebra[Polynomial[C], C] {
//   implicit override val scalar: CRing[C]
// }

// trait PolynomialOverField[@sp(Double) C]
//     extends PolynomialOverRing[C]
//     with EuclideanRing[Polynomial[C]]
//     with VectorSpace[Polynomial[C], C]
//     with FieldAssociativeAlgebra[Polynomial[C], C] { self =>
//   implicit override val scalar: Field[C]

//   override def divr(x: Polynomial[C], k: C): Polynomial[C] = x :/ k

//   def euclideanFunction(x: Polynomial[C]): BigInt = x.degree
//   def equot(x: Polynomial[C], y: Polynomial[C]): Polynomial[C] = equotmod(x, y)._1
//   def emod(x: Polynomial[C], y: Polynomial[C]): Polynomial[C] = equotmod(x, y)._2
//   override def equotmod(x: Polynomial[C], y: Polynomial[C]): (Polynomial[C], Polynomial[C]) = {
//     require(!y.isZero, "Can't divide by polynomial of zero!")
//     (x: @unchecked) match {
//       case xd: poly.PolyDense[C] => poly.PolyDense.quotmodDense(xd, y)
//       case xs: poly.PolySparse[C] =>
//         val ys = (y: @unchecked) match {
//           case yd: poly.PolyDense[C]   => poly.PolySparse.dense2sparse(yd)
//           case ys1: poly.PolySparse[C] => ys1
//         }
//         poly.PolySparse.quotmodSparse(xs, ys)
//     }
//   }
// }