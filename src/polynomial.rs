use std::{collections::{BTreeMap, btree_map}, fmt::{Debug, Display}, iter::Enumerate, ops::*, slice::Iter};

use num::{traits::{ConstOne, ConstZero}, Num, One, Zero};

use crate::poly::{const_content::ConstContent, dense_content::DenseContent, sparse_content::SparseContent, term::Term};

/// Refer to spire's <a href="https://github.com/typelevel/spire/blob/main/core/src/main/scala/spire/math/Polynomial.scala">Polynomial</a>
pub enum Polynomial<C: Num> {

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

impl<C: Num> Polynomial<C> {

    pub fn constant(c: C) -> Polynomial<C> {
        if c.is_zero() {
            Polynomial::Zero()
        } else {
            Polynomial::Constant(ConstContent(c))
        }
    }

    pub fn dense_from_vec(mut coeffs: Vec<C>) -> Polynomial<C>  {

        while !coeffs.is_empty() && coeffs.last().unwrap().is_zero() {
            coeffs.pop();
        }

        match coeffs.len() {
            0 => Polynomial::<C>::constant(C::zero()),
            1 => Polynomial::<C>::constant(coeffs.pop().unwrap()),
            _ => {
                coeffs.shrink_to_fit();
                Polynomial::Dense(DenseContent(coeffs))
            }
        }
    }

    pub fn sparse_from_map(mut coeffs: BTreeMap<usize, C>) -> Polynomial<C> {

        coeffs.retain(|_, v| !v.is_zero());

        match coeffs.len() {
            0 => Polynomial::<C>::constant(C::zero()),
            1 => {
                let mut key = 0_usize;
                for k in coeffs.keys() { key = *k; }
                let value = coeffs.remove(&key).unwrap();
                Polynomial::<C>::constant(value)
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
            Polynomial::Constant(cc) => cc.nth(n),
            Polynomial::Dense(dc) => dc.nth(n),
            Polynomial::Sparse(sc) => sc.nth(n),
        }
    }
     
    /// Return dense expression of this polynomial if this is sparse, otherwise (zero or constant) self.
    pub fn to_dense(self) -> Polynomial<C> {
        match self {
            Polynomial::Sparse(mut sc) => {
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

    /// Return sparse expression of this polynomial if this is dense, otherwise (zero or constant) self.
    pub fn to_sparse(self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(mut dc) => {
                let mut map = BTreeMap::new();
                while let Some(e) = dc.0.pop() {
                    map.insert(dc.0.len(), e);
                }
                Polynomial::<C>::sparse_from_map(map)
            },
            _ => self
        }
    }

    //********** Iteration **********/
    pub fn coeffs_iter<'a>(&'a self) -> CoeffsIter<'a, C> {
        match self {
            Polynomial::Zero() => CoeffsIter::Zero(),
            Polynomial::Constant(cc) => CoeffsIter::Constant(Some(Some(&cc.0))),
            Polynomial::Dense(dc) => CoeffsIter::Dense(dc.0.iter()),
            Polynomial::Sparse(sc) => {
                let mut map_iter = sc.0.iter();
                CoeffsIter::Sparse{ index: 0, current: map_iter.next(), map_iter}
            },
        }
    }

    pub fn into_coeffs_iter(self) -> IntoCoeffsIter<C> {
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

    // iterate coefficients of non zero terms
    pub fn non_zero_coeffs_iter<'a>(&'a self) -> NonZeroCoeffsIter<'a, C> {
        match self {
            Polynomial::Zero() => NonZeroCoeffsIter::Zero(),
            Polynomial::Constant(cc) => NonZeroCoeffsIter::Constant(Some(&cc.0)),
            Polynomial::Dense(dc) => NonZeroCoeffsIter::Dense(dc.0.iter().enumerate()),
            Polynomial::Sparse(sc) => NonZeroCoeffsIter::Sparse(sc.0.iter()),
        }
    }

    pub fn into_non_zero_coeffs_iter(self) -> IntoNonZeroCoeffsIter<C> {
        match self {
            Polynomial::Zero() => IntoNonZeroCoeffsIter::Zero(),
            Polynomial::Constant(cc) => IntoNonZeroCoeffsIter::Constant(Some(cc.0)),
            Polynomial::Dense(dc) => IntoNonZeroCoeffsIter::Dense(dc.0.into_iter().enumerate()),
            Polynomial::Sparse(sc) => IntoNonZeroCoeffsIter::Sparse(sc.0.into_iter()),
        }
    }


    // the following methods are used?
    pub fn max_order_term_coeff(&self) -> Option<&C> {
        match self {
            Polynomial::Zero() => None,
            _ => self.nth(self.degree())
        }
    }

    pub fn split(&self) -> (Vec<usize>, Vec<&C>) {
        let n = self.degree();
        let mut es = Vec::with_capacity(n);
        let mut cs = Vec::with_capacity(n);

        self.non_zero_coeffs_iter().for_each(|(e, c)| {
            es.push(e);
            cs.push(c);
        });

        (es, cs)
    }

    /// Returns a polynomial with the max term removed.
    pub fn reductum(&self) -> Polynomial<C> {
        match self {
            Polynomial::Dense(dc) => dc.reductum(),
            Polynomial::Sparse(sc) => sc.reductum(),
            _ => Polynomial::Zero(),
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

impl<C: Num + Clone> Polynomial<C> {

    /// Create dense clone of this polynomial if the self is sparse, otherwise (zero, constant or dense) return self.
    pub fn dense_clone(&self) -> Polynomial<C> {
        match self {
            Polynomial::Sparse(sc) => {
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

    /// Create sparse clone of this polynomial if self is dense, otherwise (zero, constant or sparse) return self.
    pub fn sparse_clone(&self) -> Polynomial<C> {
        match self {
            d @ Polynomial::Dense(_) => {
                let mut map = BTreeMap::new();
                for (i, e) in d.non_zero_coeffs_iter() {
                    map.insert(i, e.clone());  // note e is not zero
                }
                Polynomial::<C>::sparse_from_map(map)
            },
            _ => self.clone()
        }
    }

    pub fn terms<'a>(&'a self) -> impl Iterator<Item=Term<'a, C>> {
        self.non_zero_coeffs_iter().map(|(e, c)| Term::<C>::from_ref(e, c))
    }

    pub fn max_term<'a>(&'a self) -> Term<'a, C> {
        match self.max_order_term_coeff() {
            Some(c) => Term::<C>::from_ref(self.degree(), c),
            _ => Term::<C>::from_value(0, C::zero()),
        }
    }

    pub fn min_term<'a>(&'a self) -> Term<'a, C> {
        match self.terms().next() {
            Some(t) => t,
            _ => Term::<C>::from_value(0, C::zero()),
        }
    }

    pub fn integral(&self) -> Polynomial<C> {
        match self {
            Polynomial::Zero() => Polynomial::ZERO,
            Polynomial::Constant(cc) => Polynomial::<C>::linear_monomial(cc.0.clone()),
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
    //     foreachNonZero { (e, c) =>
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

impl<C: Num + Clone> Clone for Polynomial<C> {

    fn clone(&self) -> Self {
        match self {
            Self::Zero() => Self::Zero(),
            Self::Constant(cc) => Self::Constant(cc.clone()),
            Self::Dense(dc) => Self::Dense(dc.clone()),
            Self::Sparse(sc) => Self::Sparse(sc.clone()),
        }
    }
}
 
//********** Display and Debug **********
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
            Polynomial::Sparse(_) => "Sparse",
        };
        f.write_fmt(format_args!("Polynomial::<{}>::{}[{}]", std::any::type_name::<C>(), s, self))
    }
}

//********** Iterator **********
/// <code>next()</code> method returns <code>Some(None)</code> or <code>Some(Some(0))</code> if the coefficient is zero.
pub enum CoeffsIter<'a, C: Num> {
    Zero(),
    Constant(Option<Option<&'a C>>),
    Dense(Iter<'a, C>),
    Sparse{
        index: usize,
        current: Option<(&'a usize, &'a C)>,
        map_iter: btree_map::Iter<'a, usize, C>,
    },
}

impl<'a, C: Num> Iterator for CoeffsIter<'a, C> {

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
            CoeffsIter::Sparse{ index, current, map_iter } => 
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
}

pub enum IntoCoeffsIter<C: Num> {
    Zero(),
    Constant(Option<C>),
    Dense(std::vec::IntoIter<C>),
    Sparse{
        index: usize,
        current: Option<(usize, C)>,
        map_iter: btree_map::IntoIter<usize, C>
    },
}

impl<C: Num> Iterator for IntoCoeffsIter<C> {
    type Item = C;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoCoeffsIter::Zero() => None,
            IntoCoeffsIter::Constant(value) => value.take(),
            IntoCoeffsIter::Dense(vec_iter) =>  vec_iter.next(),
            IntoCoeffsIter::Sparse{ index, current, map_iter } => {
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
}

pub enum NonZeroCoeffsIter<'a, C: Num> {
    Zero(),
    Constant(Option<&'a C>),
    Dense(Enumerate<Iter<'a, C>>),
    Sparse(btree_map::Iter<'a, usize, C>),
}

impl<'a, C: Num> Iterator for NonZeroCoeffsIter<'a, C> {
    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            NonZeroCoeffsIter::Zero() => None,
            NonZeroCoeffsIter::Constant(value) => 
                match value.take() {
                    Some(c) => Some((0, c)),
                    None => None,
                },
            NonZeroCoeffsIter::Dense(vec_iter) => 
                loop {
                    match vec_iter.next() {
                        Some((e, c)) =>  if !c.is_zero() { return Some((e, c)); },
                        None => return None,
                    }
                },
            NonZeroCoeffsIter::Sparse(map_iter) => map_iter.next().map(|(e, c)| (*e, c)),
        }
    }
}

pub enum IntoNonZeroCoeffsIter<C: Num> {
    Zero(),
    Constant(Option<C>),
    Dense(Enumerate<std::vec::IntoIter<C>>),
    Sparse(btree_map::IntoIter<usize, C>),
}

impl<C: Num> Iterator for IntoNonZeroCoeffsIter<C> {
    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IntoNonZeroCoeffsIter::Zero() => None,
            IntoNonZeroCoeffsIter::Constant(value) => 
                match value.take() {
                    Some(c) => Some((0, c)),
                    None => None,
                },
            IntoNonZeroCoeffsIter::Dense(vec_iter) => 
                loop {
                    match vec_iter.next() {
                        Some((e, c)) =>  if !c.is_zero() { return Some((e, c)); },
                        None => return None,
                    }
                },
            IntoNonZeroCoeffsIter::Sparse(map_iter) => map_iter.next(),
        }
    }
}

//********** Eq, PartialEq, Zero and One **********
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

impl<C: Num + Eq> Eq for Polynomial<C> {}

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

impl<C: Num> One for Polynomial<C> {

    fn one() -> Self { 
        Polynomial::Constant(ConstContent(C::one()))  // raw creation
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

//********** NumOp **********
impl<C> Neg for Polynomial<C> where C: Num + Neg<Output=C>{

    type Output = Polynomial<C>;

    fn neg(self) -> Self::Output {
        match self {
            Self::Zero() => Polynomial::Zero(),
            Self::Constant(cc) => cc.neg(),
            Self::Dense(dd) => dd.neg(),
            Self::Sparse(sc) => sc.neg(),
        }
    }
}

impl<'a, C> Neg for &'a Polynomial<C> where C: Num + Clone + Neg<Output=C> {

    type Output = Polynomial<C>;

    fn neg(self) -> Self::Output {
        match self {
            Polynomial::Zero() => Polynomial::Zero(),
            Polynomial::Constant(cc) => cc.neg_ref(),
            Polynomial::Dense(dd) => dd.neg_ref(),
            Polynomial::Sparse(sc) => sc.neg_ref(),
        }
    }
}

impl<C: Num> Add for Polynomial<C> {

    type Output = Polynomial<C>;

    fn add(self, other: Self) -> Self::Output {
        match (self, other) {
            (lhs, Polynomial::Zero()) => lhs,
            (Polynomial::Zero(), rhs) => rhs,
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) => Polynomial::constant(lhs.0 + rhs.0),
            (lhs @ Polynomial::Dense(_), rhs) | 
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => add_dense(lhs, rhs),
            (lhs @ Polynomial::Sparse(_), rhs) |
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => add_sparse(lhs, rhs),
        }
    }
}

/// Return (max, min, first_arg_is_longer)
fn max_min(x: usize, y: usize) -> (usize, usize, bool) {
    if x >= y { (x, y, true) } else { (y, x, false) }
}

fn add_dense<C: Num>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> {
    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.into_coeffs_iter();
    let mut rhs_iter = rhs.into_coeffs_iter();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(x), Some(y)) => v.push(x + y),
            _ => panic!(),
        }
    }

    let rest_iter = if lhs_is_longer { lhs_iter } else { rhs_iter };
    v.extend(rest_iter);

    Polynomial::dense_from_vec(v)
}

fn add_sparse<C: Num>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> {
    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.into_non_zero_coeffs_iter();
    let mut rhs_iter = rhs.into_non_zero_coeffs_iter();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            map.insert(x.0, x.1 + y.1);

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1);
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(lhs_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1);
        map.extend(rhs_iter);
    }

    Polynomial::sparse_from_map(map)
} 

impl<'a, C: Num + Clone> Add for &'a Polynomial<C> {

    type Output = Polynomial<C>;

    fn add(self, other: Self) -> Self::Output {
        match (self, other) {
            (lhs, Polynomial::Zero()) => lhs.clone(),
            (Polynomial::Zero(), rhs) => rhs.clone(),
            (Polynomial::Constant(lhs), Polynomial::Constant(rhs)) =>
                Polynomial::constant(lhs.0.clone() + rhs.0.clone()),
            (lhs @ Polynomial::Dense(_), rhs) | 
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Dense(_)) => add_dense_ref(lhs, rhs),
            (lhs @ Polynomial::Sparse(_), rhs) |
            (lhs @ Polynomial::Constant(_), rhs @ Polynomial::Sparse(_)) => add_sparse_ref(lhs, rhs),
        }
    }
}

fn add_dense_ref<'a, 'b, C: Num + Clone>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> {
    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs_iter();
    let mut rhs_iter = rhs.coeffs_iter();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(some_x), Some(some_y)) => match (some_x, some_y) {
                (Some(x), Some(y)) => v.push(x.clone() + y.clone()),
                (Some(x), None) => v.push(x.clone()),
                (None, Some(y)) => v.push(y.clone()),
                (None, None) => v.push(C::zero()),
            },
            _ => panic!(),
        }
    }

    let rest_iter = if lhs_is_longer { lhs_iter } else { rhs_iter };
    for c in rest_iter {
        match c {
            Some(x) => v.push(x.clone()),
            None => v.push(C::zero()),
        }
    }
    // v.extend(rest_iter.map(|c| match c {
    //     Some(x) => x.clone(),
    //     None => C::zero(),
    // }));

    Polynomial::dense_from_vec(v)
}

fn add_sparse_ref<'a, 'b, C: Num + Clone>(lhs: &'a Polynomial<C>, rhs: &'a Polynomial<C>) -> Polynomial<C> {
    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.non_zero_coeffs_iter();
    let mut rhs_iter = rhs.non_zero_coeffs_iter();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            map.insert(x.0, x.1.clone() + y.1.clone());

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.clone());
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(lhs_iter.map(|(e, c)|(e, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.clone());
        map.extend(rhs_iter.map(|(e, c)|(e, c.clone())));
    }

    Polynomial::sparse_from_map(map)
} 
 
//    final private def addSparse[C: Eq: Semiring: ClassTag](lhs: PolySparse[C], rhs: PolySparse[C]): PolySparse[C] = {
//      val PolySparse(lexp, lcoeff) = lhs
//      val PolySparse(rexp, rcoeff) = rhs
 
//      val len = countSumTerms(lhs, rhs)
//      val es = new Array[Int](len)
//      val cs = new Array[C](len)
 
//      @tailrec
//      def sum(i: Int, j: Int, k: Int): PolySparse[C] =
//        if (i < lexp.length && j < rexp.length) {
//          val ei = lexp(i)
//          val ej = rexp(j)
//          if (ei == ej) {
//            es(k) = ei
//            cs(k) = lcoeff(i) + rcoeff(j)
//            sum(i + 1, j + 1, k + 1)
//          } else if (ei < ej) {
//            es(k) = ei
//            cs(k) = lcoeff(i)
//            sum(i + 1, j, k + 1)
//          } else {
//            es(k) = ej
//            cs(k) = rcoeff(j)
//            sum(i, j + 1, k + 1)
//          }
//        } else {
//          var k0 = k
//          cfor(i)(_ < lexp.length, _ + 1) { i0 =>
//            es(k0) = lexp(i0)
//            cs(k0) = lcoeff(i0)
//            k0 += 1
//          }
//          cfor(j)(_ < rexp.length, _ + 1) { j0 =>
//            es(k0) = rexp(j0)
//            cs(k0) = rcoeff(j0)
//            k0 += 1
//          }
//          PolySparse.safe(es, cs)
//        }
 
//      sum(0, 0, 0)
//    }

impl<C: Num + Neg<Output=C>> Sub for Polynomial<C>
        where C: Num, for<'a> &'a C: Neg<Output=C>{

    type Output = Polynomial<C>;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<C: Num> Mul for Polynomial<C> {

    type Output = Polynomial<C>;

    fn mul(self, other: Self) -> Self::Output {
        todo!()
    }
}

impl<C: Num> Div for Polynomial<C> {

    type Output = Polynomial<C>;

    fn div(self, other: Self) -> Self::Output {
        todo!()
    }
}

impl<C: Num> Rem for Polynomial<C> {

    type Output = Polynomial<C>;

    fn rem(self, other: Self) -> Self::Output {
        todo!()
    }
}

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