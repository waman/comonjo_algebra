use std::ops::Neg;

use num::Num;

use crate::polynomial::Polynomial;

#[derive(Clone)]
pub struct DenseContent<C: Num>(pub(crate) Vec<C>);

impl<C: Num> DenseContent<C> {

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    pub fn nth(&self, n: usize) -> Option<&C> {
        self.0.get(n)
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
}

impl<C> DenseContent<C> where C: Num + Neg<Output=C>{

    pub fn neg(self) -> Polynomial<C> {
        let v: Vec<C> = self.0.into_iter().map(|e|-e).collect();
        Polynomial::Dense(DenseContent(v))
    }
}

impl<C> DenseContent<C> where C: Num + Clone + Neg<Output=C> {

    pub fn neg_ref(&self) -> Polynomial<C> {
        let v: Vec<C> = self.0.iter().map(|e|-e.clone()).collect();
        Polynomial::Dense(DenseContent(v))
    }
}


/// Return (max, min, first_arg_is_longer)
fn max_min(x: usize, y: usize) -> (usize, usize, bool) {
    if x >= y { (x, y, true) } else { (y, x, false) }
}

pub(crate) fn add_dense<C: Num>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> {
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

pub(crate) fn add_dense_ref<'a, 'b, C: Num + Clone>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> {
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
    v.extend(rest_iter.map(|c| match c {
        Some(x) => x.clone(),
        None => C::zero(),
    }));

    Polynomial::dense_from_vec(v)
}

pub(crate) fn sub_dense<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C> 
        where C: Num + Neg<Output=C>{

    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.into_coeffs_iter();
    let mut rhs_iter = rhs.into_coeffs_iter();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(x), Some(y)) => v.push(x - y),
            _ => panic!(),
        }
    }

    if lhs_is_longer {
        v.extend(lhs_iter);
    } else {
        v.extend(rhs_iter.map(|c|-c));
    }

    Polynomial::dense_from_vec(v)
}

pub(crate) fn sub_dense_ref<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C> 
        where C: Num + Clone + Neg<Output=C>{
    let (d_max, d_min, lhs_is_longer) = max_min(lhs.degree(), rhs.degree());

    let mut v: Vec<C> = Vec::with_capacity(d_max);

    let mut lhs_iter = lhs.coeffs_iter();
    let mut rhs_iter = rhs.coeffs_iter();

    for _ in 0..=d_min {
        match (lhs_iter.next(), rhs_iter.next()) {
            (Some(some_x), Some(some_y)) => match (some_x, some_y) {
                (Some(x), Some(y)) => v.push(x.clone() - y.clone()),
                (Some(x), None) => v.push(x.clone()),
                (None, Some(y)) => v.push(y.clone()),
                (None, None) => v.push(C::zero()),
            },
            _ => panic!(),
        }
    }

    if lhs_is_longer {
        v.extend(lhs_iter.map(|c| match c {
            Some(x) => x.clone(),
            None => C::zero(),
        }));
    } else {
        v.extend(rhs_iter.map(|c| match c {
            Some(x) => -x.clone(),
            None => C::zero(),
        }));
    }

    Polynomial::dense_from_vec(v)
}

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

// def unary_-(implicit ring: Rng[C]): Polynomial[C] = {
//   val negArray = new Array[C](coeffs.length)
//   cfor(0)(_ < coeffs.length, _ + 1) { i => negArray(i) = -coeffs(i) }
//   new PolyDense(negArray)
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

// final private[math] def quotmodDense[@sp(Double) C: Field: Eq: ClassTag](lhs: PolyDense[C],
//                                                                          rhs: Polynomial[C]
// ): (Polynomial[C], Polynomial[C]) = {
//   def zipSum(lcs: Array[C], rcs: Array[C]): Array[C] =
//     (lcs + rcs).tail

//   def polyFromCoeffsLE(cs: Array[C]): Polynomial[C] =
//     Polynomial.dense(cs)

//   def polyFromCoeffsBE(cs: Array[C]): Polynomial[C] = {
//     val ncs = cs.dropWhile(_ === Field[C].zero)
//     Polynomial.dense(ncs.reverse)
//   }

//   @tailrec def eval(q: Array[C], u: Array[C], n: Int): (Polynomial[C], Polynomial[C]) = {
//     if (u.isEmpty || n < 0) {
//       (polyFromCoeffsLE(q), polyFromCoeffsBE(u))
//     } else {
//       val v0 = if (rhs.isZero) Field[C].zero else rhs.maxOrderTermCoeff
//       val q0 = u(0) / v0
//       val uprime = zipSum(u, rhs.coeffsArray.reverse.map(_ * -q0))
//       eval(Array(q0) ++ q, uprime, n - 1)
//     }
//   }

//   val cs = rhs.coeffsArray
//   if (cs.length == 0) {
//     throw new ArithmeticException("/ by zero polynomial")
//   } else if (cs.length == 1) {
//     val c = cs(0)
//     val q = Polynomial.dense(lhs.coeffs.map(_ / c))
//     val r = Polynomial.dense(new Array[C](0))
//     (q, r)
//   } else {
//     eval(new Array[C](0), lhs.coeffs.reverse, lhs.degree - rhs.degree)
//   }
// }

// }