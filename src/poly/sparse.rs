use std::collections::BTreeMap;

use crate::{algebra::{Field, Ring, Semiring}, poly::{CoeffsIterator, Polynomial}};

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

//********** Add **********/
pub(crate) fn add_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Semiring {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1 + y.1;
            if !sum.is_zero() { map.insert(x.0, sum); }

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

pub(crate) fn add_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1 + y.1;
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

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
        map.insert(x.0, x.1);
        map.extend(lhs_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.clone());
        map.extend(rhs_iter.map(|(i, c)| (i, c.clone())));
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn add_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1.ref_add(y.1);
            if !sum.is_zero() { map.insert(x.0, sum); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

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
        map.insert(x.0, x.1.clone());
        map.extend(lhs_iter.map(|(i, c)| (i, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1);
        map.extend(rhs_iter);
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn add_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'a Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let sum = x.1.ref_add(y.1);
            if !sum.is_zero() { map.insert(x.0, sum); }

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

//********** Sub **********/
pub(crate) fn sub_vv<C>(lhs: Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Ring {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1 - y.1;
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, -y.1);
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(lhs_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, -y.1);
        map.extend(rhs_iter.map(|(i, c)| (i, -c)));
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn sub_vr<'b, C>(lhs: Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1 - y.1;
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1);

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.ref_neg());
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1);
        map.extend(lhs_iter);
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.ref_neg());
        map.extend(rhs_iter.map(|(i, c)| (i, c.ref_neg())));
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn sub_rv<'a, C>(lhs: &'a Polynomial<C>, rhs: Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1.ref_sub(y.1);
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, -y.1);
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(lhs_iter.map(|(i, c)| (i, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, -y.1);
        map.extend(rhs_iter.map(|(i, c)| (i, -c)));
    }

    Polynomial::sparse_from_map(map)
} 

pub(crate) fn sub_rr<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'a Polynomial<C>) -> Polynomial<C>
        where C: Ring + Clone {

    let mut map = BTreeMap::new();

    let mut lhs_iter = lhs.nonzero_coeffs();
    let mut rhs_iter = rhs.nonzero_coeffs();

    let mut x_next = lhs_iter.next();
    let mut y_next = rhs_iter.next();

    loop {
        let x = x_next.unwrap();
        let y = y_next.unwrap();

        if x.0 == y.0 {
            let dif = x.1.ref_sub(y.1);
            if !dif.is_zero() { map.insert(x.0, dif); }

            x_next = lhs_iter.next();
            y_next = rhs_iter.next();
            if x_next.is_none() || y_next.is_none() { break; }

        } else if x.0 < y.0 {
            map.insert(x.0, x.1.clone());

            x_next = lhs_iter.next();
            y_next = Some(y);
            if x_next.is_none() { break; }

        } else {  // x.0 > y.0
            map.insert(y.0, y.1.ref_neg());
            
            x_next = Some(x);
            y_next = rhs_iter.next();
            if y_next.is_none() { break; }
        }
    }

    if let Some(x) = x_next {
        map.insert(x.0, x.1.clone());
        map.extend(lhs_iter.map(|(e, c)|(e, c.clone())));
    } else if let Some(y) = y_next {
        map.insert(y.0, y.1.ref_neg());
        map.extend(rhs_iter.map(|(e, c)|(e, c.ref_neg())));
    }

    Polynomial::sparse_from_map(map)
} 

//********** Mul **********/
pub(crate) fn mul<'a, 'b, C>(lhs: &'a Polynomial<C>, rhs: &'b Polynomial<C>) -> Polynomial<C>
        where C: Semiring + Clone {
    let mut map: BTreeMap<usize, C> = BTreeMap::new();

    for (i, x) in lhs.nonzero_coeffs() {
        for (j, y) in rhs.nonzero_coeffs() {
            let k = i + j;
            let z = x.ref_mul(y);
            match map.get_mut(&k){
                Some(c) => *c = c.ref_add(z),
                None => { map.insert(k, z); },
            }
        }
    }

    Polynomial::sparse_from_map(map)
}

//********** Div & Rem **********/
pub(crate) fn div_rem<C>(mut u: BTreeMap<usize, C>, rhs: &Polynomial<C>) -> (Polynomial<C>, Polynomial<C>) 
        where C: Field + Clone {

    let d_rhs = rhs.degree();
    let v0: &C = rhs.max_order_term_coeff().unwrap();
    let mut q: BTreeMap<usize, C> = BTreeMap::new();

    while !u.is_empty() {
        let i_last = *u.last_key_value().unwrap().0;
        if i_last < d_rhs { break; }
        let c_last = u.remove(&i_last).unwrap();

        let q0: C = c_last.ref_div(v0);
        let offset = i_last - d_rhs;  // the last of u is already removed
        for (i, c_rhs) in rhs.nonzero_coeffs() {
            if i == d_rhs { break; }
            let e = u.entry(offset + i).or_insert(C::zero());
            *e = e.ref_sub(c_rhs.ref_mul(&q0));
        }
        q.insert(i_last - d_rhs, q0);
    }

    (Polynomial::sparse_from_map(q), Polynomial::sparse_from_map(u))
}