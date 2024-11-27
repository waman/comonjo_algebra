use num::Num;

use crate::polynomial::{CoeffsIter, IntoCoeffsIter, NonZeroCoeffsIter, IntoNonZeroCoeffsIter};

#[derive(Clone)]
pub struct ConstContent<C: Num>(pub(crate) C);

impl<C: Num> ConstContent<C> {

    pub fn nth(&self, n: usize) -> Option<&C> {
        if n == 0 {  // note self.0 != 0
            Some(&self.0)
        } else {
            None
        }
    }

    pub fn coeffs_iter<'a>(&'a self) -> CoeffsIter<'a, C> {
        CoeffsIter::Constant(ConstCoeffsIter { is_visited: false, ref_to_value: &self.0 })
    }

    pub fn into_coeffs_iter(self) -> IntoCoeffsIter<C> {
        IntoCoeffsIter::Constant(ConstIntoCoeffsIter { is_visited: false, value: Some(self.0) })
    }

    pub fn non_zero_coeffs_iter<'a>(&'a self) -> NonZeroCoeffsIter<'a, C> {
        NonZeroCoeffsIter::Constant(ConstNzcIter { is_visited: false, ref_to_value: &self.0 })
    }

    pub fn into_non_zero_coeffs_iter(self) -> IntoNonZeroCoeffsIter<C> {
        IntoNonZeroCoeffsIter::Constant(ConstIntoNzcIter(Some(self.0)))
    }
}

pub struct ConstCoeffsIter<'a, C: Num> {
    is_visited: bool,
    ref_to_value: &'a C
}

impl<'a, C: Num> Iterator for ConstCoeffsIter<'a, C> {
    type Item = Option<&'a C>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_visited {
            self.is_visited = true;
            Some(Some(self.ref_to_value))
        } else {
            None
        }
    }
}

pub struct ConstIntoCoeffsIter<C: Num> {
    is_visited: bool,
    value: Option<C>
}

impl<C: Num> Iterator for ConstIntoCoeffsIter<C> {
    type Item = C;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_visited {
            self.is_visited = true;
            Some(self.value.take().unwrap())
        } else {
            None
        }
    }
}

pub struct ConstNzcIter<'a, C: Num> {
    is_visited: bool,
    ref_to_value: &'a C
}

impl<'a, C: Num> Iterator for ConstNzcIter<'a, C> {
    type Item = (usize, &'a C);

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_visited {
            self.is_visited = true;
            Some((0, self.ref_to_value))
        } else {
            None
        }
    }
}

pub struct ConstIntoNzcIter<C: Num>(Option<C>);

impl<C: Num> Iterator for ConstIntoNzcIter<C> {
    type Item = (usize, C);

    fn next(&mut self) -> Option<Self::Item> {
        if self.0.is_some() {
            Some((0, self.0.take().unwrap()))
        } else {
            None
        }
    }
}