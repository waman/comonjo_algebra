use std::ops::Neg;

use num::Num;

use crate::polynomial::Polynomial;

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
}

impl<C> ConstContent<C> where C: Num + Neg<Output=C>{

    pub fn neg(self) -> Polynomial<C> {
        Polynomial::Constant(ConstContent(-self.0))
    }
}

impl<C> ConstContent<C> where C: Num + Clone + Neg<Output=C> {

    pub fn neg_ref(&self) -> Polynomial<C> {
        Polynomial::Constant(ConstContent(-self.0.clone()))
    }
}