use crate::{algebra::algebra::Ring, polynomial::Polynomial};

#[derive(Clone)]
pub struct ConstContent<C>(pub(crate) C);

impl<C> ConstContent<C> {

    pub fn nth(&self, n: usize) -> Option<&C> {
        if n == 0 {  // note self.0 != 0
            Some(&self.0)
        } else {
            None
        }
    }
}

impl<C> ConstContent<C> where C: Ring {

    pub fn neg(self) -> Polynomial<C> {
        Polynomial::Constant(ConstContent(-self.0))
    }
}

impl<C> ConstContent<C> where C: Ring + Clone {

    pub fn neg_ref(&self) -> Polynomial<C> {
        Polynomial::Constant(ConstContent(-self.0.clone()))
    }
}