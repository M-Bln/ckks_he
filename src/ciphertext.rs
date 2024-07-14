use crate::algebra::arithmetic::{Rescale, RingMod};
use crate::algebra::big_int::BigInt;
use crate::algebra::cyclotomic_ring::CyclotomicRing;
use bnum::types::I256;
use std::ops::{Add, Mul};

//type CiphertextRing = CyclotomicRing<RingMod<I256>>;
pub type CiphertextRing<T: BigInt> = CyclotomicRing<RingMod<T>>;

struct RawCiphertext<T: BigInt>(CiphertextRing<T>, CiphertextRing<T>);

// Should I also keep a ref to public / eval key?
pub struct Ciphertext<T: BigInt> {
    raw: RawCiphertext<T>,
    level: u32,
    upper_bound_message: T,
    upper_bound_error: T,
}

impl<'a, T: BigInt> Add<&'a RawCiphertext<T>> for RawCiphertext<T> {
    type Output = RawCiphertext<T>;
    fn add(self, other: &RawCiphertext<T>) -> RawCiphertext<T> {
        RawCiphertext::<T>(self.0 + &other.0, self.1 + &other.1)
    }
}

impl<'a, T: BigInt> Mul<&'a RawCiphertext<T>> for CiphertextRing<T> {
    type Output = RawCiphertext<T>;
    fn mul(self, other: &'a RawCiphertext<T>) -> RawCiphertext<T> {
        RawCiphertext(self.clone() * &other.0, self * &other.1)
    }
}

impl<T: BigInt> Rescale<T> for RawCiphertext<T> {
    fn rescale(&mut self, scalar: T) {
        self.0.rescale(scalar);
        self.1.rescale(scalar);
    }
}

impl<T: BigInt> RawCiphertext<T> {
    pub fn mul(&self, other: &RawCiphertext<T>, evk: &RawCiphertext<T>, P: T) -> RawCiphertext<T> {
        let (d_0, d_1, d_2) = (
            self.0.clone() * &other.0,
            self.0.clone() * &other.1 + &(self.1.clone() * &other.0),
            self.1.clone() * &other.1,
        );
        let mut summand = d_2 * evk;
        summand.rescale(P);
        RawCiphertext(d_0, d_1) + &summand
    }
}