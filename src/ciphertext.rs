use crate::algebra::arithmetic::{Rescale, RingMod};
use crate::algebra::big_int::BigInt;
use crate::algebra::cyclotomic_ring::CyclotomicRing;
use bnum::types::I256;
use std::ops::{Add, Mul};

//type CiphertextRing = CyclotomicRing<RingMod<I256>>; pub type
pub type Message<T: BigInt> = CyclotomicRing<RingMod<T>>;

#[derive(Clone, Debug)]
pub struct RawCiphertext<T: BigInt>(pub Message<T>, pub Message<T>);

// Should I also keep a ref to public / eval key?
#[derive(Clone, Debug)]
pub struct Ciphertext<T: BigInt> {
    pub raw: RawCiphertext<T>,
    pub level: u32,
    pub upper_bound_message: f64,
    pub upper_bound_error: f64,
}

impl<T: BigInt> Ciphertext<T> {
    pub fn new(
        raw: RawCiphertext<T>,
        level: u32,
        upper_bound_message: f64,
        upper_bound_error: f64,
    ) -> Self {
        Ciphertext {
            raw,
            level,
            upper_bound_message,
            upper_bound_error,
        }
    }
}

impl<'a, T: BigInt> Add<&'a RawCiphertext<T>> for RawCiphertext<T> {
    type Output = RawCiphertext<T>;
    fn add(self, other: &RawCiphertext<T>) -> RawCiphertext<T> {
        RawCiphertext::<T>(self.0 + &other.0, self.1 + &other.1)
    }
}

impl<'a, T: BigInt> Mul<&'a RawCiphertext<T>> for Message<T> {
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

impl<T: BigInt> RawCiphertext<T> {
    /// Multiply the raw ciphertext by the message other. However the modulus of the ciphertext
    /// remains unchanged, even if message has a smaller modulus.
    pub fn scalar_mul_keep_modulus(&self, other: &Message<T>) -> RawCiphertext<T> {
        self.scalar_mul_integer(&other.to_integer())
    }
    pub fn scalar_mul_integer(&self, other: &CyclotomicRing<T>) -> RawCiphertext<T> {
        let other_mod = other.modulo(self.0.ref_coefficients()[0].modulus);
        RawCiphertext(self.0.clone() * &other_mod, self.1.clone() * &other_mod)
    }
}
