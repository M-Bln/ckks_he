use crate::algebra::big_int::BigInt;
use crate::ciphertext::{CiphertextRing, RawCiphertext};

#[derive(Clone, Debug)]
pub struct EvaluationKey<T: BigInt> {
    pub dimension_exponent: u32, // Ciphertexts in cyclotomic ring $R[X]/(1+X^N)$ with N = 2^dimension_exponent
    pub mul_scaling: T,          // Rescaling factor used in homomorphic multiplication
    pub q_0: T,                  // minimal modulus
    pub q: T, // modulus per level, the total modulus in coefficients ring is initialy mul_scaling * q_0 * q^{level_max}
    pub level_max: u32,
    pub variance: f64, // standard deviation of the Gaussian distribution of error sampling
    pub raw_key: RawCiphertext<T>,
}

impl<T: BigInt> EvaluationKey<T> {
    pub fn new(
        dimension_exponent: u32,
        mul_scaling: T,
        q_0: T,
        q: T,
        level_max: u32,
        variance: f64,
        raw_key: RawCiphertext<T>,
    ) -> Self {
        EvaluationKey {
            dimension_exponent,
            mul_scaling,
            q_0,
            q,
            level_max,
            variance,
            raw_key,
        }
    }
}
