use crate::algebra::big_int::{BigInt};
use crate::algebra::polynomial::{Polynomial};
use crate::ciphertext::CiphertextRing;
use crate::random_distributions::HWTDistribution;

pub struct SecretKey<T: BigInt> {
    dimension_exponent: u32, // Ciphertexts in cyclotomic ring $R[X]/(1+X^N)$ with N = 2^dimension_exponent
    hamming_weight: usize, // Number of non zero coefficients in rawKey.1
    mul_scaling: T, // Rescaling factor used in homomorphic multiplication
    q: T,
    level_max: u32,
    standard_deviation: f64, // standard deviation of the Gaussian distribution of error sampling
    raw_key: CiphertextRing<T>,
}

impl<T: BigInt> SecretKey<T> {
    pub fn new(
        dimension_exponent: u32,
        hamming_weight: usize,
        mul_scaling: T,
        q: T,
        level_max: u32,
        standard_deviation: f64,
    ) -> Self {
        let n = 1 << dimension_exponent; // 2^dimension_exponent
	let mut hwt_distribution = HWTDistribution::new(n, hamming_weight);
	let raw_key = Polynomial::new(hwt_distribution.sample::<T>()).modulo(q).to_cyclotomic(2_u32.pow(dimension_exponent));
//        let rawKey = CiphertextRing::new(secret_key_coefficients, 2_u32.pow(dimension_exponent));

        SecretKey {
            dimension_exponent,
            hamming_weight,
            mul_scaling,
            q,
            level_max,
            standard_deviation,
            raw_key,
        }
    }
}
