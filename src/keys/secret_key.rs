use crate::algebra::big_int::BigInt;
use crate::algebra::polynomial::Polynomial;
use crate::ciphertext::{Ciphertext, Message, RawCiphertext};
use crate::keys::key_generator::KeyGenerationParameters;
use crate::random_distributions::HWTDistribution;

#[derive(Clone, Debug)]
pub struct SecretKey<T: BigInt> {
    // pub dimension_exponent: u32, // Ciphertexts in cyclotomic ring $R[X]/(1+X^N)$ with N = 2^dimension_exponent
    // pub hamming_weight: usize,   // Number of non zero coefficients in rawKey.1
    // pub mul_scaling: T,          // Rescaling factor used in homomorphic multiplication
    // pub q_0: T,                  // minimal modulus
    // pub q: T, // modulus per level, the total modulus in coefficients ring is initialy mul_scaling * q_0 * q^{level_max}
    // pub level_max: u32,
    // pub variance: f64, // standard deviation of the Gaussian distribution of error sampling
    pub parameters: KeyGenerationParameters<T>,
    pub key_s: Message<T>,
}

impl<T: BigInt> SecretKey<T> {
    pub fn new(
        // dimension_exponent: u32,
        // hamming_weight: usize,
        // mul_scaling: T,
        // q_0: T,
        // q: T,
        // level_max: u32,
        parameters: KeyGenerationParameters<T>,
        //        variance: f64,
    ) -> Self {
        let n = 1 << parameters.dimension_exponent; // 2^dimension_exponent
        let modulus =
            parameters.q_0 * parameters.q.fast_exp(parameters.level_max); // the total modulus in coefficients ring is initialy mul_scaling * q_0 * q^{level_max}
        let mut hwt_distribution = HWTDistribution::new(n, parameters.hamming_weight);

        let key_coefficients = hwt_distribution.sample::<T>();
        let key_polynomial = Polynomial::new(key_coefficients);
        let key_s = key_polynomial
            .modulo(modulus)
            .to_cyclotomic(parameters.dimension_exponent);

        SecretKey {
            // dimension_exponent,
            // hamming_weight,
            // mul_scaling,
            // q_0,
            // q,
            // level_max,
            // variance,
            parameters,
            key_s,
        }
    }

    pub fn decrypt(&self, cipher: &Ciphertext<T>) -> Message<T> {
        self.decrypt_raw(&cipher.raw)
        //cipher.raw.clone().0 + &(cipher.raw.clone().1 * &self.key_s)
    }

    pub fn decrypt_raw(&self, raw: &RawCiphertext<T>) -> Message<T> {
        raw.clone().0 + &(raw.clone().1 * &self.key_s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::big_int::BigInt;
    use crate::ciphertext::RawCiphertext;
    use crate::keys::public_key::PublicKey;
    use bnum::types::I256;

    #[test]
    fn test_secret_key_constructor() {
        // let dimension_exponent = 3;
        // let hamming_weight = 4;
        // let mul_scaling = I256::new(3);
        // let q_0 = I256::new(17);
        // let q = I256::new(19);
        // let level_max = 3;
        // let variance = 3.2;

        let params = KeyGenerationParameters {
            dimension_exponent: 3,
            hamming_weight: 3,
            mul_scaling: I256::from(3),
            q_0: I256::from(7),
            q: I256::from(8),
            level_max: 2,
            standard_deviation: 3.2,
        };

        let secret_key = SecretKey::new(params);

        println!("SecretKey: {:?}", secret_key);

        assert_eq!(secret_key.parameters, params);
        // assert_eq!(secret_key.dimension_exponent, dimension_exponent);
        // assert_eq!(secret_key.hamming_weight, hamming_weight);
        // assert_eq!(secret_key.mul_scaling, mul_scaling);
        // assert_eq!(secret_key.q_0, q_0);
        // assert_eq!(secret_key.q, q);
        // assert_eq!(secret_key.level_max, level_max);
        // assert_eq!(secret_key.variance, variance);
    }
}
