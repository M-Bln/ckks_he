use crate::algebra::big_int::BigInt;
use crate::algebra::polynomial::Polynomial;
use crate::ciphertext::{Ciphertext, Message, RawCiphertext};
use crate::keys::key_generator::KeyGenerationParameters;
use crate::random_distributions::HWTDistribution;

#[derive(Clone, Debug)]
pub struct SecretKey<T: BigInt> {
    pub parameters: KeyGenerationParameters<T>,
    pub key_s: Message<T>,
}

impl<T: BigInt> SecretKey<T> {
    pub fn new(parameters: KeyGenerationParameters<T>) -> Self {
        let n = 1 << parameters.dimension_exponent; // 2^dimension_exponent

        // the initial modulus in coefficients ring is initialy mul_scaling * q_0 * q^{level_max}
        let modulus = parameters.q_0 * parameters.q.fast_exp(parameters.level_max);

        // Sample the coefficients of the secret key uniformly on vector with coefficients -1, 0 and 1,
        // and exactly hamming_weight non zero coefficients.
        let mut hwt_distribution = HWTDistribution::new(n, parameters.hamming_weight);
        let key_coefficients = hwt_distribution.sample::<T>();
        let key_polynomial = Polynomial::new(key_coefficients);
        let key_s = key_polynomial
            .modulo(modulus)
            .to_cyclotomic(parameters.dimension_exponent);

        SecretKey { parameters, key_s }
    }

    pub fn decrypt(&self, cipher: &Ciphertext<T>) -> Message<T> {
        self.decrypt_raw(&cipher.raw)
    }

    pub fn decrypt_raw(&self, raw: &RawCiphertext<T>) -> Message<T> {
        raw.clone().0 + &(raw.clone().1 * &self.key_s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bnum::types::I256;

    #[test]
    fn test_secret_key_constructor() {
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
    }
}
