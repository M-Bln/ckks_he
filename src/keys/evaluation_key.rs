use crate::algebra::arithmetic::{Rescale, RingMod};
use crate::algebra::big_int::BigInt;
use crate::ciphertext::{Ciphertext, Message, RawCiphertext};
use crate::keys::key_generator::KeyGenerationParameters;
use crate::keys::public_key::ComputationNoise;

#[derive(Clone, Debug)]
pub struct EvaluationKey<T: BigInt> {
    // pub dimension_exponent: u32, // Ciphertexts in cyclotomic ring $R[X]/(1+X^N)$ with N = 2^dimension_exponent
    // pub mul_scaling: T,          // Rescaling factor used in homomorphic multiplication
    // pub q_0: T,                  // minimal modulus
    // pub q: T, // modulus per level, the total modulus in coefficients ring is initialy mul_scaling * q_0 * q^{level_max}
    // pub level_max: u32,
    // pub variance: f64, // standard deviation of the Gaussian distribution of error sampling
    pub parameters: KeyGenerationParameters<T>,
    pub noise: ComputationNoise,
    pub raw_key: RawCiphertext<T>,
}

impl<T: BigInt> EvaluationKey<T> {
    pub fn new(
        // dimension_exponent: u32,
        // mul_scaling: T,
        // q_0: T,
        // q: T,
        // level_max: u32,
        // variance: f64,
        parameters: KeyGenerationParameters<T>,
        noise: ComputationNoise,
        raw_key: RawCiphertext<T>,
    ) -> Self {
        EvaluationKey {
            parameters,
            noise,
            raw_key,
        }
    }

    pub fn add(&self, ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Ciphertext<T> {
        let raw = ct1.raw.clone() + &ct2.raw;
        let level = std::cmp::min(ct1.level, ct2.level);
        let upper_bound_message = ct1.upper_bound_message + ct2.upper_bound_message;
        let upper_bound_error = ct1.upper_bound_error + ct2.upper_bound_error;

        let modulus = self.parameters.q.fast_exp(level) * &self.parameters.q_0;
        assert_eq!(raw.0.polynomial.ref_coefficients()[0].modulus, modulus);
        Ciphertext::<T>::new(raw, level, upper_bound_message, upper_bound_error)
    }

    // /// Devide coefficients of the ciphertext by q^level_decrement, modify level accordingly
    // pub fn rescale(&self, ct: Ciphertext<T>, level_decrement: u32) -> Ciphertext<T> {
    // 	assert!(ct.level > level_decrement);
    // 	let scaling_factor = self.parameters.q.fast_exp(level_decrement);
    // 	let raw = ct.raw.rescale(scaling_factor);
    // 	ct.level = ct.level - level_decrement;
    // 	ct.upper_bound_message = ct.upper_bound_message / scaling_factor.to_float();
    // 	ct.upper_bound_error = ct.upper_bound_error / scaling_factor.to_float() + self.noise.rescaling_noise;
    // }

    /// Devide coefficients of the ciphertext by q^level_decrement, modify level accordingly
    pub fn rescale(
        &self,
        ct: &mut Ciphertext<T>,
        level_decrement: u32,
    ) -> Result<(), OperationError> {
        if ct.level <= level_decrement {
            return Err(OperationError::LevelTooLow);
        }
        let scaling_factor = self.parameters.q.fast_exp(level_decrement);
        ct.raw.rescale(scaling_factor);
        ct.level -= level_decrement;
        ct.upper_bound_message /= scaling_factor.to_float();
        ct.upper_bound_error =
            ct.upper_bound_error / scaling_factor.to_float() + self.noise.rescaling_noise;
        Ok(())
    }

    // /// Devide coefficients of the ciphertext by q^level_decrement, modify level accordingly
    // pub fn rescale(&self, ct: &mut Ciphertext<T>, level_decrement: u32) {
    // 	assert!(ct.level > level_decrement);
    // 	let scaling_factor = self.parameters.q.fast_exp(level_decrement);
    // 	ct.raw.rescale(scaling_factor);
    // 	ct.level = ct.level - level_decrement;
    // 	ct.upper_bound_message = ct.upper_bound_message / scaling_factor.to_float();
    // 	ct.upper_bound_error = ct.upper_bound_error / scaling_factor.to_float() + self.noise.rescaling_noise;
    // }

    /// Devide coefficients of the ciphertext by factor. Beware, it does manage the level, it had to be modify beside
    pub fn rescale_factor(&self, ct: &mut Ciphertext<T>, factor: T) {
        let modulus = ct.raw.0.polynomial.ref_coefficients()[0].modulus;
        assert!(factor % modulus == T::from(0));
        ct.raw.rescale(factor);
        ct.upper_bound_message = ct.upper_bound_message / factor.to_float();
        ct.upper_bound_error =
            ct.upper_bound_error / factor.to_float() + self.noise.rescaling_noise;
    }
}

#[derive(Debug)]
pub enum OperationError {
    LevelTooLow,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::big_int::BigInt;
    use crate::algebra::polynomial::Polynomial;
    use crate::keys::key_generator::generate_keys;
    use bnum::types::I256;

    #[test]
    fn test_encrypt_add_decrypt() {
        // Define parameters for key generation
        let dimension_exponent = 4;
        let q = I256::from(1 << 16);
        let level_max = 3;

        // Generate keys using the provided helper function
        let (mut public_key, evaluation_key, secret_key) =
            generate_keys(dimension_exponent, q.clone(), level_max);

        // Create sample messages
        let message1_coefficients = vec![I256::from(1 << 13); 2_usize.pow(dimension_exponent)];
        let message2_coefficients = vec![I256::from(1 << 13); 2_usize.pow(dimension_exponent)];

        let message1 = Polynomial::new(message1_coefficients)
            .modulo(q.clone().fast_exp(level_max))
            .to_cyclotomic(dimension_exponent);

        let message2 = Polynomial::new(message2_coefficients)
            .modulo(q.clone().fast_exp(level_max))
            .to_cyclotomic(dimension_exponent);

        // Encrypt the messages
        let upper_bound_message = (1 << 13) as f64; // Example value, adjust as needed
        let ciphertext1 = public_key.encrypt(&message1, upper_bound_message);
        let ciphertext2 = public_key.encrypt(&message2, upper_bound_message);

        // Add the ciphertexts
        let added_ciphertext = evaluation_key.add(&ciphertext1, &ciphertext2);

        // Decrypt the resulting ciphertext
        let decrypted_message = secret_key.decrypt(&added_ciphertext);

        // Create the expected message sum
        let expected_message_coefficients =
            vec![I256::from(1 << 14); 2_usize.pow(dimension_exponent)];
        let expected_message = Polynomial::new(expected_message_coefficients)
            .modulo(q.clone().fast_exp(level_max))
            .to_cyclotomic(dimension_exponent);

        println!("expected error: {:?}", added_ciphertext.upper_bound_error);
        // Verify that the decrypted message is close to the sum of the original messages
        for (expected, decrypted) in expected_message
            .polynomial
            .coefficients()
            .iter()
            .zip(decrypted_message.polynomial.coefficients().iter())
        {
            let diff = *expected - *decrypted;
            println!(
                "Expected: {:?}, Decrypted: {:?}, Difference: {:?}",
                expected, decrypted, diff
            );
            assert!(diff.value < I256::from(100), "Difference too large!");
        }
    }

    #[test]
    fn test_encrypt_rescale_decrypt() {
        // Define parameters for key generation
        let dimension_exponent = 4;
        let q = I256::from(1 << 16);
        let level_max = 3;

        // Generate keys using the provided helper function
        let (mut public_key, evaluation_key, secret_key) =
            generate_keys(dimension_exponent, q.clone(), level_max);

        // Create sample messages
        let message_coefficients = vec![I256::from(1 << 29); 2_usize.pow(dimension_exponent)];

        let message = Polynomial::new(message_coefficients)
            .modulo(q.clone().fast_exp(level_max))
            .to_cyclotomic(dimension_exponent);

        // Encrypt the messages
        let upper_bound_message = (1 << 29) as f64; // Example value, adjust as needed
        let mut ciphertext = public_key.encrypt(&message, upper_bound_message);

        // Add the ciphertexts
        evaluation_key.rescale(&mut ciphertext, 1).unwrap();

        // Decrypt the resulting ciphertext
        let decrypted_message = secret_key.decrypt(&ciphertext);

        // Create the expected message sum
        let expected_message_coefficients =
            vec![I256::from(1 << 13); 2_usize.pow(dimension_exponent)];
        let expected_message = Polynomial::new(expected_message_coefficients)
            .modulo(q.clone().fast_exp(level_max - 1))
            .to_cyclotomic(dimension_exponent);

        println!("expected error: {:?}", ciphertext.upper_bound_error);
        // Verify that the decrypted message is close to the sum of the original messages
        for (expected, decrypted) in expected_message
            .polynomial
            .coefficients()
            .iter()
            .zip(decrypted_message.polynomial.coefficients().iter())
        {
            let diff = *expected - *decrypted;
            println!(
                "Expected: {:?}, Decrypted: {:?}, Difference: {:?}",
                expected, decrypted, diff
            );
            assert!(diff.value < I256::from(100), "Difference too large!");
        }

        assert!(matches!(
            evaluation_key.rescale(&mut ciphertext, level_max),
            Err(OperationError::LevelTooLow)
        ));
    }
}
