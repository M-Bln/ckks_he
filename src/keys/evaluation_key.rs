use crate::algebra::arithmetic::{Rescale, RingMod};
use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex, C64};
use crate::ciphertext::{Ciphertext, Message, RawCiphertext};
use crate::keys::key_generator::KeyGenerationParameters;
use crate::keys::public_key::ComputationNoise;

#[derive(Clone, Debug)]
/// Key necessary to perform homomorphic multiplication of ciphertexts
pub struct EvaluationKey<T: BigInt> {
    pub parameters: KeyGenerationParameters<T>,
    pub noise: ComputationNoise,
    pub raw_key: RawCiphertext<T>,
}

impl<T: BigInt> EvaluationKey<T> {
    pub fn new(
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

    /// Homomorphically multiply two ciphertext, then rescale and decreases the level by one
    pub fn mul(
        &self,
        ct1: &Ciphertext<T>,
        ct2: &Ciphertext<T>,
    ) -> Result<Ciphertext<T>, OperationError> {
        let mut result = self.pure_mul(ct1, ct2);
        self.rescale(&mut result, 1)?;
        Ok(result)
    }

    /// Homomorphically multiply two ciphertexts without rescaling the result
    pub fn pure_mul(&self, ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Ciphertext<T> {
        let raw = self.raw_mul(&ct1.raw, &ct2.raw);
        let level = std::cmp::min(ct1.level, ct2.level);
        let upper_bound_message = ct1.upper_bound_message * ct2.upper_bound_message;

        let modulus_f = raw.0.polynomial.ref_coefficients()[0].modulus.to_float();
        let mul_scaling_f = self.parameters.mul_scaling.to_float();
        let mul_error = ct1.upper_bound_message * ct2.upper_bound_error
            + ct1.upper_bound_error * ct2.upper_bound_message
            + ct1.upper_bound_error * ct2.upper_bound_error;
        let key_switch_error = self.noise.key_switch_noise * modulus_f / mul_scaling_f;
        let new_error = key_switch_error + mul_error + self.noise.rescaling_noise;
        Ciphertext::<T>::new(raw, level, upper_bound_message, new_error)
    }

    pub fn raw_mul(&self, rct1: &RawCiphertext<T>, rct2: &RawCiphertext<T>) -> RawCiphertext<T> {
        let (d0, d1, d2) = (
            rct1.0.clone() * &rct2.0,
            rct1.1.clone() * &rct2.0 + &(rct1.0.clone() * &rct2.1),
            rct1.1.clone() * &rct2.1,
        );
        let mut summand = self.raw_key.scalar_mul_keep_modulus(&d2);
        println!(
            "summand modulus: {:?}",
            summand.0.polynomial.ref_coefficients()[0].modulus
        );
        summand.rescale(self.parameters.mul_scaling);
        println!(
            "rescaled summand modulus: {:?}",
            summand.0.polynomial.ref_coefficients()[0].modulus
        );
        RawCiphertext(d0 + &summand.0, d1 + &summand.1)
    }

    /// For swk an encryption of a private key s' under the private key s, key_switch transforms
    /// a message m encrypted with s into the same message encrypted with s'.
    pub fn key_switch(&self, ct: &Ciphertext<T>, swk: &RawCiphertext<T>) -> Ciphertext<T> {
        let raw = self.key_switch_raw(&ct.raw, swk);
        let q_f = self.parameters.q.to_float();
        let mul_scaling_f = self.parameters.mul_scaling.to_float();
        let new_error =
            self.noise.key_switch_noise * q_f / mul_scaling_f + self.noise.rescaling_noise;
        let upper_bound_error = ct.upper_bound_error + new_error;
        Ciphertext::new(raw, ct.level, ct.upper_bound_message, upper_bound_error)
    }

    pub fn key_switch_raw(
        &self,
        rct: &RawCiphertext<T>,
        swk: &RawCiphertext<T>,
    ) -> RawCiphertext<T> {
        let mut to_rescale = swk.scalar_mul_keep_modulus(&rct.1);
        // let mut to_rescale = RawCiphertext(
        //     rct.1.clone()*&swk.0,
        //     rct.1.clone()*&swk.1,
        // );
        to_rescale.rescale(self.parameters.mul_scaling);
        RawCiphertext(rct.0.clone() + &to_rescale.0, rct.1.clone() + &to_rescale.1)
    }

    /// Devide coefficients of the ciphertext by factor. Beware, it does manage the level, it had to be modify beside
    pub fn rescale_factor(&self, ct: &mut Ciphertext<T>, factor: T) {
        let modulus = ct.raw.0.polynomial.ref_coefficients()[0].modulus;
        assert!(factor % modulus == T::from(0));
        ct.raw.rescale(factor);
        ct.upper_bound_message = ct.upper_bound_message / factor.to_float();
        ct.upper_bound_error =
            ct.upper_bound_error / factor.to_float() + self.noise.rescaling_noise;
    }

    pub fn raise_to_powers_of_two(
        &self,
        ct: &Ciphertext<T>,
        n: usize,
    ) -> Result<Vec<Ciphertext<T>>, OperationError> {
        if n >= (ct.level as usize) {
            return Err(OperationError::LevelTooLow);
        }
        let mut result: Vec<Ciphertext<T>> = vec![ct.clone()];
        for i in 0..n {
            result.push(self.mul(&result[i], &result[i])?);
        }
        Ok(result)
    }
}

#[derive(Debug)]
pub enum OperationError {
    LevelTooLow,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::big_int::{BigInt, FromFloat};
    use crate::algebra::polynomial::Polynomial;
    use crate::keys::key_generator::generate_keys;
    use bnum::types::I256;

    #[test]
    fn test_encrypt_pure_mul_decrypt() {
        let dimension_exponent = 2;
        let q = I256::from(1 << 20);
        let level_max = 4;

        // Generate keys
        let (mut public_key, evaluation_key, secret_key) =
            generate_keys(dimension_exponent, q.clone(), level_max);

        // Create sample messages
        let mut message1_coefficients = vec![I256::from(1 << 30); 2_usize.pow(dimension_exponent)];
        // let mut message1_coefficients = vec![I256::from(19*100*10000)];
        // message1_coefficients.append(&mut vec![I256::from(0); 2_usize.pow(dimension_exponent)-1]);

        let mut message2_coefficients = vec![I256::from(1 << 30); 2_usize.pow(dimension_exponent)];
        // let mut message2_coefficients = vec![I256::from(19*100*10000)];
        // message2_coefficients.append(&mut vec![I256::from(0); 2_usize.pow(dimension_exponent)-1]);

        let mut message1 = Polynomial::new(message1_coefficients)
            .modulo(q.clone().fast_exp(level_max))
            .to_cyclotomic(dimension_exponent);

        let mut message2 = Polynomial::new(message2_coefficients)
            .modulo(q.clone().fast_exp(level_max))
            .to_cyclotomic(dimension_exponent);

        // Encrypt and rescale the messages
        let upper_bound_message = (1 << 30) as f64;
        let mut ciphertext1 = public_key.encrypt(&message1, upper_bound_message);
        evaluation_key.rescale(&mut ciphertext1, 1).unwrap();
        let mut ciphertext2 = public_key.encrypt(&message2, upper_bound_message);
        evaluation_key.rescale(&mut ciphertext2, 1).unwrap();

        // Multiply the ciphertexts
        let mut pure_mul_ciphertext = evaluation_key.pure_mul(&ciphertext1, &ciphertext2);
        //	evaluation_key.rescale(&mut pure_mul_ciphertext,1).unwrap();

        // Decrypt the resulting ciphertext
        let decrypted_message = secret_key.decrypt(&pure_mul_ciphertext);

        // Create the expected result
        message1.rescale(q);
        message2.rescale(q);
        // println!("message1: {:?}", message1);
        // println!("message2: {:?}", message2);
        let mut expected_message = message1 * &message2;
        // println!("expected message: {:?}", expected_message);

        // println!(
        //     "expected error: {:?}",
        //     pure_mul_ciphertext.upper_bound_error
        // );

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
            assert!(
                diff.value < I256::from_float(pure_mul_ciphertext.upper_bound_error)
                    && diff.value > I256::from_float(-pure_mul_ciphertext.upper_bound_error),
                "Difference too large!"
            );
        }
    }

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

        // println!("expected error: {:?}", added_ciphertext.upper_bound_error);
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
        let dimension_exponent = 4;
        let q = I256::from(1 << 16);
        let level_max = 3;

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

        // println!("expected error: {:?}", ciphertext.upper_bound_error);
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
