use std::ops::Mul;

use crate::algebra::arithmetic::{Rescale, RingMod};
use crate::algebra::big_int::{BigInt, Zero};
use crate::algebra::complex::{Complex, C64};
use crate::algebra::polynomial::{degree_from_coefs, Polynomial, ScalarMul};
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
        // println!(
        //     "summand modulus: {:?}",
        //     summand.0.polynomial.ref_coefficients()[0].modulus
        // );
        summand.rescale(self.parameters.mul_scaling);
        // println!(
        //     "rescaled summand modulus: {:?}",
        //     summand.0.polynomial.ref_coefficients()[0].modulus
        // );
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

    /// For ct an encryption of a message m, computes an encryption
    /// of q(m/q), q(m/q)^2, q(m/q)^{2^2}, ..., q(m/q)^{2^{n-1}}
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

    // pub fn apply_polynomial<'a, U>(&self, polynomial: &Polynomial<U>, ct: &Ciphertext<T>) -> Ciphertext<T>
    // where
    // 	U: Mul<&'a Message<T>, Output = Message<T>> + Clone,
    // 	T: 'a,
    // {
    // 	let degree = polynomial.degree();
    // 	if degree = -1 {
    // 	    return self.trivial_encryption_scalar(T::from(0));
    // 	}
    // 	if degree = 0 {
    // 	    return self.trivial_encryption_scalar()
    // 	}
    // 	let largest_power = largest_power_of_two_less_than();
    // 	// Your implementation here
    // }
    pub fn apply_polynomial_coefficients(
        &self,
        polynomial_coefs: &[T],
        ct: &Ciphertext<T>,
    ) -> Result<Ciphertext<T>, OperationError> {
        let degree = degree_from_coefs(polynomial_coefs);
        if degree == -1 {
            return Ok(self.trivial_encryption_scalar(T::from(0)));
        }
        if degree == 0 {
            return Ok(self.trivial_encryption_scalar(polynomial_coefs[0]));
        }
        if degree == 1 {
            let constant_term = self.trivial_encryption_scalar(polynomial_coefs[0]);
            let degree_1_term = polynomial_coefs[1].scalar_mul(ct);
            return Ok(self.add(&constant_term, &degree_1_term));
        }
        let largest_power = largest_power_of_two_less_than(degree as u32);
        let powers = self.raise_to_powers_of_two(ct, largest_power as usize + 1)?;
        self.recursive_apply_polynomial(&polynomial_coefs[..degree as usize +1], ct, &powers)
    }

    pub fn recursive_apply_polynomial(
	&self,
	polynomial_coefs: &[T],
	ct: &Ciphertext<T>,
	powers: &[Ciphertext<T>],
    ) -> Result<Ciphertext<T>, OperationError> {
	let degree = polynomial_coefs.len() - 1;
	
	if degree <= 1 {
            return self.apply_polynomial_coefficients(polynomial_coefs, ct);
	}

	let largest_power = largest_power_of_two_less_than(degree as u32);
	let cut = largest_power as usize;

	let first_part = &polynomial_coefs[..=cut];
	let second_part = &polynomial_coefs[cut+1..];

	let first_eval = self.recursive_apply_polynomial(first_part, ct, powers)?;
	let second_eval = self.recursive_apply_polynomial(second_part, ct, powers)?;
	let second_rescaled = self.mul(&second_eval, &powers[largest_power as usize])?;

	Ok(self.add(&first_eval, &second_rescaled))
    }

    
    // pub fn recursive_apply_polynomial(
    //     &self,
    //     polynomial_coefs: &[T],
    //     ct: &Ciphertext<T>,
    //     powers: &[Ciphertext<T>],
    // ) -> Result<Ciphertext<T>, OperationError> {
    //     if polynomial_coefs.len() <= 2 {
    //         return self.apply_polynomial_coefficients(polynomial_coefs, ct);
    //     }
    //     let cut = 1 << powers.len();
    //     let first_part = &polynomial_coefs[..cut];
    //     let second_part = &polynomial_coefs[cut..];
    //     let first_eval =
    //         self.recursive_apply_polynomial(first_part, ct, &powers[..powers.len() - 1])?;
    //     let second_eval =
    //         self.recursive_apply_polynomial(second_part, ct, &powers[..powers.len() - 1])?;
    //     let second_rescaled = self.mul(&second_eval, &powers[powers.len() - 1])?;
    //     Ok(self.add(&first_eval, &second_rescaled))
    // }

    pub fn apply_polynomial(
        &self,
        polynomial: &Polynomial<T>,
        ct: &Ciphertext<T>,
    ) -> Result<Ciphertext<T>, OperationError> {
        // let mut truncated = polynomial.clone();
        // truncated.degree_truncate();
        self.apply_polynomial_coefficients(polynomial.ref_coefficients(), ct)
        // Your implementation here
    }

    pub fn message_from_scalar(&self, scalar: T) -> Message<T> {
        let dimension = 1 << self.parameters.dimension_exponent;
        let mut coefficients = vec![scalar];
        coefficients.append(&mut vec![scalar.zero(); dimension - 1]);
        let polynomial = Polynomial::new(coefficients);

        let modulus = self.parameters.q.fast_exp(self.parameters.level_max) * &self.parameters.q_0;
        polynomial
            .modulo(modulus)
            .to_cyclotomic(self.parameters.dimension_exponent)
    }

    pub fn raw_trivial_encryption(&self, message: &Message<T>) -> RawCiphertext<T> {
        RawCiphertext(
            message.clone(),
            self.message_from_scalar(T::from(0))
                .modulo(message.modulus()),
        )
    }

    pub fn trivial_encryption_scalar(&self, scalar: T) -> Ciphertext<T> {
        let raw =
            self.raw_trivial_encryption(&self.message_from_scalar(scalar * self.parameters.q));
        Ciphertext::<T> {
            raw,
            level: self.parameters.level_max,
            upper_bound_message: scalar.to_float().abs(),
            upper_bound_error: 0.0,
        }
    }

    pub fn trivial_encryption_scalar_no_rescale(&self, scalar: T) -> Ciphertext<T> {
        let raw = self.raw_trivial_encryption(&self.message_from_scalar(scalar));
        Ciphertext::<T> {
            raw,
            level: self.parameters.level_max,
            upper_bound_message: scalar.to_float().abs(),
            upper_bound_error: 0.0,
        }
    }
}

pub fn largest_power_of_two_less_than(d: u32) -> u32 {
    assert!(d > 1, "call largest_power_of_two_less_than with d < 2");
    let mut power = 1;
    let mut result = 0;
    while power < d / 2 {
        power *= 2;
        result = result + 1;
    }
    result
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
