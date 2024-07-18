use crate::algebra::big_int::BigInt;
use crate::algebra::polynomial::{Polynomial, ScalarMul};
use crate::ciphertext::RawCiphertext;
use crate::keys::client_key::ClientKey;
use crate::keys::evaluation_key::EvaluationKey;
use crate::keys::public_key::{ComputationNoise, PublicKey};
use crate::keys::secret_key::SecretKey;
use crate::keys::server_key::ServerKey;
use crate::random_distributions::{sample_n, DiscreteGaussian};
use bnum::types::I1024;

use bnum::types::I512;
use rand::distributions::uniform::UniformSampler;

#[derive(PartialEq, Clone, Copy, Debug)]
pub struct KeyGenerationParameters<T: BigInt> {
    // The dimension of the cyclotomic space of messages is 2^dimension_exponent
    pub dimension_exponent: u32,

    // The number of non zero coefficients of the secret key
    pub hamming_weight: usize,

    // The rescale factor used in multiplication (and key switch),
    // it is denoted by P in the original article.
    // The modulus of key switching key is q^level_max * q_0 * mul_scaling
    pub mul_scaling: T,

    // The modulus when all the level are consumed
    pub q_0: T,

    // The base modulus
    pub q: T,

    // The initial modulus of ciphertext is q^level_max * q_0
    pub level_max: u32,

    // Standard diviation of the Gaussian distribution used in learning with error
    pub standard_deviation: f64,
}

/// Generate a pair (client_key, server_key) with no security but fast computation
///
/// # Examples
///
/// ```
/// use ckks_he::keys::key_generator::generate_pair_keys_toy;
///
/// let (client_key, server_key) = generate_pair_keys_toy();
/// assert_eq!(client_key.secret_key.parameters.dimension_exponent, 4);
/// ```
pub fn generate_pair_keys_toy() -> (ClientKey<I1024>, ServerKey<I1024>) {
    generate_pair_keys_default(4, 5)
}

pub fn generate_pair_keys_default<T: BigInt>(
    dimension_exponent: u32,
    level_max: u32,
) -> (ClientKey<T>, ServerKey<T>) {
    generate_pair_keys_all_parameters(generate_with_default_q(dimension_exponent, level_max))
}

pub fn generate_pair_keys<T: BigInt>(
    dimension_exponent: u32,
    q: T,
    level_max: u32,
) -> (ClientKey<T>, ServerKey<T>) {
    generate_pair_keys_all_parameters(generate_most_parameters(dimension_exponent, q, level_max))
}

pub fn generate_pair_keys_all_parameters<T: BigInt>(
    params: KeyGenerationParameters<T>,
) -> (ClientKey<T>, ServerKey<T>) {
    let (mut public_key, mut evaluation_key, mut secret_key) = generate_keys_all_parameters(params);
    let mut client_key = ClientKey::new(params, secret_key, public_key.clone());
    let mut server_key = ServerKey::new(params, public_key, evaluation_key);
    (client_key, server_key)
}

pub fn generate_with_default_q<T: BigInt>(
    dimension_exponent: u32,
    level_max: u32,
) -> KeyGenerationParameters<T> {
    generate_most_parameters(dimension_exponent, T::from(1 << 30), level_max)
}

/// generate all parameters from a subset
pub fn generate_most_parameters<T: BigInt>(
    dimension_exponent: u32,
    q: T,
    level_max: u32,
) -> KeyGenerationParameters<T> {
    // Should be smaller than the dimension
    let hamming_weight = if dimension_exponent <= 12 {
        2_usize.pow(dimension_exponent / 2)
    } else {
        64 // For dimension large enough, we take the value 64 suggested in the original article
    };

    // Should be of the same order as the maximum modulus of ciphertext
    let mul_scaling = q.fast_exp(level_max) / T::from(11);
    let q_0 = T::from(1);
    let standard_deviation = 3.2;
    KeyGenerationParameters {
        dimension_exponent,
        hamming_weight,
        mul_scaling,
        q_0,
        q,
        level_max,
        standard_deviation,
    }
}

pub fn generate_keys<T: BigInt>(
    dimension_exponent: u32,
    q: T,
    level_max: u32,
) -> (PublicKey<T>, EvaluationKey<T>, SecretKey<T>) {
    generate_keys_all_parameters(generate_most_parameters(dimension_exponent, q, level_max))
}

pub fn generate_keys_all_parameters<T: BigInt>(
    params: KeyGenerationParameters<T>,
) -> (PublicKey<T>, EvaluationKey<T>, SecretKey<T>) {
    let noise = ComputationNoise::new(params);
    let secret_key = generate_secret_key(params);

    let public_key = generate_public_key(params, noise, &secret_key);

    let evaluation_key = generate_evaluation_key(params, noise, &secret_key);

    (public_key, evaluation_key, secret_key)
}

fn generate_secret_key<T: BigInt>(parameters: KeyGenerationParameters<T>) -> SecretKey<T> {
    SecretKey::<T>::new(parameters)
}

fn generate_public_key<T: BigInt>(
    params: KeyGenerationParameters<T>,
    noise: ComputationNoise,
    secret_key: &SecretKey<T>,
) -> PublicKey<T> {
    // Dimension of the cyclotomic ring of messages
    let dimension = 2_usize.pow(params.dimension_exponent);
    let mut rng = rand::thread_rng();

    // Initial modulus of ciphertexts
    let modulus = params.q_0.clone() * params.q.fast_exp(params.level_max);

    // The public key is sampled uniformly on the message space
    let sampler = T::sampler(-modulus.clone() / T::from(2), modulus.clone() / T::from(2));
    let public_key_coefficients = sample_n(sampler, dimension, &mut rng);
    let public_key_a = Polynomial::<T>::new(public_key_coefficients)
        .modulo(modulus.clone())
        .to_cyclotomic(params.dimension_exponent);

    let mut gaussian_sampler = DiscreteGaussian::new(0.0, params.standard_deviation);
    let error_coefficients = gaussian_sampler.sample_n(dimension);
    let error = Polynomial::<T>::new(error_coefficients)
        .modulo(modulus.clone())
        .to_cyclotomic(params.dimension_exponent);
    // The public key is an encryption of zero
    let raw_public_key = RawCiphertext::<T>(
        error - &(public_key_a.clone() * &secret_key.key_s),
        public_key_a,
    );

    PublicKey::<T>::new(params, noise, raw_public_key)
}

fn generate_evaluation_key<T: BigInt>(
    params: KeyGenerationParameters<T>,
    noise: ComputationNoise,
    secret_key: &SecretKey<T>,
) -> EvaluationKey<T> {
    // Dimension of the cyclotomic ring of messages
    let dimension = 2_usize.pow(params.dimension_exponent);
    let mut rng = rand::thread_rng();

    // Initial modulus of ciphertexts
    let modulus = params.q_0.clone() * params.q.fast_exp(params.level_max);

    // Modulus of the evaluation key (roughly a key switching key for secret_key^2)
    let modulus_eval = modulus.clone() * &params.mul_scaling;

    // Coefficients of the evaluation key are sampled uniformly
    let sampler_eval = T::sampler(T::from(0), modulus_eval.clone());
    let eval_key_coefficients = sample_n(sampler_eval, dimension, &mut rng);
    let eval_key_a = Polynomial::<T>::new(eval_key_coefficients)
        .modulo(modulus_eval.clone())
        .to_cyclotomic(params.dimension_exponent);

    // Raise the modulus of the secret key from modulus to modulus_eval
    let secret_key_modulo_eval = secret_key.key_s.to_integer().modulo(modulus_eval);
    let secret_key_squared = secret_key_modulo_eval.clone() * &secret_key_modulo_eval;

    let mut gaussian_sampler = DiscreteGaussian::new(0.0, params.standard_deviation);
    let eval_error_coefficients = gaussian_sampler.sample_n(dimension);
    let eval_error = Polynomial::<T>::new(eval_error_coefficients)
        .modulo(modulus_eval.clone())
        .to_cyclotomic(params.dimension_exponent);
    let raw_eval_key = RawCiphertext::<T>(
        (eval_error - &(secret_key_modulo_eval * &eval_key_a))
            + &params.mul_scaling.scalar_mul(secret_key_squared),
        eval_key_a.clone(),
    );

    EvaluationKey::<T>::new(params, noise, raw_eval_key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::big_int::{BigInt, FromFloat};
    use crate::algebra::conversion_rounding::f64_to_i256;
    use crate::ciphertext::Ciphertext;
    use bnum::types::I256;

    #[test]
    fn test_generate_keys_all_parameters() {
        let params = KeyGenerationParameters {
            dimension_exponent: 3,
            hamming_weight: 3,
            mul_scaling: I256::from(3),
            q_0: I256::from(7),
            q: I256::from(8),
            level_max: 2,
            standard_deviation: 3.2,
        };

        let (public_key, evaluation_key, secret_key) = generate_keys_all_parameters(params);

        // println!("PublicKey: {:?}", public_key);
        // println!("EvaluationKey: {:?}", evaluation_key);
        // println!("SecretKey: {:?}", secret_key);

        assert_eq!(secret_key.parameters, params);
        assert_eq!(public_key.parameters, params);
        assert_eq!(evaluation_key.parameters, params);
    }

    #[test]
    fn test_encrypt_decrypt() {
        // Define parameters for key generation
        let dimension_exponent = 10;
        let q = I256::from(16);
        let level_max = 5;
        let modulus = q.clone() * q.clone().fast_exp(level_max);

        // Generate keys using the provided helper function
        let (mut public_key, _evaluation_key, secret_key) =
            generate_keys(dimension_exponent, q.clone(), level_max);

        // Create a sample message
        let message_coefficients = vec![I256::from(1209); 2_usize.pow(dimension_exponent)];
        let message = Polynomial::new(message_coefficients)
            .modulo(modulus)
            .to_cyclotomic(dimension_exponent);

        // Encrypt the message
        let upper_bound_message = 1209.0; // Example value, adjust as needed
        let ciphertext = public_key.encrypt(&message, upper_bound_message);

        // Decrypt the ciphertext
        let decrypted_message = secret_key.decrypt(&ciphertext);

        // Verify that the decrypted message is close to the original message
        for (original, decrypted) in message
            .polynomial
            .coefficients()
            .iter()
            .zip(decrypted_message.polynomial.coefficients().iter())
        {
            let diff = (original.value - decrypted.value).remainder(&modulus);
            println!(
                "Original: {:?}, Decrypted: {:?}, Difference: {}",
                original, decrypted, diff
            );
            println!(
                "theoretical encryption error: {:?} \n {:?} ",
                f64_to_i256(public_key.noise.clean_noise),
                public_key.noise.clean_noise
            );
            assert!(
                diff < f64_to_i256(public_key.noise.clean_noise),
                "Difference too large!"
            );
            assert!(
                diff > f64_to_i256(-public_key.noise.clean_noise),
                "Difference too large!"
            );
        }
    }

    #[test]
    /// Check that the raw eval key is indeed an encryption of the rescaled squared secret key
    fn test_decrypt_raw_eval_key() {
        // Define parameters for key generation
        let dimension_exponent = 10;
        let q = I256::from(19);
        let level_max = 5;

        let (_public_key, evaluation_key, secret_key) =
            generate_keys(dimension_exponent, q.clone(), level_max);

        // Decrypt the raw evaluation key
        let decrypted_eval_key = secret_key.decrypt(&Ciphertext::new(
            evaluation_key.raw_key.clone(),
            level_max,
            (1 << dimension_exponent) as f64,
            evaluation_key.noise.clean_noise,
        ));

        // Compute the expected value: mul_scaling * key_s * key_s
        let key_s_squared = secret_key.key_s.clone() * &secret_key.key_s;
        let expected_value = secret_key.parameters.mul_scaling.scalar_mul(key_s_squared);

        // Verify that the decrypted evaluation key is close to the expected value
        for (expected, decrypted) in expected_value
            .polynomial
            .coefficients()
            .iter()
            .zip(decrypted_eval_key.polynomial.coefficients().iter())
        {
            let diff = *expected - decrypted;
            println!(
                "Expected: {:?}, Decrypted: {:?}, Difference: {:?}",
                expected, decrypted, diff
            );
            assert!(
                diff.value < I256::from_float(evaluation_key.noise.clean_noise),
                "Difference too large!"
            );
        }
    }
}
