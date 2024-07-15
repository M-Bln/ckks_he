use crate::algebra::big_int::BigInt;
use crate::algebra::polynomial::{Polynomial, ScalarMul};
use crate::ciphertext::RawCiphertext;
use crate::keys::evaluation_key::EvaluationKey;
use crate::keys::public_key::{ComputationNoise, PublicKey};
use crate::keys::secret_key::SecretKey;
use crate::random_distributions::{sample_n, DiscreteGaussian};
use rand::distributions::uniform::UniformSampler;

#[derive(PartialEq, Clone, Copy, Debug)]
pub struct KeyGenerationParameters<T: BigInt> {
    pub dimension_exponent: u32,
    pub hamming_weight: usize,
    pub mul_scaling: T,
    pub q_0: T,
    pub q: T,
    pub level_max: u32,
    pub standard_deviation: f64,
}

pub fn generate_most_parameters<T: BigInt>(
    dimension_exponent: u32,
    q: T,
    level_max: u32,
) -> KeyGenerationParameters<T> {
    let hamming_weight = 2_usize.pow(dimension_exponent / 2);
    let mul_scaling = q.fast_exp(level_max) / T::from(3);
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
    SecretKey::<T>::new(
        parameters, // params.dimension_exponent,
                   // params.hamming_weight,
                   // params.mul_scaling.clone(),
                   // params.q_0.clone(),
                   // params.q.clone(),
                   // params.level_max,
                   // params.standard_deviation,
    )
}

fn generate_public_key<T: BigInt>(
    params: KeyGenerationParameters<T>,
    noise: ComputationNoise,
    secret_key: &SecretKey<T>,
) -> PublicKey<T> {
    let dimension = 2_usize.pow(params.dimension_exponent);
    let mut rng = rand::thread_rng();
    let modulus = params.q_0.clone() * params.q.fast_exp(params.level_max);
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
    let raw_public_key = RawCiphertext::<T>(
        error - &(public_key_a.clone() * &secret_key.key_s),
        public_key_a,
    );

    PublicKey::<T>::new(
        // params.dimension_exponent,
        // params.hamming_weight,
        // params.mul_scaling.clone(),
        // params.q_0.clone(),
        // params.q.clone(),
        // params.level_max,
        // params.variance,
        params,
        noise,
        raw_public_key,
    )
}

fn generate_evaluation_key<T: BigInt>(
    params: KeyGenerationParameters<T>,
    noise: ComputationNoise,
    secret_key: &SecretKey<T>,
) -> EvaluationKey<T> {
    let dimension = 2_usize.pow(params.dimension_exponent);
    let mut rng = rand::thread_rng();
    let modulus = params.q_0.clone() * params.q.fast_exp(params.level_max);
    let modulus_eval = modulus.clone() * &params.mul_scaling;
    let sampler_eval = T::sampler(T::from(0), modulus_eval.clone());
    let eval_key_coefficients = sample_n(sampler_eval, dimension, &mut rng);
    let eval_key_a = Polynomial::<T>::new(eval_key_coefficients)
        .modulo(modulus_eval.clone())
        .to_cyclotomic(params.dimension_exponent);
    //    println!("eval_key_a: {:?}", eval_key_a);
    //    println!("secret_key: {:?}", secret_key);
    let secret_key_modulo_eval = secret_key.key_s.to_integer().modulo(modulus_eval);
    //    println!("secret_key_modulo_eval: {:?}", secret_key_modulo_eval);
    let secret_key_squared = secret_key_modulo_eval.clone() * &secret_key_modulo_eval;
    //    println!("secret_key_squared: {:?}", secret_key_squared);
    let mut gaussian_sampler = DiscreteGaussian::new(0.0, params.standard_deviation);
    let eval_error_coefficients = gaussian_sampler.sample_n(dimension);
    let eval_error = Polynomial::<T>::new(eval_error_coefficients)
        .modulo(modulus_eval.clone())
        .to_cyclotomic(params.dimension_exponent);
    //    println!("eval_error: {:?}", eval_error);
    let raw_eval_key = RawCiphertext::<T>(
        (eval_error - &(secret_key_modulo_eval * &eval_key_a))
            + &params.mul_scaling.scalar_mul(secret_key_squared),
        //             &(params.mul_scaling.scalar_mul(secret_key.key_s.clone()) * &secret_key.key_s),
        eval_key_a.clone(),
    );
    //    println!("raw_eval_key: {:?}", raw_eval_key);

    EvaluationKey::<T>::new(
        // params.dimension_exponent,
        // params.mul_scaling.clone(),
        // params.q_0.clone(),
        // params.q.clone(),
        // params.level_max,
        // params.variance,
        params,
        noise,
        raw_eval_key,
    )
}

// pub fn generate_keys_all_parameters<T: BigInt>(
//     dimension_exponent: u32,
//     hamming_weight: usize,
//     mul_scaling: T,
//     q_0: T,
//     q: T,
//     level_max: u32,
//     variance: f64,
// ) -> (PublicKey<T>, EvaluationKey<T>, SecretKey<T>) {
//     let secret_key = SecretKey::<T>::new(
//         dimension_exponent,
//         hamming_weight,
//         mul_scaling,
//         q_0,
//         q,
//         level_max,
//         variance,
//     );
//     let dimension = 2_usize.pow(dimension_exponent);
//     let mut rng = rand::thread_rng();
//     let modulus = q_0 * q.fast_exp(level_max); // the total modulus in coefficients ring is initialy mul_scaling * q_0 * q^{level_max}
//     let sampler = T::sampler(T::from(0), modulus);
//     let public_key_coefficients = sample_n(sampler, dimension, &mut rng);
//     let public_key_a = Polynomial::<T>::new(public_key_coefficients)
//         .modulo(modulus)
//         .to_cyclotomic(dimension_exponent);

//     let mut gaussian_sampler = DiscreteGaussian::new(0.0, variance);
//     let error_coefficients = gaussian_sampler.sample_n(dimension);
//     let error = Polynomial::<T>::new(error_coefficients)
//         .modulo(modulus)
//         .to_cyclotomic(dimension_exponent);
//     let raw_public_key = RawCiphertext::<T>(
//         public_key_a.clone(),
//         error - &(public_key_a * &secret_key.key_s),
//     );

//     let public_key = PublicKey::<T>::new(
//         dimension_exponent,
//         mul_scaling,
//         q_0,
//         q,
//         level_max,
//         variance,
//         raw_public_key,
//     );

//     let modulus_eval = modulus * &mul_scaling;
//     let sampler_eval = T::sampler(T::from(0), modulus_eval);
//     let eval_key_coefficients = sample_n(sampler_eval, dimension, &mut rng);
//     let eval_key_a = Polynomial::<T>::new(eval_key_coefficients)
//         .modulo(modulus_eval)
//         .to_cyclotomic(dimension_exponent);

//     let eval_error_coefficients = gaussian_sampler.sample_n(dimension);
//     let eval_error = Polynomial::<T>::new(eval_error_coefficients)
//         .modulo(modulus_eval)
//         .to_cyclotomic(dimension_exponent);
//     let raw_eval_key = RawCiphertext::<T>(
//         eval_key_a.clone(),
//         (eval_error - &eval_key_a)
//             + &(mul_scaling.scalar_mul(secret_key.key_s.clone()) * &secret_key.key_s),
//     );

//     let evaluation_key = EvaluationKey::<T>::new(
//         dimension_exponent,
//         mul_scaling,
//         q_0,
//         q,
//         level_max,
//         variance,
//         raw_eval_key,
//     );
//     (public_key, evaluation_key, secret_key)
// }

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

        println!("PublicKey: {:?}", public_key);
        println!("EvaluationKey: {:?}", evaluation_key);
        println!("SecretKey: {:?}", secret_key);

        // Verify that the keys are generated correctly
        // assert_eq!(secret_key.dimension_exponent, params.dimension_exponent);
        // assert_eq!(secret_key.hamming_weight, params.hamming_weight);
        // assert_eq!(secret_key.mul_scaling, params.mul_scaling);
        // assert_eq!(secret_key.q_0, params.q_0);
        // assert_eq!(secret_key.q, params.q);
        // assert_eq!(secret_key.level_max, params.level_max);
        // assert_eq!(secret_key.variance, params.variance);
        assert_eq!(secret_key.parameters, params);

        // assert_eq!(public_key.dimension_exponent, params.dimension_exponent);
        // assert_eq!(public_key.mul_scaling, params.mul_scaling);
        // assert_eq!(public_key.q_0, params.q_0);
        // assert_eq!(public_key.q, params.q);
        // assert_eq!(public_key.level_max, params.level_max);
        // assert_eq!(public_key.variance, params.variance);
        assert_eq!(public_key.parameters, params);

        // assert_eq!(evaluation_key.dimension_exponent, params.dimension_exponent);
        // assert_eq!(evaluation_key.mul_scaling, params.mul_scaling);
        // assert_eq!(evaluation_key.q_0, params.q_0);
        // assert_eq!(evaluation_key.q, params.q);
        // assert_eq!(evaluation_key.level_max, params.level_max);
        // assert_eq!(evaluation_key.variance, params.variance);
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
    fn test_decrypt_raw_eval_key() {
        // Define parameters for key generation
        let dimension_exponent = 10;
        let q = I256::from(19);
        let level_max = 5;

        // Generate keys using the provided helper function
        let (public_key, evaluation_key, secret_key) =
            generate_keys(dimension_exponent, q.clone(), level_max);

        // Decrypt the raw evaluation key
        let decrypted_eval_key = secret_key.decrypt(&Ciphertext::new(
            evaluation_key.raw_key.clone(),
            level_max,
            (1 << dimension_exponent) as f64, // Not used in this context
            evaluation_key.noise.clean_noise, // Not used in this context
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

    // #[test]
    // fn test_generate_keys_all_parameters() {
    //     let dimension_exponent = 2;
    //     let hamming_weight = 2;
    //     let mul_scaling = I256::from(3);
    //     let q_0 = I256::from(5);
    //     let q = I256::from(4);
    //     let level_max = 2;
    //     let variance = 9.0;

    //     let (public_key, evaluation_key, secret_key) = generate_keys_all_parameters(
    //         dimension_exponent,
    //         hamming_weight,
    //         mul_scaling.clone(),
    //         q_0.clone(),
    //         q.clone(),
    //         level_max,
    //         variance,
    //     );

    //     println!("PublicKey: {:?}", public_key);
    //     println!("EvaluationKey: {:?}", evaluation_key);
    //     println!("SecretKey: {:?}", secret_key);

    //     // Verify that the keys are generated correctly
    //     assert_eq!(secret_key.dimension_exponent, dimension_exponent);
    //     assert_eq!(secret_key.hamming_weight, hamming_weight);
    //     assert_eq!(secret_key.mul_scaling, mul_scaling);
    //     assert_eq!(secret_key.q_0, q_0);
    //     assert_eq!(secret_key.q, q);
    //     assert_eq!(secret_key.level_max, level_max);
    //     assert_eq!(secret_key.variance, variance);

    //     assert_eq!(public_key.dimension_exponent, dimension_exponent);
    //     assert_eq!(public_key.mul_scaling, mul_scaling);
    //     assert_eq!(public_key.q_0, q_0);
    //     assert_eq!(public_key.q, q);
    //     assert_eq!(public_key.level_max, level_max);
    //     assert_eq!(public_key.variance, variance);

    //     assert_eq!(evaluation_key.dimension_exponent, dimension_exponent);
    //     assert_eq!(evaluation_key.mul_scaling, mul_scaling);
    //     assert_eq!(evaluation_key.q_0, q_0);
    //     assert_eq!(evaluation_key.q, q);
    //     assert_eq!(evaluation_key.level_max, level_max);
    //     assert_eq!(evaluation_key.variance, variance);
    // }
}
