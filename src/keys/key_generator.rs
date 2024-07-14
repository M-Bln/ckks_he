use crate::algebra::big_int::BigInt;
use crate::algebra::polynomial::{Polynomial, ScalarMul};
use crate::ciphertext::RawCiphertext;
use crate::keys::evaluation_key::EvaluationKey;
use crate::keys::public_key::PublicKey;
use crate::keys::secret_key::SecretKey;
use crate::random_distributions::{sample_n, DiscreteGaussian};
use rand::distributions::uniform::UniformSampler;

pub fn generate_keys_all_parameters<T: BigInt>(
    dimension_exponent: u32,
    hamming_weight: usize,
    mul_scaling: T,
    q_0: T,
    q: T,
    level_max: u32,
    variance: f64,
) -> (PublicKey<T>, EvaluationKey<T>, SecretKey<T>) {
    let secret_key = SecretKey::<T>::new(
        dimension_exponent,
        hamming_weight,
        mul_scaling,
        q_0,
        q,
        level_max,
        variance,
    );
    let dimension = 2_usize.pow(dimension_exponent);
    let mut rng = rand::thread_rng();
    let modulus = q_0 * q.fast_exp(level_max); // the total modulus in coefficients ring is initialy mul_scaling * q_0 * q^{level_max}
    let sampler = T::sampler(T::from(0), modulus);
    let public_key_coefficients = sample_n(sampler, dimension, &mut rng);
    let public_key_a = Polynomial::<T>::new(public_key_coefficients)
        .modulo(modulus)
        .to_cyclotomic(dimension_exponent);

    let mut gaussian_sampler = DiscreteGaussian::new(0.0, variance);
    let error_coefficients = gaussian_sampler.sample_n(dimension);
    let error = Polynomial::<T>::new(error_coefficients)
        .modulo(modulus)
        .to_cyclotomic(dimension_exponent);
    let raw_public_key = RawCiphertext::<T>(
        public_key_a.clone(),
        error - &(public_key_a * &secret_key.key_s),
    );

    let public_key = PublicKey::<T>::new(
        dimension_exponent,
        mul_scaling,
        q_0,
        q,
        level_max,
        variance,
        raw_public_key,
    );

    let modulus_eval = modulus * &mul_scaling;
    let sampler_eval = T::sampler(T::from(0), modulus_eval);
    let eval_key_coefficients = sample_n(sampler_eval, dimension, &mut rng);
    let eval_key_a = Polynomial::<T>::new(eval_key_coefficients)
        .modulo(modulus_eval)
        .to_cyclotomic(dimension_exponent);

    let eval_error_coefficients = gaussian_sampler.sample_n(dimension);
    let eval_error = Polynomial::<T>::new(eval_error_coefficients)
        .modulo(modulus_eval)
        .to_cyclotomic(dimension_exponent);
    let raw_eval_key = RawCiphertext::<T>(
        eval_key_a.clone(),
        (eval_error - &eval_key_a)
            + &(mul_scaling.scalar_mul(secret_key.key_s.clone()) * &secret_key.key_s),
    );

    let evaluation_key = EvaluationKey::<T>::new(
        dimension_exponent,
        mul_scaling,
        q_0,
        q,
        level_max,
        variance,
        raw_eval_key,
    );
    (public_key, evaluation_key, secret_key)
}
