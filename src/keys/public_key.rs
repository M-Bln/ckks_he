use crate::algebra::big_int::BigInt;
use crate::algebra::polynomial::Polynomial;
use crate::ciphertext::{Ciphertext, Message, RawCiphertext};
use crate::keys::key_generator::KeyGenerationParameters;
use crate::random_distributions::{DiscreteGaussian, ZODistribution};

#[derive(Clone, Debug)]
pub struct PublicKey<T: BigInt> {
    //    pub dimension_exponent: u32, // Ciphertexts in cyclotomic ring $R[X]/(1+X^N)$ with N = 2^dimension_exponent
    //    pub hamming_weight: usize,   // Number of non zero coefficients in rawKey.1
    //    pub mul_scaling: T,          // Rescaling factor used in homomorphic multiplication
    //    pub q_0: T,                  // minimal modulus
    //    pub q: T, // modulus per level, the total modulus in coefficients ring is initialy mul_scaling * q_0 * q^{level_max}
    //    pub level_max: u32,
    //    pub variance: f64, // standard deviation of the Gaussian distribution of error sampling
    //    pub standard_deviation: f64,
    pub raw_key: RawCiphertext<T>,
    pub parameters: KeyGenerationParameters<T>,
    pub noise: ComputationNoise,
    gaussian: DiscreteGaussian,
    zod: ZODistribution,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComputationNoise {
    pub clean_noise: f64,
    pub rescaling_noise: f64,
    pub key_switch_noise: f64,
    pub inv_mul_scaling: f64,
}

impl ComputationNoise {
    pub fn new<T: BigInt>(parameters: KeyGenerationParameters<T>) -> Self {
        let dimension = (1 << parameters.dimension_exponent) as f64;
        let h = parameters.hamming_weight as f64;
        let clean_noise = parameters.standard_deviation
            * (8.0 * dimension * (2.0 as f64).sqrt()
                + 6.0 * dimension.sqrt()
                + 16.0 * (h * dimension).sqrt());
        let rescaling_noise = (dimension / 3.0).sqrt() * (3.0 + 8.0 * h.sqrt());
        let key_switch_noise =
            8.0 * parameters.standard_deviation * dimension / (3.0 as f64).sqrt();
        let inv_mul_scaling = 1.0 / parameters.mul_scaling.to_float();

        ComputationNoise {
            clean_noise,
            rescaling_noise,
            key_switch_noise,
            inv_mul_scaling,
        }
    }
    // pub fn mul_noise<T: BigInt>(&self, level: T) -> f64 {
    //     self.inv_mul_scaling * self.key_switch_noise * level.to_float() + self.rescaling_noise
    // }
}

impl<T: BigInt> PublicKey<T> {
    pub fn new(
        // dimension_exponent: u32,
        // hamming_weight: usize,
        // mul_scaling: T,
        // q_0: T,
        // q: T,
        // level_max: u32,
        // variance: f64,
        parameters: KeyGenerationParameters<T>,
        noise: ComputationNoise,
        raw_key: RawCiphertext<T>,
    ) -> Self {
        // let standard_deviation = variance.sqrt();
        let dimension = 2_usize.pow(parameters.dimension_exponent);
        let _dimension_f = dimension as f64;
        let _hamming_f = parameters.hamming_weight as f64;
        let gaussian = DiscreteGaussian::new(
            0.0,
            parameters.standard_deviation * parameters.standard_deviation,
        );
        let zod = ZODistribution::new(0.5);
        PublicKey {
            // dimension_exponent,
            // hamming_weight,
            // mul_scaling,
            // q_0,
            // q,
            // level_max,
            // variance,
            // standard_deviation,
            parameters,
            raw_key,
            noise,
            gaussian,
            zod,
        }
    }

    pub fn encrypt(&mut self, message: &Message<T>, upper_bound_message: f64) -> Ciphertext<T> {
        let raw = self.encrypt_raw(message);
        let level = self.parameters.level_max;
        let upper_bound_error = self.noise.clean_noise;
        Ciphertext::new(raw, level, upper_bound_message, upper_bound_error)
    }

    pub fn encrypt_raw(&mut self, message: &Message<T>) -> RawCiphertext<T> {
        let modulus = self.parameters.q.fast_exp(self.parameters.level_max) * self.parameters.q_0;
        let dimension = 2_usize.pow(self.parameters.dimension_exponent);

        let small_factor_coefficients = self.zod.sample_n(dimension);
        let small_factor = Polynomial::new(small_factor_coefficients)
            .modulo(modulus)
            .to_cyclotomic(self.parameters.dimension_exponent);

        let e0_coefficients = self.gaussian.sample_n(dimension);
        let e0 = Polynomial::new(e0_coefficients)
            .modulo(modulus)
            .to_cyclotomic(self.parameters.dimension_exponent);

        let e1_coefficients = self.gaussian.sample_n::<T>(dimension);
        let e1 = Polynomial::new(e1_coefficients)
            .modulo(modulus)
            .to_cyclotomic(self.parameters.dimension_exponent);
        RawCiphertext(
            small_factor.clone() * &self.raw_key.0 + &e0 + message,
            small_factor * &self.raw_key.1 + &e1,
        )
    }
}
