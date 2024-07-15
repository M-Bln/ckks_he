use crate::algebra::big_int::BigInt;
use crate::algebra::polynomial::Polynomial;
use crate::ciphertext::{Ciphertext, CiphertextRing, RawCiphertext};
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
    pub encryption_error: f64,
    pub parameters: KeyGenerationParameters<T>,
    rng: rand::rngs::ThreadRng,
    gaussian: DiscreteGaussian,
    zod: ZODistribution,
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
        raw_key: RawCiphertext<T>,
    ) -> Self {
        // let standard_deviation = variance.sqrt();
        let dimension = 2_usize.pow(parameters.dimension_exponent);
        let dimension_f = dimension as f64;
        let hamming_f = parameters.hamming_weight as f64;
        let encryption_error = parameters.standard_deviation
            * (dimension_f * 8.0 * (2.0 as f64).sqrt()
                + 6.0 * dimension_f.sqrt()
                + 16.0 * (hamming_f * dimension_f).sqrt());
        let mut rng = rand::thread_rng();
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
            encryption_error,
            rng,
            gaussian,
            zod,
        }
    }

    pub fn encrypt(
        &mut self,
        message: &CiphertextRing<T>,
        upper_bound_message: f64,
    ) -> Ciphertext<T> {
        let raw = self.encrypt_raw(message);
        let level = self.parameters.level_max;
        let upper_bound_error = self.encryption_error;
        Ciphertext::new(raw, level, upper_bound_message, upper_bound_error)
    }

    pub fn encrypt_raw(&mut self, message: &CiphertextRing<T>) -> RawCiphertext<T> {
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
