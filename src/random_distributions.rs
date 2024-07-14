use crate::algebra::big_int::BigInt;
use bnum::random::UniformInt as UniformBInt;
use bnum::types::I256;
use rand::distributions::uniform::{UniformInt, UniformSampler};
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;
use rand::Rng;
use rand_distr::{Distribution, Normal};

pub trait UniformSamplable: Sized + PartialOrd {
    type Sampler: UniformSampler<X = Self>;

    fn sampler(low: Self, high: Self) -> Self::Sampler;
}

impl UniformSamplable for I256 {
    type Sampler = UniformBInt<I256>;
    fn sampler(low: Self, high: Self) -> Self::Sampler {
        Self::Sampler::new(low, high)
    }
}

impl UniformSamplable for i64 {
    type Sampler = UniformInt<i64>;
    fn sampler(low: Self, high: Self) -> Self::Sampler {
        Self::Sampler::new(low, high)
    }
}

// impl UniformSampler for Sampler<I256> {
//     type Item = I256;
//     fn sample(&mut self, min: I256, max: I256) -> I256 {

//     }
// }

// impl<T> UniformSampler<T>
// where
//     T: PartialOrd + Copy,
// {
//     pub fn sample(&mut self, min: T, max: T) -> T {
//         let distribution = Uniform::new_inclusive(min, max);
//         distribution.sample(&mut self.rng)
//     }
// }

pub struct DiscreteGaussian {
    mean: f64,
    std_dev: f64,
    normal: Normal<f64>,
    rng: ThreadRng,
}

impl DiscreteGaussian {
    /// Creates a new DiscreteGaussian with the given mean and variance.
    pub fn new(mean: f64, variance: f64) -> Self {
        let std_dev = variance.sqrt();
        let normal = Normal::new(mean, std_dev).unwrap();
        let rng = rand::thread_rng();
        DiscreteGaussian {
            mean,
            std_dev,
            normal,
            rng,
        }
    }

    /// Generates a random integer according to the discrete Gaussian distribution.
    pub fn sample<T: BigInt>(&mut self) -> T {
        let sample = self.normal.sample(&mut self.rng);
        T::from(sample.round() as i64)
    }

    /// Generates a vector of random integers according to the discrete Gaussian distribution.
    pub fn sample_n<T: BigInt>(&mut self, n: usize) -> Vec<T> {
        let mut result = Vec::<T>::with_capacity(n);
        for _ in 0..n {
            result.push(self.sample());
        }
        result
    }
}

pub struct ZODistribution {
    rho: f64,
    rng: ThreadRng,
}

impl ZODistribution {
    /// Creates a new ZODistribution with the given rho.
    pub fn new(rho: f64) -> Self {
        let rng = rand::thread_rng();
        ZODistribution { rho, rng }
    }

    /// Generates a random BigInt according to the custom distribution.
    pub fn sample<T: BigInt>(&mut self) -> T {
        let prob: f64 = self.rng.gen();
        let value = if prob < self.rho / 2.0 {
            -1
        } else if prob < self.rho {
            1
        } else {
            0
        };
        T::from(value)
    }

    /// Generates a vector of random BigInts according to the custom distribution.
    pub fn sample_n<T: BigInt>(&mut self, n: usize) -> Vec<T> {
        (0..n).map(|_| self.sample()).collect()
    }
}

pub struct HWTDistribution {
    n: usize,
    h: usize,
    rng: ThreadRng,
}

impl HWTDistribution {
    /// Creates a new HWTDistribution with given length `n` and exactly `h` non-zero coefficients.
    pub fn new(n: usize, h: usize) -> Self {
        let rng = rand::thread_rng();
        HWTDistribution { n, h, rng }
    }

    /// Generates a random vector of length `n` with exactly `h` non-zero BigInt coefficients.
    pub fn sample<T: From<i64> + Clone>(&mut self) -> Vec<T> {
        let mut vec = vec![T::from(0); self.n];
        let mut indices: Vec<usize> = (0..self.n).collect();
        indices.shuffle(&mut self.rng);

        for &index in &indices[..self.h] {
            vec[index] = if self.rng.gen() {
                T::from(1)
            } else {
                T::from(-1)
            };
        }

        vec
    }
}

// pub struct HWTDistribution {
//     weight: usize,
//     rng: ThreadRng,
// }

// impl HWTDistribution {
//     pub fn new(weight: usize) -> Self {
// 	let rng = rand::thread_rng();
// 	HWTDistribution{
// 	    weight, rng,
// 	}
//     }

//     pub fn sample<T: BigInt>(&mut self) -> T {
// 	if self.rng.gen() {
// 	    T::from(1 as i64)
// 	} else {
// 	    T::from(-1 as i64)
// 	}
//     }

//     pub fn sample_n<T: BigInt>(&mut self, n: usize) -> Vec<T> {
// 	assert!(n >= weight)
//     }
// }

// /// Generates a random f64 number according to a Gaussian distribution with given mean and variance.
// pub fn random_gaussian(mean: f64, variance: f64) -> f64 {
//     let std_dev = variance.sqrt();
//     let normal = Normal::new(mean, std_dev).unwrap();
//     let mut rng = rand::thread_rng();
//     normal.sample(&mut rng)
// }

// /// Generates a random integer according to a Gaussian distribution with given mean and variance.
// pub fn random_gaussian_integer(mean: f64, variance: f64) -> i64 {
//     let std_dev = variance.sqrt();
//     let normal = Normal::new(mean, std_dev).unwrap();
//     let mut rng = rand::thread_rng();
//     let sample = normal.sample(&mut rng);
//     sample.round() as i64
// }

#[cfg(test)]
mod tests {
    use super::*;
    use bnum::types::I256;

    #[test]
    fn test_random_gaussian_bigint() {
        let mean = 0.0;
        let variance = 3.2 * 3.2;
        let mut discrete_gaussian = DiscreteGaussian::new(mean, variance);

        for _ in 0..10 {
            let random_integer: I256 = discrete_gaussian.sample();
            println!(
                "Random BigInt from Gaussian distribution: {}",
                random_integer
            );
        }
    }

    #[test]
    fn test_sample_n() {
        let mean = 0.0;
        let variance = 3.2 * 3.2;
        let mut discrete_gaussian = DiscreteGaussian::new(mean, variance);

        let random_integers: Vec<I256> = discrete_gaussian.sample_n(10);
        println!(
            "Random BigInt vector from Gaussian distribution: {:?}",
            random_integers
        );
        assert_eq!(random_integers.len(), 10);
    }

    #[test]
    fn test_zo_distribution() {
        let rho = 0.5;
        let mut zo_distribution = ZODistribution::new(rho);

        for _ in 0..30 {
            let random_integer: I256 = zo_distribution.sample();
            println!("Random BigInt from custom distribution: {}", random_integer);
        }
    }

    #[test]
    fn test_hwt_distribution() {
        let n = 30;
        let h = 8;
        let mut hwt_distribution = HWTDistribution::new(n, h);

        let result: Vec<I256> = hwt_distribution.sample();
        println!("Random vector: {:?}", result);

        let non_zero_count = result.iter().filter(|&&x| x != I256::from(0)).count();
        assert_eq!(non_zero_count, h);
    }

    #[test]
    fn test_uniform_sampler_i64() {
        let min = -1;
        let max = 1;
        let sampler = i64::sampler(min, max);
        let mut rng = rand::thread_rng();
        for _ in 0..30 {
            let sample = sampler.sample(&mut rng);
            println!("Uniform sample (i64): {}", sample);
            assert!(sample >= min && sample <= max);
        }
    }

    #[test]
    fn test_uniform_sampler_i256() {
        let min = I256::from(-1);
        let max = I256::from(1);
        let sampler = I256::sampler(min.clone(), max.clone());
        let mut rng = rand::thread_rng();
        for _ in 0..30 {
            let sample = sampler.sample(&mut rng);
            println!("Uniform sample (I256): {}", sample);
            assert!(sample >= min && sample <= max);
        }
    }
}
