use crate::algebra::arithmetic::RingMod;
use crate::algebra::big_int::Zero;
use crate::algebra::complex::{Complex, C64};
use crate::algebra::cyclotomic_ring::CyclotomicRing;
use bnum::types::I256;
use std::ops::{Add, Mul, Sub};

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Polynomial<T> {
    coefficients: Vec<T>,
}

impl<T> Polynomial<T> {
    pub fn new(coefficients: Vec<T>) -> Self {
        Polynomial { coefficients }
    }

    pub fn coefficients(self) -> Vec<T> {
        self.coefficients
    }

    pub fn ref_coefficients(&self) -> &[T] {
        &self.coefficients[..]
    }

    pub fn mut_coefficients(&mut self) -> &mut Vec<T> {
        &mut self.coefficients
    }
}

impl<T> Polynomial<T>
where
    T: Add<Output = T> + Mul<T, Output = T> + Clone + Zero,
{
    pub fn eval(&self, x: T) -> T
    where
        T: Add<Output = T> + Mul<Output = T> + Clone,
    {
        let mut result = self
            .coefficients
            .last()
            .cloned()
            .unwrap_or_else(|| x.zero());
        for coeff in self.coefficients.iter().rev().skip(1) {
            result = coeff.clone() + result * x.clone();
        }
        result
    }
}

impl<T> Add<&Polynomial<T>> for Polynomial<T>
where
    T: Add<Output = T> + Clone,
{
    type Output = Self;

    fn add(self, other: &Polynomial<T>) -> Self::Output {
        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut result_coeffs = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let coeff1 = self.coefficients.get(i).cloned();
            let coeff2 = other.coefficients.get(i).cloned();

            let sum = match (coeff1, coeff2) {
                (Some(c1), Some(c2)) => c1 + c2,
                (Some(c1), None) => c1,
                (None, Some(c2)) => c2,
                (None, None) => unreachable!(),
            };

            result_coeffs.push(sum);
        }

        Polynomial::new(result_coeffs)
    }
}

impl<T> Sub<&Polynomial<T>> for Polynomial<T>
where
    T: Sub<Output = T> + Clone + Zero,
{
    type Output = Self;

    fn sub(self, other: &Polynomial<T>) -> Self::Output {
        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut result_coeffs = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let coeff1 = self.coefficients.get(i).cloned();
            let coeff2 = other.coefficients.get(i).cloned();

            let sum = match (coeff1, coeff2) {
                (Some(c1), Some(c2)) => c1 - c2,
                (Some(c1), None) => c1,
                (None, Some(c2)) => c2.zero() - c2,
                (None, None) => unreachable!(),
            };

            result_coeffs.push(sum);
        }

        Polynomial::new(result_coeffs)
    }
}

impl<T> Mul<&Polynomial<T>> for Polynomial<T>
where
    T: Mul<Output = T> + Add<Output = T> + Clone + Zero,
{
    type Output = Self;

    fn mul(self, other: &Polynomial<T>) -> Self::Output {
        if self.coefficients.is_empty() | other.coefficients.is_empty() {
            return Polynomial::new(vec![]);
        }

        let result_len = self.coefficients.len() + other.coefficients.len() - 1;
        let mut result_coeffs = vec![self.coefficients[0].zero(); result_len];

        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                if i + j < result_len {
                    result_coeffs[i + j] = result_coeffs[i + j].clone()
                        + (self.coefficients[i].clone() * other.coefficients[j].clone());
                }
            }
        }

        Polynomial::new(result_coeffs)
    }
}

impl Polynomial<I256> {
    // pub fn to_f64(&self) -> Polynomial<f64> {
    //     let coefficients = self
    //         .coefficients
    //         .iter()
    //         .map(|&coeff| i256_to_f64(coeff))
    //         .collect();
    //     Polynomial { coefficients }
    // }

    pub fn modulo(&self, modulus: I256) -> Polynomial<RingMod<I256>> {
        let coefficients = self
            .coefficients
            .iter()
            .map(|&coeff| RingMod::<I256>::new(coeff, modulus))
            .collect();
        Polynomial { coefficients }
    }
}

// impl Polynomial<RingMod<I256>> {
//     pub fn to_cyclotomic(self, dimension_exponent: u32) -> CyclotomicRing<RingMod<I256>> {
//         CyclotomicRing::new(self.coefficients(), 2_usize.pow(dimension_exponent))
//     }
// }

// impl Polynomial<f64> {
//     pub fn to_i256(&self) -> Polynomial<I256> {
//         let coefficients = self
//             .coefficients
//             .iter()
//             .map(|&coeff| f64_to_i256(coeff))
//             .collect();
//         Polynomial { coefficients }
//     }
// }

// impl Polynomial<C64> {
//     pub fn to_i256(&self) -> Polynomial<I256> {
//         let coefficients = self
//             .coefficients
//             .iter()
//             .map(|&coeff| f64_to_i256(coeff.real()))
//             .collect();
//         Polynomial { coefficients }
//     }
// }

// impl Polynomial<RingMod<I256>> {
//     pub fn to_c64(&self) -> Polynomial<C64> {
//         let coefficients = self
//             .coefficients
//             .iter()
//             .map(|coeff| coeff.to_c64())
//             .collect();
//         Polynomial { coefficients }
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::arithmetic::RingMod;
    use bnum::types::I256;

    #[test]
    fn test_polynomial_creation() {
        let coefficients = vec![
            RingMod::new(10, 13),
            RingMod::new(5, 13),
            RingMod::new(3, 13),
        ];
        let poly = Polynomial::new(coefficients.clone());
        assert_eq!(poly.coefficients, coefficients);
    }

    #[test]
    fn test_polynomial_addition() {
        let modulus = 13;
        let poly1 = Polynomial::new(vec![
            RingMod::new(1, modulus),
            RingMod::new(2, modulus),
            RingMod::new(3, modulus),
        ]);
        let poly2 = Polynomial::new(vec![
            RingMod::new(4, modulus),
            RingMod::new(5, modulus),
            RingMod::new(6, modulus),
        ]);

        let expected_sum = Polynomial::new(vec![
            RingMod::new(5, modulus),
            RingMod::new(7, modulus),
            RingMod::new(9, modulus),
        ]);

        let result_sum = poly1.clone() + &poly2;
        assert_eq!(result_sum, expected_sum);
    }

    #[test]
    fn test_polynomial_multiplication() {
        let modulus = 13;
        let poly1 = Polynomial::new(vec![
            RingMod::new(1, modulus),
            RingMod::new(2, modulus),
            RingMod::new(3, modulus),
        ]);
        let poly2 = Polynomial::new(vec![
            RingMod::new(4, modulus),
            RingMod::new(5, modulus),
            RingMod::new(6, modulus),
        ]);

        let expected_product = Polynomial::new(vec![
            RingMod::new(4, modulus),
            RingMod::new(13, modulus),
            RingMod::new(28, modulus),
            RingMod::new(27, modulus),
            RingMod::new(18, modulus),
        ]);

        let result_product = poly1.clone() * &poly2;
        assert_eq!(result_product, expected_product);
    }

    #[test]
    fn test_polynomial_evaluation() {
        let poly = Polynomial::new(vec![1, 2, 3]); // Represents 1 + 2x + 3x^2
        let value = poly.eval(2); // Evaluate at x = 2
        assert_eq!(value, 17); // 1 + 2*2 + 3*2^2 = 17

        let poly = Polynomial::new(vec![0, 1, 0, 1]); // Represents x + x^3
        let value = poly.eval(2); // Evaluate at x = 2
        assert_eq!(value, 10); // 2 + 2^3 = 10
    }

    #[test]
    fn test_polynomial_multiplication_256() {
        let modulus = I256::from(13);
        let poly1 = Polynomial::new(vec![
            RingMod::new(I256::from(1), modulus),
            RingMod::new(I256::from(2), modulus),
            RingMod::new(I256::from(3), modulus),
        ]);
        let poly2 = Polynomial::new(vec![
            RingMod::new(I256::from(4), modulus),
            RingMod::new(I256::from(5), modulus),
            RingMod::new(I256::from(6), modulus),
        ]);

        let expected_product = Polynomial::new(vec![
            RingMod::new(I256::from(4), modulus),
            RingMod::new(I256::from(13 % 13), modulus),
            RingMod::new(I256::from(28 % 13), modulus),
            RingMod::new(I256::from(27 % 13), modulus),
            RingMod::new(I256::from(18 % 13), modulus),
        ]);

        let result_product = poly1.clone() * &poly2;
        assert_eq!(result_product, expected_product);
    }

    #[test]
    fn test_polynomial_evaluation_256() {
        let modulus = I256::from(13);
        let poly = Polynomial::new(vec![
            RingMod::new(I256::from(1), modulus),
            RingMod::new(I256::from(2), modulus),
            RingMod::new(I256::from(3), modulus),
        ]);
        let x = RingMod::new(I256::from(2), modulus);
        let value = poly.eval(x); // Evaluate at x = 2
        assert_eq!(value.value(), I256::from(17 % 13)); // 1 + 2*2 + 3*2^2 = 17 % 13 = 4
    }

    // #[test]
    // fn test_polynomial_conversion() {
    //     // Create a Polynomial<I256>
    //     let poly_i256 = Polynomial {
    //         coefficients: vec![I256::from(1), I256::from(2), I256::from(3)],
    //     };

    //     // Convert to Polynomial<f64>
    //     let poly_f64 = poly_i256.to_f64();
    //     assert_eq!(poly_f64.coefficients, vec![1.0, 2.0, 3.0]);

    //     // Convert back to Polynomial<I256>
    //     let poly_i256_converted_back = poly_f64.to_i256();
    //     assert_eq!(poly_i256, poly_i256_converted_back);
    // }

    // #[test]
    // fn test_polynomial_conversion() {
    //     // Create a Polynomial<I256> with large coefficients
    //     let poly_i256 = Polynomial {
    //         coefficients: vec![
    //             I256::from(12345678901234567890u64),
    //             I256::from(3765432109876543210u64),
    //             I256::from(11223344556677889900u64),
    //         ],
    //     };

    //     println!("Original Polynomial<I256>: {:?}", poly_i256);

    //     // Convert to Polynomial<f64>
    //     let poly_f64 = poly_i256.to_f64();
    //     println!("Converted to Polynomial<f64>: {:?}", poly_f64);

    //     // Convert back to Polynomial<I256>
    //     let poly_i256_converted_back = poly_f64.to_i256();
    //     println!("Converted back to Polynomial<I256>: {:?}", poly_i256_converted_back);

    //     // Check if the original polynomial and the converted-back polynomial are equal
    //     for (a, b) in poly_i256.coefficients.iter().zip(poly_i256_converted_back.coefficients.iter()) {
    //         assert_eq!(a, b, "Mismatch: original {} != converted back {}", a, b);
    //     }
    // }
}
