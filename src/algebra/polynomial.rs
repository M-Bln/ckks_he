use crate::algebra::arithmetic::{Rescale, RingMod};
use crate::algebra::big_int::{BigInt, Zero};
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

pub fn degree_from_coefs<T>(coefs: &[T]) -> i64
where
    T: Add<Output = T> + Mul<T, Output = T> + Clone + Zero,
{
    let mut degree = coefs.len() as i64 - 1;
    while degree >= 0 {
        if coefs[degree as usize].is_zero() {
            degree -= 1
        } else {
            break;
        }
    }
    degree
}

impl<T> Polynomial<T>
where
    T: Add<Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Clone + Zero,
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

    pub fn degree(&self) -> i64 {
        degree_from_coefs(self.ref_coefficients())
    }

    pub fn degree_truncate(&mut self) -> i64 {
        let degree = degree_from_coefs(self.ref_coefficients());
        self.mut_coefficients().truncate((degree + 1) as usize);
        degree
    }

    pub fn karatsuba_mul(&self, other: &Polynomial<T>) -> Polynomial<T> {
        let n = self.coefficients.len().max(other.coefficients.len());
        if n <= 1 {
            return Polynomial::new(vec![
                self.coefficients[0].clone() * other.coefficients[0].clone(),
            ]);
        }

        let k = n / 2;
        let (low_self, high_self) = self.split(k);
        let (low_other, high_other) = other.split(k);

        let z0 = low_self.karatsuba_mul(&low_other);
        let z2 = high_self.karatsuba_mul(&high_other);
        let z1 = (low_self + &high_self).karatsuba_mul(&(low_other + &high_other))
            - &z0.clone()
            - &z2.clone();

        z0.add_shifted(z1, k).add_shifted(z2, 2 * k)
    }

    fn split(&self, at: usize) -> (Polynomial<T>, Polynomial<T>) {
        let low = Polynomial::new(self.coefficients[..at].to_vec());
        let high = Polynomial::new(self.coefficients[at..].to_vec());
        (low, high)
    }

    fn add_shifted(self, other: Polynomial<T>, shift: usize) -> Polynomial<T> {
        let mut result = self.coefficients;
        if result.len() == 0 {
            return Polynomial::new(result);
        }
        result.resize(
            result.len().max(shift + other.coefficients.len()),
            result[0].zero(),
        );
        for (i, coeff) in other.coefficients.iter().enumerate() {
            result[i + shift] = result[i + shift].clone() + coeff.clone();
        }
        Polynomial::new(result)
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
    T: Add<Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Clone + Zero,
    //    T: Mul<Output = T> + Add<Output = T> + Clone + Zero,
{
    type Output = Self;

    fn mul(self, other: &Polynomial<T>) -> Self::Output {
        self.karatsuba_mul(other)
        // if self.coefficients.is_empty() | other.coefficients.is_empty() {
        //     return Polynomial::new(vec![]);
        // }

        // let result_len = self.coefficients.len() + other.coefficients.len() - 1;
        // let mut result_coeffs = vec![self.coefficients[0].zero(); result_len];

        // for i in 0..self.coefficients.len() {
        //     for j in 0..other.coefficients.len() {
        //         if i + j < result_len {
        //             result_coeffs[i + j] = result_coeffs[i + j].clone()
        //                 + (self.coefficients[i].clone() * other.coefficients[j].clone());
        //         }
        //     }
        // }

        // Polynomial::new(result_coeffs)
    }
}

// impl<T> Mul<&Polynomial<T>> for Polynomial<T>
// where
//     T: Mul<Output = T> + Add<Output = T> + Clone + Zero,
// {
//     type Output = Self;

//     fn mul(self, other: &Polynomial<T>) -> Self::Output {
//         if self.coefficients.is_empty() | other.coefficients.is_empty() {
//             return Polynomial::new(vec![]);
//         }

//         let result_len = self.coefficients.len() + other.coefficients.len() - 1;
//         let mut result_coeffs = vec![self.coefficients[0].zero(); result_len];

//         for i in 0..self.coefficients.len() {
//             for j in 0..other.coefficients.len() {
//                 if i + j < result_len {
//                     result_coeffs[i + j] = result_coeffs[i + j].clone()
//                         + (self.coefficients[i].clone() * other.coefficients[j].clone());
//                 }
//             }
//         }

//         Polynomial::new(result_coeffs)
//     }
// }

pub trait ScalarMul<RHS> {
    type Output;

    fn scalar_mul(self, rhs: RHS) -> Self::Output;
}

impl<T: BigInt> ScalarMul<Polynomial<T>> for T {
    type Output = Polynomial<T>;

    fn scalar_mul(self, rhs: Polynomial<T>) -> Polynomial<T> {
        let coefficients = rhs.coefficients.iter().map(|c| self * *c).collect();
        Polynomial::new(coefficients)
    }
}

impl<'a, T: BigInt> ScalarMul<&'a Polynomial<T>> for T {
    type Output = Polynomial<T>;

    fn scalar_mul(self, rhs: &'a Polynomial<T>) -> Polynomial<T> {
        let coefficients = rhs.coefficients.iter().map(|c| self * *c).collect();
        Polynomial::new(coefficients)
    }
}

impl<T: BigInt> ScalarMul<Polynomial<RingMod<T>>> for T {
    type Output = Polynomial<RingMod<T>>;

    fn scalar_mul(self, rhs: Polynomial<RingMod<T>>) -> Polynomial<RingMod<T>> {
        let coefficients = rhs
            .coefficients
            .iter()
            .map(|c| RingMod::new(self * c.value, c.modulus.clone()))
            .collect();
        Polynomial::new(coefficients)
    }
}

impl<'a, T: BigInt> ScalarMul<&'a Polynomial<RingMod<T>>> for T {
    type Output = Polynomial<RingMod<T>>;

    fn scalar_mul(self, rhs: &'a Polynomial<RingMod<T>>) -> Polynomial<RingMod<T>> {
        let coefficients = rhs
            .coefficients
            .iter()
            .map(|c| RingMod::new(self * c.value, c.modulus.clone()))
            .collect();
        Polynomial::new(coefficients)
    }
}

impl<'a, T: BigInt> ScalarMul<Polynomial<RingMod<T>>> for &'a T {
    type Output = Polynomial<RingMod<T>>;

    fn scalar_mul(self, rhs: Polynomial<RingMod<T>>) -> Polynomial<RingMod<T>> {
        let coefficients = rhs
            .coefficients
            .iter()
            .map(|c| RingMod::new(*self * c.value, c.modulus.clone()))
            .collect();
        Polynomial::new(coefficients)
    }
}

impl<'a, 'b, T: BigInt> ScalarMul<&'a Polynomial<RingMod<T>>> for &'b T {
    type Output = Polynomial<RingMod<T>>;

    fn scalar_mul(self, rhs: &'a Polynomial<RingMod<T>>) -> Polynomial<RingMod<T>> {
        let coefficients = rhs
            .coefficients
            .iter()
            .map(|c| RingMod::new(*self * c.value, c.modulus.clone()))
            .collect();
        Polynomial::new(coefficients)
    }
}

impl<T: BigInt> Rescale<T> for Polynomial<T> {
    fn rescale(&mut self, scalar: T) {
        for coeff in &mut self.coefficients {
            *coeff = coeff.clone() / scalar.clone();
        }
    }
}

impl<T: BigInt> Rescale<T> for Polynomial<RingMod<T>> {
    fn rescale(&mut self, scalar: T) {
        for coeff in &mut self.coefficients {
            coeff.rescale(scalar)
        }
    }
}

impl<T: BigInt> Rescale<RingMod<T>> for Polynomial<RingMod<T>> {
    fn rescale(&mut self, scalar: RingMod<T>) {
        for coeff in &mut self.coefficients {
            coeff.rescale(scalar.value)
        }
    }
}

// impl Polynomial<I256> {
//     pub fn modulo(&self, modulus: I256) -> Polynomial<RingMod<I256>> {
//         let coefficients = self
//             .coefficients
//             .iter()
//             .map(|&coeff| RingMod::<I256>::new(coeff, modulus))
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

    #[test]
    fn test_rescale_polynomial() {
        let coefficients = vec![I256::new(10), I256::new(20), I256::new(30)];
        let mut polynomial = Polynomial::new(coefficients.clone());

        println!("Original Polynomial: {:?}", polynomial);

        let scalar = I256::new(2);
        polynomial.rescale(scalar.clone());

        println!("Rescaled Polynomial: {:?}", polynomial);

        for (i, coeff) in polynomial.coefficients().iter().enumerate() {
            assert_eq!(*coeff, coefficients[i].clone() / scalar.clone());
        }
    }

    #[test]
    fn test_rescale_ringmod_polynomial_by_bigint() {
        let modulus = I256::new(100);
        let value1 = I256::new(50);
        let value2 = I256::new(-25);
        let ringmod1 = RingMod::new(value1.clone(), modulus.clone());
        let ringmod2 = RingMod::new(value2.clone(), modulus.clone());
        let mut polynomial = Polynomial::new(vec![ringmod1.clone(), ringmod2.clone()]);

        println!("Original Polynomial<RingMod<I256>>: {:?}", polynomial);

        let scalar = I256::new(4);
        polynomial.rescale(scalar.clone());

        println!(
            "Rescaled Polynomial<RingMod<I256>> with scalar I256: {:?}",
            polynomial
        );

        assert_eq!(polynomial.ref_coefficients()[0].value, value1 / scalar);
        assert_eq!(polynomial.ref_coefficients()[0].modulus, modulus / scalar);
        assert_eq!(polynomial.ref_coefficients()[1].value, value2 / scalar);
        assert_eq!(polynomial.ref_coefficients()[1].modulus, modulus / scalar);
    }

    #[test]
    fn test_rescale_ringmod_polynomial() {
        let modulus = I256::new(100);
        let value1 = I256::new(50);
        let value2 = I256::new(-25);
        let ringmod1 = RingMod::new(value1.clone(), modulus.clone());
        let ringmod2 = RingMod::new(value2.clone(), modulus.clone());
        let mut polynomial = Polynomial::new(vec![ringmod1.clone(), ringmod2.clone()]);

        println!("Original Polynomial<RingMod<I256>>: {:?}", polynomial);

        let scalar_value = I256::new(5);
        let scalar_ringmod = RingMod::new(scalar_value.clone(), modulus.clone());
        polynomial.rescale(scalar_ringmod);

        println!(
            "Rescaled Polynomial<RingMod<I256>> with scalar RingMod<I256>: {:?}",
            polynomial
        );

        assert_eq!(
            polynomial.ref_coefficients()[0].value,
            value1 / scalar_value
        );
        assert_eq!(
            polynomial.ref_coefficients()[0].modulus,
            modulus / scalar_value
        );
        assert_eq!(
            polynomial.ref_coefficients()[1].value,
            value2 / scalar_value
        );
        assert_eq!(
            polynomial.ref_coefficients()[1].modulus,
            modulus / scalar_value
        );
    }

    #[test]
    fn test_polynomial_scalar_multiplication() {
        let poly = Polynomial::new(vec![I256::from(1), I256::from(2), I256::from(3)]);
        let scalar = I256::from(2);
        let result = scalar.scalar_mul(poly);
        let expected_coefficients = vec![I256::from(2), I256::from(4), I256::from(6)];

        for (i, coeff) in result.coefficients().iter().enumerate() {
            assert_eq!(
                coeff, &expected_coefficients[i],
                "Coefficient mismatch at index {}: expected {}, got {}",
                i, expected_coefficients[i], coeff
            );
        }
    }

    #[test]
    fn test_polynomial_ringmod_scalar_multiplication() {
        let poly = Polynomial::new(vec![
            RingMod::new(I256::from(1), I256::from(7)),
            RingMod::new(I256::from(2), I256::from(7)),
            RingMod::new(I256::from(3), I256::from(7)),
        ]);
        let scalar = I256::from(2);
        let result = scalar.scalar_mul(poly);
        let expected_coefficients = vec![
            RingMod::new(I256::from(2), I256::from(7)),
            RingMod::new(I256::from(4), I256::from(7)),
            RingMod::new(I256::from(6), I256::from(7)),
        ];

        for (i, coeff) in result.coefficients().iter().enumerate() {
            assert_eq!(
                coeff, &expected_coefficients[i],
                "Coefficient mismatch at index {}: expected {:?}, got {:?}",
                i, expected_coefficients[i], coeff
            );
        }
    }

    #[test]
    fn test_polynomial_ringmod_scalar_multiplication_ref() {
        let poly = Polynomial::new(vec![
            RingMod::new(I256::from(1), I256::from(7)),
            RingMod::new(I256::from(2), I256::from(7)),
            RingMod::new(I256::from(3), I256::from(7)),
        ]);
        let scalar = I256::from(2);
        let result = scalar.scalar_mul(&poly);
        let expected_coefficients = vec![
            RingMod::new(I256::from(2), I256::from(7)),
            RingMod::new(I256::from(4), I256::from(7)),
            RingMod::new(I256::from(6), I256::from(7)),
        ];

        for (i, coeff) in result.coefficients().iter().enumerate() {
            assert_eq!(
                coeff, &expected_coefficients[i],
                "Coefficient mismatch at index {}: expected {:?}, got {:?}",
                i, expected_coefficients[i], coeff
            );
        }
    }

    #[test]
    fn test_polynomial_ringmod_scalar_multiplication_ref_scalar() {
        let poly = Polynomial::new(vec![
            RingMod::new(I256::from(1), I256::from(7)),
            RingMod::new(I256::from(2), I256::from(7)),
            RingMod::new(I256::from(3), I256::from(7)),
        ]);
        let scalar = I256::from(2);
        let result = (&scalar).scalar_mul(poly);
        let expected_coefficients = vec![
            RingMod::new(I256::from(2), I256::from(7)),
            RingMod::new(I256::from(4), I256::from(7)),
            RingMod::new(I256::from(6), I256::from(7)),
        ];

        for (i, coeff) in result.coefficients().iter().enumerate() {
            assert_eq!(
                coeff, &expected_coefficients[i],
                "Coefficient mismatch at index {}: expected {:?}, got {:?}",
                i, expected_coefficients[i], coeff
            );
        }
    }

    #[test]
    fn test_polynomial_ringmod_scalar_multiplication_ref_all() {
        let poly = Polynomial::new(vec![
            RingMod::new(I256::from(1), I256::from(7)),
            RingMod::new(I256::from(2), I256::from(7)),
            RingMod::new(I256::from(3), I256::from(7)),
        ]);
        let scalar = I256::from(2);
        let result = (&scalar).scalar_mul(&poly);
        let expected_coefficients = vec![
            RingMod::new(I256::from(2), I256::from(7)),
            RingMod::new(I256::from(4), I256::from(7)),
            RingMod::new(I256::from(6), I256::from(7)),
        ];

        for (i, coeff) in result.coefficients().iter().enumerate() {
            assert_eq!(
                coeff, &expected_coefficients[i],
                "Coefficient mismatch at index {}: expected {:?}, got {:?}",
                i, expected_coefficients[i], coeff
            );
        }
    }

    #[test]
    fn test_degree() {
        let poly1 = Polynomial::new(vec![1, 2, 3]);
        assert_eq!(poly1.degree(), 2, "Degree should be 2");

        let poly2 = Polynomial::new(vec![0]);
        assert_eq!(poly2.degree(), -1, "Degree should be -1");

        let poly3 = Polynomial::new(vec![5, 0, 0, 0]);
        assert_eq!(poly3.degree(), 0, "Degree should be 3");

        let poly4: Polynomial<i64> = Polynomial::new(vec![]);
        assert_eq!(
            poly4.degree(),
            -1,
            "Degree should be -1 for an empty polynomial"
        );
    }

    #[test]
    fn test_degree_truncate() {
        let mut poly1 = Polynomial::new(vec![1, 2, 3, 0, 0]);
        poly1.degree_truncate();
        assert_eq!(
            poly1.ref_coefficients(),
            &[1, 2, 3],
            "Truncated polynomial should be [1, 2, 3]"
        );
        assert_eq!(poly1.degree(), 2, "Degree should be 2 after truncation");

        let mut poly2 = Polynomial::new(vec![0, 0, 0, 0]);
        poly2.degree_truncate();
        assert_eq!(poly2.ref_coefficients(), &[]);
        assert_eq!(poly2.degree(), -1, "Degree should be -1");

        let mut poly4: Polynomial<i64> = Polynomial::new(vec![]);
        poly4.degree_truncate();
        assert_eq!(
            poly4.ref_coefficients(),
            &[],
            "Empty polynomial should remain empty"
        );
        assert_eq!(
            poly4.degree(),
            -1,
            "Degree should be 0 for an empty polynomial"
        );
    }
}
