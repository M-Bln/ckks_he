//use crate::algebra::arithmetic::RingMod;
use crate::algebra::big_int::Zero;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::arithmetic::RingMod;

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
}
