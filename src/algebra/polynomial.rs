use crate::algebra::arithmetic::RingMod;
use crate::algebra::big_int::{BigInt, Zero};
use std::ops::{Add, Mul};

#[derive(Clone, PartialEq, Eq, Debug)]
struct Polynomial<T> {
    coefficients: Vec<T>,
}

impl<T> Polynomial<T> {
    pub fn new(coefficients: Vec<T>) -> Self {
        Polynomial { coefficients }
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
    use crate::algebra::big_int::BigInt;

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
}
