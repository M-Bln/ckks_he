use crate::algebra::arithmetic::RingMod;
use crate::algebra::big_int::BigInt;
use std::ops::Add;

// Polynomial which are always checked to remove leading zeros so that degree = coefficients.len() -1
#[derive(Clone, PartialEq, Eq, Debug)]
struct CheckedPolynomial<T> {
    coefficients: Vec<T>,
}

impl<T> CheckedPolynomial<T> {
    pub fn new(coefficients: Vec<T>) -> Self {
        CheckedPolynomial { coefficients }
    }

    // Get the degree of the polynomial, we agree that degree(0) = 0
    pub fn degree(&self) -> usize {
        if self.coefficients.len() >= 1 {
            self.coefficients.len() - 1
        } else {
            0
        }
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }
}

impl<T: BigInt> Add<&CheckedPolynomial<RingMod<T>>> for CheckedPolynomial<RingMod<T>> {
    type Output = Self;
    fn add(self, other: &CheckedPolynomial<RingMod<T>>) -> CheckedPolynomial<RingMod<T>> {
        if self.is_zero() & other.is_zero() {
            return CheckedPolynomial::<RingMod<T>>::new(vec![]);
        }
        if other.is_zero() {
            return self;
        }
        if self.is_zero() {
            return other.clone();
        }
        let modulus: T = self.coefficients[0].modulus();
        let zero = RingMod::<T>::new(T::from(0), modulus);

        let max_length = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut current_length = 0;
        let mut coefficients: Vec<RingMod<T>> = Vec::new();
        for i in 0..max_length {
            let mut coeff = RingMod::<T>::new(T::from(0), modulus);
            if i < self.coefficients.len() {
                coeff = coeff + self.coefficients[i];
            }
            if i < other.coefficients.len() {
                coeff = coeff + other.coefficients[i];
            }
            if coeff != zero {
                current_length = i + 1;
            }
            coefficients.push(coeff);
        }
        coefficients.truncate(current_length);
        CheckedPolynomial::<RingMod<T>>::new(coefficients)
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
        let poly = CheckedPolynomial::new(coefficients.clone());
        assert_eq!(poly.coefficients, coefficients);
    }

    #[test]
    fn test_polynomial_degree() {
        let coefficients = vec![
            RingMod::new(10, 13),
            RingMod::new(5, 13),
            RingMod::new(3, 13),
        ];
        let poly = CheckedPolynomial::new(coefficients);
        assert_eq!(poly.degree(), 2);

        let zero_poly: CheckedPolynomial<RingMod<i64>> = CheckedPolynomial::new(vec![]);
        assert_eq!(zero_poly.degree(), 0);
    }

    #[test]
    fn test_polynomial_is_zero() {
        let zero_poly: CheckedPolynomial<RingMod<i64>> = CheckedPolynomial::new(vec![]);
        assert!(zero_poly.is_zero());

        let non_zero_poly = CheckedPolynomial::new(vec![RingMod::new(10, 13)]);
        assert!(!non_zero_poly.is_zero());
    }

    #[test]
    fn test_polynomial_addition() {
        let modulus = 13;
        let poly1 = CheckedPolynomial::new(vec![
            RingMod::new(1, modulus),
            RingMod::new(2, modulus),
            RingMod::new(3, modulus),
        ]);
        let poly2 = CheckedPolynomial::new(vec![
            RingMod::new(4, modulus),
            RingMod::new(5, modulus),
            RingMod::new(6, modulus),
        ]);

        let expected_sum = CheckedPolynomial::new(vec![
            RingMod::new(5, modulus),
            RingMod::new(7, modulus),
            RingMod::new(9, modulus),
        ]);

        let result_sum = poly1 + &poly2;

        assert_eq!(result_sum, expected_sum);
    }

    #[test]
    fn test_polynomial_addition_distinct_degrees() {
        let modulus = 13;
        let poly1 =
            CheckedPolynomial::new(vec![RingMod::new(1, modulus), RingMod::new(2, modulus)]);
        let poly2 = CheckedPolynomial::new(vec![
            RingMod::new(4, modulus),
            RingMod::new(5, modulus),
            RingMod::new(6, modulus),
        ]);

        let expected_sum = CheckedPolynomial::new(vec![
            RingMod::new(5, modulus),
            RingMod::new(7, modulus),
            RingMod::new(6, modulus),
        ]);

        let result_sum = poly1 + &poly2;

        assert_eq!(result_sum, expected_sum);
    }

    #[test]
    fn test_polynomial_addition_smaller_degree() {
        let modulus = 13;
        let poly1 = CheckedPolynomial::new(vec![
            RingMod::new(1, modulus),
            RingMod::new(2, modulus),
            RingMod::new(7, modulus),
        ]);
        let poly2 = CheckedPolynomial::new(vec![
            RingMod::new(4, modulus),
            RingMod::new(5, modulus),
            RingMod::new(6, modulus),
        ]);

        let expected_sum =
            CheckedPolynomial::new(vec![RingMod::new(5, modulus), RingMod::new(7, modulus)]);

        let result_sum = poly1 + &poly2;

        assert_eq!(result_sum, expected_sum);
    }

    #[test]
    fn test_polynomial_addition_result_zero() {
        let modulus = 13;
        let poly1 = CheckedPolynomial::new(vec![
            RingMod::new(1, modulus),
            RingMod::new(2, modulus),
            RingMod::new(7, modulus),
        ]);
        let poly2 = CheckedPolynomial::new(vec![
            RingMod::new(12, modulus),
            RingMod::new(11, modulus),
            RingMod::new(-7, modulus),
        ]);

        let expected_sum = CheckedPolynomial::new(vec![]);

        let result_sum = poly1 + &poly2;

        assert_eq!(result_sum, expected_sum);
    }

    #[test]
    fn test_polynomial_addition_with_zero() {
        let modulus = 13;
        let poly1 = CheckedPolynomial::new(vec![
            RingMod::new(1, modulus),
            RingMod::new(2, modulus),
            RingMod::new(3, modulus),
        ]);
        let zero_poly: CheckedPolynomial<RingMod<i64>> = CheckedPolynomial::new(vec![]);

        let result_sum = poly1.clone() + &zero_poly;
        assert_eq!(result_sum, poly1);

        let result_sum = zero_poly + &poly1;
        assert_eq!(result_sum, poly1);
    }
}
