use crate::algebra::arithmetic::RingMod;
use crate::algebra::big_int::Zero;
use crate::algebra::complex::{Complex, C64};
use crate::algebra::polynomial::Polynomial;
use std::ops::{Add, Mul, Sub};

use bnum::types::I256;

#[derive(Clone, Debug, PartialEq)]
pub struct CyclotomicRing<T> {
    pub polynomial: Polynomial<T>,
    pub dimension: usize,
}

impl<T> CyclotomicRing<T>
where
    T: Add<Output = T>
        + for<'a> Add<&'a T, Output = T>
        + Sub<Output = T>
        + for<'a> Sub<&'a T, Output = T>
        + Clone
        + Zero,
{
    pub fn new(mut coefficients: Vec<T>, dimension: usize) -> Self {
        assert!(
            !coefficients.is_empty(),
            "Don't create new cyclotomic ring with empty vector"
        ); // This is because we need to know the modulus
        extend_by_zero(&mut coefficients, dimension);
        let polynomial = Polynomial::new(coefficients);
        let mut cyclo = CyclotomicRing {
            polynomial,
            dimension,
        };
        cyclo.reduce();
        cyclo
    }

    fn reduce(&mut self) {
        let n = self.dimension;
        let coefficients = self.polynomial.ref_coefficients().to_vec();

        let mut new_coefficients = coefficients[..n].to_vec();
        for (index, coeff) in coefficients[n..].iter().enumerate() {
            // beware! index == 0 corresponds to coefficient of degree n
            let shifted_index = index % n;
            let shifted_coeff = new_coefficients[shifted_index].clone();
            if ((index + n) / n) % 2 == 0 {
                new_coefficients[shifted_index] = shifted_coeff + coeff.clone();
            } else {
                new_coefficients[shifted_index] = shifted_coeff - coeff.clone();
            }
        }

        self.polynomial
            .mut_coefficients()
            .splice(.., new_coefficients);
        self.polynomial.mut_coefficients().truncate(n);
    }
}

fn extend_by_zero<T: Zero + Clone>(coefficients: &mut Vec<T>, minimal_length: usize) {
    assert!(
        !coefficients.is_empty(),
        "cannot extend empty vector by zero, we might need to know the modulus"
    );
    if coefficients.len() >= minimal_length {
        return;
    }
    let extension_length = minimal_length - coefficients.len();
    let mut extension = vec![coefficients[0].zero(); extension_length];
    coefficients.append(&mut extension)
}

impl<'a, T> Add<&'a CyclotomicRing<T>> for CyclotomicRing<T>
where
    T: Add<Output = T>
        + for<'b> Add<&'b T, Output = T>
        + Sub<Output = T>
        + for<'b> Sub<&'b T, Output = T>
        + Clone
        + Zero,
{
    type Output = Self;

    fn add(self, other: &'a CyclotomicRing<T>) -> Self::Output {
        assert_eq!(
            self.dimension, other.dimension,
            "Dimensions must be equal for addition"
        );
        let new_poly = self.polynomial + &other.polynomial;
        CyclotomicRing::new(new_poly.coefficients(), self.dimension)
    }
}

impl<'a, T> Sub<&'a CyclotomicRing<T>> for CyclotomicRing<T>
where
    T: Add<Output = T>
        + for<'b> Add<&'b T, Output = T>
        + Sub<Output = T>
        + for<'b> Sub<&'b T, Output = T>
        + Clone
        + Zero,
{
    type Output = Self;

    fn sub(self, other: &'a CyclotomicRing<T>) -> Self::Output {
        assert_eq!(
            self.dimension, other.dimension,
            "Dimensions must be equal for subtraction"
        );
        let new_poly = self.polynomial - &other.polynomial;
        CyclotomicRing::new(new_poly.coefficients(), self.dimension)
    }
}

impl<'a, T> Mul<&'a CyclotomicRing<T>> for CyclotomicRing<T>
where
    T: Add<Output = T>
        + for<'b> Add<&'b T, Output = T>
        + Sub<Output = T>
        + for<'b> Sub<&'b T, Output = T>
        + Clone
        + Zero
        + Mul<Output = T>,
{
    type Output = Self;

    fn mul(self, other: &'a CyclotomicRing<T>) -> Self::Output {
        assert_eq!(
            self.dimension, other.dimension,
            "Dimensions must be equal for multiplication"
        );
        let new_poly = self.polynomial * &other.polynomial;
        CyclotomicRing::new(new_poly.coefficients(), self.dimension)
    }
}

// impl CyclotomicRing<RingMod<I256>> {
//     pub fn to_c64(&self) -> Polynomial<C64> {
//         let coefficients = self
//             .polynomial
//             .ref_coefficients()
//             .iter()
//             .map(|coeff| coeff.to_c64())
//             .collect();
//         Polynomial::new(coefficients)
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::arithmetic::RingMod;
    use bnum::types::I256;

    #[test]
    fn test_extend_with_zeros() {
        let mut vec = vec![1, 2, 3];
        extend_by_zero(&mut vec, 5);
        assert_eq!(vec, vec![1, 2, 3, 0, 0]);
    }

    #[test]
    fn test_cyclotomic_ring_creation() {
        let coefficients = vec![1, 2, 3];
        let dimension = 5;
        let cyclo_ring = CyclotomicRing::new(coefficients, dimension);
        assert_eq!(cyclo_ring.polynomial.coefficients(), vec![1, 2, 3, 0, 0]);
    }

    #[test]
    fn test_reduce() {
        let dimension = 2;
        let mut cyclo_ring = CyclotomicRing::new(vec![1, 2, -2, 3, 4, 7, 8], dimension);
        cyclo_ring.reduce();
        let expected_coefficients = vec![-1, 6];
        assert_eq!(cyclo_ring.polynomial.coefficients(), expected_coefficients);
    }

    #[test]
    fn test_cyclotomic_ring_addition() {
        let dimension = 5;
        let poly1 = CyclotomicRing::new(vec![1, 2, 3], dimension);
        let poly2 = CyclotomicRing::new(vec![4, 5, 6], dimension);

        let expected_sum = CyclotomicRing::new(vec![5, 7, 9], dimension);

        let result_sum = poly1 + &poly2;
        assert_eq!(
            result_sum.polynomial.coefficients(),
            expected_sum.polynomial.coefficients()
        );
    }

    #[test]
    fn test_cyclotomic_ring_subtraction() {
        let dimension = 5;
        let poly1 = CyclotomicRing::new(vec![1, 2, 3], dimension);
        let poly2 = CyclotomicRing::new(vec![4, 5, 6], dimension);

        let expected_difference = CyclotomicRing::new(vec![-3, -3, -3, 0, 0], dimension);

        let result_difference = poly1 - &poly2;
        assert_eq!(result_difference.polynomial, expected_difference.polynomial);
    }

    #[test]
    fn test_cyclotomic_ring_multiplication() {
        let dimension = 2;
        let poly1 = CyclotomicRing::new(vec![1, 2, 3], dimension);
        let poly2 = CyclotomicRing::new(vec![4, 5, 6], dimension);

        let expected_product = CyclotomicRing::new(vec![-6, -14], dimension);

        let result_product = poly1 * &poly2;
        assert_eq!(result_product.polynomial, expected_product.polynomial);
    }

    #[test]
    fn test_cyclotomic_ring_subtraction_256() {
        let dimension = 5;
        let modulus = I256::from(13);
        let poly1 = CyclotomicRing::new(
            vec![
                RingMod::new(I256::from(1), modulus),
                RingMod::new(I256::from(2), modulus),
                RingMod::new(I256::from(3), modulus),
            ],
            dimension,
        );
        let poly2 = CyclotomicRing::new(
            vec![
                RingMod::new(I256::from(4), modulus),
                RingMod::new(I256::from(5), modulus),
                RingMod::new(I256::from(6), modulus),
            ],
            dimension,
        );

        let expected_difference = CyclotomicRing::new(
            vec![
                RingMod::new(I256::from(-3), modulus),
                RingMod::new(I256::from(-3), modulus),
                RingMod::new(I256::from(-3), modulus),
            ],
            dimension,
        );

        let result_difference = poly1 - &poly2;
        assert_eq!(
            result_difference.polynomial.coefficients(),
            expected_difference.polynomial.coefficients()
        );
    }

    #[test]
    fn test_cyclotomic_ring_multiplication_256() {
        let dimension = 2;
        let modulus = I256::from(13);
        let poly1 = CyclotomicRing::new(
            vec![
                RingMod::new(I256::from(1), modulus),
                RingMod::new(I256::from(2), modulus),
                RingMod::new(I256::from(3), modulus),
            ],
            dimension,
        );
        let poly2 = CyclotomicRing::new(
            vec![
                RingMod::new(I256::from(4), modulus),
                RingMod::new(I256::from(5), modulus),
                RingMod::new(I256::from(6), modulus),
            ],
            dimension,
        );

        let expected_product = CyclotomicRing::new(
            vec![
                RingMod::new(I256::from(-6), modulus),
                RingMod::new(I256::from(12), modulus),
            ],
            dimension,
        );

        let result_product = poly1 * &poly2;
        assert_eq!(
            result_product.polynomial.coefficients(),
            expected_product.polynomial.coefficients()
        );
    }
}
