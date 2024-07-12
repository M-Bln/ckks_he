use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex, C64};
use crate::algebra::polynomial::Polynomial;
use crate::algebra::cyclotomic_ring::CyclotomicRing;

use std::f64::consts::PI;

#[derive(Debug)]
pub struct Encoder<T: BigInt> {
    dimension_exponent: u32,
    modulus: T,
    sigma_inverse_matrix: Vec<Vec<C64>>,
    sigma_matrix: Vec<Vec<C64>>,
}

impl<T: BigInt> Encoder<T> {
    pub fn new(dimension_exponent: u32, modulus: T) -> Self {
        let roots = C64::all_2_to_the_h_th_roots_of_unity(dimension_exponent);
        let vandermonde_matrix = Self::generate_vandermonde_matrix(&roots, dimension_exponent);
        let sigma_inverse_matrix = Self::generate_sigma_inverse_matrix(&vandermonde_matrix, &roots);
        let sigma_matrix =
            Self::generate_sigma_matrix(&vandermonde_matrix, &roots, dimension_exponent);

        Encoder {
            dimension_exponent,
            modulus,
            sigma_inverse_matrix,
            sigma_matrix,
        }
    }

    fn generate_vandermonde_matrix(roots: &[C64], dimension_exponent: u32) -> Vec<Vec<C64>> {
        let n = 2_usize.pow(dimension_exponent);
        let mut matrix = vec![vec![C64::new(0.0, 0.0); n]; n];

        for i in 0..n {
            for j in 0..n {
                matrix[i][j] = roots[(i * j) % n].clone();
            }
        }

        matrix
    }

    fn generate_sigma_inverse_matrix(
        vandermonde_matrix: &[Vec<C64>],
        roots: &[C64],
    ) -> Vec<Vec<C64>> {
        let n = roots.len();
        let mut matrix = vandermonde_matrix.to_vec();

        for j in 0..n {
            let factor = roots[j];
            for i in 0..n {
                matrix[i][j] = matrix[i][j] * factor;
            }
        }

        matrix
    }

    fn generate_sigma_matrix(
        vandermonde_matrix: &[Vec<C64>],
        roots: &[C64],
        dimension_exponent: u32,
    ) -> Vec<Vec<C64>> {
        let n = 2_usize.pow(dimension_exponent);
        let mut matrix = vandermonde_matrix.to_vec();

        for i in 0..n {
            let factor = roots[n - i -1];
            for j in 0..n {
                matrix[i][j] = matrix[i][j] * factor;
            }
        }

        matrix
    }

    pub fn apply_sigma_inverse(&self, poly: &Polynomial<C64>) -> Vec<C64> {
        let mut result = vec![C64::new(0.0, 0.0); poly.ref_coefficients().len()];
        for (i, row) in self.sigma_inverse_matrix.iter().enumerate() {
            result[i] = row.iter()
                .zip(poly.ref_coefficients())
                .fold(C64::new(0.0, 0.0), |sum, (a, b)| sum + a.clone() * b);
        }
        result
    }
}

// #[derive(Debug)]
// pub struct Encoder<T: BigInt> {
//     dimension_exponent: u32,
//     modulus: T,
//     sigma_inverse_matrix: Vec<Vec<C64>>,
//     sigma_matrix: Vec<Vec<C64>>,
// }

// impl<T: BigInt> Encoder<T> {
//     pub fn new(dimension_exponent: u32, modulus: T) -> Self {
//         let roots = C64::all_2_to_the_h_th_roots_of_unity(dimension_exponent);
//         let vandermonde_matrix = Self::generate_vandermonde_matrix(&roots, dimension_exponent);

//         Encoder {
//             dimension_exponent,
//             modulus,
//             vandermonde_matrix,
//         }
//     }

//     fn generate_vandermonde_matrix(roots: &[C64], dimension_exponent: u32) -> Vec<Vec<C64>> {
//         let n = 2_usize.pow(dimension_exponent);
//         let mut matrix = vec![vec![C64::new(0.0, 0.0); n]; n];

//         for i in 0..n {
//             for j in 0..n {
//                 matrix[i][j] = roots[(i * j) % n].clone();
//             }
//         }

//         matrix
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::big_int::BigInt;
    use crate::algebra::complex::C64;
    use crate::algebra::arithmetic::RingMod;
    use bnum::types::I256;

    
    #[test]
    fn test_sigma_inverse_matrix_evaluation() {
        let dimension_exponent = 3;
        let modulus = I256::new(13);
        let encoder = Encoder::new(dimension_exponent, modulus.clone());

        // Example polynomial in the cyclotomic ring
        let poly_ringmod = Polynomial::new( 
            vec![
                RingMod::new(I256::new(1), modulus.clone()),
                RingMod::new(I256::new(2), modulus.clone()),
                RingMod::new(I256::new(3), modulus.clone()),
		RingMod::new(I256::new(1), modulus.clone()),
                RingMod::new(I256::new(2), modulus.clone()),
                RingMod::new(I256::new(3), modulus.clone()),
		RingMod::new(I256::new(1), modulus.clone()),
                RingMod::new(I256::new(2), modulus.clone()),
            ],
      );

        // Convert polynomial to Polynomial<C64>
        let poly_c64 = poly_ringmod.to_c64();

        // Generate primitive 2^{h+1} roots of unity
        let roots = C64::primitive_2_to_the_h_th_roots_of_unity(dimension_exponent+1);// C64::primitive_2_to_the_h_plus_1_th_roots_of_unity(dimension_exponent as usize);

        // Evaluate polynomial at these roots
        let evaluations: Vec<C64> = roots.iter().map(|root| poly_c64.eval(root.clone())).collect();

        // Apply sigma_inverse_matrix to the polynomial coefficients
        let sigma_inverse_result = encoder.apply_sigma_inverse(&poly_c64);

        // Compare results
        for (i, eval) in evaluations.iter().enumerate() {
            assert!(
                (eval.real() - sigma_inverse_result[i].real()).abs() < 1e-10,
                "Real part mismatch at index {}: expected {}, got {}",
                i,
                eval.real(),
                sigma_inverse_result[i].real()
            );
            assert!(
                (eval.imaginary() - sigma_inverse_result[i].imaginary()).abs() < 1e-10,
                "Imaginary part mismatch at index {}: expected {}, got {}",
                i,
                eval.imaginary(),
                sigma_inverse_result[i].imaginary()
            );
        }
    }
    // #[test]
    // fn test_print_vandermonde_matrix() {
    //     let dimension_exponent = 3;
    //     let modulus: i64 = 13;
    //     let encoder = Encoder::new(dimension_exponent, modulus);

    //     for row in encoder.vandermonde_matrix.iter() {
    //         for element in row.iter() {
    //             print!("({:.4}, {:.4}) ", element.real(), element.imaginary());
    //         }
    //         println!();
    //     }
    // }
}

// pub struct Encoder<T: BigInt> {
//     modulus: T,
//     dimension_exponent: usize, // h such that  M = 2^h is the degree of the cyclotomic polynomial
//     primitive_roots: Vec<C64>, // The primitive 2^h-th roots of unity
//     matrix_sigma: Vec<Vec<C64>>, // The matrix of the linear map sigma
//     matrix_sigma_inver: Vec<Vec<C64>>, // The matrix of the inverse of sigma
// }
