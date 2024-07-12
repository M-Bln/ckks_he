use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex, C64};
use crate::algebra::cyclotomic_ring::CyclotomicRing;
use crate::algebra::polynomial::Polynomial;

use std::f64::consts::PI;

#[derive(Debug)]
pub struct Encoder<T: BigInt> {
    dimension_exponent: u32, // h such that  M = 2^h is the degree of the cyclotomic polynomial
    modulus: T,
    sigma_inverse_matrix: Vec<Vec<C64>>,
    sigma_matrix: Vec<Vec<C64>>,
}

impl<T: BigInt> Encoder<T> {
    pub fn new(dimension_exponent: u32, modulus: T) -> Self {
        let mut roots = C64::all_2_to_the_h_th_roots_of_unity(dimension_exponent + 1);
        roots.truncate(2_usize.pow(dimension_exponent));
        let roots_vandermonde = C64::all_2_to_the_h_th_roots_of_unity(dimension_exponent);
        let vandermonde_matrix =
            Self::generate_vandermonde_matrix(&roots_vandermonde, dimension_exponent);
        let sigma_inverse_matrix = Self::generate_sigma_inverse_matrix(&vandermonde_matrix, &roots);


	let mut inverse_roots = C64::inverse_all_2_to_the_h_th_roots_of_unity(dimension_exponent + 1);
        inverse_roots.truncate(2_usize.pow(dimension_exponent));
        let inverse_roots_vandermonde = C64::inverse_all_2_to_the_h_th_roots_of_unity(dimension_exponent);
        let inverse_vandermonde_matrix =
            Self::generate_vandermonde_matrix(&inverse_roots_vandermonde, dimension_exponent);
	
        let sigma_matrix =
            Self::generate_sigma_matrix(&inverse_vandermonde_matrix, &inverse_roots, dimension_exponent);

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
	let n_complex = C64::new(n as f64, 0.0);
        let mut matrix = vandermonde_matrix.to_vec();

        for i in 0..n {
            let factor = roots[i];
            for j in 0..n {
                matrix[i][j] = matrix[i][j] * factor / n_complex ;
            }
        }

        matrix
    }

    pub fn apply_sigma_inverse(&self, poly: &Polynomial<C64>) -> Vec<C64> {
        let mut result = vec![C64::new(0.0, 0.0); poly.ref_coefficients().len()];
        for (i, row) in self.sigma_inverse_matrix.iter().enumerate() {
            result[i] = row
                .iter()
                .zip(poly.ref_coefficients())
                .fold(C64::new(0.0, 0.0), |sum, (a, b)| sum + a.clone() * b);
        }
        result
    }
}

fn multiply_matrices(a: &[Vec<C64>], b: &[Vec<C64>]) -> Vec<Vec<C64>> {
    let n = a.len();
    let mut result = vec![vec![C64::new(0.0, 0.0); n]; n];
    
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                result[i][j] = result[i][j] + a[i][k] * b[k][j];
            }
        }
    }

    result
}

fn identity_matrix(n: usize) -> Vec<Vec<C64>> {
    let mut matrix = vec![vec![C64::new(0.0, 0.0); n]; n];
    for i in 0..n {
        matrix[i][i] = C64::new(1.0, 0.0);
    }
    matrix
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
    use crate::algebra::arithmetic::RingMod;
    use crate::algebra::big_int::BigInt;
    use crate::algebra::complex::C64;
    use bnum::types::I256;

    #[test]
    fn test_sigma_inverse_matrix_evaluation() {
        let dimension_exponent = 3;
        let modulus = I256::new(13);
        let encoder = Encoder::new(dimension_exponent, modulus.clone());

        // Example polynomial in the cyclotomic ring
        let poly_ringmod = CyclotomicRing::new(
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
            2_usize.pow(dimension_exponent),
        );

        // Convert polynomial to Polynomial<C64>
        let poly_c64 = poly_ringmod.to_c64();
        println!("poly_c64: {:?}", poly_c64);

        // Generate primitive 2^{h+1} roots of unity
        let roots = C64::primitive_2_to_the_h_th_roots_of_unity(dimension_exponent + 1); // C64::primitive_2_to_the_h_plus_1_th_roots_of_unity(dimension_exponent as usize);
                                                                                         //let roots = C64::all_2_to_the_h_th_roots_of_unity(dimension_exponent);

        println!("root: {:?}", roots);
        // Evaluate polynomial at these roots
        let evaluations: Vec<C64> = roots
            .iter()
            .map(|root| poly_c64.eval(root.clone()))
            .collect();

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

    #[test]
    fn test_print_matrix() {
        let dimension_exponent = 3;
        let modulus: i64 = 13;
        let encoder = Encoder::new(dimension_exponent, modulus);

        for row in encoder.sigma_inverse_matrix.iter() {
            for element in row.iter() {
                print!("({:.4}, {:.4}) ", element.real(), element.imaginary());
            }
            println!();
        }
    }

    #[test]
    fn test_sigma_matrices_are_inverses() {
        let dimension_exponent = 3;
        let modulus = I256::new(13);
        let encoder = Encoder::new(dimension_exponent, modulus.clone());

        // Multiply sigma_matrix by sigma_inverse_matrix
        let product = multiply_matrices(&encoder.sigma_matrix, &encoder.sigma_inverse_matrix);

        // Generate the identity matrix
        let identity = identity_matrix(product.len());

        // Check if the product is equal to the identity matrix
        for i in 0..product.len() {
            for j in 0..product.len() {
                assert!(
                    (product[i][j].real() - identity[i][j].real()).abs() < 1e-10,
                    "Real part mismatch at ({}, {}): expected {}, got {}",
                    i, j, identity[i][j].real(), product[i][j].real()
                );
                assert!(
                    (product[i][j].imaginary() - identity[i][j].imaginary()).abs() < 1e-10,
                    "Imaginary part mismatch at ({}, {}): expected {}, got {}",
                    i, j, identity[i][j].imaginary(), product[i][j].imaginary()
                );
            }
        }
    }
}

// pub struct Encoder<T: BigInt> {
//     modulus: T,
//     dimension_exponent: usize, // h such that  M = 2^h is the degree of the cyclotomic polynomial
//     primitive_roots: Vec<C64>, // The primitive 2^h-th roots of unity
//     matrix_sigma: Vec<Vec<C64>>, // The matrix of the linear map sigma
//     matrix_sigma_inver: Vec<Vec<C64>>, // The matrix of the inverse of sigma
// }
