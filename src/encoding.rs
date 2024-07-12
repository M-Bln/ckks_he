use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex, C64};
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
            let factor = roots[n - i];
            for j in 0..n {
                matrix[i][j] = matrix[i][j] * factor;
            }
        }

        matrix
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
