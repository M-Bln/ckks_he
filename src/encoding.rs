use crate::algebra::arithmetic::{c64_to_ring_mod_256, RingMod};
use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex, C64};
use crate::algebra::cyclotomic_ring::CyclotomicRing;
use crate::algebra::linear_algebra::{apply_matrix, identity_matrix, multiply_matrices};
use crate::algebra::polynomial::Polynomial;

use bnum::types::I256;
use std::f64::consts::PI;

#[derive(Debug)]
pub struct Encoder<T: BigInt> {
    dimension_exponent: u32, // h such that  M = 2^h is the degree of the cyclotomic polynomial
    modulus: T, // The space of ciphertexts is Z/(modulus Z) [X] / (1 + X ^(2^(dimension_exponent)))
    sigma_inverse_matrix: Vec<Vec<C64>>, // Matrix to compute canonical embedding
    sigma_matrix: Vec<Vec<C64>>, // Matrix of sigma, the inverse of the canonical embedding
}

impl Encoder<I256> {
    pub fn encode(&self, plaintext: &[C64]) -> CyclotomicRing<RingMod<I256>> {
        let expected_length = 2_usize.pow(self.dimension_exponent - 1);
        assert_eq!(
            plaintext.len(),
            expected_length,
            "Plaintext length does not match half dimension of encoder"
        );

        // Step 1: Apply projection_inverse to plaintext
        let projection_inverse_result = self.projection_inverse(&plaintext);

        // Step 2: Apply sigma_inverse to the result of projection_inverse
        let sigma_inverse_result = self.sigma_inverse(&projection_inverse_result);

        // Step 3: Convert the polynomial to I256
        let integer_polynomial = sigma_inverse_result.to_i256();

        // Step 4: Reduce the polynomial modulo modulus
        let modular_polynomial = integer_polynomial.modulo(self.modulus.clone());
        println!("Modular Polynomial: {:?}", modular_polynomial);

        // Step 5: Reduce modulo the cyclotomic polynomial
        let cyclotomic_polynomial = modular_polynomial.to_cyclotomic(self.dimension_exponent);
        println!("Cyclotomic Polynomial: {:?}", cyclotomic_polynomial);

        cyclotomic_polynomial
    }

    pub fn decode(&self, message: CyclotomicRing<RingMod<I256>>) -> Vec<C64> {
        let expected_dimension = 2_usize.pow(self.dimension_exponent);
        assert_eq!(
            message.dimension, expected_dimension,
            "Message dimension does not match expected dimension"
        );

        // Step 1: Convert to complex coefficients
        let complex_polynomial = message.to_c64();
        println!("Complex Polynomial: {:?}", complex_polynomial);

        // Step 2: Apply sigma
        let canonical_embedding = self.sigma(&complex_polynomial);
        println!("Canonical Embedding: {:?}", canonical_embedding);

        // Step 3: Project to message space
        let projection_result = self.projection(&canonical_embedding);
        println!("Projection Result: {:?}", projection_result);

        projection_result
    }
}

// impl Encoder<I256> {
//     pub fn encode(&self, plaintext: &[C64]) -> CyclotomicRing<RingMod<I256>> {
//         let expected_length = 2_usize.pow(self.dimension_exponent - 1);
//         assert_eq!(
//             plaintext.len(),
//             expected_length,
//             "Plaintext length does not match half dimension of encoder"
//         );

//         let integer_polynomial = self.sigma_inverse(&self.projection_inverse(&plaintext)).to_i256();
//         let modular_polynomial = integer_polynomial.modulo(self.modulus);
//         modular_polynomial.to_cyclotomic(self.dimension_exponent)
//     }

//     pub fn decode(&self, message: CyclotomicRing<RingMod<I256>>) -> Vec<C64> {
//         let expected_dimension = 2_usize.pow(self.dimension_exponent);
//         assert_eq!(
//             message.dimension,
//             expected_dimension,
//             "Plaintext length does not match half dimension of encoder"
//         );
// 	let complex_polynomial = message.to_c64();
// 	let canonical_embedding = self.sigma(&complex_polynomial);
// 	self.projection(&canonical_embedding)
//     }
// }

impl<T: BigInt> Encoder<T> {
    /// Generate an encoder for the space of ciphertexts $\ZZ/(modulus \ZZ) [X] / (1 + X ^{2^{dimension_exponent}})$
    /// Decoding $decode: \ZZ/(modulus \ZZ) [X] / (1 + X ^{2^{dimension_exponent}}) \to \CC^{2^{dimension_exponent-1}}$
    /// is the composition $decode = \pi\inv \circ \sigma\inv$ with $\sigma$ the canonical embedding. In the encoder
    /// $\sigma$ and $\sigma\inv$ are stored as matrices of roots of unity.
    // TODO use FFT rather than matrices (no emergency, probably not the bottleneck)
    pub fn new(dimension_exponent: u32, modulus: T) -> Self {
        // We take a primitive 2^{h+1}-th root of unity \zeta = \exp(2 i \pi/ 2^{h+1})
        // Generate (1, \zeta, \zeta^2, \zeta^3, ..., \zeta^{2^{h+1}-1})
        assert!(
            dimension_exponent >= 1,
            "Dimension exponent must be greater than 1"
        );
        let mut roots = C64::all_2_to_the_h_th_roots_of_unity(dimension_exponent + 1);

        // Generate (1, \zeta^2, \zeta^4, \zeta^6, ..., \zeta^{2^{h+2}-2})
        let roots_vandermonde = even_index_elements(&roots);

        // Keep (1, \zeta, \zeta^2, ..., \zeta^{2^h-1})
        roots.truncate(2_usize.pow(dimension_exponent));

        // Compute the Vandermonde matrix with coefficient at line i and column j: (\zeta^{2(i-1)(j-1)})
        let vandermonde_matrix =
            Self::generate_vandermonde_matrix(&roots_vandermonde, dimension_exponent);
        let sigma_matrix = Self::generate_sigma_matrix(&vandermonde_matrix, &roots);

        let inverse_roots: Vec<C64> = roots.iter().map(|z| z.conjugate()).collect();
        let inverse_roots_vandermonde: Vec<C64> =
            roots_vandermonde.iter().map(|z| z.conjugate()).collect();

        let inverse_vandermonde_matrix =
            Self::generate_vandermonde_matrix(&inverse_roots_vandermonde, dimension_exponent);

        let sigma_inverse_matrix = Self::generate_sigma_inverse_matrix(
            &inverse_vandermonde_matrix,
            &inverse_roots,
            dimension_exponent,
        );

        Encoder {
            dimension_exponent,
            modulus,
            sigma_inverse_matrix,
            sigma_matrix,
        }
    }

    fn projection(&self, z: &[C64]) -> Vec<C64> {
        assert_eq!(
            z.len(),
            2_usize.pow(self.dimension_exponent),
            "Length of the argument of projection not equal to dimension"
        );
        even_index_elements(&z)
    }

    fn projection_inverse(&self, z: &[C64]) -> Vec<C64> {
        let n = 2_usize.pow(self.dimension_exponent - 1);
        assert_eq!(
            z.len(),
            n,
            "Unexpected length of the argument of projection inverse"
        );

        let mut result = Vec::with_capacity(2 * n);

        for i in 0..n {
            result.push(z[i].clone());
            result.push(z[n - i - 1].conjugate());
        }

        result
    }

    fn sigma(&self, p: &Polynomial<C64>) -> Vec<C64> {
        apply_matrix(&self.sigma_matrix, p.ref_coefficients())
    }

    fn sigma_inverse(&self, z: &[C64]) -> Polynomial<C64> {
        Polynomial::<C64>::new(apply_matrix(&self.sigma_inverse_matrix, z))
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

    fn generate_sigma_matrix(vandermonde_matrix: &[Vec<C64>], roots: &[C64]) -> Vec<Vec<C64>> {
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

    fn generate_sigma_inverse_matrix(
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
                matrix[i][j] = matrix[i][j] * factor / n_complex;
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

fn even_index_elements<T: Clone>(vec: &[T]) -> Vec<T> {
    vec.iter()
        .enumerate()
        .filter(|&(i, _)| i % 2 == 0)
        .map(|(_, x)| x.clone())
        .collect()
}

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

// fn multiply_matrices(a: &[Vec<C64>], b: &[Vec<C64>]) -> Vec<Vec<C64>> {
//     let n = a.len();
//     let mut result = vec![vec![C64::new(0.0, 0.0); n]; n];

//     for i in 0..n {
//         for j in 0..n {
//             for k in 0..n {
//                 result[i][j] = result[i][j] + a[i][k] * b[k][j];
//             }
//         }
//     }

//     result
// }

// fn identity_matrix(n: usize) -> Vec<Vec<C64>> {
//     let mut matrix = vec![vec![C64::new(0.0, 0.0); n]; n];
//     for i in 0..n {
//         matrix[i][i] = C64::new(1.0, 0.0);
//     }
//     matrix
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::arithmetic::RingMod;
    use crate::algebra::big_int::BigInt;
    use crate::algebra::complex::C64;
    use bnum::types::I256;

    // #[test]
    // fn test_encode_decode() {
    //     let dimension_exponent = 4;
    //     let modulus = I256::new(1300000);
    //     let encoder = Encoder::new(dimension_exponent, modulus.clone());

    //     // Generate a vector of large complex numbers as input
    //     let plaintext: Vec<C64> = (0..2_usize.pow(dimension_exponent - 1))
    //         .map(|i| C64::new((i*100) as f64 + 10000.0, (i * i * 100)  as f64 + 76522.1))
    //         .collect();
    // 	println!("plaintext: {:?}", plaintext);

    //     // Encode the plaintext
    //     let encoded = encoder.encode(&plaintext);
    // 	println!("encoded: {:?}", encoded);

    //     // Decode the encoded message
    //     let decoded = encoder.decode(encoded);
    // 	println!("decoded: {:?}", decoded);

    //     // Verify that the decoded message matches the original plaintext
    //     let tol = 100.0;
    //     for (a, b) in plaintext.iter().zip(decoded.iter()) {
    //         assert!(
    //             approx_eq(a.real(), b.real(), tol),
    //             "Real part mismatch: expected {}, got {}",
    //             a.real(),
    //             b.real()
    //         );
    //         assert!(
    //             approx_eq(a.imaginary(), b.imaginary(), tol),
    //             "Imaginary part mismatch: expected {}, got {}",
    //             a.imaginary(),
    //             b.imaginary()
    //         );
    //     }
    // }

    #[test]
    fn test_real_coeff() {
        let dimension_exponent = 4;
        let modulus = I256::new(13);
        let encoder = Encoder::new(dimension_exponent, modulus.clone());

        // Generate a large vector of complex numbers as input
        let plaintext: Vec<C64> = (0..2_usize.pow(dimension_exponent - 1))
            .map(|i| C64::new(i as f64 * 1e5, (i * i) as f64 * 1e5))
            .collect();

        println!("Original Plaintext: {:?}", plaintext);

        // Apply projection_inverse to plaintext
        let projection_inverse_result = encoder.projection_inverse(&plaintext);
        println!("Projection Inverse Result: {:?}", projection_inverse_result);

        // Apply sigma_inverse to the result of projection_inverse
        let sigma_inverse_result = encoder.sigma_inverse(&projection_inverse_result);
        println!("Sigma Inverse Result: {:?}", sigma_inverse_result);

        // Check that the coefficients of sigma_inverse_result are real
        let tol = 110.0;
        for coeff in sigma_inverse_result.coefficients().iter() {
            assert!(
                approx_eq(coeff.imaginary(), 0.0, tol),
                "Imaginary part is not zero: got {}",
                coeff.imaginary()
            );
        }
    }

    #[test]
    fn test_sigma_y() {
        let dimension_exponent = 4;
        let modulus = I256::new(1298736);
        let encoder = Encoder::new(dimension_exponent, modulus.clone());

        // Create a polynomial with non-imaginary coefficients
        let coefficients: Vec<C64> = (0..2_usize.pow(dimension_exponent))
            .map(|i| C64::new(i as f64, 0.0))
            .collect();
        let polynomial = Polynomial::new(coefficients);

        println!("Original Polynomial: {:?}", polynomial);

        // Apply sigma to the polynomial
        let sigma_result = encoder.sigma(&polynomial);
        println!("Sigma Result: {:?}", sigma_result);

        // Check and print each coefficient of the sigma result
        for (i, coeff) in sigma_result.iter().enumerate() {
            println!(
                "Coefficient {}: ({:.4}, {:.4})",
                i,
                coeff.real(),
                coeff.imaginary()
            );
        }
    }

    #[test]
    fn test_encode_decode() {
        let dimension_exponent = 4;
        let modulus = I256::new(1000000000000000);
        let encoder = Encoder::new(dimension_exponent, modulus.clone());

        // Generate a large vector of complex numbers as input
        let plaintext: Vec<C64> = (0..2_usize.pow(dimension_exponent - 1))
            .map(|i| C64::new(i as f64 * 1e5, (i * i) as f64 * 1e5))
            .collect();

        println!("Original Plaintext: {:?}", plaintext);

        // Encode the plaintext
        let encoded = encoder.encode(&plaintext);
        println!("Encoded Message: {:?}", encoded);

        // Decode the encoded message
        let decoded = encoder.decode(encoded);
        println!("Decoded Plaintext: {:?}", decoded);

        // Verify that the decoded message matches the original plaintext
        let tol = 2.0;
        for (a, b) in plaintext.iter().zip(decoded.iter()) {
            assert!(
                approx_eq(a.real(), b.real(), tol),
                "Real part mismatch: expected {}, got {}",
                a.real(),
                b.real()
            );
            assert!(
                approx_eq(a.imaginary(), b.imaginary(), tol),
                "Imaginary part mismatch: expected {}, got {}",
                a.imaginary(),
                b.imaginary()
            );
        }
    }

    #[test]
    fn test_sigma_inverse_matrix_evaluation() {
        // We check that sigma_inverse indeed compute canonical embedding, i.e.
        // the evaluation of the polynomial at primitive roots of unity
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
        let roots = C64::primitive_2_to_the_h_th_roots_of_unity(dimension_exponent + 1);

        println!("root: {:?}", roots);
        // Evaluate polynomial at these roots
        let evaluations: Vec<C64> = roots
            .iter()
            .map(|root| poly_c64.eval(root.clone()))
            .collect();

        // Apply sigma_inverse_matrix to the polynomial coefficients
        //        let sigma_inverse_result = encoder.apply_sigma_inverse(&poly_c64);
        let sigma_inverse_result_coeffs = encoder.sigma(&poly_c64);
        // Compare results
        for (i, eval) in evaluations.iter().enumerate() {
            assert!(
                (eval.real() - sigma_inverse_result_coeffs[i].real()).abs() < 1e-10,
                "Real part mismatch at index {}: expected {}, got {}",
                i,
                eval.real(),
                sigma_inverse_result_coeffs[i].real()
            );
            assert!(
                (eval.imaginary() - sigma_inverse_result_coeffs[i].imaginary()).abs() < 1e-10,
                "Imaginary part mismatch at index {}: expected {}, got {}",
                i,
                eval.imaginary(),
                sigma_inverse_result_coeffs[i].imaginary()
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

        for row in product.iter() {
            for element in row.iter() {
                print!("({:.4}, {:.4}) ", element.real(), element.imaginary());
            }
            println!();
        }

        // Generate the identity matrix
        let identity = identity_matrix(product.len());

        // Check if the product is equal to the identity matrix
        for i in 0..product.len() {
            for j in 0..product.len() {
                assert!(
                    (product[i][j].real() - identity[i][j].real()).abs() < 1e-10,
                    "Real part mismatch at ({}, {}): expected {}, got {}",
                    i,
                    j,
                    identity[i][j].real(),
                    product[i][j].real()
                );
                assert!(
                    (product[i][j].imaginary() - identity[i][j].imaginary()).abs() < 1e-10,
                    "Imaginary part mismatch at ({}, {}): expected {}, got {}",
                    i,
                    j,
                    identity[i][j].imaginary(),
                    product[i][j].imaginary()
                );
            }
        }
    }

    #[test]
    fn test_sigma_and_sigma_inverse() {
        let dimension_exponent = 3;
        let modulus = I256::new(13);
        let encoder = Encoder::new(dimension_exponent, modulus.clone());

        // Generate a random polynomial of type Polynomial<C64>
        let coefficients = vec![
            C64::new(1.0, 0.5),
            C64::new(-1.5, 2.0),
            C64::new(3.0, -3.5),
            C64::new(0.0, 1.0),
            C64::new(-2.5, -0.5),
            C64::new(4.5, -2.0),
            C64::new(-3.5, 3.5),
            C64::new(2.0, -1.0),
        ];
        let polynomial = Polynomial::new(coefficients);

        // Apply sigma to the polynomial
        let sigma_result = encoder.sigma(&polynomial);

        // Apply sigma_inverse to the resulting vector
        let sigma_inverse_result = encoder.sigma_inverse(&sigma_result);

        // Verify that the original polynomial and the resulting polynomial are approximately equal
        let tol = 1e-10;
        for (a, b) in polynomial
            .coefficients()
            .iter()
            .zip(sigma_inverse_result.coefficients().iter())
        {
            assert!(
                approx_eq(a.real(), b.real(), tol),
                "Real part mismatch: expected {}, got {}",
                a.real(),
                b.real()
            );
            assert!(
                approx_eq(a.imaginary(), b.imaginary(), tol),
                "Imaginary part mismatch: expected {}, got {}",
                a.imaginary(),
                b.imaginary()
            );
        }
    }

    #[test]
    fn test_projection_inverse() {
        let dimension_exponent = 3;
        let modulus = I256::new(13);
        let encoder = Encoder::new(dimension_exponent, modulus.clone());

        // Example input vector z
        let z = vec![
            C64::new(1.0, 1.0),
            C64::new(2.0, 2.0),
            C64::new(3.0, 3.0),
            C64::new(4.0, 4.0),
        ];

        // Expected output
        let expected = vec![
            C64::new(1.0, 1.0),
            C64::new(4.0, -4.0),
            C64::new(2.0, 2.0),
            C64::new(3.0, -3.0),
            C64::new(3.0, 3.0),
            C64::new(2.0, -2.0),
            C64::new(4.0, 4.0),
            C64::new(1.0, -1.0),
        ];

        // Compute the result
        let result = encoder.projection_inverse(&z);

        // Check if the result matches the expected output
        let tol = 1e-10;
        for (r, e) in result.iter().zip(expected.iter()) {
            assert!(
                approx_eq(r.real(), e.real(), tol),
                "Real part mismatch: expected {}, got {}",
                e.real(),
                r.real()
            );
            assert!(
                approx_eq(r.imaginary(), e.imaginary(), tol),
                "Imaginary part mismatch: expected {}, got {}",
                e.imaginary(),
                r.imaginary()
            );
        }
    }

    #[test]
    fn test_projection_and_inverse() {
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
        let roots = C64::primitive_2_to_the_h_th_roots_of_unity(dimension_exponent + 1);

        println!("root: {:?}", roots);
        // Evaluate polynomial at these roots
        let z: Vec<C64> = roots
            .iter()
            .map(|root| poly_c64.eval(root.clone()))
            .collect();

        // Apply projection followed by projection_inverse
        let projected = encoder.projection(&z);
        let projected_inverse = encoder.projection_inverse(&projected);

        // Verify projection_inverse(projection(z)) == z
        let tol = 1e-10;
        for (a, b) in z.iter().zip(projected_inverse.iter()) {
            assert!(
                approx_eq(a.real(), b.real(), tol),
                "Real part mismatch: expected {}, got {}",
                a.real(),
                b.real()
            );
            assert!(
                approx_eq(a.imaginary(), b.imaginary(), tol),
                "Imaginary part mismatch: expected {}, got {}",
                a.imaginary(),
                b.imaginary()
            );
        }

        // Apply projection_inverse followed by projection
        let projected_inverse_again = encoder.projection_inverse(&projected);
        let projected_again = encoder.projection(&projected_inverse_again);

        // Verify projection(projection_inverse(z)) == z
        for (a, b) in projected.iter().zip(projected_again.iter()) {
            assert!(
                approx_eq(a.real(), b.real(), tol),
                "Real part mismatch: expected {}, got {}",
                a.real(),
                b.real()
            );
            assert!(
                approx_eq(a.imaginary(), b.imaginary(), tol),
                "Imaginary part mismatch: expected {}, got {}",
                a.imaginary(),
                b.imaginary()
            );
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
