use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex, C64};
use crate::ciphertext::Ciphertext;
use crate::encoding::Encoder;
use crate::keys::key_generator::KeyGenerationParameters;
use crate::keys::public_key::PublicKey;
use crate::keys::secret_key::SecretKey;

type Plaintext = Vec<C64>;

pub struct ClientKey<T: BigInt> {
    pub encoder: Encoder<T>,
    pub secret_key: SecretKey<T>,
    pub public_key: PublicKey<T>,
    pub plaintext_dimension: usize,
}

#[derive(Debug)]
pub enum EncryptionError {
    UnexpectedPlaintextDimension,
    ClearLargerThanUpperBound,
}

/// Encapsulate the data necessary to the client side i.e. the secret key tue public key and an encoder
impl<T: BigInt> ClientKey<T> {
    pub fn new(
        parameters: KeyGenerationParameters<T>,
        secret_key: SecretKey<T>,
        public_key: PublicKey<T>,
    ) -> Self {
        assert!(parameters.dimension_exponent > 0);
        let plaintext_dimension = 1 << (parameters.dimension_exponent - 1); // 2^(dimension_exponent-1)
        let modulus = parameters.q.fast_exp(parameters.level_max) * &parameters.mul_scaling;
        let encoder = Encoder::new(
            parameters.dimension_exponent,
            modulus,
            parameters.q.to_float(),
        );
        ClientKey {
            encoder,
            secret_key,
            public_key,
            plaintext_dimension,
        }
    }

    pub fn encrypt(
        &mut self,
        clear: &Plaintext,
        upper_bound: f64,
    ) -> Result<Ciphertext<T>, EncryptionError> {
        if clear.len() != self.plaintext_dimension {
            return Err(EncryptionError::UnexpectedPlaintextDimension);
        }
        if !bounded_by(&clear, upper_bound) {
            return Err(EncryptionError::ClearLargerThanUpperBound);
        }
        let encoded = self.encoder.encode(&clear);
        let encrypted = self
            .public_key
            .encrypt(&encoded, upper_bound * self.encoder.scaling_factor);
        Ok(encrypted)
    }

    pub fn decrypt(&mut self, ct: &Ciphertext<T>) -> Plaintext {
        let encoded = self.secret_key.decrypt(ct);
        self.encoder.decode(&encoded)
    }

    pub fn rescaled_error(&self, ct: &Ciphertext<T>) -> f64 {
        ct.upper_bound_error / self.encoder.scaling_factor
    }

    pub fn rescaled_upperbound_message(&self, ct: &Ciphertext<T>) -> f64 {
        ct.upper_bound_message / self.encoder.scaling_factor
    }
}

pub fn bounded_by(z: &[C64], bound: f64) -> bool {
    z.iter().all(|&c| c.magnitude() <= bound)
}

pub fn to_plaintext(real: &[f64]) -> Vec<C64> {
    real.iter().map(|&r| C64::new(r, 0.0)).collect()
}

pub fn multiply_plaintexts(plaintext1: &[C64], plaintext2: &[C64]) -> Plaintext {
    assert_eq!(
        plaintext1.len(),
        plaintext2.len(),
        "Plaintexts must have the same length to be multiplied"
    );

    plaintext1
        .iter()
        .zip(plaintext2.iter())
        .map(|(&a, &b)| a * b)
        .collect()
}

pub fn scalar_mul_plaintext(scalar: C64, plaintext: &[C64]) -> Plaintext {
    plaintext.iter().map(|&a| scalar * a).collect()
}

pub fn add_plaintexts(plaintext1: &[C64], plaintext2: &[C64]) -> Plaintext {
    assert_eq!(
        plaintext1.len(),
        plaintext2.len(),
        "Plaintexts must have the same length to be added"
    );

    plaintext1
        .iter()
        .zip(plaintext2.iter())
        .map(|(&a, &b)| a + b)
        .collect()
}

#[cfg(test)]
mod tests {
    use bnum::types::{I1024, I256, I512};

    use super::*;
    use crate::algebra::big_int::{BigInt, ToFloat};
    use crate::algebra::complex::{raise_to_powers_of_two, raise_to_powers_of_two_rescale, C64};
    use crate::algebra::polynomial::Polynomial;
    use crate::keys::key_generator::{generate_pair_keys, generate_pair_keys_default};
    use crate::random_distributions::generate_random_vector;

    #[test]
    fn test_encrypt_decrypt_plaintext() {
        // Define parameters for key generation
        let dimension_exponent = 4;
        let q = I256::from(1 << 13);
        let level_max = 4;

        // Generate pair of keys
        let (mut client_key, _server_key) =
            generate_pair_keys_default::<I512>(dimension_exponent, level_max);

        // Create a sample message as a vector of f64
        let message_real = vec![12.0, 0.6, 31.0, 32.0, 44.0, 21.0, 432.0, 34.0];
        let message_plaintext = to_plaintext(&message_real);

        // Encrypt the message
        let ciphertext = client_key.encrypt(&message_plaintext, 50000.0).unwrap();

        // Decrypt the message
        let decrypted_plaintext = client_key.decrypt(&ciphertext);

        let error_max = ciphertext.upper_bound_error;
        println!("upperbound error: {:?}", error_max);
        // Verify that the decrypted message is close to the original message
        for (original, decrypted) in message_plaintext.iter().zip(decrypted_plaintext.iter()) {
            let diff = (*original - decrypted).magnitude();
            println!(
                "Original: {}, Decrypted: {}, Difference: {}",
                original, decrypted, diff
            );
            assert!(diff < error_max, "Difference in real part too large!");
        }
    }

    #[test]
    fn test_encrypt_add_decrypt_plaintext() {
        // Define parameters for key generation
        let dimension_exponent = 4;
        let q = I256::from(1 << 13);
        let level_max = 4;

        // Generate pair of keys
        let (mut client_key, mut server_key) =
            generate_pair_keys(dimension_exponent, q.clone(), level_max);

        // Create a sample message as a vector of f64
        let message_real1 = vec![
            12343.0, 23344.0, 32431.0, 43432.0, 3435.0, 3321.0, 432.0, 34532.0,
        ];
        let message_plaintext1 = to_plaintext(&message_real1);

        let message_real2 = vec![
            143.0, 234.0, 3231.0, 432.0, 34435.0, 33221.0, 43432.0, 34532.0,
        ];
        let message_plaintext2 = to_plaintext(&message_real2);

        // Encrypt the message
        let ct1 = client_key.encrypt(&message_plaintext1, 50000.0).unwrap();
        let ct2 = client_key.encrypt(&message_plaintext2, 50000.0).unwrap();

        let result = server_key.add(&ct1, &ct2);
        // Decrypt the message
        let clear_result = client_key.decrypt(&result);

        let expected_result = add_plaintexts(&message_plaintext1, &message_plaintext2);

        let error_max = result.upper_bound_error;
        println!("upperbound error: {:?}", error_max);
        // Verify that the decrypted message is close to the original message
        for (expected, decrypted) in expected_result.iter().zip(clear_result.iter()) {
            let diff = (*expected - decrypted).magnitude();
            println!(
                "Original: {}, Decrypted: {}, Difference: {}",
                expected, decrypted, diff
            );
            assert!(diff < error_max, "Difference in real part too large!");
        }
    }

    #[test]
    fn test_encrypt_pure_mul_decrypt_plaintext() {
        // Define parameters for key generation
        let dimension_exponent = 4;
        //let q = I256::from(1 << 13);
        let level_max = 5;

        // Generate pair of keys
        let (mut client_key, mut server_key) =
            generate_pair_keys_default::<I512>(dimension_exponent, level_max);

        // Create a sample message as a vector of f64
        let message_real1 = vec![0.0003, 0.01, 0.0003, 0.33, 0.00001, 0.1, 0.2, 0.01];
        let message_plaintext1 = to_plaintext(&message_real1);

        let message_real2 = vec![
            0.143, 0.234, 0.3231, 0.432, 0.34435, 0.33221, 0.43432, 0.34532,
        ];
        let message_plaintext2 = to_plaintext(&message_real2);

        // Encrypt the message
        let ct1 = client_key.encrypt(&message_plaintext1, 50000.0).unwrap();
        let ct2 = client_key.encrypt(&message_plaintext2, 50000.0).unwrap();

        let result = server_key.pure_mul(&ct1, &ct2);
        // Decrypt the message
        let clear_result = client_key.decrypt(&result);

        let expected_result = multiply_plaintexts(&message_plaintext1, &message_plaintext2);

        let error_max = result.upper_bound_error;
        println!("upperbound error: {:?}", error_max);
        // Verify that the decrypted message is close to the original message
        for (expected, decrypted) in expected_result.iter().zip(clear_result.iter()) {
            let diff = (*expected - decrypted).magnitude();
            println!(
                "Original: {}, Decrypted: {}, Difference: {}",
                expected, decrypted, diff
            );
            assert!(diff < error_max, "Difference in real part too large!");
        }
    }

    #[test]
    fn test_encrypt_mul_decrypt_plaintext() {
        // Define parameters for key generation
        let dimension_exponent = 5;
        // let q = I256::from(1 << 14);
        // let q_sqrt = (1 << 7) as f64;
        // let q_inverse = 1.0 / (1 << 14) as f6
        4;
        let level_max = 4;

        // Generate pair of keys
        let (mut client_key, mut server_key) =
            generate_pair_keys_default::<I512>(dimension_exponent, level_max);

        // Create a sample message as a vector of f64
        let message_real1 = vec![
            60.0, 70.0, 50.0, 42.0, 45.0, 32.0, 42.0, 72.0, 60.0, 70.0, 50.0, 42.0, 45.0, 32.0,
            42.0, 72.0,
        ];
        let message_plaintext1 = to_plaintext(&message_real1);
        // let message_plaintext1 =
        //     scalar_mul_plaintext(C64::new(q_sqrt, 0.0), &to_plaintext(&message_real1));

        let message_real2 = vec![
            60.0, 70.0, 50.0, 43.0, 45.0, 32.0, 42.0, 73.0, 60.0, 70.0, 50.0, 42.0, 45.0, 32.0,
            42.0, 72.0,
        ];
        let message_plaintext2 = to_plaintext(&message_real2);
        // let message_plaintext2 =
        //     scalar_mul_plaintext(C64::new(q_sqrt, 0.0), &to_plaintext(&message_real2));

        // Encrypt the message
        let ct1 = client_key.encrypt(&message_plaintext1, 73.0).unwrap();
        let ct2 = client_key.encrypt(&message_plaintext2, 73.0).unwrap();

        let result = server_key.mul(&ct1, &ct2).unwrap();
        // Decrypt the message
        let clear_result = client_key.decrypt(&result);

        let expected_result = multiply_plaintexts(&message_plaintext1, &message_plaintext2);
        // let expected_result = scalar_mul_plaintext(
        //     C64::new(q_inverse, 0.0),
        //     &multiply_plaintexts(&message_plaintext1, &message_plaintext2),
        // );

        let error_max = client_key.rescaled_error(&result);
        println!("upperbound error: {:?}", error_max);
        // Verify that the decrypted message is close to the original message
        for (expected, decrypted) in expected_result.iter().zip(clear_result.iter()) {
            let diff = (*expected - decrypted).magnitude();
            println!(
                "Original: {}, Decrypted: {}, Difference: {}",
                expected, decrypted, diff
            );
            assert!(diff < error_max, "Difference in real part too large!");
        }
    }

    #[test]
    fn test_raise_to_powers_of_two() {
        // Define parameters for key generation
        let dimension_exponent = 7;
        // let q = I256::from(1 << 14);
        // let qf = (1 << 14) as f64;
        // let q_sqrt = (1 << 7) as f64;
        // let q_inverse = 1.0 / (1 << 14) as f64;
        let level_max = 6;
        let n = 4;

        // Generate pair of keys
        let (mut client_key, mut server_key) =
            generate_pair_keys_default::<I1024>(dimension_exponent, level_max);

        let message_real = generate_random_vector(1 << (dimension_exponent - 1), -1.5, 1.5);
        // let message_real = vec![
        //     1.27, 1.01, 0.79, 1.49, 0.73, 1.06, 0.64, 1.29, 0.80, 0.69, 1.48, 0.70, 1.27, 1.39,
        //     1.18, 4.0,
        // ];

        // Create a sample message as a vector of f64
        // let message_real = vec![
        //     1.0, .0, 1.5, 0.1, 2.1, 1.0, 0.1, 0.06, 1.10, 1.0, 0.1, 2.0, 0.45, 0.32,
        //     2.0, 0.1,
        // ];
        let message_plaintext = to_plaintext(&message_real);
        //            scalar_mul_plaintext(C64::new(q_sqrt, 0.0), &to_plaintext(&message_real));
        let expected_result: Vec<Vec<C64>> = message_plaintext
            .iter()
            .map(|c| raise_to_powers_of_two(*c, n))
            .collect();

        let ciphertext = client_key.encrypt(&message_plaintext, 4.0).unwrap();
        let result = server_key
            .evaluation_key
            .raise_to_powers_of_two(&ciphertext, n)
            .unwrap();
        let clear_result: Vec<Plaintext> = result.iter().map(|c| client_key.decrypt(c)).collect();
        for i in 0..expected_result.len() {
            for j in 0..result.len() {
                let expected = expected_result[i][j];
                let obtained = clear_result[j][i];
                let error = (expected - obtained).magnitude();
                let expected_error = client_key.rescaled_error(&result[j]);
                // println!("indices: i:{} j:{}", i, j);
                // println!("expected: {}", expected);
                // println!("optained: {}", obtained);
                // println!("relative error: {}", error / expected.magnitude());
                // println!("expected error: {}", expected_error);
                assert!(error < expected_error, "error too big!");
            }
        }
    }

    #[test]
    fn test_apply_polynomial() {
        let dimension_exponent = 4;
        let level_max = 5;
        let n = 4;

        // Generate pair of keys
        let (mut client_key, mut server_key) =
            generate_pair_keys_default::<I1024>(dimension_exponent, level_max);

        let message_real = generate_random_vector(1 << (dimension_exponent - 1), -1.5, 1.5);

        let message_plaintext = to_plaintext(&message_real);
        //	let polynomial = Polynomial::<I256>::new(vec![I256::from(1), I256::from(2), I256::from(3), I256::from(4)]);
        let polynomial = Polynomial::<I1024>::new(vec![
            I1024::from(1),
            I1024::from(2),
            I1024::from(1),
            I1024::from(1),
        ]);
        let complex_polynomial = polynomial.to_c64();

        // Create a sample message as a vector of f64
        //let message_real = vec![1.0, 70.0, 50.0, 42.0, 45.0, 32.0, 42.0, 72.0];
        // let message_real = vec![
        //     1.0, 70.0, 50.0, 42.0, 45.0, 32.0, 42.0, 72.0, 60.0, 70.0, 50.0, 42.0, 45.0, 32.0,
        //     42.0, 72.0,
        // ];
        // let message_plaintext =
        //     scalar_mul_plaintext(C64::new(q_sqrt, 0.0), &to_plaintext(&message_real));
        let expected_result: Vec<C64> = message_plaintext
            .iter()
            .map(|c| {
                complex_polynomial.eval(*c)
                // + complex_polynomial.ref_coefficients()[0] * C64::new(1.0 - qf, 0.0)
            })
            .collect();

        let ciphertext = client_key.encrypt(&message_plaintext, 1.5).unwrap();
        let result = server_key
            .apply_polynomial(&polynomial, &ciphertext)
            .unwrap();
        let clear_result = client_key.decrypt(&result);
        //        let clear_result: Vec<Plaintext> = result.iter().map(|c| client_key.decrypt(c)).collect();
        for i in 0..clear_result.len() {
            let expected = expected_result[i];
            let obtained = clear_result[i];
            let error = (expected - obtained).magnitude();
            let expected_error = result.upper_bound_error;
            println!("index: i:{} ", i);
            println!("expected: {}", expected);
            println!("optained: {}", obtained);
            println!("relative error: {}", error / expected.magnitude());
            println!("expected error: {}", expected_error);
            assert!(error <= expected_error, "error too big!");
        }
    }

    #[test]
    fn test_trivial_encryption_scalar() {
        let dimension_exponent = 13;
        let level_max = 4;
        let (mut client_key, mut server_key) =
            generate_pair_keys_default::<I512>(dimension_exponent, level_max);

        let scalar = I256::from(1234);
        let scalar_f = scalar.to_float();
        let trivial_encryption = server_key
            .evaluation_key
            .trivial_encryption_scalar(I512::from(1234));
        let result = client_key.decrypt(&trivial_encryption);
        println!("expected_error: {}", trivial_encryption.upper_bound_error);
        println!("expected_result: {}", scalar_f);
        for r in result {
            println!("result: {}", r);
            assert!((r - C64::new(1234.0, 0.0)).magnitude() < 1.0);

            // // Define parameters for key generation

            // let q = I256::from(1 << 14);
            // let qf = (1 << 14) as f64;
            // let q_sqrt = (1 << 7) as f64;
            // let q_inverse = 1.0 / (1 << 14) as f64;
            // let level_max = 4;
            // let n = 3;

            // // Generate pair of keys
            // let (mut client_key, mut server_key) =
            //     generate_pair_keys(dimension_exponent, q.clone(), level_max);
            // let scalar = I256::from(1234);
            // let scalar_f = scalar.to_float();
            // let trivial_encryption = server_key.trivial_encryption_scalar(I256::from(1234));
            // let result = client_key.decrypt(&trivial_encryption);
            // println!("expected_error: {}", trivial_encryption.upper_bound_error);
            // println!("expected_result: {}", scalar_f);
            // for r in result {
            //     println!("result: {}", r);
            //     assert!((r - C64::new(1234.0, 0.0)).magnitude() < 1.0);
        }
    }
}
