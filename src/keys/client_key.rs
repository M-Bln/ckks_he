use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex,C64};
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
	let plaintext_dimension = 1 << (parameters.dimension_exponent-1); // 2^(dimension_exponent-1)
        let modulus = parameters.q.fast_exp(parameters.level_max) * &parameters.mul_scaling;
        let encoder = Encoder::new(parameters.dimension_exponent, modulus);
        ClientKey {
            encoder,
            secret_key,
            public_key,
	    plaintext_dimension,
        }
    }

    pub fn encrypt(&mut self, clear: &Plaintext, upper_bound: f64) -> Result<Ciphertext<T>,EncryptionError>{
	if clear.len() != self.plaintext_dimension {
	    return Err(EncryptionError::UnexpectedPlaintextDimension);
	}
	if !bounded_by(&clear, upper_bound) {
	    return Err(EncryptionError::ClearLargerThanUpperBound);
	}
	let encoded = self.encoder.encode(&clear);
	let encrypted = self.public_key.encrypt(&encoded, upper_bound);
	Ok(encrypted)
    }

    pub fn decrypt(&mut self, ct: &Ciphertext<T>) -> Plaintext{
	let encoded = self.secret_key.decrypt(ct);
	self.encoder.decode(&encoded)
    }
}

pub fn bounded_by(z: &[C64], bound: f64) -> bool {
    z.iter().all(|&c| c.magnitude() <= bound)    
}

pub fn to_plaintext(real: &[f64]) -> Vec<C64> {
    real.iter().map(|&r| C64::new(r, 0.0)).collect()
}

pub fn multiply_plaintexts(plaintext1: &[C64], plaintext2: &[C64]) -> Vec<C64> {
    assert_eq!(plaintext1.len(), plaintext2.len(), "Plaintexts must have the same length to be multiplied");

    plaintext1.iter().zip(plaintext2.iter()).map(|(&a, &b)| a * b).collect()
}

pub fn add_plaintexts(plaintext1: &[C64], plaintext2: &[C64]) -> Vec<C64> {
    assert_eq!(plaintext1.len(), plaintext2.len(), "Plaintexts must have the same length to be added");

    plaintext1.iter().zip(plaintext2.iter()).map(|(&a, &b)| a + b).collect()
}

#[cfg(test)]
mod tests {
    use bnum::types::I256;
    
    use super::*;
    use crate::algebra::big_int::BigInt;
    use crate::algebra::complex::C64;
    use crate::keys::key_generator::{generate_pair_keys};
    
    #[test]
    fn test_encrypt_decrypt_plaintext() {
        // Define parameters for key generation
        let dimension_exponent = 4;
        let q = I256::from(1 << 13);
        let level_max = 4;

        // Generate pair of keys
        let (mut client_key, _server_key) = generate_pair_keys(dimension_exponent, q.clone(), level_max);

        // Create a sample message as a vector of f64
        let message_real = vec![12343.0, 23344.0, 32431.0, 43432.0, 34435.0, 33221.0, 43432.0, 34532.0];
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
            println!("Original: {}, Decrypted: {}, Difference: {}", original, decrypted, diff);
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
        let (mut client_key, mut server_key) = generate_pair_keys(dimension_exponent, q.clone(), level_max);

        // Create a sample message as a vector of f64
        let message_real1 = vec![12343.0, 23344.0, 32431.0, 43432.0, 3435.0, 3321.0, 432.0, 34532.0];
        let message_plaintext1 = to_plaintext(&message_real1);

        let message_real2 = vec![143.0, 234.0, 3231.0, 432.0, 34435.0, 33221.0, 43432.0, 34532.0];
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
            println!("Original: {}, Decrypted: {}, Difference: {}", expected, decrypted, diff);
            assert!(diff < error_max, "Difference in real part too large!");
        }
    }
}
