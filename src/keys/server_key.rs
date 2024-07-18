use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex, C64};
use crate::algebra::polynomial::Polynomial;
use crate::ciphertext::Ciphertext;
use crate::encoding::Encoder;
use crate::keys::evaluation_key::{EvaluationKey, OperationError};
use crate::keys::key_generator::KeyGenerationParameters;
use crate::keys::public_key::PublicKey;
use crate::keys::secret_key::SecretKey;

/// Encapsulates the data necessary for the server side, including the public key, evaluation key, and an encoder.
///
/// The `ServerKey` struct is responsible for managing the homomorphic operations that can be performed on the ciphertexts.
/// It holds the public key for encryption, the evaluation key for homomorphic operations, and an encoder for encoding and decoding
/// plaintext messages. This struct provides methods for homomorphic addition, multiplication, rescaling, and evaluation of polynomials
/// on ciphertexts.
///
/// # Examples
///
/// Basic usage:
///
/// ```
/// use ckks_he::keys::key_generator::generate_pair_keys_toy;
/// use ckks_he::random_distributions::generate_random_vector;
/// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};
///
/// let (mut client_key, server_key) = generate_pair_keys_toy();
/// let plaintext_dimension = client_key.plaintext_dimension();
/// let plaintext_bound = 5.0;
///
/// let real_plaintext1 = generate_random_vector(
///     plaintext_dimension,
///     -plaintext_bound,
///     plaintext_bound,
/// );
/// let plaintext1 = to_plaintext(&real_plaintext1);
/// let ciphertext1 = client_key.encrypt(&plaintext1, plaintext_bound).unwrap();
///
/// let real_plaintext2 = generate_random_vector(
///     plaintext_dimension,
///     -plaintext_bound,
///     plaintext_bound,
/// );
/// let plaintext2 = to_plaintext(&real_plaintext2);
/// let ciphertext2 = client_key.encrypt(&plaintext2, plaintext_bound).unwrap();
/// let result = server_key.add(&ciphertext1, &ciphertext2);
/// let decrypted = client_key.decrypt(&result);
///
/// let expected: Plaintext = plaintext1.iter().zip(plaintext2.iter()).map(|(c1, c2)| *c1+c2).collect();
/// let error = calculate_error(&expected, &decrypted);
/// assert!(error < client_key.rescaled_error(&result), "Error larger than expected");
/// ```
pub struct ServerKey<T: BigInt> {
    pub public_key: PublicKey<T>,
    pub evaluation_key: EvaluationKey<T>,
    pub encoder: Encoder<T>,
}

impl<T: BigInt> ServerKey<T> {
    /// Creates a new `ServerKey` with the given parameters, public key, and evaluation key.
    ///
    /// This function initializes the `ServerKey` by setting the public key, evaluation key, and encoder.
    ///
    /// # Arguments
    ///
    /// * `parameters` - The key generation parameters.
    /// * `public_key` - The public key used for encryption.
    /// * `evaluation_key` - The evaluation key used for homomorphic operations.
    ///
    /// # Panics
    ///
    /// Panics if `parameters.dimension_exponent` is less than or equal to zero.
    pub fn new(
        parameters: KeyGenerationParameters<T>,
        public_key: PublicKey<T>,
        evaluation_key: EvaluationKey<T>,
    ) -> Self {
        assert!(parameters.dimension_exponent > 0);
        let plaintext_dimension = 1 << (parameters.dimension_exponent - 1); // 2^(dimension_exponent-1)
        let modulus = parameters.q.fast_exp(parameters.level_max) * &parameters.mul_scaling;
        let encoder = Encoder::new(
            parameters.dimension_exponent,
            modulus,
            parameters.q.to_float(),
        );
        ServerKey {
            public_key,
            evaluation_key,
            encoder,
        }
    }

    /// Returns an upper bound of the error expected after decryption.
    pub fn rescaled_error(&self, ct: &Ciphertext<T>) -> f64 {
        ct.upper_bound_error / self.encoder.scaling_factor
    }

    /// Returns an upper bound of the message expected after decryption.
    pub fn rescaled_upperbound_message(&self, ct: &Ciphertext<T>) -> f64 {
        ct.upper_bound_message / self.encoder.scaling_factor
    }
}

/// Macro to implement homomorphic operations on the `ServerKey` by delegating to `evaluation_key`.
macro_rules! delegate_to_eval_key {
    (
        $(
            $(#[$meta:meta])*
            $method:ident($($arg:ident: $arg_type:ty),*) -> $ret:ty
        ),*
    ) => {
        $(
            $(#[$meta])*
            pub fn $method(&self, $($arg: $arg_type),*) -> $ret {
                self.evaluation_key.$method($($arg),*)
            }
        )*
    };
}

impl<T: BigInt> ServerKey<T> {
    delegate_to_eval_key!(
    /// Adds two ciphertexts.
    ///
    /// # Arguments
    ///
    /// * `ct1` - The first ciphertext encrypting a message `m1`.
    /// * `ct2` - The second ciphertext encrypting a message `m2`.
    ///
    /// # Returns
    ///
    /// A new ciphertext representing the sum of the two inputs, encrypting `m1 + m2`.
    ///
    /// # Examples
    ///
    /// ```
    /// use ckks_he::keys::key_generator::generate_pair_keys_toy;
    /// use ckks_he::random_distributions::generate_random_vector;
    /// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};
    ///
    /// let (mut client_key, server_key) = generate_pair_keys_toy();
    /// let plaintext_dimension = client_key.plaintext_dimension();
    /// let plaintext_bound = 5.0;
    ///
    /// let real_plaintext1 = generate_random_vector(
    ///     plaintext_dimension,
    ///     -plaintext_bound,
    ///     plaintext_bound,
    /// );
    /// let plaintext1 = to_plaintext(&real_plaintext1);
    /// let ciphertext1 = client_key.encrypt(&plaintext1, plaintext_bound).unwrap();
    ///
    /// let real_plaintext2 = generate_random_vector(
    ///     plaintext_dimension,
    ///     -plaintext_bound,
    ///     plaintext_bound,
    /// );
    /// let plaintext2 = to_plaintext(&real_plaintext2);
    /// let ciphertext2 = client_key.encrypt(&plaintext2, plaintext_bound).unwrap();
    /// let result = server_key.add(&ciphertext1, &ciphertext2);
    /// let decrypted = client_key.decrypt(&result);
    ///
    /// let expected: Plaintext = plaintext1.iter().zip(plaintext2.iter()).map(|(c1, c2)| *c1 + *c2).collect();
    /// let error = calculate_error(&expected, &decrypted);
    /// assert!(error < client_key.rescaled_error(&result), "Error larger than expected");
    /// ```
        add(ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Ciphertext<T>,

        /// Rescales the ciphertext by the given level decrement.
        ///
        /// # Arguments
        ///
        /// * `ct` - The ciphertext to be rescaled.
        /// * `level_decrement` - The amount by which to decrement the level.
        ///
        /// # Returns
        ///
        /// A result indicating success or failure.
        rescale(ct: &mut Ciphertext<T>, level_decrement: u32) -> Result<(), OperationError>,

    /// Homomorphically multiplies two ciphertexts.
    ///
    /// # Arguments
    ///
    /// * `ct1` - The first ciphertext encrypting a message `m1`.
    /// * `ct2` - The second ciphertext encrypting a message `m2`.
    ///
    /// # Returns
    ///
    /// A new ciphertext encrypting `m1 * m2`.
    ///
    /// # Examples
    ///
    /// ```
    /// use ckks_he::keys::key_generator::generate_pair_keys_toy;
    /// use ckks_he::random_distributions::generate_random_vector;
    /// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};
    ///
    /// let (mut client_key, server_key) = generate_pair_keys_toy();
    /// let plaintext_dimension = client_key.plaintext_dimension();
    /// let plaintext_bound = 5.0;
    ///
    /// let real_plaintext1 = generate_random_vector(
    ///     plaintext_dimension,
    ///     -plaintext_bound,
    ///     plaintext_bound,
    /// );
    /// let plaintext1 = to_plaintext(&real_plaintext1);
    /// let ciphertext1 = client_key.encrypt(&plaintext1, plaintext_bound).unwrap();
    ///
    /// let real_plaintext2 = generate_random_vector(
    ///     plaintext_dimension,
    ///     -plaintext_bound,
    ///     plaintext_bound,
    /// );
    /// let plaintext2 = to_plaintext(&real_plaintext2);
    /// let ciphertext2 = client_key.encrypt(&plaintext2, plaintext_bound).unwrap();
    /// let result = server_key.mul(&ciphertext1, &ciphertext2).unwrap();
    /// let decrypted = client_key.decrypt(&result);
    ///
    /// let expected: Plaintext = plaintext1.iter().zip(plaintext2.iter()).map(|(c1, c2)| *c1 * *c2).collect();
    /// let error = calculate_error(&expected, &decrypted);
    /// assert!(error < client_key.rescaled_error(&result), "Error larger than expected");
    /// ```
        mul(ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Result<Ciphertext<T>, OperationError>,

    /// Raises a ciphertext to powers of two up to `2^n`.
    ///
    /// This function returns a vector of ciphertexts where each element is the input ciphertext raised
    /// to the power of `2^i` for `i` in `0..n`.
    ///
    /// # Arguments
    ///
    /// * `ct` - The ciphertext to be raised to powers of two.
    /// * `n` - The number of powers to compute.
    ///
    /// # Returns
    ///
    /// A vector of ciphertexts where each element is the input ciphertext raised to the power of `2^i`.
        raise_to_powers_of_two(ct: &Ciphertext<T>, n: usize) -> Result<Vec<Ciphertext<T>>, OperationError>,

        /// Trivially encrypts a scalar value.
        ///
        /// # Arguments
        ///
        /// * `scalar` - The scalar value to be encrypted.
        ///
        /// # Returns
        ///
        /// A ciphertext representing the encrypted scalar.
        ///
        /// # Examples
        ///
        /// ```
    /// use ckks_he::algebra::complex::{Complex, C64};
        /// use ckks_he::keys::key_generator::generate_pair_keys_toy;
        /// use ckks_he::random_distributions::generate_random_vector;
        /// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};
    /// use bnum::types::I1024;
        ///
        /// let (mut client_key, server_key) = generate_pair_keys_toy();
        /// let scalar = I1024::from(10);
        /// let ciphertext = server_key.trivial_encryption_scalar(scalar);
        /// let decrypted = client_key.decrypt(&ciphertext);
        ///
        /// let error = (C64::new(10.0, 0.0) - &decrypted[0]).magnitude();
        /// assert!(error <= client_key.rescaled_error(&ciphertext), "Error larger than expected");
        /// ```
        trivial_encryption_scalar(scalar: T) -> Ciphertext<T>,

        /// Evaluates a polynomial on a ciphertext.
        ///
        /// # Arguments
        ///
        /// * `polynomial` - The polynomial to be evaluated.
        /// * `ct` - The ciphertext on which to evaluate the polynomial.
        ///
        /// # Returns
        ///
        /// A ciphertext representing the polynomial evaluation.
        ///
        /// # Examples
        ///
        /// ```
        /// use ckks_he::keys::key_generator::generate_pair_keys_toy;
        /// use ckks_he::random_distributions::generate_random_vector;
        /// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};
        /// use ckks_he::algebra::polynomial::Polynomial;
    /// use bnum::types::I1024;
        ///
        /// let (mut client_key, server_key) = generate_pair_keys_toy();
        /// let plaintext_dimension = client_key.plaintext_dimension();
        /// let plaintext_bound = 5.0;
        ///
        /// let real_plaintext = generate_random_vector(
        ///     plaintext_dimension,
        ///     -plaintext_bound,
        ///     plaintext_bound,
        /// );
        /// let plaintext = to_plaintext(&real_plaintext);
        /// let ciphertext = client_key.encrypt(&plaintext, plaintext_bound).unwrap();
        ///
        /// let polynomial = Polynomial::new(vec![
        ///     I1024::from(2),
        ///     I1024::from(1),
        ///     I1024::from(2),
        ///     I1024::from(1),
        /// ]);
        /// let result = server_key.apply_polynomial(&polynomial, &ciphertext).unwrap();
        /// let decrypted = client_key.decrypt(&result);
        ///
    /// let complex_polynomial = polynomial.to_c64();
        /// let expected: Plaintext = plaintext.iter().map(|c| complex_polynomial.eval(*c)).collect();
        /// let error = calculate_error(&expected, &decrypted);
        /// assert!(error < client_key.rescaled_error(&result), "Error larger than expected");
        /// ```
        apply_polynomial(polynomial: &Polynomial<T>, ct: &Ciphertext<T>) -> Result<Ciphertext<T>, OperationError>
    );
}

// Adds two ciphertexts.

// # Arguments

// * `ct1` - The first ciphertext.
// * `ct2` - The second ciphertext.

// # Returns

// A new ciphertext representing the sum of the two inputs.

// # Examples

// ```
// use ckks_he::keys::key_generator::generate_pair_keys_toy;
// use ckks_he::random_distributions::generate_random_vector;
// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};

// let (mut client_key, server_key) = generate_pair_keys_toy();
// let plaintext_dimension = client_key.plaintext_dimension();
// let plaintext_bound = 5.0;

// let real_plaintext1 = generate_random_vector(
//     plaintext_dimension,
//     -plaintext_bound,
//     plaintext_bound,
// );
// let plaintext1 = to_plaintext(&real_plaintext1);
// let ciphertext1 = client_key.encrypt(&plaintext1, plaintext_bound).unwrap();

// let real_plaintext2 = generate_random_vector(
//     plaintext_dimension,
//     -plaintext_bound,
//     plaintext_bound,
// );
// let plaintext2 = to_plaintext(&real_plaintext2);
// let ciphertext2 = client_key.encrypt(&plaintext2, plaintext_bound).unwrap();
// let result = server_key.add(&ciphertext1, &ciphertext2);
// let decrypted = client_key.decrypt(&result);

// let expected: Plaintext = plaintext1.iter().zip(plaintext2.iter()).map(|(c1, c2)| *c1 + *c2).collect();
// let error = calculate_error(&expected, &decrypted);
// assert!(error < client_key.rescaled_error(&result), "Error larger than expected");
// ```

// Homomorphically multiplies two ciphertexts without rescaling.

// # Arguments

// * `ct1` - The first ciphertext.
// * `ct2` - The second ciphertext.

// # Returns

// A new ciphertext representing the product of the two inputs without rescaling.

// # Examples

// ```
// use ckks_he::keys::key_generator::generate_pair_keys_toy;
// use ckks_he::random_distributions::generate_random_vector;
// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};

// let (mut client_key, server_key) = generate_pair_keys_toy();
// let plaintext_dimension = client_key.plaintext_dimension();
// let plaintext_bound = 5.0;

// let real_plaintext1 = generate_random_vector(
//     plaintext_dimension,
//     -plaintext_bound,
//     plaintext_bound,
// );
// let plaintext1 = to_plaintext(&real_plaintext1);
// let ciphertext1 = client_key.encrypt(&plaintext1, plaintext_bound).unwrap();

// let real_plaintext2 = generate_random_vector(
//     plaintext_dimension,
//     -plaintext_bound,
//     plaintext_bound,
// );
// let plaintext2 = to_plaintext(&real_plaintext2);
// let ciphertext2 = client_key.encrypt(&plaintext2, plaintext_bound).unwrap();
// let result = server_key.pure_mul(&ciphertext1, &ciphertext2);
// let decrypted = client_key.decrypt(&result);

// let expected: Plaintext = plaintext1.iter().zip(plaintext2.iter()).map(|(c1, c2)| *c1 * *c2).collect();
// let error = calculate_error(&expected, &decrypted);
// assert!(error < client_key.rescaled_error(&result), "Error larger than expected");
// ```

// Raises the ciphertext to powers of two.

// # Arguments

// * `ct` - The ciphertext to be raised.
// * `n` - The maximum power of two.

// # Returns

// A vector of ciphertexts representing the powers of two.

// # Examples

// ```
// use ckks_he::algebra::complex::raise_to_powers_of_two;
// use ckks_he::keys::key_generator::generate_pair_keys_toy;
// use ckks_he::random_distributions::generate_random_vector;
// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};

// let (mut client_key, server_key) = generate_pair_keys_toy();
// let plaintext_dimension = client_key.plaintext_dimension();
// let plaintext_bound = 5.0;
// let n = 3;

// let real_plaintext = generate_random_vector(
//     plaintext_dimension,
//     -plaintext_bound,
//     plaintext_bound,
// );
// let plaintext = to_plaintext(&real_plaintext);
// let ciphertext = client_key.encrypt(&plaintext, plaintext_bound).unwrap();
// let powers = server_key.raise_to_powers_of_two(&ciphertext, n).unwrap();
// let clear_result = powers.iter().map(|ct| client_key.decrypt(&ct)).collect::<Vec<_>>();
// let expected_result: Vec<Vec<C64>> = message_plaintext
//    .iter()
//    .map(|c| raise_to_powers_of_two(*c, n))
//    .collect();
// for i in 0..expected_result.len() {
//    for j in 0..result.len() {
//        let expected = expected_result[i][j];
//        let obtained = clear_result[j][i];
//        let error = (expected - obtained).magnitude();
//        let expected_error = client_key.rescaled_error(&result[j]);
//        assert!(error <= expected_error, "error too big!");
//    }
// }
// ```

// /// Macro to implement homomorphic operations on the ServerKey by delegating to evaluation_key
// macro_rules! delegate_to_eval_key {
//     ($($method:ident($($arg:ident: $arg_type:ty),*) -> $ret:ty),*) => {
//         $(
//             pub fn $method(&self, $($arg: $arg_type),*) -> $ret {
//                 self.evaluation_key.$method($($arg),*)
//             }
//         )*
//     };
// }

// impl<T: BigInt> ServerKey<T> {
//     delegate_to_eval_key!(
//         add(ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Ciphertext<T>,
//         rescale(ct: &mut Ciphertext<T>, level_decrement: u32) -> Result<(), OperationError>,
//         mul(ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Result<Ciphertext<T>, OperationError>,
//         pure_mul(ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Ciphertext<T>,
// 	raise_to_powers_of_two(ct: &Ciphertext<T>, n: usize) -> Result<Vec<Ciphertext<T>>,OperationError>,
// 	trivial_encryption_scalar(scalar: T) -> Ciphertext<T>,
// 	apply_polynomial(polynomial: &Polynomial<T>, ct: &Ciphertext<T>) -> Result<Ciphertext<T>,OperationError>
//     );
// }

// # Examples

// ```
// use ckks_he::keys::key_generator::generate_pair_keys_toy;
// use ckks_he::random_distributions::generate_random_vector;
// use ckks_he::keys::client_key::{Plaintext, calculate_error, to_plaintext};

// let (mut client_key, server_key) = generate_pair_keys_toy();
// let plaintext_dimension = client_key.plaintext_dimension();
// let plaintext_bound = 5.0;

// let real_plaintext1 = generate_random_vector(
//     plaintext_dimension,
//     -plaintext_bound,
//     plaintext_bound,
// );
// let plaintext1 = to_plaintext(&real_plaintext1);
// let ciphertext1 = client_key.encrypt(&plaintext1, plaintext_bound).unwrap();

// let real_plaintext2 = generate_random_vector(
//     plaintext_dimension,
//     -plaintext_bound,
//     plaintext_bound,
// );
// let plaintext2 = to_plaintext(&real_plaintext2);
// let ciphertext2 = client_key.encrypt(&plaintext2, plaintext_bound).unwrap();
// let result = server_key.mul(&ciphertext1, &ciphertext2).unwrap();
// let decrypted = client_key.decrypt(&result);

// let expected: Plaintext = plaintext1.iter().zip(plaintext2.iter()).map(|(c1, c2)| *c1 * *c2).collect();
// let error = calculate_error(&expected, &decrypted);
// assert!(error < client_key.rescaled_error(&result), "Error larger than expected");
// ```
