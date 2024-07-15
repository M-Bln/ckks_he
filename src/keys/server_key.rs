use crate::algebra::big_int::BigInt;
use crate::algebra::complex::{Complex,C64};
use crate::ciphertext::Ciphertext;
use crate::encoding::Encoder;
use crate::keys::key_generator::KeyGenerationParameters;
use crate::keys::evaluation_key::{EvaluationKey, OperationError};
use crate::keys::public_key::PublicKey;
use crate::keys::secret_key::SecretKey;


pub struct ServerKey<T: BigInt> {
    pub public_key: PublicKey<T>,
    pub evaluation_key: EvaluationKey<T>,
    pub encoder: Encoder<T>,
}


/// Encapsulate the data necessary to the server side i.e. the public key, the evaluation key and an encoder
impl<T: BigInt> ServerKey<T> {
    pub fn new(
        parameters: KeyGenerationParameters<T>,
        public_key: PublicKey<T>,
	evaluation_key: EvaluationKey<T>,
    ) -> Self {
	assert!(parameters.dimension_exponent > 0);
	let plaintext_dimension = 1 << (parameters.dimension_exponent-1); // 2^(dimension_exponent-1)
        let modulus = parameters.q.fast_exp(parameters.level_max) * &parameters.mul_scaling;
	let encoder = Encoder::new(parameters.dimension_exponent, modulus);
	ServerKey{
	    public_key,
	    evaluation_key,
	    encoder,
	}
    }

    
}

/// Macro to implement homomorphic operations on the ServerKey by delegating to evaluation_key
macro_rules! delegate_to_eval_key {
    ($($method:ident($($arg:ident: $arg_type:ty),*) -> $ret:ty),*) => {
        $(
            pub fn $method(&self, $($arg: $arg_type),*) -> $ret {
                self.evaluation_key.$method($($arg),*)
            }
        )*
    };
}

impl<T: BigInt> ServerKey<T> {
    delegate_to_eval_key!(
        add(ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Ciphertext<T>,
        rescale(ct: &mut Ciphertext<T>, level_decrement: u32) -> Result<(), OperationError>,
        mul(ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Result<Ciphertext<T>, OperationError>,
        pure_mul(ct1: &Ciphertext<T>, ct2: &Ciphertext<T>) -> Ciphertext<T>
    );
}
