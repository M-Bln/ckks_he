use crate::algebra::big_int::BigInt;
use crate::encoding::Encoder;
use crate::keys::key_generator::KeyGenerationParameters;
use crate::keys::public_key::PublicKey;
use crate::keys::secret_key::SecretKey;

pub struct ClientKey<T: BigInt> {
    pub encoder: Encoder<T>,
    pub secret_key: SecretKey<T>,
    pub public_key: PublicKey<T>,
}

impl<T: BigInt> ClientKey<T> {
    pub fn new(
        parameters: KeyGenerationParameters<T>,
        secret_key: SecretKey<T>,
        public_key: PublicKey<T>,
    ) -> Self {
        let modulus = parameters.q.fast_exp(parameters.level_max) * &parameters.mul_scaling;
        let encoder = Encoder::new(parameters.dimension_exponent, modulus);
        ClientKey {
            encoder,
            secret_key,
            public_key,
        }
    }
}
