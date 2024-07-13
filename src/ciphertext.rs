use std::ops::Add;
use bnum::types::I256;
use crate::algebra::cyclotomic_ring::CyclotomicRing;
use crate::algebra::arithmetic::RingMod;

type CiphertextRing = CyclotomicRing<RingMod<I256>>;

struct RawCiphertext(CiphertextRing, CiphertextRing);

impl<'a> Add<&'a RawCiphertext> for RawCiphertext {
    type Output = RawCiphertext;
    fn add(self, other: &RawCiphertext) -> RawCiphertext {
	RawCiphertext(self.0 + &other.0, self.1 + &other.1)
    }
}

impl RawCiphertext {
    pub fn mul(&self, other: &RawCiphertext, evk: &RawCiphertext, P: I256) -> RawCiphertext {
	let (d_0, d_1, d_2) = (
	    self.0.clone() * &other.0,
	    self.0.clone() * &other.1 + &(self.1.clone() * &other.0),
	    self.1.clone() * &other.1,
	);
	RawCiphertext(
	    d_0 ,
	    d_1,
	)
    }
}
