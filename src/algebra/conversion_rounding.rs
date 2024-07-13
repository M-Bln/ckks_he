use crate::algebra::arithmetic::RingMod;
use crate::algebra::complex::{Complex, C64};
use crate::algebra::cyclotomic_ring::CyclotomicRing;
use crate::algebra::polynomial::Polynomial;
use bnum::types::I256;

/// Convert a big integer to a float, loses precision
pub fn i256_to_f64(n: I256) -> f64 {
    n.to_string().parse::<f64>().unwrap()
}

/// Convert a float to a large integer, loses the decimals
pub fn f64_to_i256(n: f64) -> I256 {
    I256::from_str_radix(&format!("{:.0}", n), 10).unwrap()
}

impl RingMod<I256> {
    pub fn to_c64(&self) -> C64 {
        let real = i256_to_f64(self.value.clone());
        C64::new(real, 0.0)
    }
}

pub fn c64_to_ring_mod_256(complex: &C64, modulus: I256) -> RingMod<I256> {
    RingMod::new(f64_to_i256(complex.real()), modulus)
}

impl CyclotomicRing<RingMod<I256>> {
    pub fn to_c64(&self) -> Polynomial<C64> {
        let coefficients = self
            .polynomial
            .ref_coefficients()
            .iter()
            .map(|coeff| coeff.to_c64())
            .collect();
        Polynomial::new(coefficients)
    }
}

impl Polynomial<I256> {
    pub fn to_f64(&self) -> Polynomial<f64> {
        let coefficients = self
            .ref_coefficients()
            .iter()
            .map(|coeff| i256_to_f64(coeff.clone()))
            .collect();
        Polynomial::new(coefficients)
    }
}

impl Polynomial<RingMod<I256>> {
    pub fn to_cyclotomic(self, dimension_exponent: u32) -> CyclotomicRing<RingMod<I256>> {
        CyclotomicRing::new(self.coefficients(), 2_usize.pow(dimension_exponent))
    }
}

impl Polynomial<f64> {
    pub fn to_i256(&self) -> Polynomial<I256> {
        let coefficients = self
            .ref_coefficients()
            .iter()
            .map(|coeff| f64_to_i256(*coeff))
            .collect();
        Polynomial::new(coefficients)
    }
}

impl Polynomial<C64> {
    pub fn to_i256(&self) -> Polynomial<I256> {
        let coefficients = self
            .ref_coefficients()
            .iter()
            .map(|coeff| f64_to_i256(coeff.real()))
            .collect();
        Polynomial::new(coefficients)
    }
}

impl Polynomial<RingMod<I256>> {
    pub fn to_c64(&self) -> Polynomial<C64> {
        let coefficients = self
            .ref_coefficients()
            .iter()
            .map(|coeff| coeff.to_c64())
            .collect();
        Polynomial::new(coefficients)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bnum::types::I256;
    #[test]
    fn test_i256_to_f64() {
        // Test with a positive number
        let big_int = I256::from(123456789012345678901234567890i128);
        let float_value = i256_to_f64(big_int);
        assert!((float_value - 1.2345678901234568e29).abs() < 1e15);

        // Test with a negative number
        let big_int = I256::from(-123456789012345678901234567890i128);
        let float_value = i256_to_f64(big_int);
        assert!((float_value + 1.2345678901234568e29).abs() < 1e15);

        // Test with zero
        let big_int = I256::from(0);
        let float_value = i256_to_f64(big_int);
        assert_eq!(float_value, 0.0);
    }

    const EPSILON: f64 = 1e-10;

    #[test]
    fn test_f64_to_i256() {
        // Test with a positive number
        let float_value = 1.2345678901234568e29;
        println!("Original f64: {}", float_value);
        let big_int = f64_to_i256(float_value);
        println!("Converted to I256: {}", big_int);
        let converted_back = i256_to_f64(big_int);
        println!("Converted back to f64: {}", converted_back);
        let expected_value = I256::from(123456789012345678901234567890i128);
        println!("Expected I256: {}", expected_value);
        assert!((converted_back - float_value).abs() < EPSILON);

        // Test with a negative number
        let float_value = -1.2345678901234568e29;
        println!("Original f64: {}", float_value);
        let big_int = f64_to_i256(float_value);
        println!("Converted to I256: {}", big_int);
        let converted_back = i256_to_f64(big_int);
        println!("Converted back to f64: {}", converted_back);
        let expected_value = I256::from(-123456789012345678901234567890i128);
        println!("Expected I256: {}", expected_value);
        assert!((converted_back - float_value).abs() < EPSILON);

        // Test with zero
        let float_value = 0.0;
        println!("Original f64: {}", float_value);
        let big_int = f64_to_i256(float_value);
        println!("Converted to I256: {}", big_int);
        let converted_back = i256_to_f64(big_int);
        println!("Converted back to f64: {}", converted_back);
        let expected_value = I256::from(0);
        println!("Expected I256: {}", expected_value);
        assert!((converted_back - float_value).abs() < EPSILON);
    }
}
