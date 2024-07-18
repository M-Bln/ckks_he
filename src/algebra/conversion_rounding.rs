use crate::algebra::arithmetic::RingMod;
use crate::algebra::big_int::{BigInt, ToFloat, Zero};
use crate::algebra::complex::{Complex, C64};
use crate::algebra::cyclotomic_ring::CyclotomicRing;
use crate::algebra::polynomial::Polynomial;
use bnum::types::{I1024, I256, I512};
use std::ops::{Add, Sub};

impl<T> Polynomial<T>
where
    T: Add<Output = T>
        + for<'a> Add<&'a T, Output = T>
        + Sub<Output = T>
        + for<'a> Sub<&'a T, Output = T>
        + Clone
        + Zero,
{
    pub fn to_cyclotomic(self, dimension_exponent: u32) -> CyclotomicRing<T> {
        CyclotomicRing::new(self.coefficients(), 2_usize.pow(dimension_exponent))
    }
}

impl<T: BigInt> Polynomial<T> {
    /// Reduces the coefficients of a polynomial modulo modulus
    pub fn modulo(&self, modulus: T) -> Polynomial<RingMod<T>> {
        let coefficients = self
            .ref_coefficients()
            .iter()
            .map(|c| c.modulo(modulus))
            .collect();
        Polynomial::new(coefficients)
    }
}

// Conversion involving floats, thus  with loss of precision
/// Convert a big integer to a float, loses precision
pub fn i256_to_f64(n: I256) -> f64 {
    n.to_string().parse::<f64>().unwrap()
}

pub fn i512_to_f64(n: I512) -> f64 {
    n.to_string().parse::<f64>().unwrap()
}

pub fn i1024_to_f64(n: I1024) -> f64 {
    n.to_string().parse::<f64>().unwrap()
}

/// Convert a float to a large integer, loses the decimals
pub fn f64_to_i256(n: f64) -> I256 {
    I256::from_str_radix(&format!("{:.0}", n), 10).unwrap()
}

pub fn f64_to_i512(n: f64) -> I512 {
    I512::from_str_radix(&format!("{:.0}", n), 10).unwrap()
}

pub fn f64_to_i1024(n: f64) -> I1024 {
    I1024::from_str_radix(&format!("{:.0}", n), 10).unwrap()
}

/// Convert the integer coefficients into complex coefficients
impl RingMod<I256> {
    #[allow(dead_code)]
    pub fn to_c64(&self) -> C64 {
        let real = i256_to_f64(self.value.clone());
        C64::new(real, 0.0)
    }
}

// /// Round the real part
// pub fn c64_to_ring_mod_256(complex: &C64, modulus: I256) -> RingMod<I256> {
//     RingMod::new(f64_to_i256(complex.real()), modulus)
// }

/// Round the real part of the complex coefficients
impl CyclotomicRing<RingMod<I256>> {
    #[allow(dead_code)]
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

impl<T> CyclotomicRing<T>
where
    T: Add<Output = T>
        + for<'a> Add<&'a T, Output = T>
        + Sub<Output = T>
        + for<'a> Sub<&'a T, Output = T>
        + Clone
        + Zero
        + ToFloat,
{
    /// Converts the integer representatives of the coefficients into complex coefficients
    pub fn to_c64(&self) -> Polynomial<C64> {
        let coefficients = self
            .ref_coefficients()
            .iter()
            .map(|coeff| C64::new(coeff.to_float(), 0.0))
            .collect();
        Polynomial::new(coefficients)
    }
}

impl<T: BigInt> CyclotomicRing<RingMod<T>> {
    /// Pick integer representatives of coefficients in the range ]-modulus/2, modulus/2]
    pub fn to_integer(&self) -> CyclotomicRing<T> {
        let coefficients = self
            .ref_coefficients()
            .iter()
            .map(|coeff| coeff.value)
            .collect();
        CyclotomicRing::new(coefficients, self.dimension)
    }
}

impl<T: BigInt> Polynomial<T> {
    pub fn to_c64(&self) -> Polynomial<C64> {
        Polynomial::new(
            self.ref_coefficients()
                .iter()
                .map(|c| C64::new(c.to_float(), 0.0))
                .collect(),
        )
    }
}

impl Polynomial<C64> {
    pub fn to_integer<T: BigInt>(&self) -> Polynomial<T> {
        let coefficients = self
            .ref_coefficients()
            .iter()
            .map(|coeff| T::from_float(coeff.real()))
            .collect();
        Polynomial::new(coefficients)
    }
}

// impl Polynomial<C64> {
//     pub fn to_i256(&self) -> Polynomial<I256> {
//         let coefficients = self
//             .ref_coefficients()
//             .iter()
//             .map(|coeff| f64_to_i256(coeff.real()))
//             .collect();
//         Polynomial::new(coefficients)
//     }
// }

// impl Polynomial<RingMod<I256>> {
//     pub fn to_c64(&self) -> Polynomial<C64> {
//         let coefficients = self
//             .ref_coefficients()
//             .iter()
//             .map(|coeff| coeff.to_c64())
//             .collect();
//         Polynomial::new(coefficients)
//     }
// }

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
