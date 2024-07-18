use bnum::types::{I1024, I256, I512};

use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use crate::algebra::arithmetic::RingMod;
use crate::algebra::conversion_rounding::{
    f64_to_i1024, f64_to_i256, f64_to_i512, i1024_to_f64, i256_to_f64, i512_to_f64,
};
use crate::random_distributions::UniformSamplable;

pub trait Zero: Sized {
    fn zero(&self) -> Self;
    fn is_zero(&self) -> bool;
}

// Use local trait rather than From<f64> in order to implement BigInt for non local type i64
pub trait FromFloat {
    fn from_float(value: f64) -> Self;
}

pub trait ToFloat {
    fn to_float(&self) -> f64;
}

/// A trait representing large integers with common arithmetic operations.
pub trait BigInt:
    Sized
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> Div<&'a Self, Output = Self>
    + for<'a> Rem<&'a Self, Output = Self>
    + Copy
    + Default
    + fmt::Display
    + fmt::Debug
    + From<i64>
    + FromFloat
    + ToFloat
    + Eq
    + PartialOrd
    + Ord
    + UniformSamplable
    + Zero
{
    fn new(value: i64) -> Self {
        Self::from(value)
    }

    /// Reduce an integer modulo a modulus
    fn modulo(&self, modulus: Self) -> RingMod<Self> {
        RingMod::new(self.clone() % modulus, modulus)
    }

    /// Performs fast exponentiation of the value to the given exponent.
    ///
    /// # Examples
    ///
    /// ```
    /// use ckks_he::algebra::big_int::BigInt;
    /// use bnum::types::I1024;
    ///
    /// let exponent = 9;
    /// let big_int_value = I1024::from(2 as i64);
    ///
    /// let result = big_int_value.fast_exp(exponent);
    /// assert_eq!(result, I1024::from(512 as i64)); // 2^9
    /// ```
    fn fast_exp(&self, exponent: u32) -> Self {
        let mut base = *self;
        let mut exp = exponent;
        let mut result = Self::from(1);

        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base;
            }
            base = base * base;
            exp /= 2;
        }
        result
    }

    /// Picks a representative of `self` modulo `other` in the range `]-other/2, other/2]`.
    ///
    /// Assumes `other` is greater than 0.
    ///
    /// # Examples
    ///
    /// ```
    /// use ckks_he::algebra::big_int::BigInt;
    /// use bnum::types::I1024;
    ///
    /// let big_int_value = I1024::from(38);
    /// let big_int_other = I1024::from(10);
    ///
    /// let remainder = big_int_value.remainder(&big_int_other);
    /// assert_eq!(remainder, I1024::from(-2)); // 38 % 10 == -2
    /// ```
    fn remainder(&self, other: &Self) -> Self {
        let mut remainder = *self % other;
        if Self::from(2) * remainder > *other {
            remainder = remainder - other;
        } else if Self::from(-2) * remainder >= *other {
            remainder = remainder + other;
        }
        remainder
    }
}

impl<T: BigInt> Zero for T {
    fn zero(&self) -> Self {
        T::from(0)
    }
    fn is_zero(&self) -> bool {
        *self == T::from(0)
    }
}

impl BigInt for i64 {}

impl FromFloat for i64 {
    fn from_float(value: f64) -> i64 {
        value as i64
    }
}

impl ToFloat for i64 {
    fn to_float(&self) -> f64 {
        *self as f64
    }
}

impl BigInt for I256 {}

impl FromFloat for I256 {
    fn from_float(value: f64) -> I256 {
        f64_to_i256(value)
    }
}

impl ToFloat for I256 {
    fn to_float(&self) -> f64 {
        i256_to_f64(*self)
    }
}

impl BigInt for I512 {}

impl FromFloat for I512 {
    fn from_float(value: f64) -> I512 {
        f64_to_i512(value)
    }
}

impl ToFloat for I512 {
    fn to_float(&self) -> f64 {
        i512_to_f64(*self)
    }
}

impl BigInt for I1024 {}

impl FromFloat for I1024 {
    fn from_float(value: f64) -> I1024 {
        f64_to_i1024(value)
    }
}

impl ToFloat for I1024 {
    fn to_float(&self) -> f64 {
        i1024_to_f64(*self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bnum::types::I256;

    #[test]
    fn test_new_i64() {
        let big_int = i64::new(42);
        assert_eq!(big_int, 42);
    }

    #[test]
    fn test_add_i64() {
        let a = i64::new(10);
        let b = i64::new(20);
        let result = a + b;
        assert_eq!(result, 30);
    }

    #[test]
    fn test_sub_i64() {
        let a = i64::new(30);
        let b = i64::new(10);
        let result = a - b;
        assert_eq!(result, 20);
    }

    #[test]
    fn test_mul_i64() {
        let a = i64::new(5);
        let b = i64::new(6);
        let result = a * b;
        assert_eq!(result, 30);
    }

    #[test]
    fn test_div_i64() {
        let a = i64::new(50);
        let b = i64::new(5);
        let result = a / b;
        assert_eq!(result, 10);
    }

    #[test]
    fn test_rem_i64() {
        let a = i64::new(17);
        let b = i64::new(5);
        let result = a % b;
        assert_eq!(result, 2);
    }

    #[test]
    fn test_add_ref_i64() {
        let a = i64::new(15);
        let b = i64::new(25);
        let result = a + &b;
        assert_eq!(result, 40);
    }

    #[test]
    fn test_sub_ref_i64() {
        let a = i64::new(100);
        let b = i64::new(30);
        let result = a - &b;
        assert_eq!(result, 70);
    }

    #[test]
    fn test_mul_ref_i64() {
        let a = i64::new(7);
        let b = i64::new(8);
        let result = a * &b;
        assert_eq!(result, 56);
    }

    #[test]
    fn test_div_ref_i64() {
        let a = i64::new(80);
        let b = i64::new(4);
        let result = a / &b;
        assert_eq!(result, 20);
    }

    #[test]
    fn test_rem_ref_i64() {
        let a = i64::new(20);
        let b = i64::new(4);
        let result = a % &b;
        assert_eq!(result, 0);
    }

    #[test]
    fn test_new_i256() {
        let big_int = I256::new(42);
        assert_eq!(big_int, I256::from(42));
    }

    #[test]
    fn test_add_i256() {
        let a = I256::new(10);
        let b = I256::new(20);
        let result = a + b;
        assert_eq!(result, I256::from(30));
    }

    #[test]
    fn test_sub_i256() {
        let a = I256::new(30);
        let b = I256::new(10);
        let result = a - b;
        assert_eq!(result, I256::from(20));
    }

    #[test]
    fn test_mul_i256() {
        let a = I256::new(5);
        let b = I256::new(6);
        let result = a * b;
        assert_eq!(result, I256::from(30));
    }

    #[test]
    fn test_div_i256() {
        let a = I256::new(50);
        let b = I256::new(5);
        let result = a / b;
        assert_eq!(result, I256::from(10));
    }

    #[test]
    fn test_rem_i256() {
        let a = I256::new(17);
        let b = I256::new(5);
        let result = a % b;
        assert_eq!(result, I256::from(2));
    }

    #[test]
    fn test_add_ref_i256() {
        let a = I256::new(15);
        let b = I256::new(25);
        let result = a + &b;
        assert_eq!(result, I256::from(40));
    }

    #[test]
    fn test_sub_ref_i256() {
        let a = I256::new(100);
        let b = I256::new(30);
        let result = a - &b;
        assert_eq!(result, I256::from(70));
    }

    #[test]
    fn test_mul_ref_i256() {
        let a = I256::new(7);
        let b = I256::new(8);
        let result = a * &b;
        assert_eq!(result, I256::from(56));
    }

    #[test]
    fn test_div_ref_i256() {
        let a = I256::new(80);
        let b = I256::new(4);
        let result = a / &b;
        assert_eq!(result, I256::from(20));
    }

    #[test]
    fn test_rem_ref_i256() {
        let a = I256::new(20);
        let b = I256::new(4);
        let result = a % &b;
        assert_eq!(result, I256::from(0));
    }

    #[test]
    fn test_fast_exp() {
        let base = I256::from(2);
        let exponent = 10;
        let result = base.fast_exp(exponent);
        let expected = I256::from(1024);
        assert_eq!(result, expected);

        let base = I256::from(5);
        let exponent = 3;
        let result = base.fast_exp(exponent);
        let expected = I256::from(125);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_remainder() {
        let modulus = I256::from(10);

        let a = I256::from(25);
        let rem = a.remainder(&modulus);
        println!("a = {}, modulus = {}, remainder = {}", a, modulus, rem);
        assert_eq!(rem, I256::from(5));

        let b = I256::from(-25);
        let rem = b.remainder(&modulus);
        println!("b = {}, modulus = {}, remainder = {}", b, modulus, rem);
        assert_eq!(rem, I256::from(5));

        let f = I256::from(7);
        let rem = f.remainder(&modulus);
        println!("f = {}, modulus = {}, remainder = {}", f, modulus, rem);
        assert_eq!(rem, I256::from(-3));

        let e = I256::from(-9);
        let rem = e.remainder(&modulus);
        println!("e = {}, modulus = {}, remainder = {}", e, modulus, rem);
        assert_eq!(rem, I256::from(1));

        let f = I256::from(-19);
        let rem = f.remainder(&modulus);
        println!("f = {}, modulus = {}, remainder = {}", f, modulus, rem);
        assert_eq!(rem, I256::from(1));

        let f = I256::from(5);
        let rem = f.remainder(&modulus);
        println!("f = {}, modulus = {}, remainder = {}", f, modulus, rem);
        assert_eq!(rem, I256::from(5));

        let f = I256::from(-5);
        let rem = f.remainder(&modulus);
        println!("f = {}, modulus = {}, remainder = {}", f, modulus, rem);
        assert_eq!(rem, I256::from(5));
    }
}
