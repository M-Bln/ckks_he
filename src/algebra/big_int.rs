use std::fmt;
use std::ops::{Add, Div, Mul, Rem, Sub};
use bnum::types::I256;

pub trait Zero: Sized {
    fn zero(&self) -> Self;
}

pub trait BigInt:
    Sized
    + Add<Output = Self>
    + Sub<Output = Self>
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
    + Eq
    + Zero
{
    fn new(value: i64) -> Self {
	Self::from(value)
    }
}

impl<T: BigInt> Zero for T {
    fn zero(&self) -> Self {
        T::from(0)
    }
}

impl BigInt for i64 {}

impl BigInt for I256 {}


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
}


