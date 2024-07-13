use bnum::types::I256;
use std::fmt;
use std::ops::{Add, Div, Mul, Rem, Sub};

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

// pub fn i256_to_f64(n: I256) -> f64 {
//     n.to_string().parse::<f64>().unwrap()
// }

// pub fn f64_to_i256(n: f64) -> I256 {
//     I256::from_str_radix(&format!("{:.0}", n), 10).unwrap()
// }

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

    // #[test]
    // fn test_i256_to_f64() {
    //     // Test with a positive number
    //     let big_int = I256::from(123456789012345678901234567890i128);
    //     let float_value = i256_to_f64(big_int);
    //     assert!((float_value - 1.2345678901234568e29).abs() < 1e15);

    //     // Test with a negative number
    //     let big_int = I256::from(-123456789012345678901234567890i128);
    //     let float_value = i256_to_f64(big_int);
    //     assert!((float_value + 1.2345678901234568e29).abs() < 1e15);

    //     // Test with zero
    //     let big_int = I256::from(0);
    //     let float_value = i256_to_f64(big_int);
    //     assert_eq!(float_value, 0.0);
    // }

    // const EPSILON: f64 = 1e-10;

    // #[test]
    // fn test_f64_to_i256() {
    //     // Test with a positive number
    //     let float_value = 1.2345678901234568e29;
    //     println!("Original f64: {}", float_value);
    //     let big_int = f64_to_i256(float_value);
    //     println!("Converted to I256: {}", big_int);
    //     let converted_back = i256_to_f64(big_int);
    //     println!("Converted back to f64: {}", converted_back);
    //     let expected_value = I256::from(123456789012345678901234567890i128);
    //     println!("Expected I256: {}", expected_value);
    //     assert!((converted_back - float_value).abs() < EPSILON);

    //     // Test with a negative number
    //     let float_value = -1.2345678901234568e29;
    //     println!("Original f64: {}", float_value);
    //     let big_int = f64_to_i256(float_value);
    //     println!("Converted to I256: {}", big_int);
    //     let converted_back = i256_to_f64(big_int);
    //     println!("Converted back to f64: {}", converted_back);
    //     let expected_value = I256::from(-123456789012345678901234567890i128);
    //     println!("Expected I256: {}", expected_value);
    //     assert!((converted_back - float_value).abs() < EPSILON);

    //     // Test with zero
    //     let float_value = 0.0;
    //     println!("Original f64: {}", float_value);
    //     let big_int = f64_to_i256(float_value);
    //     println!("Converted to I256: {}", big_int);
    //     let converted_back = i256_to_f64(big_int);
    //     println!("Converted back to f64: {}", converted_back);
    //     let expected_value = I256::from(0);
    //     println!("Expected I256: {}", expected_value);
    //     assert!((converted_back - float_value).abs() < EPSILON);
    // }
}
