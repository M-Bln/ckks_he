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
    fn new(value: i64) -> Self;
}

impl<T: BigInt> Zero for T {
    fn zero(&self) -> Self {
        T::from(0)
    }
}

impl BigInt for i64 {
    fn new(value: i64) -> Self {
        value
    }
}

#[cfg(test)]
mod tests {
    use crate::algebra::big_int::BigInt;

    #[test]
    fn test_new() {
        let big_int = i64::new(42);
        assert_eq!(big_int, 42);
    }

    #[test]
    fn test_add() {
        let a = i64::new(10);
        let b = i64::new(20);
        let result = a + b;
        assert_eq!(result, 30);
    }

    #[test]
    fn test_sub() {
        let a = i64::new(30);
        let b = i64::new(10);
        let result = a - b;
        assert_eq!(result, 20);
    }

    #[test]
    fn test_mul() {
        let a = i64::new(5);
        let b = i64::new(6);
        let result = a * b;
        assert_eq!(result, 30);
    }

    #[test]
    fn test_div() {
        let a = i64::new(50);
        let b = i64::new(5);
        let result = a / b;
        assert_eq!(result, 10);
    }

    #[test]
    fn test_rem() {
        let a = i64::new(17);
        let b = i64::new(5);
        let result = a % b;
        assert_eq!(result, 2);
    }

    #[test]
    fn test_add_ref() {
        let a = i64::new(15);
        let b = i64::new(25);
        let result = a + &b;
        assert_eq!(result, 40);
    }

    #[test]
    fn test_sub_ref() {
        let a = i64::new(100);
        let b = i64::new(30);
        let result = a - &b;
        assert_eq!(result, 70);
    }

    #[test]
    fn test_mul_ref() {
        let a = i64::new(7);
        let b = i64::new(8);
        let result = a * &b;
        assert_eq!(result, 56);
    }

    #[test]
    fn test_div_ref() {
        let a = i64::new(80);
        let b = i64::new(4);
        let result = a / &b;
        assert_eq!(result, 20);
    }

    #[test]
    fn test_mod_ref() {
        let a = i64::new(20);
        let b = i64::new(4);
        let result = a % &b;
        assert_eq!(result, 0);
    }
}
