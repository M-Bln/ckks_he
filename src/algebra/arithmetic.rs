use crate::algebra::big_int::{BigInt, Zero};
use std::cmp::{Eq, PartialEq};
use std::ops::{Add, Mul, Sub};

#[derive(Copy, Clone, Debug)]
pub struct RingMod<T: BigInt> {
    value: T,
    modulus: T,
}

impl<T: BigInt> Zero for RingMod<T> {
    fn zero(&self) -> Self {
        Self::new(T::from(0), self.modulus)
    }
}
impl<T: BigInt> RingMod<T> {
    pub fn new(value: T, modulus: T) -> Self {
        RingMod {
            modulus,
            value: value % &modulus,
        }
    }

    // Getter for the value
    pub fn value(&self) -> T {
        self.value
    }

    // Getter for the modulus
    pub fn modulus(&self) -> T {
        self.modulus
    }
}

// Addition of RingMod
impl<T: BigInt> Add for RingMod<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        assert!(
            self.modulus == other.modulus,
            "Moduli must be equal for addition"
        );
        let value = (self.value + other.value) % self.modulus;
        RingMod {
            modulus: self.modulus,
            value,
        }
    }
}

// Addition of RingMod with reference
impl<'a, T: BigInt> Add<&'a RingMod<T>> for RingMod<T> {
    type Output = Self;

    fn add(self, other: &'a RingMod<T>) -> Self {
        assert!(
            self.modulus == other.modulus,
            "Moduli must be equal for addition"
        );
        let value = (self.value + other.value) % self.modulus;
        RingMod {
            modulus: self.modulus,
            value,
        }
    }
}

impl<T: BigInt> Sub for RingMod<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        assert!(
            self.modulus == other.modulus,
            "Moduli must be equal for addition"
        );
        let value = (self.value - other.value) % self.modulus;
        RingMod {
            modulus: self.modulus,
            value,
        }
    }
}

// Addition of RingMod with reference
impl<'a, T: BigInt> Sub<&'a RingMod<T>> for RingMod<T> {
    type Output = Self;

    fn sub(self, other: &'a RingMod<T>) -> Self {
        assert!(
            self.modulus == other.modulus,
            "Moduli must be equal for addition"
        );
        let value = (self.value - other.value) % self.modulus;
        RingMod {
            modulus: self.modulus,
            value,
        }
    }
}

// Multiplication of RingMod
impl<T: BigInt> Mul for RingMod<T> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        assert!(
            self.modulus == other.modulus,
            "Moduli must be equal for multiplication"
        );
        let value = (self.value * other.value) % self.modulus;
        RingMod {
            modulus: self.modulus,
            value,
        }
    }
}

// Multiplication of RingMod with reference
impl<'a, T: BigInt> Mul<&'a RingMod<T>> for RingMod<T> {
    type Output = Self;

    fn mul(self, other: &'a RingMod<T>) -> Self {
        assert!(
            self.modulus == other.modulus,
            "Moduli must be equal for multiplication"
        );
        let value = (self.value * other.value) % self.modulus;
        RingMod {
            modulus: self.modulus,
            value,
        }
    }
}

impl<T: BigInt> PartialEq for RingMod<T> {
    fn eq(&self, other: &Self) -> bool {
        self.modulus == other.modulus && (self.value - other.value) % self.modulus == T::new(0)
    }
}

impl<T: BigInt> Eq for RingMod<T> {}

#[cfg(test)]
mod tests {
    use super::*;
    use bnum::types::I256;

    #[test]
    fn test_ring_mod_addition_i256() {
        let modulus = I256::from(10);
        let value1 = I256::from(7);
        let value2 = I256::from(5);

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        let result_add = ring_mod1 + ring_mod2;
        assert_eq!(result_add.value(), I256::from(2)); // (7 + 5) % 10 = 2
    }

    #[test]
    fn test_ring_mod_multiplication_i256() {
        let modulus = I256::from(10);
        let value1 = I256::from(7);
        let value2 = I256::from(5);

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        let result_mul = ring_mod1 * ring_mod2;
        assert_eq!(result_mul.value(), I256::from(5)); // (7 * 5) % 10 = 5
    }

    #[test]
    fn test_ring_mod_equality_i256() {
        let modulus = I256::from(10);
        let value1 = I256::from(7);
        let value2 = I256::from(7);

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        assert_eq!(ring_mod1, ring_mod2);
    }

    #[test]
    fn test_ring_mod_inequality_i256() {
        let modulus = I256::from(10);
        let value1 = I256::from(7);
        let value2 = I256::from(5);

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        assert_ne!(ring_mod1, ring_mod2);
    }

    #[test]
    fn test_ring_mod_subtraction_i256() {
        let modulus = I256::from(10);
        let value1 = I256::from(7);
        let value2 = I256::from(5);

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        let result_sub = ring_mod1 - ring_mod2;
        assert_eq!(result_sub.value(), I256::from(2)); // (7 - 5) % 10 = 2
    }

    #[test]
    fn test_ring_mod_addition() {
        let modulus = 10;
        let value1 = 7;
        let value2 = 5;

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        let result_add = ring_mod1 + ring_mod2;
        assert_eq!(result_add.value(), 2); // (7 + 5) % 10 = 2
    }

    #[test]
    fn test_ring_mod_multiplication() {
        let modulus = 10;
        let value1 = 7;
        let value2 = 5;

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        let result_mul = ring_mod1 * ring_mod2;
        assert_eq!(result_mul.value(), 5); // (7 * 5) % 10 = 5
    }

    #[test]
    fn test_ring_mod_equality() {
        let modulus = 10;
        let value1 = 7;
        let value2 = 7;

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        assert_eq!(ring_mod1, ring_mod2);
    }

    #[test]
    fn test_ring_mod_inequality() {
        let modulus = 10;
        let value1 = 7;
        let value2 = 5;

        let ring_mod1 = RingMod::new(value1, modulus);
        let ring_mod2 = RingMod::new(value2, modulus);

        assert_ne!(ring_mod1, ring_mod2);
    }
}
