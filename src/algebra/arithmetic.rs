use crate::algebra::big_int::BigInt;
use std::cmp::{Eq, PartialEq};
use std::ops::{Add, Mul};

#[derive(Copy, Clone, Debug)]
struct RingMod<T: BigInt> {
    modulus: T,
    value: T,
}

impl<T: BigInt> RingMod<T> {
    pub fn new(modulus: T, value: T) -> Self {
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

    #[test]
    fn test_ring_mod_addition() {
        let modulus = 10;
        let value1 = 7;
        let value2 = 5;

        let ring_mod1 = RingMod::new(modulus, value1);
        let ring_mod2 = RingMod::new(modulus, value2);

        let result_add = ring_mod1 + ring_mod2;
        assert_eq!(result_add.value(), 2); // (7 + 5) % 10 = 2
    }

    #[test]
    fn test_ring_mod_multiplication() {
        let modulus = 10;
        let value1 = 7;
        let value2 = 5;

        let ring_mod1 = RingMod::new(modulus, value1);
        let ring_mod2 = RingMod::new(modulus, value2);

        let result_mul = ring_mod1 * ring_mod2;
        assert_eq!(result_mul.value(), 5); // (7 * 5) % 10 = 5
    }

    #[test]
    fn test_ring_mod_equality() {
        let modulus = 10;
        let value1 = 7;
        let value2 = 7;

        let ring_mod1 = RingMod::new(modulus, value1);
        let ring_mod2 = RingMod::new(modulus, value2);

        assert_eq!(ring_mod1, ring_mod2);
    }

    #[test]
    fn test_ring_mod_inequality() {
        let modulus = 10;
        let value1 = 7;
        let value2 = 5;

        let ring_mod1 = RingMod::new(modulus, value1);
        let ring_mod2 = RingMod::new(modulus, value2);

        assert_ne!(ring_mod1, ring_mod2);
    }
}