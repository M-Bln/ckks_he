use crate::algebra::big_int::{BigInt, Zero};
use crate::algebra::complex::{Complex, C64};
use crate::algebra::conversion_rounding::{f64_to_i256, i256_to_f64};
use std::cmp::{Eq, PartialEq};
use std::ops::{Add, Mul, Sub};

use bnum::types::I256;

#[derive(Copy, Clone, Debug)]
pub struct RingMod<T: BigInt> {
    pub value: T,
    pub modulus: T,
}

impl<T: BigInt> Zero for RingMod<T> {
    fn zero(&self) -> Self {
        Self::new(T::from(0), self.modulus)
    }
    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
}
impl<T: BigInt> RingMod<T> {
    pub fn new(value: T, modulus: T) -> Self {
        RingMod {
            modulus,
            value: value.remainder(&modulus),
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

    pub fn pow(&self, n: u64) -> Self {
        let mut base = *self;
        let mut exp = n;
        let mut result = RingMod::new(T::from(1), self.modulus);

        while exp > 0 {
            if exp % 2 == 1 {
                result = result.mul(base);
            }
            base = base.mul(base);
            exp /= 2;
        }

        result
    }
}

pub trait Rescale<T> {
    fn rescale(&mut self, scalar: T);
}

impl<T: BigInt> Rescale<T> for RingMod<T> {
    fn rescale(&mut self, scalar: T) {
        self.value = self.value / scalar;
        self.modulus = self.modulus / scalar;
    }
}

impl<T: BigInt> Rescale<RingMod<T>> for RingMod<T> {
    fn rescale(&mut self, scalar: RingMod<T>) {
        self.value = self.value / scalar.value;
        self.modulus = self.modulus / scalar.value;
    }
}

// /// Addition is performed modulo the smallest of each modulus
// impl<T: BigInt> Add for RingMod<T> {
//     type Output = Self;

//     fn add(self, other: Self) -> Self {
//         if self.modulus < other.modulus {
//             let value = self.value + other.value
//             RingMod {
//                 modulus: self.modulus,
//                 value,
//             }
//         } else {
//             let value =
//                 ((self.value % other.modulus.clone()) + other.value) % other.modulus.clone();
//             RingMod {
//                 modulus: other.modulus,
//                 value,
//             }
//         }
//     }
// }

/// Addition is performed modulo the smallest of each modulus
impl<T: BigInt> Add for RingMod<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let modulus = std::cmp::min(self.modulus, other.modulus);
        let value = self.value + other.value;
        Self::new(value, modulus)
    }
}

/// Addition of RingMod with reference, it is performed modulo the smallest of each modulus
impl<'a, T: BigInt> Add<&'a RingMod<T>> for RingMod<T> {
    type Output = Self;

    fn add(self, other: &Self) -> Self {
        let modulus = std::cmp::min(self.modulus.clone(), other.modulus.clone());
        let value = self.value + other.value.clone();
        Self::new(value, modulus)
    }
}

/// Subtraction is performed modulo the smallest of each modulus
impl<T: BigInt> Sub for RingMod<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let modulus = std::cmp::min(self.modulus.clone(), other.modulus.clone());
        let value = self.value - other.value;
        Self::new(value, modulus)
    }
}

/// Subtraction of RingMod with reference
impl<'a, T: BigInt> Sub<&'a RingMod<T>> for RingMod<T> {
    type Output = Self;

    fn sub(self, other: &'a RingMod<T>) -> Self {
        let modulus = std::cmp::min(self.modulus.clone(), other.modulus.clone());
        let value = self.value - other.value.clone();
        Self::new(value, modulus)
    }
}

/// Multiplication is performed modulo the smallest of each modulus
impl<T: BigInt> Mul for RingMod<T> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let modulus = std::cmp::min(self.modulus.clone(), other.modulus.clone());
        let value = self.value * other.value;
        Self::new(value, modulus)
    }
}

/// Multiplication of RingMod with reference
impl<'a, T: BigInt> Mul<&'a RingMod<T>> for RingMod<T> {
    type Output = Self;

    fn mul(self, other: &'a RingMod<T>) -> Self {
        let modulus = std::cmp::min(self.modulus.clone(), other.modulus.clone());
        let value = self.value * other.value.clone();
        Self::new(value, modulus)
    }
}

// /// Addition of RingMod with reference, it is performed modulo the smallest of each modulus
// impl<'a, T: BigInt> Add<&'a RingMod<T>> for RingMod<T> {
//     type Output = Self;

//     fn add(self, other: &Self) -> Self {
//         if self.modulus < other.modulus {
//             let value = (self.value + (other.value % self.modulus.clone())) % self.modulus.clone();
//             RingMod {
//                 modulus: self.modulus,
//                 value,
//             }
//         } else {
//             let value =
//                 ((self.value % other.modulus.clone()) + other.value) % other.modulus.clone();
//             RingMod {
//                 modulus: other.modulus,
//                 value,
//             }
//         }
//     }
// }

// /// Subtraction is performed modulo the smallest of each modulus
// impl<T: BigInt> Sub for RingMod<T> {
//     type Output = Self;

//     fn sub(self, other: Self) -> Self {
//         if self.modulus < other.modulus {
//             let value = (self.value - (other.value % self.modulus.clone())) % self.modulus.clone();
//             RingMod {
//                 modulus: self.modulus,
//                 value,
//             }
//         } else {
//             let value =
//                 ((self.value % other.modulus.clone()) - other.value) % other.modulus.clone();
//             RingMod {
//                 modulus: other.modulus,
//                 value,
//             }
//         }
//     }
// }

// /// Subtraction of RingMod with reference
// impl<'a, T: BigInt> Sub<&'a RingMod<T>> for RingMod<T> {
//     type Output = Self;

//     fn sub(self, other: &'a RingMod<T>) -> Self {
//         if self.modulus < other.modulus {
//             let value = (self.value - (other.value % self.modulus.clone())) % self.modulus.clone();
//             RingMod {
//                 modulus: self.modulus,
//                 value,
//             }
//         } else {
//             let value =
//                 ((self.value % other.modulus.clone()) - other.value) % other.modulus.clone();
//             RingMod {
//                 modulus: other.modulus,
//                 value,
//             }
//         }
//     }
// }

// /// Multiplication is performed modulo the smallest of each modulus
// impl<T: BigInt> Mul for RingMod<T> {
//     type Output = Self;

//     fn mul(self, other: Self) -> Self {
//         if self.modulus < other.modulus {
//             let value = (self.value * (other.value % self.modulus.clone())) % self.modulus.clone();
//             RingMod {
//                 modulus: self.modulus,
//                 value,
//             }
//         } else {
//             let value =
//                 ((self.value % other.modulus.clone()) * other.value) % other.modulus.clone();
//             RingMod {
//                 modulus: other.modulus,
//                 value,
//             }
//         }
//     }
// }

// /// Multiplication of RingMod with reference
// impl<'a, T: BigInt> Mul<&'a RingMod<T>> for RingMod<T> {
//     type Output = Self;

//     fn mul(self, other: &'a RingMod<T>) -> Self {
//         if self.modulus < other.modulus {
//             let value = (self.value * (other.value % self.modulus.clone())) % self.modulus.clone();
//             RingMod {
//                 modulus: self.modulus,
//                 value,
//             }
//         } else {
//             let value =
//                 ((self.value % other.modulus.clone()) * other.value) % other.modulus.clone();
//             RingMod {
//                 modulus: other.modulus,
//                 value,
//             }
//         }
//     }
// }

// impl<T: BigInt> Sub for RingMod<T> {
//     type Output = Self;

//     fn sub(self, other: Self) -> Self {
//         assert!(
//             self.modulus == other.modulus,
//             "Moduli must be equal for addition"
//         );
//         let value = (self.value - other.value) % self.modulus;
//         RingMod {
//             modulus: self.modulus,
//             value,
//         }
//     }
// }

// // Addition of RingMod with reference
// impl<'a, T: BigInt> Sub<&'a RingMod<T>> for RingMod<T> {
//     type Output = Self;

//     fn sub(self, other: &'a RingMod<T>) -> Self {
//         assert!(
//             self.modulus == other.modulus,
//             "Moduli must be equal for addition"
//         );
//         let value = (self.value - other.value) % self.modulus;
//         RingMod {
//             modulus: self.modulus,
//             value,
//         }
//     }
// }

// // Multiplication of RingMod
// impl<T: BigInt> Mul for RingMod<T> {
//     type Output = Self;

//     fn mul(self, other: Self) -> Self {
//         assert!(
//             self.modulus == other.modulus,
//             "Moduli must be equal for multiplication"
//         );
//         let value = (self.value * other.value) % self.modulus;
//         RingMod {
//             modulus: self.modulus,
//             value,
//         }
//     }
// }

// // Multiplication of RingMod with reference
// impl<'a, T: BigInt> Mul<&'a RingMod<T>> for RingMod<T> {
//     type Output = Self;

//     fn mul(self, other: &'a RingMod<T>) -> Self {
//         assert!(
//             self.modulus == other.modulus,
//             "Moduli must be equal for multiplication"
//         );
//         let value = (self.value * other.value) % self.modulus;
//         RingMod {
//             modulus: self.modulus,
//             value,
//         }
//     }
// }

impl<T: BigInt> PartialEq for RingMod<T> {
    fn eq(&self, other: &Self) -> bool {
        self.modulus == other.modulus && (self.value - other.value) % self.modulus == T::new(0)
    }
}

impl<T: BigInt> Eq for RingMod<T> {}

// impl RingMod<I256> {
//     pub fn to_c64(&self) -> C64 {
//         let real = i256_to_f64(self.value.clone());
//         C64::new(real, 0.0)
//     }
// }

// pub fn c64_to_ring_mod_256(complex: &C64, modulus: I256) -> RingMod<I256> {
//     RingMod::new(f64_to_i256(complex.real()), modulus)
// }

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

    #[test]
    fn test_ring_mod_rescale() {
        let modulus = I256::new(100);
        let value = I256::new(50);
        let mut ringmod = RingMod::new(value.clone(), modulus.clone());

        println!("Original RingMod: {:?}", ringmod);

        let scalar = I256::new(2);
        ringmod.rescale(scalar.clone());

        println!("Rescaled RingMod: {:?}", ringmod);

        assert_eq!(ringmod.value, value / scalar.clone());
        assert_eq!(ringmod.modulus, modulus / scalar);
    }

    #[test]
    fn test_add_ringmod_different_moduli() {
        let modulus1 = I256::new(100);
        let value1 = I256::new(50);
        let ringmod1 = RingMod {
            modulus: modulus1,
            value: value1,
        };

        let modulus2 = I256::new(200);
        let value2 = I256::new(75);
        let ringmod2 = RingMod {
            modulus: modulus2,
            value: value2,
        };

        // Addition should be done in the smaller ring
        let result = ringmod1 + ringmod2;
        println!("Result of addition with different moduli: {:?}", result);

        assert_eq!(result.modulus, I256::new(100));
        assert_eq!(result.value, I256::new((50 + (75 % 100)) % 100));
    }

    #[test]
    fn test_sub_ringmod_different_moduli() {
        let modulus1 = I256::new(100);
        let value1 = I256::new(50);
        let ringmod1 = RingMod {
            modulus: modulus1,
            value: value1,
        };

        let modulus2 = I256::new(200);
        let value2 = I256::new(75);
        let ringmod2 = RingMod {
            modulus: modulus2,
            value: value2,
        };

        // Subtraction should be done in the smaller ring
        let result = ringmod1 - ringmod2;
        println!("Result of subtraction with different moduli: {:?}", result);

        assert_eq!(result.modulus, I256::new(100));
        assert_eq!(result.value, I256::new((50 - (75 % 100)) % 100));
    }

    #[test]
    fn test_mul_ringmod_different_moduli() {
        let modulus1 = I256::new(100);
        let value1 = I256::new(50);
        let ringmod1 = RingMod {
            modulus: modulus1,
            value: value1,
        };

        let modulus2 = I256::new(200);
        let value2 = I256::new(75);
        let ringmod2 = RingMod {
            modulus: modulus2,
            value: value2,
        };

        // Multiplication should be done in the smaller ring
        let result = ringmod1 * ringmod2;
        println!(
            "Result of multiplication with different moduli: {:?}",
            result
        );

        assert_eq!(result.modulus, I256::new(100));
        assert_eq!(result.value, I256::new((50 * (75 % 100)) % 100));
    }
}
