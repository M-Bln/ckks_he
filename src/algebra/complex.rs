use std::f64::consts::PI;
use std::fmt;
use std::ops::{Add, Neg, Div, Mul, Sub};

use crate::algebra::big_int::Zero;

pub trait Complex:
    Sized
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> Div<&'a Self, Output = Self>
    + Copy
    + Default
    + fmt::Display
    + fmt::Debug
    + Zero
{
    type Real: From<f64> + Neg<Output = Self::Real>;
    fn new(real: Self::Real, imaginary: Self::Real) -> Self;
    fn real(&self) -> Self::Real;
    fn imaginary(&self) -> Self::Real;
    fn magnitude(&self) -> Self::Real;
    fn phase(&self) -> Self::Real;

    fn conjugate(&self) -> Self {
	Self::new(self.real(), -self.imaginary())
    }

    fn primitive_root_of_unity(n: u32) -> Self {
        let angle = 2.0 * PI / (n as f64);
        Self::new(Self::Real::from(angle.cos()), Self::Real::from(angle.sin()))
    }

    fn primitive_2_to_the_h_th_roots_of_unity(h: u32) -> Vec<Self> {
        if h == 0 {
            return vec![Self::new(Self::Real::from(1.0), Self::Real::from(0.0))];
        }
        let mut result = Vec::with_capacity(2_usize.pow(h - 1));
        let angle_0 = 2.0 * PI / (2_u64.pow(h) as f64);
        result.push(Self::new(
            Self::Real::from(angle_0.cos()),
            Self::Real::from(angle_0.sin()),
        ));
        for k in 1..2_u64.pow(h - 1) {
            let angle = angle_0 * ((2 * k) as f64 + 1.0);
            result.push(Self::new(
                Self::Real::from(angle.cos()),
                Self::Real::from(angle.sin()),
            ));
        }
        result
    }

    fn all_2_to_the_h_th_roots_of_unity(h: u32) -> Vec<Self> {
        if h == 0 {
            return vec![Self::new(Self::Real::from(1.0), Self::Real::from(0.0))];
        }
        let mut result = Vec::with_capacity(2_usize.pow(h));
        let angle_0 = 2.0 * PI / (2_u64.pow(h) as f64);
        result.push(Self::new(Self::Real::from(1.0), Self::Real::from(0.0)));
        for k in 1..2_u64.pow(h) {
            let angle = angle_0 * (k as f64);
            result.push(Self::new(
                Self::Real::from(angle.cos()),
                Self::Real::from(angle.sin()),
            ));
        }
        result
    }
}

#[derive(Clone, Copy)]
pub struct C64 {
    real: f64,
    imaginary: f64,
}

impl Complex for C64 {
    type Real = f64;

    fn new(real: Self::Real, imaginary: Self::Real) -> Self {
        C64 { real, imaginary }
    }

    fn real(&self) -> Self::Real {
        self.real
    }

    fn imaginary(&self) -> Self::Real {
        self.imaginary
    }

    fn magnitude(&self) -> f64 {
        (self.real * self.real + self.imaginary * self.imaginary).sqrt()
    }

    fn phase(&self) -> f64 {
        self.imaginary.atan2(self.real)
    }
}

impl Zero for C64 {
    fn zero(&self) -> Self {
        Self::new(0.0, 0.0)
    }
}

impl Add for C64 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        C64 {
            real: self.real + other.real,
            imaginary: self.imaginary + other.imaginary,
        }
    }
}

impl<'a> Add<&'a C64> for C64 {
    type Output = Self;

    fn add(self, other: &'a C64) -> Self {
        C64 {
            real: self.real + other.real,
            imaginary: self.imaginary + other.imaginary,
        }
    }
}

impl Sub for C64 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        C64 {
            real: self.real - other.real,
            imaginary: self.imaginary - other.imaginary,
        }
    }
}

impl<'a> Sub<&'a C64> for C64 {
    type Output = Self;

    fn sub(self, other: &'a C64) -> Self {
        C64 {
            real: self.real - other.real,
            imaginary: self.imaginary - other.imaginary,
        }
    }
}

impl Mul for C64 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        C64 {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real,
        }
    }
}

impl<'a> Mul<&'a C64> for C64 {
    type Output = Self;

    fn mul(self, other: &'a C64) -> Self {
        C64 {
            real: self.real * other.real - self.imaginary * other.imaginary,
            imaginary: self.real * other.imaginary + self.imaginary * other.real,
        }
    }
}

impl Div for C64 {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        let denominator = other.real * other.real + other.imaginary * other.imaginary;
        C64 {
            real: (self.real * other.real + self.imaginary * other.imaginary) / denominator,
            imaginary: (self.imaginary * other.real - self.real * other.imaginary) / denominator,
        }
    }
}

impl<'a> Div<&'a C64> for C64 {
    type Output = Self;

    fn div(self, other: &'a C64) -> Self {
        let denominator = other.real * other.real + other.imaginary * other.imaginary;
        C64 {
            real: (self.real * other.real + self.imaginary * other.imaginary) / denominator,
            imaginary: (self.imaginary * other.real - self.real * other.imaginary) / denominator,
        }
    }
}

impl fmt::Display for C64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.imaginary >= 0.0 {
            write!(f, "{} + {}i", self.real, self.imaginary)
        } else {
            write!(f, "{} - {}i", self.real, -self.imaginary)
        }
    }
}

impl fmt::Debug for C64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}i)", self.real, self.imaginary)
    }
}

impl Default for C64 {
    fn default() -> Self {
        Self::new(0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_creation_and_accessors() {
        let complex = C64::new(3.0, 4.0);
        assert_eq!(complex.real(), 3.0);
        assert_eq!(complex.imaginary(), 4.0);
    }

    #[test]
    fn test_addition() {
        let a = C64::new(1.0, 2.0);
        let b = C64::new(3.0, 4.0);
        let result = a + b;
        assert_eq!(result.real(), 4.0);
        assert_eq!(result.imaginary(), 6.0);
    }

    #[test]
    fn test_subtraction() {
        let a = C64::new(5.0, 6.0);
        let b = C64::new(1.0, 2.0);
        let result = a - b;
        assert_eq!(result.real(), 4.0);
        assert_eq!(result.imaginary(), 4.0);
    }

    #[test]
    fn test_multiplication() {
        let a = C64::new(2.0, 3.0);
        let b = C64::new(4.0, 5.0);
        let result = a * b;
        assert_eq!(result.real(), -7.0);
        assert_eq!(result.imaginary(), 22.0);
    }

    #[test]
    fn test_division() {
        let a = C64::new(5.0, 6.0);
        let b = C64::new(2.0, 1.0);
        let result = a / b;
        assert_eq!(result.real(), 3.2);
        assert_eq!(result.imaginary(), 1.4);
    }

    #[test]
    fn test_magnitude() {
        let complex = C64::new(3.0, 4.0);
        assert_eq!(complex.magnitude(), 5.0);
    }

    #[test]
    fn test_phase() {
        let complex = C64::new(1.0, 1.0);
        assert_eq!(complex.phase(), std::f64::consts::FRAC_PI_4);
    }

    #[test]
    fn test_primitive_root_of_unity() {
        let root = C64::primitive_root_of_unity(3);
        let expected_real = -0.5;
        let expected_imaginary = (3.0_f64).sqrt() / 2.0;
        assert!((root.real() - expected_real).abs() < 1e-10);
        assert!((root.imaginary() - expected_imaginary).abs() < 1e-10);
    }

    #[test]
    fn test_primitive_2_to_the_h_th_roots_of_unity() {
        // Test for h = 0
        let roots_h0 = C64::primitive_2_to_the_h_th_roots_of_unity(0);
        assert_eq!(roots_h0.len(), 1);
        assert!((roots_h0[0].real() - 1.0).abs() < 1e-10);
        assert!((roots_h0[0].imaginary() - 0.0).abs() < 1e-10);

        // Test for h = 1
        let roots_h1 = C64::primitive_2_to_the_h_th_roots_of_unity(1);
        assert_eq!(roots_h1.len(), 1);
        assert!((roots_h1[0].real() + 1.0).abs() < 1e-10); // cos(pi)
        assert!((roots_h1[0].imaginary() - 0.0).abs() < 1e-10); // sin(pi)

        // Test for h = 2
        let roots_h2 = C64::primitive_2_to_the_h_th_roots_of_unity(2);
        assert_eq!(roots_h2.len(), 2);
        let expected_h2 = vec![
            (0.0, 1.0),  // cos(pi/2), sin(pi/2)
            (0.0, -1.0), // cos(3pi/2), sin(3pi/2)
        ];
        for (i, root) in roots_h2.iter().enumerate() {
            assert!((root.real() - expected_h2[i].0).abs() < 1e-10);
            assert!((root.imaginary() - expected_h2[i].1).abs() < 1e-10);
        }

        // Test for h = 3
        let roots_h3 = C64::primitive_2_to_the_h_th_roots_of_unity(3);
        assert_eq!(roots_h3.len(), 4);
        let expected_h3 = vec![
            (0.70710678118, 0.70710678118),   // cos(pi/4), sin(pi/4)
            (-0.70710678118, 0.70710678118),  // cos(3pi/4), sin(3pi/4)
            (-0.70710678118, -0.70710678118), // cos(5pi/4), sin(5pi/4)
            (0.70710678118, -0.70710678118),  // cos(7pi/4), sin(7pi/4)
        ];
        for (i, root) in roots_h3.iter().enumerate() {
            if (root.real() - expected_h3[i].0).abs() >= 1e-10
                || (root.imaginary() - expected_h3[i].1).abs() >= 1e-10
            {
                println!(
                    "Mismatch at index {}: expected ({}, {}), got ({}, {})",
                    i,
                    expected_h3[i].0,
                    expected_h3[i].1,
                    root.real(),
                    root.imaginary()
                );
            }
            assert!((root.real() - expected_h3[i].0).abs() < 1e-10);
            assert!((root.imaginary() - expected_h3[i].1).abs() < 1e-10);
        }
        // for (i, root) in roots_h3.iter().enumerate() {
        //     assert!((root.real() - expected_h3[i].0).abs() < 1e-10);
        //     assert!((root.imaginary() - expected_h3[i].1).abs() < 1e-10);
        // }
    }

    #[test]
    fn test_all_2_to_the_h_th_roots_of_unity() {
        // Test for h = 0
        let roots_h0 = C64::all_2_to_the_h_th_roots_of_unity(0);
        assert_eq!(roots_h0.len(), 1);
        assert!((roots_h0[0].real() - 1.0).abs() < 1e-10);
        assert!((roots_h0[0].imaginary() - 0.0).abs() < 1e-10);

        // Test for h = 1
        let roots_h1 = C64::all_2_to_the_h_th_roots_of_unity(1);
        assert_eq!(roots_h1.len(), 2);
        let expected_h1 = vec![
            (1.0, 0.0),  // cos(0), sin(0)
            (-1.0, 0.0), // cos(pi), sin(pi)
        ];
        for (i, root) in roots_h1.iter().enumerate() {
            assert!((root.real() - expected_h1[i].0).abs() < 1e-10);
            assert!((root.imaginary() - expected_h1[i].1).abs() < 1e-10);
        }

        // Test for h = 2
        let roots_h2 = C64::all_2_to_the_h_th_roots_of_unity(2);
        assert_eq!(roots_h2.len(), 4);
        let expected_h2 = vec![
            (1.0, 0.0),  // cos(0), sin(0)
            (0.0, 1.0),  // cos(pi/2), sin(pi/2)
            (-1.0, 0.0), // cos(pi), sin(pi)
            (0.0, -1.0), // cos(3pi/2), sin(3pi/2)
        ];
        for (i, root) in roots_h2.iter().enumerate() {
            assert!((root.real() - expected_h2[i].0).abs() < 1e-10);
            assert!((root.imaginary() - expected_h2[i].1).abs() < 1e-10);
        }

        // Test for h = 3
        let roots_h3 = C64::all_2_to_the_h_th_roots_of_unity(3);
        assert_eq!(roots_h3.len(), 8);
        let expected_h3 = vec![
            (1.0, 0.0),                       // cos(0), sin(0)
            (0.70710678118, 0.70710678118),   // cos(pi/4), sin(pi/4)
            (0.0, 1.0),                       // cos(pi/2), sin(pi/2)
            (-0.70710678118, 0.70710678118),  // cos(3pi/4), sin(3pi/4)
            (-1.0, 0.0),                      // cos(pi), sin(pi)
            (-0.70710678118, -0.70710678118), // cos(5pi/4), sin(5pi/4)
            (0.0, -1.0),                      // cos(3pi/2), sin(3pi/2)
            (0.70710678118, -0.70710678118),  // cos(7pi/4), sin(7pi/4)
        ];
        for (i, root) in roots_h3.iter().enumerate() {
            if (root.real() - expected_h3[i].0).abs() >= 1e-10
                || (root.imaginary() - expected_h3[i].1).abs() >= 1e-10
            {
                println!(
                    "Mismatch at index {}: expected ({}, {}), got ({}, {})",
                    i,
                    expected_h3[i].0,
                    expected_h3[i].1,
                    root.real(),
                    root.imaginary()
                );
            }
            assert!((root.real() - expected_h3[i].0).abs() < 1e-10);
            assert!((root.imaginary() - expected_h3[i].1).abs() < 1e-10);
        }
    }
}
