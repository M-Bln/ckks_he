use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

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
{
    type Real;
    fn new(real: Self::Real, imaginary: Self::Real) -> Self;
    fn real(&self) -> Self::Real;
    fn imaginary(&self) -> Self::Real;
    fn magnitude(&self) -> Self::Real;
    fn phase(&self) -> Self::Real;
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
}
