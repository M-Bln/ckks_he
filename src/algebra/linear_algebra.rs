use crate::algebra::complex::{Complex, C64};

pub fn identity_matrix(n: usize) -> Vec<Vec<C64>> {
    let mut matrix = vec![vec![C64::new(0.0, 0.0); n]; n];
    for i in 0..n {
        matrix[i][i] = C64::new(1.0, 0.0);
    }
    matrix
}

pub fn multiply_matrices(a: &[Vec<C64>], b: &[Vec<C64>]) -> Vec<Vec<C64>> {
    let n = a.len();
    let mut result = vec![vec![C64::new(0.0, 0.0); n]; n];

    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                result[i][j] = result[i][j] + a[i][k] * b[k][j];
            }
        }
    }

    result
}

pub fn apply_matrix(a: &[Vec<C64>], b: &[C64]) -> Vec<C64> {
    let n = a.len();
    let mut result = vec![C64::new(0.0, 0.0); n];

    for i in 0..n {
        for j in 0..n {
	    result[i] = result[i] + a[i][j] * b[j];
        }
    }

    result
}
