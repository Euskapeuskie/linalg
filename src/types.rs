use std::{fmt::Display, iter::Sum, ops::{Add, Mul}};
use num_traits;


// general 2d-matrix
#[derive(Debug, Clone)]
pub struct Matrix<T, const M: usize, const N: usize> {
    // M rows, N columns
    pub data: [[T; N]; M],
}

// column vector
pub type Vector<T, const N: usize> = Matrix<T, N, 1>;


impl<T, const M: usize, const N: usize> Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    /// create a matrix of size MxN with just 0s
    pub fn zeros() -> Self {
        let data = [[T::zero(); N]; M];
        Matrix { data: data }
    }
    /// create a matrix of size MxN with just 1s
    pub fn ones() -> Self {
        let data = [[T::one(); N]; M];
        Matrix { data: data }
    }

    pub fn transpose(&self) -> Matrix<T, N, M> {
        todo!()
    }
}


// Matrix + Matrix addition
impl<T, const M: usize, const N: usize> Add for &Matrix<T, M, N>
where 
    T: Copy + Add<Output = T> + num_traits::Num
{
    type Output = Matrix<T, M, N>;

    fn add(self, rhs: Self) -> Self::Output {

        let mut ans = Matrix::zeros();
        // for each row
        for m in 0..M {
            // for each column
            for n in 0..N {
                ans.data[m][n] = self.data[m][n] + rhs.data[m][n];
            }
        }
        ans
    }
}


// Matrix x Matrix multiplication
impl<T, U, const M: usize, const N: usize, const P: usize> Mul<&Matrix<U, N, P>> for &Matrix<T, M, N>
where 
    T: Copy + From<U> + num_traits::Num + Sum,
    U: Copy,
{
    type Output = Matrix<T, M, P>;

    fn mul(self, rhs: &Matrix<U, N, P>) -> Self::Output {
        let mut ans = Matrix::<T, M, P>::zeros();
        for i in 0..M {
            for j in 0..P {
                ans.data[i][j] = (0..N)
                    .map(|k| self.data[i][k] * T::from(rhs.data[k][j]))
                    .sum();
            }
        }
        ans
    }
}


// Matrix x Scalar multiplication
impl<U, T, const M: usize, const N: usize> Mul<U> for &Matrix<T, M, N>
where 
    T: Copy + Mul<T, Output = T> + From<U>,
    U: Copy + num_traits::Num
{
    type Output = Matrix<T, M, N>;

    fn mul(self, rhs: U) -> Self::Output {
        let mut ans = self.clone();
        for i in 0..M {
            for j in 0..N {
                ans.data[i][j] = T::from(rhs) * ans.data[i][j];
            }
        }
        ans
    }
}


// Display implementation (println!)
impl<T, const M: usize, const N: usize> Display for Matrix<T, M, N>
where
    T: Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for row in self.data.iter() {
            write!(f, "[")?;
            for (i, col) in row.iter().enumerate() {
                if i == row.len()-1 {
                    write!(f, "{col}")?;
                }
                else {
                    write!(f, "{col}, ")?;
                }
            }
            write!(f, "]\n")?;
        }
        Ok(())
    }
}




// TESTS
mod test {
    use crate::types::Matrix;
    use crate::types::Vector;

    #[test]
    fn create_matrix() {
        // zeros
        let a: Matrix<f64, 3, 2> = Matrix::zeros();
        // ones
        let b: Matrix<f64, 3, 2> = Matrix::ones();
        // for array
    }

    #[test]
    fn print_matrix() {
        let a: Matrix<f64, 3, 2> = Matrix::ones();
        println!("{a}");
        let b: Matrix<usize, 10, 10> = Matrix::zeros();
        println!("{b}");
    }

    #[test]
    fn add_matrices() {
        let a: Matrix<f64, 3, 2> = Matrix::ones();
        let b: Matrix<usize, 10, 10> = Matrix::zeros();
    }

    #[test]
    fn sub_matrices() {
        todo!();
    }

    #[test]
    fn mul_matrices() {
        let a: Matrix<f64, 3, 3> = Matrix::ones();
        let b: Matrix<f64, 3, 3> = Matrix::ones();
        let d: Matrix<u32, 3, 3> = Matrix::ones();
        let c = &(&a * &b) * &d;
        let e = &a + &b;

        let a: Matrix<f64, 3, 7> = Matrix::ones();
        let b: Matrix<f32, 7, 3> = Matrix::ones();
        let c = &a * &b;
    }

    #[test]
    fn mul_matrix_scalar() {
        // Matrix
        let a: Matrix<f64, 3, 2> = Matrix::ones();
        let s: f32 = 64.3;
        let b = &a * s;

        // Vector
        let v = Vector::<f64, 5>::ones();
        let w = &v * s;
    }
}