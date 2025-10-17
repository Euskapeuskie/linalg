use core::num;
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


// General rectangular matrices
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

    /// transposing is just swapping rows and columns, so
    /// AT_ji = A_ij
    pub fn transpose(&self) -> Matrix<T, N, M> {
        let mut ans = Matrix::<T, N, M>::zeros();
        // row
        for i in 0..M {
            // col
            for j in 0..N {
                ans.data[j][i] = self.data[i][j];
            }
        }
        ans
    }
}


// Matrix from array
impl<T, const M: usize, const N: usize> From<[[T; N]; M]> for Matrix<T, M, N> {
    fn from(value: [[T; N]; M]) -> Self {
        Self { data: value }
    }
}


// Square matrices
impl<T, const N: usize> Matrix<T, N, N>
where
    T: Copy + num_traits::Num
{
    /// Create an Identity matrix of size NxN
    pub fn identity() -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = T::one();
        }
        Self { data: data }
    }


    /// Performs the LU factorization on the matrix and returns a tuple
    /// 
    /// (L, U)
    ///  
    /// This only works for floating point type matrices.
    pub fn lu_decomposition(&self) -> Matrix<T, N, N>
    where T: num_traits::Float + PartialOrd
    {
        let mut U = self.clone();
        let mut L = Self::identity();

        // 1) Permutation - put the biggest element by absolute value in the upper most position as pivot (numeric stability)
        // 2) Elimination - row by row
        for i in 0..N {
            // find the row index where that column has the maximum value
            let maybe_row_col_max = self.data.iter().enumerate().skip(i) // only look at data starting at row i
                .map(|(row_i, row)| (row_i, row[i])) // get row index and element in data[row][i]
                .fold(None, |acc, x| match acc {
                    None => Some(x),
                    Some(y) => match x.partial_cmp(&y) {
                        Some(std::cmp::Ordering::Greater) => Some(x),
                        Some(_) => Some(y),
                        None => None,
                    },
                }).map(|x| x.0); // we're only interested in the row index
            println!("{maybe_row_col_max:?}, {i}");
            
            // swap rows if necessary
            if let Some(row_col_max) = maybe_row_col_max {
                // if the largest element of column i isn't in row i, then swap the rows
                if row_col_max != i {
                    // permutation matrix
                    U.data.swap(row_col_max, i);
                }
            }

            // eliminate rows
            for j in i+1..N {
                let (top, bottom) = U.data.split_at_mut(j);
                let row_i = &top[i];
                let row_j = &mut bottom[0];

                let factor = row_j[0] / row_i[0];
                for k in i..N {
                    println!("k is {k}, j is {j}");
                    row_j[k] = row_j[k] - factor*row_i[k];
                }
            }    
        }
        U
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
        // identity
        let i = Matrix::<i64, 7, 7>::identity();
        // from array
        let arr = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        let m = Matrix::from(arr);
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

    #[test]
    fn transpose() {
        let a: Matrix<f64, 5, 2> = Matrix::ones();
        let b = a.transpose();
    }

    #[test]
    fn lu_factors() {
        let arr = [[1., 2.], [4., 7.]];
        let a = Matrix::from(arr);
        let u = a.lu_decomposition();
        println!("{a}");
        println!("{u}");
    }
}