use std::{fmt::Display, ops::{Add, Sub, Index, IndexMut, Mul, RangeTo}};
use crate::types::Matrix;



// Matrix from array
impl<T, const M: usize, const N: usize> From<[[T; N]; M]> for Matrix<T, M, N> {
    fn from(value: [[T; N]; M]) -> Self {
        Self { data: value }
    }
}


// Matrix from slice
impl<T, const M: usize, const N: usize> From<&[T]> for Matrix<T, M, N> 
where 
    T: Copy + num_traits::Num
{
    fn from(value: &[T]) -> Self {
        assert_eq!(value.len(), N*M, "cannot build a matrix of size {M}x{N} out of a {} slice", value.len());

        let mut ans = [[T::zero();  N]; M];
        for i in 0..value.len() {
            let row_index = i / N;
            let col_index = i % N;
            ans[row_index][col_index] = value[i];
        }
        let ans = Matrix::from(ans);
        ans
    }
}


// Matrix from lambda function
impl<T, F, const M: usize, const N: usize> From<F> for Matrix<T, M, N>
where
    F: Fn(usize) -> [T; N],
    T: Copy + num_traits::Num,
{
    fn from(value: F) -> Self {
        let mut ans = Self::zeros();
        for i in 0..M {
            ans[i] = value(i);
        }
        ans
    }
}


// Matrix from function and array with specified points
impl<T, F, const M: usize, const N: usize> From<(F, [T; M])> for Matrix<T, M, N>
where
    F: Fn(T) -> [T; N],
    T: Copy + num_traits::Num, 
{
    fn from(value: (F, [T; M])) -> Self {
        let(f, ts ) = value;
        let mut ans = Self::zeros();

        for i in 0..M {
            ans[i] = f(ts[i]);
        }
        ans
    }
}


// Matrix from function and step_size
impl<T, F, const M: usize, const N: usize> From<(F, T)> for Matrix<T, M, N>
where
    F: Fn(T) -> [T; N],
    T: Copy + num_traits::Float,
{
    fn from(value: (F, T)) -> Self {
        let (f, stepsize) = value;
        let mut ans = Self::zeros();

        for i in 0..M {
            let x = T::from(i).unwrap() * stepsize;
            ans[i] = f(x);
        }
        ans
    }
}


// Matrix + Matrix addition
impl<T, const M: usize, const N: usize> Add for Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    type Output = Matrix<T, M, N>;

    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

// &Matrix + Matrix
impl<T, const M: usize, const N: usize> Add<Matrix<T, M, N>> for &Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    type Output = Matrix<T, M, N>;

    fn add(self, rhs: Matrix<T, M, N>) -> Self::Output {
        self + &rhs
    }
}


// Matrix + &Matrix
impl<T, const M: usize, const N: usize> Add<&Matrix<T, M, N>> for Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    type Output = Matrix<T, M, N>;

    fn add(self, rhs: &Matrix<T, M, N>) -> Self::Output {
        &self + rhs
    }
}


// &Matrix + &Matrix
impl<'a, 'b, T, const M: usize, const N: usize> Add<&'b Matrix<T, M, N>> for &'a Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    type Output = Matrix<T, M, N>;

    fn add(self, rhs: &'b Matrix<T, M, N>) -> Self::Output {

        let mut ans = Matrix::zeros();
        // for each row
        for i in 0..M {
            // for each column
            for j in 0..N {
                ans[i][j] = self[i][j] + rhs[i][j];
            }
        }
        ans
    }
}


// Matrix + Matrix subtraction
impl<T, const M: usize, const N: usize> Sub for &Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    type Output = Matrix<T, M, N>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut ans = Matrix::zeros();
        for i in 0..M {
            for j in 0..N {
                ans[i][j] = self[i][j] - rhs[i][j];
            }
        }
        ans
    }
}


// Matrix x Matrix multiplication rectangular
impl<T, const M: usize, const N: usize, const P: usize> Mul<&Matrix<T, N, P>> for &Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    type Output = Matrix<T, M, P>;

    fn mul(self, rhs: &Matrix<T, N, P>) -> Self::Output {
        let mut ans = Matrix::<T, M, P>::zeros();
        for i in 0..M {
            for j in 0..P {
                ans.data[i][j] = (0..N)
                    .map(|k| self.data[i][k] * rhs.data[k][j])
                    .fold(T::zero(), |acc, x| acc+x);
            }
        }
        ans
    }
}


// Matrix x Scalar multiplication
impl<T, const M: usize, const N: usize> Mul<T> for &Matrix<T, M, N>
where 
    T: Copy + num_traits::Num,
{
    type Output = Matrix<T, M, N>;

    fn mul(self, rhs: T) -> Self::Output {
        let mut ans = self.clone();
        for i in 0..M {
            for j in 0..N {
                ans.data[i][j] = rhs * ans.data[i][j];
            }
        }
        ans
    }
}


// Index single row
impl<T, const M: usize, const N: usize> Index<usize> for Matrix<T, M, N>
{
    type Output = [T; N];
    
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}
// IndexMut single row
impl<T, const M: usize, const N: usize> IndexMut<usize> for Matrix<T, M, N>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}


// Index multiple rows
// Rust guarantees that an array is contiguous in memory
impl<T, const M: usize, const N: usize> Index<RangeTo<usize>> for Matrix<T, M, N>
{
    type Output = [T];

    fn index(&self, index: RangeTo<usize>) -> &Self::Output {
        // construct &[T]: start_address, n_elements
        // start_address = start of first row
        // n_elements = index.end*n_columns

        let start_address = self.data.as_ptr() as *const T; // start address
        let n_elements = index.end * N; // n_elements

        assert!(n_elements <= M*N, "range index bigger than size of array"); // bounds check
        // build and return the pointer
        unsafe {
            std::slice::from_raw_parts(start_address, n_elements)
        }
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