use std::{fmt::Display, ops::{Add, Index, IndexMut, Mul}};
use num_traits;


// general 2d-matrix
#[derive(Debug, Clone)]
pub struct Matrix<T, const M: usize, const N: usize> {
    // M rows, N columns
    pub data: [[T; N]; M],
}


// column vector
pub type Vector<T, const M: usize> = Matrix<T, M, 1>;


// General rectangular matrices
impl<T, const M: usize, const N: usize> Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
{
    /// number of rows
    pub fn n_rows(&self) -> usize {
        self.data.len()
    }

    /// number of columns
    pub fn n_cols(&self) -> usize {
        let maybe_n_cols = self.data.first();
        if let Some(n_cols) = maybe_n_cols {
            return n_cols.len();
        }
        0
    }

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

    /// creates a matrix of size MxN with polynomial column values, starting at 0, so e.g. for poly = 2
    /// [0, 0, 1]
    /// [1, 1, 1]
    /// [4, 2, 1]
    /// [9, 3, 1]
    /// ...
    pub fn poly() -> Self
    where 
        T: From<u32>
    {
        let mut data = Self::ones();
        for i in 0..M {
            for j in 0..N {
                let x = i.pow(j as u32) as u32;
                data[i][j] = T::from(x);
            }
        }
        data
    }

    /// transposing is just swapping rows and columns, so
    /// AT_ji = A_ij
    pub fn transpose(&self) -> Matrix<T, N, M> {
        let mut ans = Matrix::<T, N, M>::zeros();
        // row
        for i in 0..M {
            // col
            for j in 0..N {
                ans[j][i] = self[i][j];
            }
        }
        ans
    }

    /// Add a column to the right end of the matrix
    /// A_mn + B_m1 = C_mn+1
    pub fn append_col(&self, col: &Vector<T, M>) -> Matrix<T, M, {N+1}> {
        let mut ans = Matrix::<T, M, {N+1}>::zeros();
        for i in 0..M {
            for j in 0..N {
                ans[i][j] = self[i][j];
            }
            ans[i][N] = col[i][0];
        }
        ans
    }


    /// Bring a matrix to reduced row echolon form
    pub fn rref(&self) -> Self {
        // Basically the same as U in LU decomposition but 0s are allowed in pivot positions
        let mut ans = Matrix::<T, N, M>::zeros();
        for i in 0..M {

        }
        todo!();
    }


    pub fn pinv(&self) -> Self {
        todo!()
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

    /// Performs the LU factorization on the matrix A to satisfy the equation PA = LU and returns a tuple
    /// 
    /// (P, L, U)
    ///  
    /// - P: Permutation matrix
    /// - L: Lower triangular matrix
    /// - U: Upper traingular matrix
    /// This only works for floating point type matrices.
    pub fn lu_decomposition(&self) -> Result<(Self, Self, Self), String>
    where T: num_traits::Float,
    {
        let mut p: Matrix<T, N, N> = Self::identity();
        let mut l = Self::zeros();
        let mut u = self.clone();
        

        // 1) Permutation - put the biggest element by absolute value in the upper most position as pivot (numeric stability)
        // 2) Elimination - row by row
        for i in 0..N {
            // find the row index where that column has the maximum value
            let maybe_row_max_index = u.data.iter().enumerate().skip(i) // only look at rows from i downwards
                .map(|(row_i, row)| (row_i, row[i])) // look at the i_th element in each row
                .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).expect("Unable to order rows"))
                .map(|(index, _)| index);
            
            // swap rows if necessary
            if let Some(row_max_index) = maybe_row_max_index {
                // if the largest element of column i isn't in row i, then swap the rows
                u.data.swap(row_max_index, i);
                l.data.swap(row_max_index, i);
                // calculate new permutation matrix
                let mut p_i: Matrix<T, N, N> = Self::identity();
                p_i.data.swap(row_max_index, i);
                p = &p_i * &p;
            }

            // eliminate rows
            let (top, bottom) = u.data.split_at_mut(i+1);
            let pivot_row = &top[i];
            let pivot_factor = pivot_row[i];

            for j in 0..bottom.len() {
                // Check if division by 0
                if pivot_factor == T::from(0).unwrap() {
                    return Err(String::from("Singular matrix cannot be factorized"));
                }
                let row_j = &mut bottom[j];
                let factor = row_j[i] / pivot_factor;
                // update L matrix: stores information
                // "how did I eliminate element [i][j]"
                l.data[i+1+j][i] = factor;
                // for every column in row j
                for k in i..N {
                    row_j[k] = row_j[k] - factor*pivot_row[k]
                }
            }

        }
        // add the ones for the main diagonal of L
        l = &l + &Self::identity();
        Ok((p, l, u))
    }


    pub fn inv(&self) -> Result<Self, String>
    where
        T: num_traits::Float,
    {
        // Use LU decomposition to calculate the inverse
        // A * A.inv() = I | with A = LU
        // LU * A.inv() = I | with y = U * A.inv()
        // Ly = I --> solve for each row to get y
        let mut y = Vector::<T, N>::zeros();
        let (p, l, u) = self.lu_decomposition()?;
        todo!();
    }
}



// Matrix from array
impl<T, const M: usize, const N: usize> From<[[T; N]; M]> for Matrix<T, M, N> {
    fn from(value: [[T; N]; M]) -> Self {
        Self { data: value }
    }
}



// Matrix + Matrix addition
impl<T, const M: usize, const N: usize> Add for &Matrix<T, M, N>
where 
    T: Copy + num_traits::Num
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


// Matrix x Matrix multiplication rectangular
impl<T, U, const M: usize, const N: usize, const P: usize> Mul<&Matrix<U, N, P>> for &Matrix<T, M, N>
where 
    T: Copy + num_traits::Num + From<U>,
    U: Copy + num_traits::Num,
{
    type Output = Matrix<T, M, P>;

    fn mul(self, rhs: &Matrix<U, N, P>) -> Self::Output {
        let mut ans = Matrix::<T, M, P>::zeros();
        for i in 0..M {
            for j in 0..P {
                ans.data[i][j] = (0..N)
                    .map(|k| self.data[i][k] * T::from(rhs.data[k][j]))
                    .fold(T::zero(), |acc, x| acc+x);
            }
        }
        ans
    }
}


// Matrix x Scalar multiplication
impl<T, U, const M: usize, const N: usize> Mul<U> for &Matrix<T, M, N>
where 
    T: Copy + From<U> + num_traits::Num,
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


// Index implementation
impl<T, const M: usize, const N: usize> Index<usize> for Matrix<T, M, N>
{
    type Output = [T; N];
    
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}
// IndexMut implementation
impl<T, const M: usize, const N: usize> IndexMut<usize> for Matrix<T, M, N>
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
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
        let _a: Matrix<f64, 3, 2> = Matrix::zeros();
        // ones
        let _b: Matrix<f64, 3, 2> = Matrix::ones();
        // identity
        let _i = Matrix::<i64, 7, 7>::identity();
        // from array
        let arr = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
        let _m = Matrix::from(arr);
    }

    #[test]
    fn print_matrix() {
        let a: Matrix<f64, 3, 2> = Matrix::ones();
        println!("{a}");
    }

    #[test]
    fn add_matrices() {
        let a: Matrix<usize, 5, 2> = Matrix::ones();
        let b: Matrix<usize, 5, 2> = Matrix::zeros();
        let _c = &a + &b;
    }

    #[test]
    fn mul_matrices() {
        let a: Matrix<f64, 3, 3> = Matrix::ones();
        let b: Matrix<f64, 3, 3> = Matrix::ones();
        let d: Matrix<f64, 3, 3> = Matrix::ones();
        let _c = &(&a * &b) * &d;
        let _e = &a + &b;

        let a: Matrix<f64, 3, 7> = Matrix::ones();
        let b: Matrix<f64, 7, 3> = Matrix::ones();
        let _c = &a * &b;
    }

    #[test]
    fn mul_matrix_scalar() {
        // Matrix
        let a: Matrix<f64, 3, 2> = Matrix::ones();
        let s: f32 = 64.3;
        let _b = &a * s;

        // Vector
        let v = Vector::<f64, 5>::ones();
        let _w = &v * s;
    }

    #[test]
    fn transpose() {
        let a: Matrix<f64, 5, 2> = Matrix::ones();
        let _b = a.transpose();
    }

    #[test]
    fn lu_factors() {
        let arr = [[1., 0., 7.4], [0., 1., 0.], [0., 0., 0.]];
        let a = Matrix::from(arr);
        let (_p, _l, _u) = a.lu_decomposition().unwrap();
    }

    #[test]
    fn index() {
        let mut a = Matrix::<f64, 5, 5>::ones();
        let _x = a[3][4]; // index borrow
        let y = &mut a[3][4]; // index mut
        *y = 4.0;
        println!("{a}"); // index borrow again
    }

    #[test]
    fn append_col() {
        let a = Matrix::<f64, 5, 5>::ones();
        let v = Vector::<f64, 5>::ones();
        let b = a.append_col(&v);
    }
}