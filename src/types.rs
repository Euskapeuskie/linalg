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


    /// As slice
    pub fn as_slice(&self) -> &[T] {
        &self[..M]
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


    /// Bring a matrix to row echelon form
    pub fn row_echelon(&self) -> Self 
    where
        T: num_traits::Float

    {
        let mut u = self.clone();
        let max_iter = M.min(N);

        // 1) Permutation - put the biggest element by absolute value in the upper most position as pivot (numeric stability)
        // 2) Elimination - row by row
        for i in 0..max_iter {
            // find the row index where that column has the maximum value
            let maybe_row_max_index = u.data.iter().enumerate().skip(i) // only look at rows from i downwards
                .map(|(row_i, row)| (row_i, row[i])) // look at the i_th element in each row
                .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).expect("Unable to order rows"))
                .map(|(index, _)| index);

            // swap rows if necessary
            if let Some(row_max_index) = maybe_row_max_index {
                // if the largest element of column i isn't in row i, then swap the rows
                u.data.swap(row_max_index, i);
            }

            // eliminate rows
            let (top, bottom) = u.data.split_at_mut(i+1);
            let pivot_row = &top[i];
            let pivot_factor = pivot_row[i];

            // Check if pivot is a 0 -> if yes go to the next row
            if pivot_factor.is_zero() {
                continue
            }

            for j in 0..bottom.len() {
                let row_j = &mut bottom[j];
                let factor = row_j[i] / pivot_factor;
                for k in i..N {
                    row_j[k] = row_j[k] - factor*pivot_row[k]
                }
            }
        }
        u
    }


    /// Compute the QR decomposition of the matrix that satisfies A = QR using Gram-Schmidt
    /// - Q: Orthonomormal basis of the column space of A
    /// - R: Upper triangular matrix
    /// Returns (q, r)
    /// 
    /// Internal workings:
    /// 1) choose the first non-zero column as my initial basis
    /// 2) calculate for each following column...
    pub fn qr_decomposition(&self) -> (Matrix<T, N, M>, Self)
    where 
        T: num_traits::Float
    {

        // transpose columns into rows just more convenient to work with
        let mut q = Self::zeros().transpose();
        let mut r = Self::zeros();

        // Build the Q-Matrix
        for i in 0..self.n_cols() {
            q[i] = {
                // b is the vector we want to orthonormalize to q
                let mut b = Vector::from(self.transpose()[i].map(|x| [x]));

                // subtract the components of b in the other orthonormal directions to get only the new component of b
                for j in 0..i {
                    let a = Vector::from(q[j].map(|x| [x]));
                    b  = &b - &(&a * (((&a.transpose() * &b)[0][0]) / ((&a.transpose() * &a)[0][0])));
                }
                // normalize b to unit length
                b = &b * (T::one()/b.magnitude());
                b.to_array()
            };
        }

        // Build the R-Matrix
        for i in 0..self.n_cols() {
            for j in 0..self.n_rows() {
                // r_ij = (a_j)T*q_i --> column vector j of A times row vector i of Q
                let a = Vector::from(self.transpose()[j].map(|x| [x]));
                let q = Vector::from(q[i].map(|x| [x]));
                r[i][j] = (&a.transpose() * &q)[0][0];
            }
        }

        (q, r)
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
            // Check if division by 0
            if pivot_factor.is_zero() {
                return Err(String::from("Singular matrix cannot be factorized"));
            }

            for j in 0..bottom.len() {
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


    /// Returns the inverse of the matrix or an error if the matrix is singular
    /// uses LU-factorization for calculating the inverse
    /// A * A.inv() = I | with PA = LU
    /// LU * A.inv() = P*I | with Y = U * A.inv()
    /// LY = P -> solve for Y
    /// Resubstitution: U*A.inv() = Y -> solve for A.inv()
    pub fn inv(&self) -> Result<Self, String>
    where
        T: num_traits::Float,
    {
        let (p, l, u) = self.lu_decomposition()?;

        // 1) Solve for LY = P
        // L is a lower triangular matrix with only 1s on the main diagonal!
        let mut y = Self::zeros();
        // i: row we're acting on of Y
        for i in 0..N {
            // j: column of Y
            for j in 0..N {
                let mut sum = T::zero();
                for k in 0..i {
                    sum = sum + l[i][k] * y[k][j];
                }
                y[i][j] = (p[i][j] - sum) / l[i][i] // l[i][i] is always 1
            }
        }

        // Solve for U*A.inv() = Y
        // U is a upper triangular matrix, so start with last row
        let mut a_inv = Self::zeros();
        // start with last row because U = [0, 0, ...., x]
        for i in (0..N).rev() {
            for j in 0..N {
                let mut ans = T::zero();
                for k in i+1..N {
                    ans = ans + u[i][k] * a_inv[k][j];
                }
                a_inv[i][j] = (y[i][j] - ans) / u[i][i];
            }
        }
        Ok(a_inv)
    }


    /// Calculates the determinant of a square matrix
    /// 1) Bring matrix in upper triangular form (keep track of row exchanges)
    /// 2) Determinant = (-1)^n_row_changes * product of pivots
    /// TODO: (-1)^n_row_changes!
    pub fn det(&self) -> T 
    where
        T: num_traits::Float
    {
        todo!();
        let u = self.row_echelon();
        let mut ans = T::one();
        for i in 0..N {
            ans = ans * u[i][i];
        }
        ans
    }
}



// Vectors
impl<T, const M: usize> Vector<T, M>
where
    T: Copy + num_traits::Num
{

    /// Returns the magnitude or length of a vector:
    /// sqrt(x_transpose * x)
    pub fn magnitude(&self) -> T
    where
        T: num_traits::Float,
    {
        let x = &(self.transpose()) * self;
        x[0][0].sqrt()
    }


    /// Vector into Array
    pub fn to_array(self) -> [T; M] {
        self.transpose()[0]
    }
}




// TESTS
#[cfg(test)]
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

        let _f = Matrix::<_, 10, _>::from(|x: usize| [1, x, x.pow(2)]);
        println!("{_f}");
        let ts = [0.1, 0.2, 0.34];
        let _ft = Matrix::<f64, _, _>::from((|x: f64| [1.0, x, x.powf(0.3)], ts));
        println!("{_ft}");

        let _fstep = Matrix::<f64, 10, _>::from((|x: f64| [x.powi(2)], 0.1));
        println!("{_fstep}");
    }

    #[test]
    fn print_matrix() {
        let a: Matrix<f64, 3, 2> = Matrix::ones();
        println!("{a}");
    }

    #[test]
    fn add_sub_matrices() {
        let a: Matrix<usize, 5, 2> = Matrix::ones();
        let b: Matrix<usize, 5, 2> = Matrix::zeros();
        let _c = &a + &b;
        let _d = &a - &b;
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
        let arr = [[1., 0., 7.4], [0., 1., 0.], [1., 0., 0.]];
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
    fn inverse() {
        let a = [
            [1., 1., 1.],
            [2., 3., 3.],
            [3., 4., 5.],
        ];
        let a = Matrix::from(a);
        let inv = a.inv().unwrap();
        println!("Inverse:");
        println!("{inv}");
    }

    #[test]
    fn row_echelon() {
        let a = Matrix::<f64, 4, 4>::identity();
        let t = a.row_echelon();
        println!("{t}");
        let det = a.det();
        println!("{det}");
    }

    #[test]
    fn qr_decomposition() {
        let a = [
            [1., 2., 4.],
            [0., 0., 5.],
            [0., 3., 6.],
        ];
        let a = Matrix::from(a);
        let (q, r) = a.qr_decomposition();
        println!("{q}");
        println!("{r}");
    }


    #[test]
    fn vector() {
        let a = [1, 2, 3, 4];
        let a = Vector::from(a.map(|x| [x]));
        let a = &a[..4];
        println!("{a:?}");
    }
}