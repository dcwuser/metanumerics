using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents the LU decomposition of a square matrix.
    /// </summary>
    /// <remarks><para>An LU decomposition is a representation of a matrix M as the product of a lower-left-triagular matrix L and
    /// and an upper-right-triangular matrix U. To reduce numerical instabilities, we actually decompose a row-wise
    /// permutation of a matrix, so that we have P A = L U, where P is a permutation matrix.</para>
    /// <para>For example, here is an LU decomposition of a permutation of a simple 3 X 3 matrix:</para>
    /// <img src="../images/LUDecomposition.png" />
    /// <para>Given an LU decomposition of a permutation A, we can solve systems of linear equations or compute the determinant or inverse of A.
    /// LU decomposition is the fastest way to solve an arbitrary system of linear equations. It is much faster,
    /// and less subject to rounding errors, to solve Ax=b by LU decomposition than than by inverting A and multiplying A<sup>-1</sup>b.</para>
    /// <para>You can use the <see cref="SquareMatrix.LUDecomposition"/> method to obtain the LU decomposition of any non-singular
    /// square matrix.</para>
    /// </remarks>
    /// <seealso cref="SquareMatrix"/>
    public sealed class LUDecomposition {

        internal LUDecomposition (double[] luStore, int[] permutation, int parity, int dimension) {
            Debug.Assert(luStore != null);
            Debug.Assert(permutation != null);
            Debug.Assert(Math.Abs(parity) == 1);
            this.luStore = luStore;
            this.permutation = permutation;
            this.parity = parity;
            this.dimension = dimension;
        }

        private readonly double[] luStore;
        private readonly int[] permutation;
        private readonly int parity;
        private readonly int dimension;

        /// <summary>
        /// Gets the dimension of the system.
        /// </summary>
        public int Dimension {
            get {
                return (dimension);
            }
        }

        /// <summary>
        /// Computes the determinant of the original matrix.
        /// </summary>
        /// <returns>The determinant of the original matrix.</returns>
        public double Determinant () {
            
            double det = parity;
            for (int i = 0; i < dimension; i++) {
                det = det * luStore[dimension * i + i];
            }
            return (det);
        }


        /// <summary>
        /// Solves A x = b.
        /// </summary>
        /// <param name="rhs">The right-hand side vector b.</param>
        /// <returns>The solution vector x.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="rhs"/> is <c>null</c>.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="rhs"/> is not the same as the dimension of the matrix.</exception>
        public ColumnVector Solve (IReadOnlyList<double> rhs) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (rhs.Count != dimension) throw new DimensionMismatchException();

            // copy rhs into an array, reshuffling according to permutations
            double[] x = new double[dimension];
            for (int i = 0; i < x.Length; i++) {
                x[i] = rhs[permutation[i]];
            }

            // solve Ly=b and Ux=y
            Blas2.dTrsv(false, true, luStore, 0, 1, dimension, x, 0, 1, dimension);
            //SquareMatrixAlgorithms.SolveLowerLeftTriangular(luStore, x, 0, dimension);
            Blas2.dTrsv(true, false, luStore, 0, 1, dimension, x, 0, 1, dimension);
            //SquareMatrixAlgorithms.SolveUpperRightTriangular(luStore, x, 0, dimension);

            return (new ColumnVector(x, dimension));

        }

        /// <summary>
        /// Computes the the inverse of the original matrix.
        /// </summary>
        /// <returns>The independent inverse of the original matrix.</returns>
        public SquareMatrix Inverse () {

            // set up permuted unit matrix
            double[] aiStore = new double[dimension * dimension];
            for (int i = 0; i < dimension; i++) {
                aiStore[dimension * permutation[i] + i] = 1.0;
            }

            // solve each column
            for (int c = 0; c < dimension; c++) {
                int i0 = dimension * c;
                SquareMatrixAlgorithms.SolveLowerLeftTriangular(luStore, aiStore, i0, dimension);
                SquareMatrixAlgorithms.SolveUpperRightTriangular(luStore, aiStore, i0, dimension);

            }

            return (new SquareMatrix(aiStore, dimension));
            
            /*
            // this is basically just inserting unit vectors into Solve,
            // but we reproduce that logic so as to take some advantage
            // of the fact that most of the components are zero

            int n = Dimension;

            SquareMatrix MI = new SquareMatrix(n);

            // iterate over the columns
            for (int c = 0; c < n; c++) {

                // we are dealing with the k'th unit vector
                // the c'th unit vector gets mapped to the k'th unit vector, where perm[k] = c

                // solve L y = e_k
                // components below the k'th are zero; this loop determines k as it zeros lower components
                int k = 0;
                while (perm[k] != c) {
                    MI.SetEntry(k, c, 0.0);
                    k++;
                }
                // the k'th component is one
                MI[k, c] = 1.0;
                // higher components are non-zero...
                for (int i = k + 1; i < n; i++) {
                    double t = 0.0;
                    // ...but loop to compute them need only go over j for which MI[j, c] != 0, i.e. j >= k
                    for (int j = k; j < i; j++) {
                        t -= lu.GetEntry(i, j) * MI.GetEntry(j, c);
                    }
                    MI.SetEntry(i, c, t);
                }

                // solve U x = y
                // at this stage, the first few components of y are zero and component k is one
                // but this doesn't let us limit the loop any further, since we already only loop over j > i 
                for (int i = n - 1; i >= 0; i--) {
                    double t = MI.GetEntry(i, c);
                    for (int j = n - 1; j > i; j--) {
                        t -= lu.GetEntry(i, j) * MI.GetEntry(j, c);
                    }
                    MI.SetEntry(i, c, t / lu.GetEntry(i, i));
                }

            }

            return (MI);
            */
        }

        /// <summary>
        /// Gets the L factor.
        /// </summary>
        /// <returns>The lower-left triangular factor L of the LU decomposition.</returns>
        /// <remarks>
        /// <para>The pivoted LU decomposition algorithm guarantees that the diagonal entries of this matrix are all one, and
        /// that the magnitudes of the sub-diagonal entries are all less than or equal to one.</para>
        /// </remarks>
        public SquareMatrix LMatrix () {
            double[] lStore = new double[dimension * dimension];
            for (int c = 0; c < dimension; c++) {
                int i0 = dimension * c + c;
                lStore[i0] = 1.0;
                Blas1.dCopy(luStore, i0 + 1, 1, lStore, i0 + 1, 1, dimension - c - 1);
            }
            return (new SquareMatrix(lStore, dimension));
        }

        /// <summary>
        /// Gets the U factor.
        /// </summary>
        /// <returns>The upper-right triangular factor U of the LU decomposition.</returns>
        public SquareMatrix UMatrix () {
            double[] uStore = new double[dimension * dimension];
            for (int c = 0; c < dimension; c++) {
                Blas1.dCopy(luStore, dimension * c, 1, uStore, dimension * c, 1, c + 1);
            }
            return (new SquareMatrix(uStore, dimension));
        }

        /// <summary>
        /// Gets the permutation matrix.
        /// </summary>
        /// <returns>The permutation matrix P in the PA = LU decomposition.</returns>
        /// <remarks>
        /// <para>A permutation matrix is just a "scrambled" identity matrix: 1 appears exactly once in each row and column, but
        /// not necessarily in the diagonal position.</para>
        /// </remarks>
        public SquareMatrix PMatrix () {
            double[] pStore = new double[dimension * dimension];
            for (int i = 0; i < dimension; i++) {
                pStore[dimension * permutation[i] + i] = 1.0;
            }
            return (new SquareMatrix(pStore, dimension));
        }

    }

}
