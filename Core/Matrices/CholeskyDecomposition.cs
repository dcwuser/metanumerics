using System;
using System.Collections.Generic;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents the Cholesky Decomposition of a symmetric, positive definite matrix. 
    /// </summary>
    /// <remarks>
    /// <para>A Cholesky decomposition represents a matrix as the product of a lower-left triangular matrix and its transpose. For example:</para>
    /// <img src="../images/CholeskyDecomposition.png" />
    /// <para>The Choleksy decomposition of a symmetric, positive definite matrix can be obtained using the
    /// <see cref="SymmetricMatrix.CholeskyDecomposition"/> method of the <see cref="SymmetricMatrix"/> class.</para>
    /// </remarks>
    /// <example>
    /// <para>Here is an example that uses a Cholesky decomposition to solve a linear algebra problem.</para>
    /// <code lang="cs">
    /// // Solve Ax = b via Cholesky decomposition
    /// CholeskyDecomposition CD = A.CholsekyDecomposition();
    /// ColumnVector b = new ColumnVector(1.0, 2.0, 3.0);
    /// ColumnVector x CD.Solve(b);
    /// </code>
    /// </example>
    /// <seealso cref="SymmetricMatrix.CholeskyDecomposition"/>
    public sealed class CholeskyDecomposition {

        internal SymmetricMatrix sqrtM;

        /// <summary>
        /// Gets the dimension of the system.
        /// </summary>
        public int Dimension {
            get {
                return (sqrtM.Dimension);
            }
        }

        /// <summary>
        /// Returns the Cholesky square root matrix.
        /// </summary>
        /// <returns>A lower-left triangular matrix A, such that A A<sup>T</sup> = M.</returns>
        public SquareMatrix SquareRootMatrix () {
            SquareMatrix A = new SquareMatrix(Dimension);
            for (int r = 0; r < Dimension; r++) {
                for (int c = 0; c <= r; c++) {
                    A[r, c] = sqrtM[r, c];
                }
            }
            return (A);
        }

        /*
        public virtual LowerTriangularMatrix LeftFactor {
            get {
                LowerTriangularMatrix L = new LowerTriangularMatrix(Dimension);
                for (int r = 0; r < Dimension; r++) {
                    for (int c = 0; c <= r; c++) {
                        L[r, c] = sqrtM[r, c];
                    }
                }
                return (L);
            }
        }

        public virtual UpperTriangularMatrix RightFactor {
            get {
                UpperTriangularMatrix R = new UpperTriangularMatrix(Dimension);
                for (int c = 0; c < Dimension; c++) {
                    for (int r = 0; r <= c; r++) {
                        R[r, c] = sqrtM[r, c];
                    }
                }
                return (R);
            }
        }
        */

        /// <summary>
        /// Computes the solution vector that, when multiplied by the original matrix, produces the given left-hand side vector.
        /// </summary>
        /// <param name="rhs">The right-hand-side vector.</param>
        /// <returns>The left-hand-side (solution) vector.</returns>
        public ColumnVector Solve (IList<double> rhs) {
            if (rhs == null) throw new ArgumentNullException("rhs");
            if (rhs.Count != Dimension) throw new DimensionMismatchException();

            // Determine Ly = x
            double[] y = new double[Dimension];
            for (int i = 0; i < Dimension; i++) {
                y[i] = rhs[i];
                for (int j = 0; j < i; j++) {
                    y[i] -= sqrtM[j, i] * y[j];
                }
                y[i] = y[i] / sqrtM[i, i];
            }

            // Determine L^{T} z = y
            double[] z = new double[Dimension];

            // Re-use y for z, since we don't need a value after it's set
            for (int i = (Dimension - 1); i >= 0; i--) {
                z[i] = y[i];
                for (int j = (Dimension - 1); j > i; j--) {
                    z[i] -= sqrtM[i, j] * z[j];
                }
                z[i] = z[i] / sqrtM[i, i];
            }
            return (new ColumnVector(z, Dimension));
        }

        /// <summary>
        /// Computes the inverse of the original matrix.
        /// </summary>
        /// <returns>M<sup>-1</sup></returns>
        public SymmetricMatrix Inverse () {
            SymmetricMatrix MI = new SymmetricMatrix(Dimension);

            // do each column as a RHS
            for (int c = 0; c < Dimension; c++) {
                MI[c, c] = 1.0 / sqrtM[c, c];
                for (int r = c + 1; r < Dimension; r++) {
                    MI[r, c] = 0.0;
                    for (int i = c; i < r; i++) {
                        MI[r, c] -= sqrtM[r, i] * MI[i, c];
                    }
                    MI[r, c] = MI[r, c] / sqrtM[r, r];
                }
            }

            // unnecessary?
            for (int c = 0; c < Dimension; c++) {
                for (int r = (Dimension - 1); r >= c; r--) {
                    for (int i = r + 1; i < Dimension; i++) {
                        MI[r, c] -= sqrtM[r, i] * MI[i, c];
                    }
                    MI[r, c] = MI[r, c] / sqrtM[r, r];
                }
            }
            // end unnecessary?

            return (MI);
        }

        /// <summary>
        /// Computes the determinant of the original matrix.
        /// </summary>
        /// <returns>det M</returns>
        public double Determinant () {
            double lnDet = 0.0;
            for (int i = 0; i < Dimension; i++) {
                lnDet += 2.0 * Math.Log(sqrtM[i, i]);
            }
            return (Math.Exp(lnDet));
        }

        internal CholeskyDecomposition (SymmetricMatrix sqrtM) {
            this.sqrtM = sqrtM;
        }

    }

}