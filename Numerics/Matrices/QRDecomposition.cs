using System;
using System.Collections.Generic;
using System.Diagnostics;


namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a QR decomposition of a matrix.
    /// </summary>
    /// <remarks>
    /// <para>A QR decomposition represents a rectangular matrix as the product of a square, orthogonal matrix Q and a rectangular,
    /// right-upper-triangular matrix R. For example:</para>
    /// <img src="../images/QRDecomposition.png" />
    /// <para>The method <see cref="RectangularMatrix.QRDecomposition"/> of the <see cref="RectangularMatrix"/> class can be used to QR decompose a
    /// rectangular matrix.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/QR_decomposition"/>
    public sealed class QRDecomposition {

        private readonly double[] qtStore;
        private readonly double[] rStore;
        private readonly int rows, cols;

        internal QRDecomposition (double[] qtStore, double[] rStore, int rows, int columns) {
            Debug.Assert(qtStore != null);
            Debug.Assert(rStore != null);
            Debug.Assert(rows > 0);
            Debug.Assert(columns > 0);
            this.qtStore = qtStore;
            this.rStore = rStore;
            this.rows = rows;
            this.cols = columns;
        }

        /// <summary>
        /// The orthogonal matrix Q.
        /// </summary>
        /// <value>The orthogonal matrix Q in the decomposition M = Q R.</value>
        /// <remarks>
        /// <para>The returned matrix is read-only. If you need to make changes to it, you can call <see cref="SquareMatrix.Copy"/> to obtain a writable copy.</para>
        /// </remarks> 
        public SquareMatrix QMatrix {
            get {
                return (new SquareMatrix(qtStore, 0, rows, 1, rows, true));
            }
        }

        /// <summary>
        /// The upper-right triangular matrix R.
        /// </summary>
        /// <value>The upper-right triangular matrix R in the decomposition M = Q R.</value>
        /// <remarks>
        /// <para>The returned matrix is read-only. If you need to make changes to it, you can call <see cref="RectangularMatrix.Copy"/> to obtain a writable copy.</para>
        /// </remarks> 
        public RectangularMatrix RMatrix {
            get {
                return (new RectangularMatrix(rStore, rows, cols, true));
            }
        }


        /// <summary>
        /// Solve the system A x = b.
        /// </summary>
        /// <param name="rhs">The right-hand-side b.</param>
        /// <returns>The column vector x for which A x is closest to b.</returns>
        public ColumnVector Solve (IReadOnlyList<double> rhs) {

            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (rhs.Count != rows) throw new DimensionMismatchException();

            // BLAS requires an array, but doesn't modify it, so if given one, use it directly
            double[] x = rhs as double[];
            if (x == null) {
                x = new double[rows];
                rhs.CopyTo(x, 0);
            }

            // Q^T x is a row-length vector, but we only need the first cols entries
            // so truncate Q^T to cols X rows, so that Q^T x is only of length cols
            double[] y = new double[cols];
            Blas2.dGemv(qtStore, 0, 1, rows, x, 0, 1, y, 0, 1, cols, rows);
            MatrixAlgorithms.SolveUpperRightTriangular(rStore, rows, cols, y, 0);

            return (new ColumnVector(y, cols));
        }
        

        /// <summary>
        /// Get the number of rows in the original matrix.
        /// </summary>
        public int RowCount {
            get {
                return (rows);
            }
        }

        /// <summary>
        /// Gets the number of columns in the original matrix.
        /// </summary>
        public int ColumnCount {
            get {
                return (cols);
            }
        }

        // Solve a linear system using QR decomposition. This isn't really the obvious place for this algorithm, but since
        // it's used by both BivariableSample and UncertainMeasurementSample, it needs to be centralized somewhere and here
        // is as good a place as any. Probably it should be moved in to RectangularMatrixAlgorithms and it should
        // work directly on the underlying array storage for increased performance.

        internal static void SolveLinearSystem (RectangularMatrix A, ColumnVector b, out ColumnVector a, out SymmetricMatrix C) {

            Debug.Assert(A.RowCount == b.Dimension);
            Debug.Assert(A.ColumnCount <= A.RowCount);

            // Our problem is to "solve" A a = b, where a is the vector of coefficients. QR gives the least squares
            // solution that minimizes | A a - y |.
            QRDecomposition QR = A.QRDecomposition();
            a = QR.Solve(b);

            // The covariance matrix C = (A^T A)^{-1} = (R^T R)^{-1}. This is a matrix multiplication plus an inversion.

            // A faster way is to use the direct solution
            //   for k = n ... 1
            //     C_{kk} = R_{kk}^{-1} \left[ R_{kk}^{-1} - \sum_{j=k+1}^{n} R_{kj} C_{kj} \right]
            //     for i = k-1 ... 1
            //       C_{ik} = -R_{ii}^{-1} \left[ \sum_{j=i+1}^{k} R_{ij} C_{jk} + \sum_{j=k+1}^{n} R_{ij} C_{kj} \right]
            //     end
            //   end
            // This is detailed in Ake Bjorck, "Numerical Methods for Least Squares Problems", pp. 118-120

            C = new SymmetricMatrix(a.Dimension);
            RectangularMatrix R = QR.RMatrix;
            double c; // used for storage so we don't call bounds-checking accessors each time
            for (int k = C.Dimension - 1; k >= 0; k--) {
                c = 1.0 / R[k, k]; ;
                for (int j = k + 1; j < C.Dimension; j++) {
                    c -= R[k, j] * C[k, j];
                }
                C[k, k] = c / R[k, k];
                for (int i = k - 1; i >= 0; i--) {
                    c = 0.0;
                    for (int j = i + 1; j <= k; j++) {
                        c += R[i, j] * C[j, k];
                    }
                    for (int j = k + 1; j < C.Dimension; j++) {
                        c += R[i, j] * C[k, j];
                    }
                    C[i, k] = -c / R[i, i];
                }
            }
        }

    }

}
