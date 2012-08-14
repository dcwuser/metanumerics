using System;
using System.Collections.Generic;
using System.Text;

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
    public sealed class QRDecomposition {

        private double[] qtStore;
        private double[] rStore;
        private int rows, cols;

        internal QRDecomposition (double[] qtStore, double[] rStore, int rows, int columns) {
            this.qtStore = qtStore;
            this.rStore = rStore;
            this.rows = rows;
            this.cols = columns;
        }

        /// <summary>
        /// The orthogonal matrix Q.
        /// </summary>
        /// <returns>The orthogonal matrix Q.</returns>
        public SquareMatrix QMatrix () {
            double[] qStore = new double[rows * rows];
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < rows; c++) {
                    qStore[rows * c + r] = qtStore[rows * r + c];
                }
            }
            return (new SquareMatrix(qStore, rows));
        }

        /// <summary>
        /// The upper-right triangular matrix R.
        /// </summary>
        /// <returns>The upper-right triangular matrix R.</returns>
        public RectangularMatrix RMatrix () {
            double[] store = new double[rows * cols];
            Array.Copy(rStore, store, rStore.Length);
            return (new RectangularMatrix(store, rows, cols));
        }

        
        /// <summary>
        /// Solve the system Q x = b.
        /// </summary>
        /// <param name="rhs">The right-hand-side b.</param>
        /// <returns>The column vector x for which Q x is closest to b.</returns>
        public ColumnVector Solve (IList<double> rhs) {

            if (rhs == null) throw new ArgumentNullException("rhs");
            if (rhs.Count != rows) throw new DimensionMismatchException();

            // copy rhs into an array, if necessary
            double[] x;
            if (rhs is double[]) {
                x = (double[]) rhs;
            } else {
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

    }

}
