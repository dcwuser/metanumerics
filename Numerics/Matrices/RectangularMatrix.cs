using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.CompilerServices;


namespace Meta.Numerics.Matrices {

    /// <summary>
    /// A rectangular matrix of real numbers.
    /// </summary>
    public sealed class RectangularMatrix : AnyRectangularMatrix {

        // initialization and storage

        /// <summary>
        /// Initializes a rectangular matrix with the given dimensions.
        /// </summary>
        /// <param name="rowCount">The number of rows.</param>
        /// <param name="columnCount">The number of columns.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="rowCount"/> or <paramref name="columnCount"/>
        /// is less than one.</exception>
        public RectangularMatrix (int rowCount, int columnCount) {
            if (rowCount < 1) throw new ArgumentOutOfRangeException(nameof(rowCount));
            if (columnCount < 1) throw new ArgumentOutOfRangeException(nameof(columnCount));
            rows = rowCount;
            cols = columnCount;
            store = MatrixAlgorithms.AllocateStorage(rows, cols, ref offset, ref rowStride, ref colStride);
        }

        /// <summary>
        /// Initializes a rectangular matrix from the given 2D array.
        /// </summary>
        /// <param name="source">The source 2D array.</param>
        /// <exception cref="ArgumentNullException"><paramref name="source"/> is <see langword="null"/>.</exception>
        public RectangularMatrix (double[,] source) {
            if (source == null) throw new ArgumentNullException(nameof(source));
            rows = source.GetLength(0);
            cols = source.GetLength(1);
            store = MatrixAlgorithms.AllocateStorage(rows, cols, ref offset, ref rowStride, ref colStride);
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < cols; c++) {
                    store[MatrixAlgorithms.GetIndex(rows, cols, r, c)] = source[r, c];
                }
            }
        }

        internal RectangularMatrix (double[] store, int offset, int rowStride, int colStride, int rows, int cols, bool isReadOnly) : base(isReadOnly) {
            Debug.Assert(store != null);
            Debug.Assert(rows > 0);
            Debug.Assert(cols > 0);
            this.store = store;
            this.offset = offset;
            this.rowStride = rowStride;
            this.colStride = colStride;
            this.rows = rows;
            this.cols = cols;
        }

        internal RectangularMatrix (double[] store, int rows, int cols, bool isReadOnly) : this(store, 0, 1, rows, rows, cols, isReadOnly) { }

        internal RectangularMatrix (double[] store, int rows, int cols) : this(store, rows, cols, false) { }

        private readonly double[] store;
        private readonly int offset, rowStride, colStride;
        private readonly int rows, cols;

        // required implementations of abstract methods

        /// <inheritdoc />
        public override double this[int r, int c] {
            get {
                CheckBounds(r, c);
                return (store[MatrixAlgorithms.GetIndex(offset, rowStride, colStride, r, c)]);
            }
            set {
                CheckBounds(r, c);
                if (IsReadOnly) throw new InvalidOperationException();
                store[MatrixAlgorithms.GetIndex(offset, rowStride, colStride, r, c)] = value;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void CheckBounds(int r, int c) {
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException(nameof(r));
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException(nameof(c));
        }

        /// <inheritdoc />
        public override int RowCount {
            get { return (rows); }
        }

        /// <inheritdoc />
        public override int ColumnCount {
            get { return (cols); }
        }

        // optional, optimized replacements for RectangularMatrixBase operations

        /// <inheritdoc />
        public override double OneNorm () {
            return (MatrixAlgorithms.OneNorm(store, offset, rowStride, colStride, rows, cols));
        }

        /// <inheritdoc />
        public override double InfinityNorm () {
            return (MatrixAlgorithms.InfinityNorm(store, offset, rowStride, colStride, rows, cols));
        }

        /// <inheritdoc />
        public override ColumnVector Column (int c) {
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException(nameof(c));
            double[] cStore = new double[rows];
            Blas1.dCopy(store, offset + colStride * c, rowStride, cStore, 0, 1, rows);
            return (new ColumnVector(cStore, rows));
        }

        /// <inheritdoc />
        public override RowVector Row (int r) {
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException(nameof(r));
            double[] rStore = new double[cols];
            Blas1.dCopy(store, offset + rowStride * r, colStride, rStore, 0, 1, cols);
            return (new RowVector(rStore, cols));
        }

        /// <summary>
        /// Adds two real, rectangular matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The sum matrix <paramref name="A"/> + <paramref name="B"/>.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="B"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="A"/> and <paramref name="B"/> do not have the same dimensions.</exception>
        public static RectangularMatrix operator + (RectangularMatrix A, RectangularMatrix B) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.rows != B.rows) throw new DimensionMismatchException();
            if (A.cols != B.cols) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Add(
                A.store, A.offset, A.rowStride, A.colStride,
                B.store, B.offset, B.rowStride, B.colStride,
                A.rows, A.cols
            );
            return (new RectangularMatrix(abStore, A.rows, A.cols));
        }

        /// <summary>
        /// Subtracts two real, rectangular matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The sum matrix <paramref name="A"/> - <paramref name="B"/>.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="B"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="A"/> and <paramref name="B"/> do not have the same dimensions.</exception>
        public static RectangularMatrix operator - (RectangularMatrix A, RectangularMatrix B) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.rows != B.rows) throw new DimensionMismatchException();
            if (A.cols != B.cols) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Subtract(
                A.store, A.offset, A.rowStride, A.colStride,
                B.store, B.offset, B.rowStride, B.colStride,
                A.rows, A.cols
            );
            return (new RectangularMatrix(abStore, A.rows, A.cols));
        }

        /// <summary>
        /// Multiplies two real, rectangular matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The product matrix AB.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="B"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The column count of <paramref name="A"/> does not equal the row count of <paramref name="B"/>.</exception>
        public static RectangularMatrix operator * (RectangularMatrix A, RectangularMatrix B) {
            // This is faster than the base operator, because it knows about the underlying storage structure
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.cols != B.rows) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Multiply(
                A.store, A.offset, A.rowStride, A.colStride,
                B.store, B.offset, B.rowStride, B.colStride,
                A.rows, A.cols, B.cols
            );
            return (new RectangularMatrix(abStore, A.rows, B.cols));
        }

        /// <summary>
        /// Multiply a real, rectangular matrix by a real constant.
        /// </summary>
        /// <param name="alpha">The constant.</param>
        /// <param name="A">The matrix.</param>
        /// <returns>The product &#x3B1;A.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> is <see langword="null"/>.</exception>
        public static RectangularMatrix operator * (double alpha, RectangularMatrix A) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            double[] store = MatrixAlgorithms.Multiply(alpha, A.store, A.offset, A.rowStride, A.colStride, A.rows, A.cols);
            return (new RectangularMatrix(store, A.rows, A.cols));
        }

        /// <summary>
        /// Divides a real, rectangular matrix by a real constant.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="alpha">The constant.</param>
        /// <returns>The quotient A/&#x3B1;.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> is <see langword="null"/>.</exception>
        public static RectangularMatrix operator / (RectangularMatrix A, double alpha) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            double[] store = MatrixAlgorithms.Multiply(1.0 / alpha, A.store, A.offset, A.rowStride, A.colStride, A.rows, A.cols);
            return (new RectangularMatrix(store, A.rows, A.cols));
        }

        /// <summary>
        /// Negates a real, rectangular matrix.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <returns>The matrix -A.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> is <see langword="null"/>.</exception>
        public static RectangularMatrix operator - (RectangularMatrix A) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            double[] store = MatrixAlgorithms.Multiply(-1.0, A.store, A.offset, A.rowStride, A.colStride, A.rows, A.cols);
            return (new RectangularMatrix(store, A.rows, A.cols));
        }

        // simple specific operations

        /// <summary>
        /// Copies the matrix.
        /// </summary>
        /// <returns>An independent copy of the matrix.</returns>
        public RectangularMatrix Copy () {
            double[] cStore = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, rows, cols);
            return (new RectangularMatrix(cStore, rows, cols));
        }

        /// <summary>
        /// Gets the transpose of the matrix.
        /// </summary>
        /// <value>M<sup>T</sup></value>
        /// <remarks>
        /// <para>The returned transpose matrix is not independent of the original matrix.
        /// Instead, it is a read-only view into the same storage as the original matrix with row an column indices reversed.
        /// This has the advantage that is can be produced with almost no time and memory cost.
        /// It has the disadvantage that any subsequent changes to the original matrix will be reflected in the returned transpose,
        /// which might not be what you expected. If you want an independent, write-able transpose matrix, call <see cref="Copy"/>
        /// on the returned matrix.
        /// </para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Transpose"/>
        public RectangularMatrix Transpose {
            get {
                // Just switch rows <-> columns on same storage
                return (new RectangularMatrix(store, offset, colStride, rowStride, cols, rows, true));
            }
            //double[] transpose = MatrixAlgorithms.Copy(store, offset, colStride, rowStride, cols, rows);
            //return (new RectangularMatrix(transpose, cols, rows));
        }



        // complicated specific operations

        /// <summary>
        /// Computes the QR decomposition of the matrix.
        /// </summary>
        /// <returns>The QR decomposition of the matrix.</returns>
        /// <remarks>
        /// <para>Only matrices with a number of rows greater than or equal to the number of columns can be QR decomposed. If your
        /// matrix has more columns than rows, you can QR decompose its transpose.</para>
        /// </remarks>
        /// <seealso cref="QRDecomposition"/>
        public QRDecomposition QRDecomposition () {

            if (rows < cols) throw new InvalidOperationException();

            double[] rStore = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, rows, cols);

            double[] qtStore = SquareMatrixAlgorithms.CreateUnitMatrix(rows);

            MatrixAlgorithms.QRDecompose(rStore, qtStore, rows, cols);

            return (new QRDecomposition(qtStore, rStore, rows, cols));

        }

        /// <summary>
        /// Computes the singular value decomposition of the matrix.
        /// </summary>
        /// <returns>The SVD of the matrix.</returns>
        public SingularValueDecomposition SingularValueDecomposition () {

            if (rows >= cols) {

                // copy the matrix so as not to disturb the original
                double[] copy = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, rows, cols);

                // bi-diagonalize it
                double[] a, b;
                MatrixAlgorithms.Bidiagonalize(copy, rows, cols, out a, out b);

                // form the U and V matrices
                double[] left = MatrixAlgorithms.AccumulateBidiagonalU(copy, rows, cols);
                double[] right = MatrixAlgorithms.AccumulateBidiagonalV(copy, rows, cols);

                // find the singular values of the bi-diagonal matrix
                MatrixAlgorithms.ExtractSingularValues(a, b, left, right, rows, cols);

                // sort them
                MatrixAlgorithms.SortValues(a, left, right, rows, cols);

                // package it all up
                return (new SingularValueDecomposition(left, a, right, rows, cols));

            } else {

                double[] scratch = MatrixAlgorithms.Transpose(store, rows, cols);

                double[] a, b;
                MatrixAlgorithms.Bidiagonalize(scratch, cols, rows, out a, out b);

                double[] left = MatrixAlgorithms.AccumulateBidiagonalU(scratch, cols, rows);
                double[] right = MatrixAlgorithms.AccumulateBidiagonalV(scratch, cols, rows);

                MatrixAlgorithms.ExtractSingularValues(a, b, left, right, cols, rows);

                MatrixAlgorithms.SortValues(a, left, right, cols, rows);

                left = MatrixAlgorithms.Transpose(left, cols, cols);
                right = MatrixAlgorithms.Transpose(right, rows, rows);

                return (new SingularValueDecomposition(right, a, left, rows, cols));


            }

        }

        /// <summary>
        /// Computes the product of a rectangular matrix and a column vector.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The column vector Av.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="v"/> is null.</exception>
        /// <exception cref="DimensionMismatchException">The column count of <paramref name="A"/> is not the same as the dimension of <paramref name="v"/>.</exception>
        public static ColumnVector operator * (RectangularMatrix A, ColumnVector v) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (v == null) throw new ArgumentNullException(nameof(v));
            if (A.cols != v.dimension) throw new DimensionMismatchException();
            double[] avStore = new double[A.rows];
            Blas2.dGemv(A.store, A.offset, A.rowStride, A.colStride, v.store, v.offset, v.stride, avStore, 0, 1, A.rows, A.cols);
            return (new ColumnVector(avStore, A.rows));
        }

        // should above method be in ColumnVector type so that all references to vector fields are in VectorBase, RowVector, and ColumnVector?

        /// <summary>
        /// Casts a rectangular matrix to a square matrix.
        /// </summary>
        /// <param name="A">The matrix to cast, which must have an equal number of rows and columns.</param>
        /// <returns>A square matrix, not independent of the original matrix.</returns>
        /// <remarks>
        /// <para>It can occur that the mode of construction of a RectangularMatrix guarantees that it is
        /// actually square. For example, if you multiply an N X M rectangular matrix by an M X N rectangular matrix,
        /// the result is an N X N square matrix. However, when determining the type of the product, the .NET
        /// compiler considers only the types of the multiplicands. Since a RectangularMatrix times a RectangularMatrix
        /// yields a RectangularMatrix, it will consider the type of the product to be RetangularMatrix, even though
        /// its rows and column dimensions will be equal. You can use this explicit cast to obtain a SquareMatrix type.
        /// </para>
        /// <para>Note that the output of the cast is not independent of the original matrix. This makes the cast operation
        /// fast, but changes to the resulting SquareMatrix will also change the original RectangularMatrix. To obtain
        /// an independent matrix, use the <see cref="Copy"/> method.</para>
        /// </remarks>
        /// <exception cref="InvalidCastException">The row and column dimensions of the matrix are not equal.</exception>
        public static explicit operator SquareMatrix (RectangularMatrix A) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (A.rows != A.cols) throw new InvalidCastException();
            return (new SquareMatrix(A.store, A.offset, A.rowStride, A.colStride, A.rows, false));
        }

    }

}
