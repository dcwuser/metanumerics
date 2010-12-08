﻿using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// A rectangular matrix of real numbers.
    /// </summary>
    public sealed class RectangularMatrix : RectangularMatrixBase {

        // initialization and storage

        /// <summary>
        /// Initializes a rectangular matrix with the given dimensions.
        /// </summary>
        /// <param name="rowCount">The number of rows.</param>
        /// <param name="columnCount">The number of columns.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="rowCount"/> or <paramref name="columnCount"/>
        /// is less than one.</exception>
        public RectangularMatrix (int rowCount, int columnCount) {
            if (rowCount < 1) throw new ArgumentOutOfRangeException("rowCount");
            if (columnCount < 1) throw new ArgumentOutOfRangeException("columnCount");
            rows = rowCount;
            cols = columnCount;
            store = MatrixAlgorithms.AllocateStorage(rows, cols);
        }

        /// <summary>
        /// Initializes a rectangular matrix from the given 2D array.
        /// </summary>
        /// <param name="source">The source 2D array.</param>
        public RectangularMatrix (double[,] source) {
            if (source == null) throw new ArgumentNullException("source");
            rows = source.GetLength(0);
            cols = source.GetLength(1);
            store = MatrixAlgorithms.AllocateStorage(rows, cols);
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < cols; c++) {
                    store[MatrixAlgorithms.GetIndex(rows, cols, r, c)] = source[r, c];
                }
            }
        }

        internal RectangularMatrix (double[] store, int rows, int cols) {
            this.store = store;
            this.rows = rows;
            this.cols = cols;
        }

        private double[] store;
        private int rows, cols;

        // required implementations of abstract methods

        /// <inheritdoc />
        public override double this[int r, int c] {
            get {
                if ((r < 0) || (r > rows)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c > cols)) throw new ArgumentOutOfRangeException("r");
                return (store[MatrixAlgorithms.GetIndex(rows, cols, r, c)]);
            }
            set {
                if ((r < 0) || (r > rows)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c > cols)) throw new ArgumentOutOfRangeException("r");
                store[MatrixAlgorithms.GetIndex(rows, cols, r, c)] = value;
            }
        }

        /// <inheritdoc />
        public override int RowCount {
            get { return (rows); }
        }

        /// <inheritdoc />
        public override int ColumnCount {
            get { return (cols); }
        }

        // optional, optomized replacements for RectangularMatrixBase operations

        /// <inheritdoc />
        public override double OneNorm () {
            return (MatrixAlgorithms.OneNorm(store, rows, cols));
        }

        /// <inheritdoc />
        public override double InfinityNorm () {
            return (MatrixAlgorithms.InfinityNorm(store, rows, cols));
        }

        /// <inheritdoc />
        public override ColumnVector Column (int c) {
            if ((c < 0) || (c >= cols)) throw new ArgumentOutOfRangeException("c");
            double[] cStore = new double[rows];
            Blas1.dCopy(store, rows * c, 1, cStore, 0, 1, rows);
            return (new ColumnVector(cStore, rows));
        }

        /// <inheritdoc />
        public override RowVector Row (int r) {
            if ((r < 0) || (r >= rows)) throw new ArgumentOutOfRangeException("r");
            double[] rStore = new double[cols];
            Blas1.dCopy(store, r, rows, rStore, 0, 1, cols);
            return (new RowVector(rStore, cols));
        }

        /// <summary>
        /// Adds two real, rectangular matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The sum matrix <paramref name="A"/> + <paramref name="B"/>.</returns>
        public static RectangularMatrix operator + (RectangularMatrix A, RectangularMatrix B) {
            if (A == null) throw new ArgumentNullException("A");
            if (B == null) throw new ArgumentNullException("B");
            if (A.rows != B.rows) throw new DimensionMismatchException();
            if (A.cols != B.cols) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Add(A.store, B.store, A.rows, A.cols);
            return (new RectangularMatrix(abStore, A.rows, A.cols));
        }

        /// <summary>
        /// Multiplies two real, rectangular matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The product matrix AB.</returns>
        public static RectangularMatrix operator * (RectangularMatrix A, RectangularMatrix B) {
            // this is faster than the base operator, because it knows about the underlying structure
            if (A == null) throw new ArgumentNullException("A");
            if (B == null) throw new ArgumentNullException("B");
            if (A.cols != B.rows) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Multiply(A.store, A.rows, A.cols, B.store, B.rows, B.cols);
            return (new RectangularMatrix(abStore, A.rows, B.cols));
        }

        // simple specific operations

        /// <summary>
        /// Creates a clone of the matrix.
        /// </summary>
        /// <returns>A clone of the matrix.</returns>
        public RectangularMatrix Clone () {
            double[] cStore = MatrixAlgorithms.Copy(store, rows, cols);
            return (new RectangularMatrix(cStore, rows, cols));
        }

        /// <summary>
        /// Returns the transpose of the matrix.
        /// </summary>
        /// <returns>M<sup>T</sup></returns>
        public RectangularMatrix Transpose () {
            double[] tStore = MatrixAlgorithms.Transpose(store, rows, cols);
            return (new RectangularMatrix(tStore, cols, rows));
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
        public QRDecomposition QRDecompose () {

            if (rows < cols) throw new InvalidOperationException();

            double[] rStore = MatrixAlgorithms.Copy(store, rows, cols);
            //double[] rStore = new double[store.Length];
            //Array.Copy(store, rStore, store.Length);

            //double[] qtStore = new double[rows * rows];
            //for (int i = 0; i < rows; i++) { qtStore[rows * i + i] = 1.0; }
            double[] qtStore = SquareMatrixAlgorithms.CreateUnitMatrix(rows);

            MatrixAlgorithms.QRDecompose(rStore, qtStore, rows, cols);

            return (new QRDecomposition(qtStore, rStore, rows, cols));

        }

    }

}
