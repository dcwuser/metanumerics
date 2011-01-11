using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace Meta.Numerics.Matrices {

    // The new design idea is as follows:
    // 1. The Matrix class is public. It is basically a pure shim over the storage format. To actually do all the operations it makes public,
    //    it calls into the MatrixAlgorithm class.
    // 2. The MatrixAlgorithm class is internal. It contains all the algorithmic logic for doing matrix operations on the storage format. It
    //    calls into BLAS to perform optomized low-level operations.
    // 3. The Blas classes have optomized low-level operations.
    // (1) and (2) should be repeated for each kind of matrix: square, symmetric, tridiagonal, etc.

#if FUTURE

    /// <summary>
    /// Represents a rectangular matrix of real numbers.
    /// </summary>
    public class Matrix : IMatrix {

        private double[] store;
        private int rows;
        private int columns;

        /// <summary>
        /// Initializes a rectangular matrix with the given dimensions.
        /// </summary>
        /// <param name="rows">The number of rows.</param>
        /// <param name="columns">The number of columns.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="rows"/> or <paramref name="columns"/>
        /// is less than one.</exception>
        public Matrix (int rows, int columns) {
            if (rows < 1) throw new ArgumentOutOfRangeException("rows");
            if (columns < 1) throw new ArgumentOutOfRangeException("columns");
            this.store = new double[rows * columns];
            this.rows = rows;
            this.columns = columns;
        }

        /// <summary>
        /// Initializes a rectangular matrix from the given CLR matrix.
        /// </summary>
        /// <param name="matrix">A CLR matrix of initialization values.</param>
        /// <remarks><para>The matrix created is independent of the CLR matrix used to initialize it.</para></remarks>
        /// <exception cref="ArgumentNullException"><paramref name="matrix"/> is null</exception>
        public Matrix (double[,] matrix) {
            if (matrix == null) throw new ArgumentNullException("matrix");
            rows = matrix.GetLength(0);
            columns = matrix.GetLength(1);
            store = new double[rows * columns];
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < columns; c++) {
                    store[rows * c + r] = matrix[r, c];
                }
            }
        }

        internal Matrix (double[] store, int rows, int columns) {
            this.store = store;
            this.rows = rows;
            this.columns = columns;
        }

        /// <summary>
        /// Gets the number of matrix rows.
        /// </summary>
        public int RowCount {
            get {
                return (rows);
            }
        }

        /// <summary>
        /// Gets the number of matrix columns.
        /// </summary>
        public int ColumnCount {
            get {
                return (columns);
            }
        }

        /// <summary>
        /// Gets or sets the value of a matrix entry.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The value of the <paramref name="r"/>,<paramref name="c"/> matrix entry.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> or <paramref name="c"/> is
        /// outside the valid range.</exception>
        public double this[int r, int c] {
            get {
                if ((r < 0) || (r >= RowCount)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= ColumnCount)) throw new ArgumentOutOfRangeException("c");
                return (store[rows * c + r]);
            }
            set {
                if ((r < 0) || (r >= RowCount)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= ColumnCount)) throw new ArgumentOutOfRangeException("c");
                store[rows * c + r] = value;
            }
        }

        /// <summary>
        /// Gets a copy of a matrix row.
        /// </summary>
        /// <param name="r">The (zero-based) row number.</param>
        /// <returns>The <paramref name="r"/>th row of the matrix.</returns>
        /// <remarks><para>The returned row vector is independent of the matrix; changing its values
        /// will not change the entries of the original matrix.</para></remarks>
        /// <exception cref="ArgumentOutOfRangeException">The specified row does not exist.</exception>
        public RowVector Row (int r) {
            if ((r < 0) || (r >= RowCount)) throw new ArgumentOutOfRangeException("r");
            RowVector v = new RowVector(ColumnCount);
            for (int c = 0; c < ColumnCount; c++) {
                v[c] = this[r, c];
            }
            return (v);
        }

        /// <summary>
        /// Gets a copy of a matrix column.
        /// </summary>
        /// <param name="c">The (zero-based) column number.</param>
        /// <returns>The <paramref name="c"/>th column of the matrix.</returns>
        /// <remarks><para>The returned column vector is independent of the matrix; changing its values
        /// will not change the entries of the original matrix.</para></remarks>
        /// <exception cref="ArgumentOutOfRangeException">The specified column does not exist.</exception>
        public ColumnVector Column (int c) {
            if ((c < 0) || (c >= ColumnCount)) throw new ArgumentOutOfRangeException("c");
            ColumnVector v = new ColumnVector(RowCount);
            for (int r = 0; r < RowCount; r++) {
                v[r] = this[r, c];
            }
            return (v);
        }

        /// <summary>
        /// Creates a clone of the matrix.
        /// </summary>
        /// <returns>A clone of the matrix.</returns>
        public Matrix Clone () {
            double[] cStore = MatrixAlgorithms.Clone(store, rows, columns);
            return (new Matrix(cStore, rows, columns));
        }

        IMatrix IMatrix.Clone () {
            return (Clone());
        }

        /// <summary>
        /// Generates the transpose of the matrix.
        /// </summary>
        /// <returns>The matrix M<sup>T</sup>.</returns>
        public Matrix Transpose () {
            double[] tStore = MatrixAlgorithms.Transpose(store, rows, columns);
            return (new Matrix(tStore, columns, rows));
        }

        IMatrix IMatrix.Transpose () {
            return (Transpose());
        }

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

            if (rows < columns) throw new InvalidOperationException();

            double[] rStore = new double[store.Length];
            Array.Copy(store, rStore, store.Length);

            double[] qtStore = new double[rows * rows];
            for (int i = 0; i < rows; i++) { qtStore[rows * i + i] = 1.0; }

            MatrixAlgorithms.QRDecompose(rStore, qtStore, rows, columns);

            return (new QRDecomposition(qtStore, rStore, rows, columns));

        }

        // Operators

        // equality

        internal static bool Equals (IMatrix M1, IMatrix M2) {
            if (((object) M1) == null) {
                if (((object) M2) == null) {
                    // both are null
                    return (true);
                } else {
                    // only M1 is null
                    return (false);
                }
            } else {
                if (((object) M2) == null) {
                    // only M2 is null
                    return (false);
                } else {
                    // neither is null
                    if (M1.RowCount != M2.RowCount) throw new DimensionMismatchException();
                    if (M1.ColumnCount != M2.ColumnCount) throw new DimensionMismatchException();
                    for (int r = 0; r < M1.RowCount; r++) {
                        for (int c = 0; c < M1.ColumnCount; c++) {
                            if (M1[r, c] != M2[r, c]) return (false);
                        }
                    }
                    return (true);
                }
            }
        }

        /// <summary>
        /// Determines whether two matrices are equal.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>True if <paramref name="M1"/> and <paramref name="M2"/> are equal, otherwise false.</returns>
        public static bool operator == (Matrix M1, Matrix M2) {
            return (Equals(M1, M2));
        }

        /// <summary>
        /// Determines whether two matrices are not equal.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>False if <paramref name="M1"/> and <paramref name="M2"/> are equal, otherwise true.</returns>
        public static bool operator != (Matrix M1, Matrix M2) {
            return (!Equals(M1, M2));
        }

        /// <summary>
        /// Determines whether the given object is an equal matrix.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if <paramref name="obj"/> is an equal matrix, otherwise false.</returns>
        public override bool Equals (object obj) {
            return (Matrix.Equals(this, obj as IMatrix));
        }

        /// <inheritdoc />
        public override int GetHashCode () {
            return base.GetHashCode();
        }

        // matrix arithmetic
        /*
        internal static Matrix Add (IMatrix M1, IMatrix M2) {
            if (M1.RowCount != M2.RowCount) throw new DimensionMismatchException();
            if (M1.ColumnCount != M2.ColumnCount) throw new DimensionMismatchException();
            Matrix N = new Matrix(M1.RowCount, M1.ColumnCount);
            for (int r = 0; r < M1.RowCount; r++) {
                for (int c = 0; c < M1.ColumnCount; c++) {
                    N[r, c] = M1[r, c] + M2[r, c];
                }
            }
            return (N);
        }
        */

        /// <summary>
        /// Computes the sum of two matrices.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>The sum <paramref name="M1"/> + <paramref name="M2"/>.</returns>
        /// <remarks>
        /// <para>Matrix addition is an O(N<sup>2</sup>) process.</para>
        /// </remarks>
        public static Matrix operator + (Matrix M1, Matrix M2) {
            if (M1.rows != M2.rows) throw new DimensionMismatchException();
            if (M1.columns != M2.columns) throw new DimensionMismatchException();
            double[] sStore = MatrixAlgorithms.Add(M1.store, M2.store, M1.rows, M1.columns);
            return(new Matrix(sStore, M1.rows, M1.columns));
        }

        internal static Matrix Subtract (IMatrix M1, IMatrix M2) {
            if (M1.RowCount != M2.RowCount) throw new DimensionMismatchException();
            if (M1.ColumnCount != M2.ColumnCount) throw new DimensionMismatchException();
            Matrix N = new Matrix(M1.RowCount, M1.ColumnCount);
            for (int r = 0; r < M1.RowCount; r++) {
                for (int c = 0; c < M1.ColumnCount; c++) {
                    N[r, c] = M1[r, c] - M2[r, c];
                }
            }
            return (N);
        }

        /// <summary>
        /// Computes the difference of two matrices.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>The difference <paramref name="M1"/> - <paramref name="M2"/>.</returns>
        /// <remarks>
        /// <para>Matrix subtraction is an O(N<sup>2</sup>) process.</para>
        /// </remarks>
        public static Matrix operator - (Matrix M1, Matrix M2) {
            return(Subtract(M1,M2));
        }

        internal static Matrix Multiply (IMatrix M1, IMatrix M2) {
            if (M1.ColumnCount != M2.RowCount) throw new DimensionMismatchException();
            Matrix N = new Matrix(M1.RowCount, M2.ColumnCount);
            for (int r = 0; r < N.RowCount; r++) {
                for (int c = 0; c < N.ColumnCount; c++) {
                    N[r, c] = 0.0;
                    for (int i = 0; i < M1.ColumnCount; i++) {
                        N[r, c] += M1[r, i] * M2[i, c];
                    }
                }
            }
            return (N);
        }

        /// <summary>
        /// Computes the product of two matrices.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>The product <paramref name="M1"/> * <paramref name="M2"/>.</returns>
        /// <remarks>
        /// <para>Note that matrix multiplication is not commutative, i.e. M1*M2 is generally not the same as M2*M1.</para>
        /// <para>Matrix multiplication is an O(N<sup>3</sup>) process.</para>
        /// </remarks>        
        public static Matrix operator * (Matrix M1, Matrix M2) {
            return (Multiply(M1, M2));
        }

        public static Matrix operator * (Matrix M1, SquareMatrix M2) {
            return (Multiply(M1, M2));
        }

        public static Matrix operator * (SquareMatrix M1, Matrix M2) {
            return (Multiply(M1, M2));
        }

        public static Matrix operator * (Matrix M1, SymmetricMatrix M2) {
            return (Multiply(M1, M2));
        }

        public static Matrix operator * (SymmetricMatrix M1, Matrix M2) {
            return (Multiply(M1, M2));
        }

        // arithmetic with doubles

        internal static Matrix Multiply (double x, IMatrix M) {
            Matrix N = new Matrix(M.RowCount, M.ColumnCount);
            for (int r=0; r<N.RowCount; r++) {
                for (int c = 0; c < N.ColumnCount; c++) {
                    N[r, c] = x * M[r, c];
                }
            }
            return (N);
        }

        /// <summary>
        /// Computes the product of a real number and a matrix.
        /// </summary>
        /// <param name="x">The real number.</param>
        /// <param name="M">The matrix.</param>
        /// <returns>The product of <paramref name="x"/> <paramref name="M"/>.</returns>
        public static Matrix operator * (double x, Matrix M) {
            return (Multiply(x, M));
        }

        /// <summary>
        /// Computes the the quotient of a matrix and a real number.
        /// </summary>
        /// <param name="M">The matrix.</param>
        /// <param name="x">The real number.</param>
        /// <returns>The quotient <paramref name="M"/>/<paramref name="x"/>.</returns>
        public static Matrix operator / (Matrix M, double x) {
            return (Multiply(1.0 / x, M));
        }

        /*
        public static Matrix operator * (Matrix M, double x) {
            return (Multiply(x, M));
        }
        */

        internal static void WriteMatrix (IMatrix matrix, TextWriter writer) {

            for (int r = 0; r < matrix.RowCount; r++) {
                writer.Write("{ ");
                for (int c = 0; c < matrix.ColumnCount; c++) {
                    writer.Write("{0,15:g12} ", matrix[r, c]);
                }
                writer.WriteLine("}");
            }

        }

#if SHO
        /// <summary>
        /// Produces the representation of the matrix for the Python interactive console.
        /// </summary>
        /// <returns>A string representation of the matrix.</returns>
        public string __repr__ () {
            StringWriter writer = new StringWriter();
            Matrix.WriteMatrix(this, writer);
            return (writer.ToString());
        }
#endif

    }

#endif

    internal static class MatrixAlgorithms {

        public static double[] AllocateStorage (int nRows, int nCols) {
            return (new double[nRows * nCols]);
        }

        // basic access: O(1) operations

        public static int GetIndex (int nRows, int nCols, int r, int c) {
            return (nRows * c + r);
        }

        public static void SetEntry (double[] store, int nRows, int nCols, int r, int c, double value) {
            store[nRows * c + r] = value;
        }

        public static double GetEntry (double[] store, int nRows, int nCols, int r, int c) {
            return (store[nRows * c + r]);
        }

        // transformations and simple arithmetic: O(N^2) operations

        public static double[] Copy (double[] store, int nRows, int nColumns) {
            double[] cloneStore = new double[store.Length];
            store.CopyTo(cloneStore, 0);
            return (cloneStore);
        }

        public static double[] Transpose (double[] store, int nRows, int nColumns) {
            double[] tStore = new double[store.Length];
            for (int r = 0; r < nRows; r++) {
                for (int c = 0; c < nColumns; c++) {
                    tStore[nColumns * r + c] = store[nRows * c + r];
                }
            }
            return (tStore);
        }

        public static double OneNorm (double[] store, int nRows, int nColumns) {
            double norm = 0.0;
            for (int c = 0; c < nColumns; c++) {
                double csum = Blas1.dNrm1(store, c * nRows, 1, nRows);
                if (csum > norm) norm = csum;
            }
            return (norm);
        }

        public static double InfinityNorm (double[] store, int nRows, int nColumns) {
            double norm = 0.0;
            for (int r = 0; r < nRows; r++) {
                double rsum = Blas1.dNrm1(store, r, nRows, nColumns);
                if (rsum > norm) norm = rsum;
            }
            return (norm);
        }

        public static double[] Add (double[] xStore, double[] yStore, int nRows, int nCols) {
            double[] store = new double[nRows * nCols];
            for (int i = 0; i < store.Length; i++) {
                store[i] = xStore[i] + yStore[i];
            }
            return (store);
        }

        public static double[] Subtract (double[] xStore, double[] yStore, int nRows, int nCols) {
            double[] store = new double[nRows * nCols];
            for (int i = 0; i < store.Length; i++) {
                store[i] = xStore[i] - yStore[i];
            }
            return (store);
        }

        public static double[] Multiply (double alpha, double[] store, int nRows, int nCols) {
            double[] result = new double[nRows * nCols];
            for (int i = 0; i < store.Length; i++) {
                result[i] = alpha * store[i];
            }
            return (result);
        }

        // multiplication: an O(N^3) algorithm

        public static double[] Multiply (double[] aStore, int aRows, int aCols, double[] bStore, int bRows, int bCols) {

            double[] abStore = new double[aRows * bCols];
            for (int j = 0; j < bCols; j++) {
                int aRows_j = aRows * j;
                int bRows_j = bRows * j;
                for (int k = 0; k < aCols; k++) {
                    double t = bStore[bRows_j + k];
                    if (t != 0.0) {
                        int aRows_k = aRows * k;
                        for (int i = 0; i < aRows; i++) {
                            abStore[aRows_j + i] += t * aStore[aRows_k + i];
                        }
                    }
                }
            }

            return (abStore);

        }

        // more complicated operations: also O(N^3)

        // QR decomposition write A = Q R, where Q is orthogonal (Q^T = Q^{-1}) and R is upper-right triangular.

        // Given a QR decomposition, it is easy to solve Ax=b by writing (QR)x=b, Q(Rx)=b, Rx=y and Qy=b.
        // Then y = Q^T b (an N^2 multiplication) and R x = y is a triangular system (also N^2 to solve).

        // To form a QR decomposition, we apply householder reflections to zero each column in turn.

        // XXXXX    ABCDE    XXXXX    XXXXX    XXXXX    XXXXX
        // XXXXX    0BCDE    0ABCD    0XXXX    0XXXX    0XXXX
        // XXXXX    0BCDE    00BCD    00ABC    00XXX    00XXX
        // XXXXX -> 0BCDE -> 00BCD -> 000BC -> 000AB -> 000XX
        // XXXXX    0BCDE    00BCD    000BC    0000B    0000A    
        // XXXXX    0BCDE    00BCD    000BC    0000B    00000
        // XXXXX    0BCDE    00BCD    000BC    0000B    00000

        public static void QRDecompose (double[] store, double[] qtStore, int rows, int cols) {

            double[] vStore = new double[rows];

            // loop over columns
            for (int k = 0; k < cols; k++) {

                int offset = rows * k + k;
                int length = rows - k;
                double a;
                SquareMatrixAlgorithms.GenerateHouseholderReflection(store, offset, 1, length, out a);

                // apply P to the other columns of A
                for (int c = k + 1; c < cols; c++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, offset, 1, store, rows * c + k, 1, length);
                }

                // apply P to Q to accumulate the transform
                // since Q is rows X rows, its column count is just rows
                for (int c = 0; c < rows; c++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, offset, 1, qtStore, rows * c + k, 1, length);
                }

                // we are done with the Householder vector now, so we can zero the first column
                store[offset] = a;
                for (int i = 1; i < length; i++) {
                    store[offset + i] = 0.0;
                }

            }

        }

        // ** SVD **

        public static void ExtractSingularValues (double[] a, double[] b, double[] uStore, double[] vStore, int rows, int cols) {

            int p = 0;
            int n = cols - 1;

            int count = 0;
            while (count < 32) {

                // move the lower boundary up as far as we can
                while ((n > 0) && (Math.Abs(b[n - 1]) <= (2.0E-16) * Math.Abs(a[n - 1] + a[n]))) {
                    b[n - 1] = 0.0;
                    if (a[n] < 0.0) {
                        a[n] = -a[n];
                        if (vStore != null) Blas1.dScal(-1.0, vStore, n * cols, 1, cols);
                    }
                    n--;
                    count = 0;
                }

                if (n == 0) return;

                GolubKahn(a, b, p, n, uStore, vStore, rows, cols);
                count++;

            }

            throw new NonconvergenceException();

        }

        // A Golub-Kahn step in the SVD is like a Francis step in a QR eigendecomposition. After each step, the superdiagonal
        // elements of the bidiagonal matrix become smaller, making it converge to diagonal form.

        // Finding the singular values of A is equivilent to finding the eigenvalues of A^T A. For a bidiagonal matrix,
        // B^T B is tridiagonal. A QR step on B^T B would find the eigenvalues of its lower right 2 X 2 part, determine
        // the transformation that would bring (B^T B - e I)

        // |XX   |    |AA   |    |ABC  |    |XA0  |    |XX   |    |XX   |    |XX   |    |XX   |    |XX   |
        // | XX  |    |BBX  |    |0BC  |    | BB  |    | ABC |    | XA0 |    | XX  |    | XX  |    | XX  |
        // |  XX | -> |  XX | -> |  XX | -> | CCX | -> | 0BC | -> |  BB | -> |  ABC| -> |  XA0| -> |  XX |
        // |   XX|    |   XX|    |   XX|    |   XX|    |   XX|    |  CCX|    |  0BC|    |   BB|    |   AB|
        // |    X|    |    X|    |    X|    |    X|    |    X|    |    X|    |    X|    |   CC|    |   0B|

        private static void GolubKahn (double[] a, double[] b, int p, int n, double[] uStore, double[] vStore, int rows, int cols) {

            //WriteBidiagonal(a, b, 0.0, n);

            int m = n - 1;

            // find the eigenvalues of the bottom 2 X 2 of the tridiagonal matrix B^T * B
            // store the eigenvalue closer to the bottom value
            double T_mm = MoreMath.Pow(a[m], 2);
            if (m > p) T_mm += MoreMath.Pow(b[m - 1], 2);
            double T_nn = MoreMath.Pow(a[n], 2) + MoreMath.Pow(b[m], 2);
            double T_mn = a[m] * b[m];
            double dd = Math.Sqrt(Math.Pow(T_mm - T_nn, 2) + 4.0 * Math.Pow(T_mn, 2));
            double e1 = (T_mm + T_nn + dd) / 2.0;
            double e2 = (T_mm + T_nn - dd) / 2.0;
            //Console.WriteLine("e = {0} {1}", e1, e2);
            double e;
            if (Math.Abs(e1 - T_nn) < Math.Abs(e2 - T_nn)) {
                e = e1;
            } else {
                e = e2;
            }
            //Console.WriteLine(e);

            // find the givens rotation that would zero the first column of (B^T B - e I)

            // variables for cosine and sine of rotation angles
            double c, s;

            double u = MoreMath.Pow(a[p], 2) - e;
            double v = a[p] * b[p];
            GenerateGivensRotation(ref u, ref v, out c, out s);

            // apply that rotation directly as B * G

            // a temporary variable to store the one non-biadiagonal element
            double t = 0.0;

            ApplyGivensRotation(c, s, ref a[p], ref b[p]);
            ApplyGivensRotation(c, s, ref t, ref a[p + 1]);

            if (vStore != null) {
                for (int j = 0; j < cols; j++) {
                    ApplyGivensRotation(c, s, ref vStore[p * cols + j], ref vStore[(p + 1) * cols + j]);
                }
            }

            //WriteBidiagonal(a, b, t, n);

            // "chase the bulge" down B
            for (int i = p; i < m; i++) {

                //Console.WriteLine("i = {0}", i);

                // G * B to zero sub-diagonal t
                GenerateGivensRotation(ref a[i], ref t, out c, out s);
                ApplyGivensRotation(c, s, ref b[i], ref a[i + 1]);
                ApplyGivensRotation(c, s, ref t, ref b[i + 1]);

                if (uStore != null) {
                    for (int j = 0; j < rows; j++) {
                        ApplyGivensRotation(c, s, ref uStore[j * rows + i], ref uStore[j * rows + i + 1]);
                    }
                }

                // B * G to zero (super+1)-diagonal t
                GenerateGivensRotation(ref b[i], ref t, out c, out s);
                ApplyGivensRotation(c, s, ref a[i + 1], ref b[i + 1]);
                ApplyGivensRotation(c, s, ref t, ref a[i + 2]);

                //WriteBidiagonal(a, b, t, n);

                if (vStore != null) {
                    for (int j = 0; j < cols; j++) {
                        ApplyGivensRotation(c, s, ref vStore[(i+1) * cols + j], ref vStore[(i + 2) * cols + j]);
                    }
                }

            }

            // G * B to remove one last sub-diagonal t
            GenerateGivensRotation(ref a[m], ref t, out c, out s);
            ApplyGivensRotation(c, s, ref b[m], ref a[n]);

            if (uStore != null) {
                for (int j = 0; j < rows; j++) {
                    ApplyGivensRotation(c, s, ref uStore[j * rows + m], ref uStore[j * rows + n]);
                }
            }

            //WriteBidiagonal(a, b, t, n);

        }

        // [  c  s ] [ u ] = [ u' ]
        // [ -s  c ] [ v ] = [ v' ]

        private static void ApplyGivensRotation (double c, double s, ref double u, ref double v) {
            double t = u;
            u = c * t + s * v;
            v = c * v - s * t;
        }

        // A Givens rotation is a simple 2D rotation to make a vector horizonal, i.e. to zero one component.

        // [  c  s ] [ a ] = [ r ]
        // [ -s  c ] [ b ] = [ 0 ]

        // This routine replaces a and b with r and 0, and returns c and s.

        private static void GenerateGivensRotation (ref double a, ref double b, out double c, out double s) {
            if (b == 0.0) {
                c = 1.0;
                s = 0.0;
            } else {
                double r = MoreMath.Hypot(a, b);
                c = a / r;
                s = b / r;
                a = r;
                b = 0.0;
            }
        }


        // SVD

        // Bidiagonalization uses seperate left and right Householder reflections to bring a matrix into a form of bandwidth two.

        // XXXXX    ABCDE    XA000    XX000    XX000
        // XXXXX    0BCDE    0BBBB    0ABCD    0XA00
        // XXXXX    0BCDE    0CCCC    00BCD    00BBB
        // XXXXX -> 0BCDE -> 0DDDD -> 00BCD -> 00CCC -> etc
        // XXXXX    0BCDE    0EEEE    00BCD    00DDD
        // XXXXX    0BCDE    0FFFF    00BCD    00EEE
        // XXXXX    0BCDE    0GGGG    00BCD    00FFF

        // The end result is B = U A V, where U and V are different (so this isn't a similiarity transform) orthogonal matrices.

        // Note that we can't get fully diagonal with this approach, because if the right transform tried to zero the first superdiagonal
        // element, it would interfere with the zeros already created.

        // The reslting diagonal and first superdiagonal elements are returned in a and b. The matrix itself is gradually
        // overwritten with the elements of the Householder reflections used. When we are done
        // the matrix is replaced by Householder vector components as follows:

        // APPPP
        // ABQQQ
        // ABCRR
        // ABCDS
        // ABCDE
        // ABCDE

        // where A, B, C, D, E are the Householder reflections applied from the left to zero columns, and P, Q, R, S are
        // the Householder reflections applied from the right to zero rows.

        public static void Bidiagonalize (double[] store, int rows, int cols, out double[] a, out double[] b) {

            Debug.Assert(rows >= cols);
            a = new double[cols];
            b = new double[cols - 1];


            for (int k = 0; k < cols; k++) {

                // generate a householder reflection which, when applied from the left, will zero sub-diagonal elements in the kth column
                // use those elements to store the reflection vector; store the resulting diagonal element seperately
                SquareMatrixAlgorithms.GenerateHouseholderReflection(store, rows * k + k, 1, rows - k, out a[k]);

                // apply that reflection to all the columns; only subsequent ones are affected since preceeding ones
                // contain only zeros in the subdiagonal rows
                for (int c = k + 1; c < cols; c++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, rows * k + k, 1, store, rows * c + k, 1, rows - k);
                }

                if ((k + 1) < cols) {
                    // generate a Householder reflection which, when applied from the right, will zero (super+1)-diagonal elements in the kth row
                    // again, store the elements of the reflection vector in the zeroed matrix elements and store the resulting superdiagonal element seperately
                    SquareMatrixAlgorithms.GenerateHouseholderReflection(store, rows * (k + 1) + k, rows, cols - (k + 1), out b[k]);

                    // apply the reflection to all the rows; only subsequent ones are affected since the preceeding ones contain only zeros in the
                    // affected columns; the already-zeroed column elements are not distrubed because those columns are not affected
                    // this restriction is why we cannot fully diagonalize using this transform
                    for (int r = k + 1; r < rows; r++) {
                        SquareMatrixAlgorithms.ApplyHouseholderReflection(store, rows * (k + 1) + k, rows, store, rows * (k + 1) + r, rows, cols - (k + 1));
                    }
                }


            }

        }

        public static double[] AccumulateBidiagonalV (double[] store, int rows, int cols) {
            //  Q_1 * ... (Q_{n-1} * (Q_n * 1))
            double[] result = SquareMatrixAlgorithms.CreateUnitMatrix(cols);
            for (int k = cols - 2; k >= 0; k--) {
                // apply Householder reflection to each column from the left
                for (int j = k + 1; j < cols; j++) {
                    //Console.WriteLine("k={0}", k);
                    //Console.WriteLine("offset={0} stride={1} offset={2} stride={3} count={4}", (k + 1) * rows + k, rows, (k + 1) * rows + (k + 1), rows, cols - k - 1);
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, (k + 1) * rows + k, rows, result, j * cols + (k + 1), 1, cols - k - 1);
                }
            }
            return (result);

        }

        public static double[] AccumulateBidiagonalU (double[] store, int rows, int cols) {
            // ((1 *  Q_n) * Q_{n-1}) ... * Q_1
            double[] result = SquareMatrixAlgorithms.CreateUnitMatrix(rows);

            // iterate over Householder reflections
            for (int k = cols - 1; k >= 0; k--) {
                // apply Householder reflection to each row from the right
                for (int j = k; j < rows; j++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, k * rows + k, 1, result, k * rows + j, rows, rows - k);
                }
            }

            return (result);
        }

        public static void SortValues (double[] values, double[] utStore, double[] vStore, int rows, int cols) {

            // this is a selection sort, an O(N^2) sort which requires fewer swaps than an insertion sort or an O(N ln N) sort

            // loop over ranks
            for (int i = 0; i < cols; i++) {

                // find the next largest value
                int j = i;
                double t = values[i];
                for (int k = i + 1; k < cols; k++) {
                    if (values[k] > t) {
                        j = k;
                        t = values[k];
                    }
                }

                // if necessary swap it with the current element
                if (j != i) {
                    Global.Swap(ref values[i], ref values[j]);
                    if (utStore != null) Blas1.dSwap(utStore, i, rows, utStore, j, rows, rows);
                    if (vStore != null) Blas1.dSwap(vStore, i * cols, 1, vStore, j * cols, 1, cols);
                }

            }

        }

    }

}
