using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a rectangular matrix of real numbers.
    /// </summary>
    public class Matrix : IMatrix {

        private double[,] M;

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
            M = new double[rows,columns];
        }

        /// <summary>
        /// Initializes a rectangular matrix from the given CLR matrix.
        /// </summary>
        /// <param name="matrix">A CLR matrix of initialization values.</param>
        /// <remarks><para>The matrix created is independent of the CLR matrix used to initialize it.</para></remarks>
        /// <exception cref="ArgumentNullException"><paramref name="matrix"/> is null</exception>
        public Matrix (double[,] matrix) {
            if (matrix == null) throw new ArgumentNullException("matrix");
            M = new double[matrix.GetLength(0), matrix.GetLength(1)];
            for (int r = 0; r < M.GetLength(0); r++) {
                for (int c = 0; c < M.GetLength(1); c++) {
                    this.M[r, c] = matrix[r, c];
                }
            }
        }

        /// <summary>
        /// Gets the number of matrix rows.
        /// </summary>
        public int RowCount {
            get {
                return(M.GetLength(0));
            }
        }

        /// <summary>
        /// Gets the number of matrix columns.
        /// </summary>
        public int ColumnCount {
            get {
                return (M.GetLength(1));
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
                return (M[r, c]);
            }
            set {
                if ((r < 0) || (r >= RowCount)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= ColumnCount)) throw new ArgumentOutOfRangeException("c");
                M[r, c] = value;
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
            Matrix N = new Matrix(RowCount, ColumnCount);
            for (int r = 0; r < RowCount; r++) {
                for (int c = 0; c < ColumnCount; c++) {
                    N[r, c] = this[r, c];
                }
            }
            return (N);
        }

        IMatrix IMatrix.Clone () {
            return (Clone());
        }

        /// <summary>
        /// Generates the transpose of the matrix.
        /// </summary>
        /// <returns>The matrix M<sup>T</sup>.</returns>
        public Matrix Transpose () {
            Matrix N = new Matrix(ColumnCount, RowCount);
            for (int r = 0; r < RowCount; r++) {
                for (int c = 0; c < ColumnCount; c++) {
                    N[c, r] = this[r, c];
                }
            }
            return (N);
        }

        IMatrix IMatrix.Transpose () {
            return (Transpose());
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

        // matrix arithmetic

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
            return (Add(M1, M2));
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

    }

}
