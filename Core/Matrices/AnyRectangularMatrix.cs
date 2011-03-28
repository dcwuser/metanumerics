using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Describes the form of all real matrices.
    /// </summary>
    /// <remarks>
    /// <para>This is an abstract class that describes any real matrix. If you wish to create a concrete
    /// instance of a real, non-square matrix, use the <see cref="RectangularMatrix"/> class. If, on the
    /// other hand, you wish to write a function that can operate on any real matrix, it's probably a good
    /// idea to accept a <see cref="AnyRectangularMatrix"/>, so that any concrete implementation
    /// can also be passed into your function.</para>
    /// </remarks>
    public abstract class AnyRectangularMatrix : AnyMatrix<double> {

#if TEMP

        /// <summary>
        /// Gets or sets the value of a matrix entry.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The value of the <paramref name="r"/>,<paramref name="c"/> matrix entry.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> or <paramref name="c"/> is
        /// outside the valid range.</exception>
        public abstract double this[int r, int c] { get; set; }

        /// <summary>
        /// Gets the number of matrix rows.
        /// </summary>
        public abstract int RowCount { get; }

        /// <summary>
        /// Gets the number of matrix columns.
        /// </summary>
        public abstract int ColumnCount { get; }

#endif

        // we an implement some operations, but many will be slow because they do not have access to the underlying storage
        // we will override them with faster implementations in implementing classes, but these ensure that the operations
        // are always defined, including, for binary operations, between different matrix types

        /// <summary>
        /// Computes the 1-norm of the matrix.
        /// </summary>
        /// <returns>||M||<sub>1</sub></returns>
        /// <remarks>
        /// <para>The 1-norm of a matrix is the largest column sum.</para>
        /// </remarks>
        public virtual double OneNorm () {
            // one-norm is maximum column sum
            double norm = 0.0;
            for (int c = 0; c < ColumnCount; c++) {
                double csum = 0.0;
                for (int r = 0; r < RowCount; r++) {
                    csum += Math.Abs(this[r, c]);
                }
                if (csum > norm) norm = csum;
            }
            return (norm);
        }

        /// <summary>
        /// Computes the &#x221E;-norm of the matrix.
        /// </summary>
        /// <returns>||M||<sub>&#x221E;</sub></returns>
        /// <remarks>
        /// <para>The &#x221E;-norm of a matrix is the largest row sum.</para>
        /// </remarks>
        public virtual double InfinityNorm () {
            // infinity-norm is maximum row sum
            double norm = 0.0;
            for (int r = 0; r < RowCount; r++) {
                double rsum = 0.0;
                for (int c = 0; c < ColumnCount; c++) {
                    rsum += Math.Abs(this[r, c]);
                }
                if (rsum > norm) norm = rsum;
            }
            return (norm);
        }

        /// <summary>
        /// Computes the Frobenius-norm of the matrix.
        /// </summary>
        /// <returns>||M||<sub>F</sub></returns>
        /// <remarks>
        /// <para>The Frobenius-norm of a matrix the square root of the sum of the squares
        /// of all the elements. In the case of a row or column vector, this reduces
        /// to the Euclidean vector norm.</para>
        /// </remarks>
        public virtual double FrobeniusNorm () {
            double norm = 0.0;
            for (int r = 0; r < RowCount; r++) {
                for (int c = 0; c < ColumnCount; c++) {
                    norm += MoreMath.Pow2(this[r, c]);
                }
            }
            return (Math.Sqrt(norm));
        }

        /// <summary>
        /// Gets a copy of the specified column.
        /// </summary>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>An independent copy of the specified column.</returns>
        public virtual ColumnVector Column (int c) {
            if ((c < 0) || (c >= ColumnCount)) throw new ArgumentOutOfRangeException("c");
            ColumnVector v = new ColumnVector(RowCount);
            for (int r = 0; r < RowCount; r++) {
                v[r] = this[r, c];
            }
            return (v);
        }

        /// <summary>
        /// Gets a copy of the specified row.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <returns>An independent copy of the specified row.</returns>
        public virtual RowVector Row (int r) {
            if ((r < 0) || (r >= RowCount)) throw new ArgumentOutOfRangeException("r");
            RowVector v = new RowVector(ColumnCount);
            for (int c = 0; c < ColumnCount; c++) {
                v[c] = this[r, c];
            }
            return (v);
        }

        /// <summary>
        /// Adds any two real, rectangular matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The sum matrix A + B.</returns>
        /// <remarks>
        /// <para>Matrix addition is an O(N<sup>2</sup>) process.</para>
        /// </remarks>
        public static RectangularMatrix operator + (AnyRectangularMatrix A, AnyRectangularMatrix B) {
            if (A == null) throw new ArgumentNullException("A");
            if (B == null) throw new ArgumentNullException("B");
            if (A.RowCount != B.RowCount) throw new DimensionMismatchException();
            if (A.ColumnCount != B.ColumnCount) throw new DimensionMismatchException();
            RectangularMatrix M = new RectangularMatrix(A.RowCount, A.ColumnCount);
            for (int i = 0; i < M.RowCount; i++) {
                for (int j = 0; j < M.ColumnCount; j++) {
                    M[i, j] = A[i, j] + B[i, j];
                }
            }
            return (M);
        }

        /// <summary>
        /// Subtracts any two real, rectangular matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The difference matrix A - B.</returns>
        /// <remarks>
        /// <para>Matrix subtraction is an O(N<sup>2</sup>) process.</para>
        /// </remarks>
        public static RectangularMatrix operator - (AnyRectangularMatrix A, AnyRectangularMatrix B) {
            if (A == null) throw new ArgumentNullException("A");
            if (B == null) throw new ArgumentNullException("B");
            if (A.RowCount != B.RowCount) throw new DimensionMismatchException();
            if (A.ColumnCount != B.ColumnCount) throw new DimensionMismatchException();
            RectangularMatrix M = new RectangularMatrix(A.RowCount, A.ColumnCount);
            for (int i = 0; i < M.RowCount; i++) {
                for (int j = 0; j < M.ColumnCount; j++) {
                    M[i, j] = A[i, j] - B[i, j];
                }
            }
            return (M);
        }

        /// <summary>
        /// Multiplies any two real, rectangular matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The product matrix AB.</returns>
        /// <remarks>
        /// <para>For matrix multiplication, the column count of the first matrix must equal the row count of the second
        /// matrix.</para>
        /// <para>Matrix multiplication is an O(N<sup>3</sup>) process.</para>
        /// </remarks>
        public static RectangularMatrix operator * (AnyRectangularMatrix A, AnyRectangularMatrix B) {
            if (A == null) throw new ArgumentNullException("A");
            if (B == null) throw new ArgumentNullException("B");
            if (A.ColumnCount != B.RowCount) throw new DimensionMismatchException();
            RectangularMatrix M = new RectangularMatrix(A.RowCount, B.ColumnCount);
            for (int i = 0; i < A.RowCount; i++) {
                for (int j = 0; j < B.ColumnCount; j++) {
                    for (int k = 0; k < B.RowCount; k++) {
                        M[i, j] += A[i, k] * B[k, j];
                    }
                }
            }
            return (M);
        }

        /*
        /// <summary>
        /// Multiplies any real, rectangular matrix by a column vector.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="x">The column vector.</param>
        /// <returns>The column vector Ax.</returns>
        /// <remarks>
        /// <para>Matrix multiplication is an O(N<sup>2</sup>) process.</para>
        /// </remarks>
        public static ColumnVector operator * (RectangularMatrixBase A, ColumnVector x) {
            if (A == null) throw new ArgumentNullException("A");
            if (x == null) throw new ArgumentNullException("x");
            if (A.RowCount != x.Dimension) throw new DimensionMismatchException();
            ColumnVector Ax = new ColumnVector(A.ColumnCount);
            for (int i = 0; i < A.ColumnCount; i++) {
                for (int j = 0; j < A.RowCount; j++) {
                    Ax[i] += A[i, j] * x[j];
                }
            }
            return (Ax);
        }
        */

        /// <summary>
        /// Multiplies any real, rectangular matrix by a real constant.
        /// </summary>
        /// <param name="alpha">The constant.</param>
        /// <param name="A">The matrix.</param>
        /// <returns>The product matrix.</returns>
        public static RectangularMatrix operator * (double alpha, AnyRectangularMatrix A) {
            if (A == null) throw new ArgumentNullException("A");
            RectangularMatrix B = new RectangularMatrix(A.RowCount, A.ColumnCount);
            for (int i = 0; i < A.RowCount; i++) {
                for (int j = 0; j < A.ColumnCount; j++) {
                    B[i, j] = alpha * A[i, j];
                }
            }
            return (B);
        }

        /// <summary>
        /// Multiplies any real, rectangular matrix with a real column vector.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The product column vector.</returns>
        public static ColumnVector operator * (AnyRectangularMatrix A, ColumnVector v) {
            if (A == null) throw new ArgumentNullException("A");
            if (v == null) throw new ArgumentNullException("v");
            if (A.ColumnCount != v.Dimension) throw new DimensionMismatchException();
            ColumnVector Av = new ColumnVector(A.RowCount);

            for (int r = 0; r < A.RowCount; r++) {
                for (int c = 0; c < A.ColumnCount; c++) {
                    Av[r] += A[r, c] * v[c];
                }
            }

            return (Av);
        }

        // equality operations
        // do not re-implement for concrete implementations, because people shouldn't be doing matrix equality comparisons anyway

        /// <summary>
        /// Determines whether two matrices are equal.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>True if <paramref name="A"/> and <paramref name="B"/> are equal, otherwise false.</returns>
        public static bool operator == (AnyRectangularMatrix A, AnyRectangularMatrix B) {
            if ((object) A == null) {
                if ((object) B == null) {
                    return (true);
                } else {
                    return (false);
                }
            } else {
                if ((object) B == null) {
                    return (false);
                } else {
                    if (A.RowCount != B.RowCount) return(false);
                    if (A.ColumnCount != B.ColumnCount) return(false);
                    for (int r = 0; r < A.RowCount; r++) {
                        for (int c = 0; c < A.ColumnCount; c++) {
                            if (A[r, c] != B[r, c]) return (false);
                        }
                    }
                    return (true);
                }
            }
        }

        /// <summary>
        /// Determines whether two matrices are not equal.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>False if <paramref name="A"/> and <paramref name="B"/> are equal, otherwise true.</returns>
        public static bool operator != (AnyRectangularMatrix A, AnyRectangularMatrix B) {
            return (!(A == B));
        }

        /// <summary>
        /// Determines whether the given object is an equal matrix.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if <paramref name="obj"/> is an equal matrix, otherwise false.</returns>
        public override bool Equals (object obj) {
            AnyRectangularMatrix M = obj as RectangularMatrix;
            if (obj == null) {
                return (false);
            } else {
                return ((this == M));
            }
        }

        /// <summary>
        /// Not a valid operation.
        /// </summary>
        /// <returns>Throws an <see cref="InvalidOperationException"/>.</returns>
        /// <remarks>
        /// <para>The <see cref="Object.GetHashCode"/> method is used to provide a quick equality test when an object
        /// is used as a key in a <see cref="System.Collections.Generic.Dictionary{TKey,TValue}"/> or <see cref="System.Collections.Hashtable"/>.
        /// Since a useful hash code of a matrix would need to involve all its elements, it is not possible to make this a fast operation.
        /// Also, since matrices are not immutable, they should not be used as hash keys. (A matrix might be changed after it
        /// had already been used as a key.) For these reasons, requesting a hash code for a matrix is not supported.
        /// </para>
        /// </remarks>
        /// <exception cref="NotSupportedException">This method always throws this exception.</exception>
        public override int GetHashCode () {
            throw new NotSupportedException();
        }

#if SHO
        /// <summary>
        /// Produces the representation of the matrix for the Python interactive console.
        /// </summary>
        /// <returns>A string representation of the matrix.</returns>
        public string __repr__ () {
            using (StringWriter writer = new StringWriter(CultureInfo.CurrentCulture)) {
                Write(writer);
                return (writer.ToString());
            }
        }
#endif

        internal void Write (TextWriter writer) {

            for (int r = 0; r < RowCount; r++) {
                writer.Write("{ ");
                for (int c = 0; c < ColumnCount; c++) {
                    writer.Write("{0,15:g12} ", this[r, c]);
                }
                writer.WriteLine("}");
            }

        }

    }


}
