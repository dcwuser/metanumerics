using System;
using System.Collections.Generic;


namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Describes the form of all real, square matrices.
    /// </summary>
    /// <remarks>
    /// <para>This is an abstract class that describes any real, square matrix. If you wish to create a concrete
    /// instance of a real, non-square matrix, use the <see cref="SquareMatrix"/> class. If, on the
    /// other hand, you wish to write a function that can operate on any real, square matrix, it's probably a good
    /// idea to accept a <see cref="AnySquareMatrix"/>, so that your function could operate on any concrete implementation.</para>
    /// </remarks>
    public abstract class AnySquareMatrix : AnyRectangularMatrix {

        /// <summary>
        /// Initializes a new instance of the AnySquareMatrix class.
        /// </summary>
        protected AnySquareMatrix () : base() { }

        internal AnySquareMatrix (bool isReadOnly) : base(isReadOnly) { }

        /// <summary>
        /// Gets or sets the dimension of the square matrix.
        /// </summary>
        public abstract int Dimension { get; }

        /// <inheritdoc />
        public override int RowCount {
            get { return (Dimension); }
        }

        /// <inheritdoc />
        public override int ColumnCount {
            get { return (Dimension); }
        }

        /// <summary>
        /// Computes the trace of the square matrix.
        /// </summary>
        /// <returns>tr(M)</returns>
        /// <remarks>
        /// <para>The trace of a square matrix is the sum of its diagonal elements.</para>
        /// </remarks>
        public virtual double Trace () {
            double tr = 0.0;
            for (int i = 0; i < Dimension; i++) {
                tr += this[i, i];
            }
            return (tr);
        }

        /// <summary>
        /// Adds any two real, square matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The sum matrix A + B.</returns>
        /// <remarks>
        /// <para>Matrix addition is an O(N<sup>2</sup>) process.</para>
        /// </remarks>
        public static SquareMatrix operator + (AnySquareMatrix A, AnySquareMatrix B) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.Dimension != B.Dimension) throw new DimensionMismatchException();
            SquareMatrix M = new SquareMatrix(A.Dimension);
            for (int r = 0; r < M.Dimension; r++) {
                for (int c = 0; c < M.Dimension; c++) {
                    M[r, c] = A[r, c] + B[r, c];
                }
            }
            return (M);
        }

        /// <summary>
        /// Subtracts any two real, square matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The difference matrix A - B.</returns>
        /// <remarks>
        /// <para>Matrix addition is an O(N<sup>2</sup>) process.</para>
        /// </remarks>
        public static SquareMatrix operator - (AnySquareMatrix A, AnySquareMatrix B) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.Dimension != B.Dimension) throw new DimensionMismatchException();
            SquareMatrix M = new SquareMatrix(A.Dimension);
            for (int r = 0; r < M.Dimension; r++) {
                for (int c = 0; c < M.Dimension; c++) {
                    M[r, c] = A[r, c] - B[r, c];
                }
            }
            return (M);
        }

        /// <summary>
        /// Multiplies any two real, square matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The product matrix AB.</returns>
        /// <remarks>
        /// <para>Matrix multiplication is an O(N<sup>3</sup>) process.</para>
        /// </remarks>
        public static SquareMatrix operator * (AnySquareMatrix A, AnySquareMatrix B) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.Dimension != B.Dimension) throw new DimensionMismatchException();
            SquareMatrix M = new SquareMatrix(A.Dimension);
            for (int r = 0; r < M.Dimension; r++) {
                for (int c = 0; c < M.Dimension; c++) {
                    for (int k = 0; k < M.Dimension; k++) {
                        M[r, c] += A[r, k] * B[k, c];
                    }
                }
            }
            return (M);
        }

    }

}
