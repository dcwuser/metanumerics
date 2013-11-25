using System;
using System.Collections.Generic;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// A row vector of real numbers.
    /// </summary>
    /// <remarks>
    /// <para>An N-dimensional row vector is an 1 X N dimensional matrix.</para>
    /// </remarks>
    public sealed class RowVector : VectorBase {

        /// <summary>
        /// Initializes a new row vector with the given dimension.
        /// </summary>
        /// <param name="dimension">The dimension of the vector, which must be positive.</param>
        public RowVector (int dimension)
            : base(dimension) {
        }

        /// <summary>
        /// Initializes a new row vector from the given component list.
        /// </summary>
        /// <param name="list">A list of vector components.</param>
        public RowVector (IList<double> list)
            : base(list) {
        }

        /// <summary>
        /// Initializes a new row vector with the given components.
        /// </summary>
        /// <param name="list">A list of vector components.</param>
        public RowVector (params double[] list)
            : base(list) {
        }

        internal RowVector (double[] store, int dim)
            : base(store, dim) {
        }

        /// <inheritdoc />
        public override int RowCount {
            get { return (1); }
        }

        /// <inheritdoc />
        public override int ColumnCount {
            get { return (dimension); }
        }

        /// <inheritdoc />
        public override double this[int r, int c] {
            get {
                if (r != 0) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
                return (store[c]);
            }
            set {
                if (r != 0) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
                store[c] = value;
            }
        }

        /// <summary>
        /// Returns the transpose of the row vector.
        /// </summary>
        /// <returns>An independent column vector with the same components as the row vector.</returns>
        public ColumnVector Transpose () {
            double[] copy = VectorAlgorithms.Copy(store, dimension);
            return (new ColumnVector(copy, dimension));
        }

        /// <summary>
        /// Returns a copy of the row vector.
        /// </summary>
        /// <returns>An independent copy of the row vector.</returns>
        public RowVector Copy () {
            double[] copy = VectorAlgorithms.Copy(store, dimension);
            return (new RowVector(copy, dimension));
        }

        /// <summary>
        /// Computes the sum of two row vectors.
        /// </summary>
        /// <param name="v1">The first row vector.</param>
        /// <param name="v2">The second row vector.</param>
        /// <returns>The sum <paramref name="v1"/> + <paramref name="v2"/>.</returns>
        public static RowVector operator + (RowVector v1, RowVector v2) {
            if (v1 == null) throw new ArgumentNullException("v1");
            if (v2 == null) throw new ArgumentNullException("v2");
            if (v1.dimension != v2.dimension) throw new DimensionMismatchException();
            double[] store = VectorAlgorithms.Add(v1.store, v2.store, v1.dimension);
            return (new RowVector(store, v1.dimension));
        }

        /// <summary>
        /// Computes the difference of two column vectors.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>The difference <paramref name="v1"/> - <paramref name="v2"/>.</returns>
        public static RowVector operator - (RowVector v1, RowVector v2) {
            if (v1 == null) throw new ArgumentNullException("v1");
            if (v2 == null) throw new ArgumentNullException("v2");
            if (v1.dimension != v2.dimension) throw new DimensionMismatchException();
            double[] store = VectorAlgorithms.Subtract(v1.store, v2.store, v1.dimension);
            return (new RowVector(store, v1.dimension));
        }

        /// <summary>
        /// Multiplies a row vector by a real, scalar constant.
        /// </summary>
        /// <param name="alpha">The real, scalar constant.</param>
        /// <param name="v">The row vector.</param>
        /// <returns>The product.</returns>
        public static RowVector operator * (double alpha, RowVector v) {
            if (v == null) throw new ArgumentNullException("v");
            double[] store = VectorAlgorithms.Multiply(alpha, v.store, v.dimension);
            return (new RowVector(store, v.dimension));
        }

        /// <summary>
        /// Divides a row vector by a real, scalar constant.
        /// </summary>
        /// <param name="alpha">The real, scalar constant.</param>
        /// <param name="v">The row vector.</param>
        /// <returns>The result.</returns>
        public static RowVector operator / (RowVector v, double alpha) {
            if (v == null) throw new ArgumentNullException("v");
            double[] store = VectorAlgorithms.Multiply(1.0 / alpha, v.store, v.dimension);
            return (new RowVector(store, v.dimension));
        }

        /// <summary>
        /// Negates a row vector.
        /// </summary>
        /// <param name="v">The row vector.</param>
        /// <returns>-v</returns>
        public static RowVector operator - (RowVector v) {
            return (-1.0 * v);
        }

        /// <summary>
        /// Multiplies any real, rectangular matrix by a row vector.
        /// </summary>
        /// <param name="v">The row vector.</param>
        /// <param name="A">The matrix.</param>
        /// <returns>The product row vector.</returns>
        public static RowVector operator * (RowVector v, AnyRectangularMatrix A) {
            if (v == null) throw new ArgumentNullException("v");
            if (A == null) throw new ArgumentNullException("A");
            if (v.Dimension != A.RowCount) throw new DimensionMismatchException();
            RowVector vA = new RowVector(A.ColumnCount);

            for (int r = 0; r < A.RowCount; r++) {
                for (int c = 0; c < A.ColumnCount; c++) {
                    vA[c] += v[r] * A[r, c];
                }
            }

            return (vA);
        }

        /// <summary>
        /// Computes the inner (scalar or dot) product of a row and a column vector.
        /// </summary>
        /// <param name="v">The row vector.</param>
        /// <param name="u">The column vector.</param>
        /// <returns>The value of the scalar product.</returns>
        public static double operator * (RowVector v, ColumnVector u) {
            if (v == null) throw new ArgumentNullException("v");
            if (u == null) throw new ArgumentNullException("u");
            if (v.dimension != u.dimension) throw new DimensionMismatchException();
            double p = Blas1.dDot(v.store, 0, 1, u.store, 0, 1, v.dimension);
            return (p);
        }

    }

}
