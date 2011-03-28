using System;
using System.Collections.Generic;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// A column vector of real numbers.
    /// </summary>
    /// <remarks>
    /// <para>An N-dimensional column vector is an N X 1 dimensional matrix.</para>
    /// </remarks>
    public sealed class ColumnVector : VectorBase {

        /// <summary>
        /// Initializes a new column vector with the given dimension.
        /// </summary>
        /// <param name="dimension">The dimension of the vector, which must be positive.</param>
        public ColumnVector (int dimension)
            : base(dimension) {
        }

        /// <summary>
        /// Initializes a new column vector from the given component list.
        /// </summary>
        /// <param name="list">A list of vector components.</param>
        public ColumnVector (IList<double> list)
            : base(list) {
        }

        /// <summary>
        /// Initializes a new column vector with the given components.
        /// </summary>
        /// <param name="list">A list of vector components.</param>
        public ColumnVector (params double[] list)
            : base(list) {
        }

        internal ColumnVector (double[] store, int dimension)
            : base(store, dimension) {
        }


        /// <inheritdoc />
        public override int RowCount {
            get { return (dimension); }
        }

        /// <inheritdoc />
        public override int ColumnCount {
            get { return (1); }
        }

        /// <inheritdoc />
        public override double this[int r, int c] {
            get {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
                if (c != 0) throw new ArgumentOutOfRangeException("c");
                return (store[r]);
            }
            set {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
                if (c != 0) throw new ArgumentOutOfRangeException("c");
                store[r] = value;
            }
        }

        /// <summary>
        /// Returns the transpose of the column vector.
        /// </summary>
        /// <returns>An independent row vector with the same components as the column vector.</returns>
        public RowVector Transpose () {
            double[] copy = VectorAlgorithms.Copy(store, dimension);
            return (new RowVector(copy, dimension));
        }

        /// <summary>
        /// Returns a copy of the column vector.
        /// </summary>
        /// <returns>An independent copy of the column vector.</returns>
        public ColumnVector Copy () {
            double[] copy = VectorAlgorithms.Copy(store, dimension);
            return (new ColumnVector(copy, dimension));
        }

        /// <summary>
        /// Computes the sum of two column vectors.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>The sum <paramref name="v1"/> + <paramref name="v2"/>.</returns>
        public static ColumnVector operator + (ColumnVector v1, ColumnVector v2) {
            if (v1 == null) throw new ArgumentNullException("v1");
            if (v2 == null) throw new ArgumentNullException("v2");
            if (v1.dimension != v2.dimension) throw new DimensionMismatchException();
            double[] store = VectorAlgorithms.Add(v1.store, v2.store, v1.dimension);
            return (new ColumnVector(store, v1.dimension));
        }

        /// <summary>
        /// Computes the difference of two column vectors.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>The difference <paramref name="v1"/> - <paramref name="v2"/>.</returns>
        public static ColumnVector operator - (ColumnVector v1, ColumnVector v2) {
            if (v1 == null) throw new ArgumentNullException("v1");
            if (v2 == null) throw new ArgumentNullException("v2");
            if (v1.dimension != v2.dimension) throw new DimensionMismatchException();
            double[] store = VectorAlgorithms.Subtract(v1.store, v2.store, v1.dimension);
            return (new ColumnVector(store, v1.dimension));
        }

        /// <summary>
        /// Multiplies a column vector by a real, scalar constant.
        /// </summary>
        /// <param name="alpha">The real, scalar constant.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The product.</returns>
        public static ColumnVector operator * (double alpha, ColumnVector v) {
            if (v == null) throw new ArgumentNullException("v");
            double[] store = VectorAlgorithms.Multiply(alpha, v.store, v.dimension);
            return (new ColumnVector(store, v.dimension));
        }

        /// <summary>
        /// Divides a column vector by a real, scalar constant.
        /// </summary>
        /// <param name="alpha">The real, scalar constant.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The product.</returns>
        public static ColumnVector operator / (ColumnVector v, double alpha) {
            if (v == null) throw new ArgumentNullException("v");
            double[] store = VectorAlgorithms.Multiply(1.0 / alpha, v.store, v.dimension);
            return (new ColumnVector(store, v.dimension));
        }

        /// <summary>
        /// Negates a column vector.
        /// </summary>
        /// <param name="v">The column vector.</param>
        /// <returns>-v</returns>
        public static ColumnVector operator - (ColumnVector v) {
            return (-1.0 * v);
        }

        /*
        public static RectangularMatrix operator * (ColumnVector u, RowVector v) {
            RectangularMatrix M = new RectangularMatrix(u.Dimension, v.Dimension);
            for (int r = 0; r < u.Dimension; r++) {
                for (int c = 0; c < v.Dimension; c++) {
                    M[r, c] = u[r] * v[c];
                }
            }
            return (M);
        }
        */

    }

}
