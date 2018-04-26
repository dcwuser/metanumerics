using System;
using System.Collections.Generic;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// A column vector of real numbers.
    /// </summary>
    /// <remarks>
    /// <para>An N-dimensional column vector is an N X 1 dimensional matrix.</para>
    /// </remarks>
    public sealed class ColumnVector : AnyVector {

        /// <summary>
        /// Initializes a new column vector with the given dimension.
        /// </summary>
        /// <param name="dimension">The dimension of the vector, which must be positive.</param>
        public ColumnVector (int dimension) : base(dimension) { }

        /// <summary>
        /// Initializes a new column vector from the given component list.
        /// </summary>
        /// <param name="list">A list of vector components.</param>
        public ColumnVector (IReadOnlyList<double> list) : base(list) { }

        /// <summary>
        /// Initializes a new column vector with the given components.
        /// </summary>
        /// <param name="list">A list of vector components.</param>
        public ColumnVector (params double[] list) : base(list) { }

        internal ColumnVector (double[] store, int offset, int stride, int dimension, bool isReadOnly) : base(store, offset, stride, dimension, isReadOnly) { }

        internal ColumnVector (double[] store, int dimension) : base(store, dimension) { }


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
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
                if (c != 0) throw new ArgumentOutOfRangeException(nameof(c));
                return (base[r]);
            }
            set {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
                if (c != 0) throw new ArgumentOutOfRangeException(nameof(c));
                base[r] = value;
            }
        }

        /// <summary>
        /// Returns the transpose of the column vector.
        /// </summary>
        /// <value>A row vector with the same components as the column vector.</value>
        public RowVector Transpose {
            get {
                return (new RowVector(store, offset, stride, dimension, true));
            }
        }

        /// <summary>
        /// Returns a copy of the column vector.
        /// </summary>
        /// <returns>An independent copy of the column vector.</returns>
        public ColumnVector Copy () {
            double[] copy = VectorAlgorithms.Copy(store, offset, stride, dimension);
            return (new ColumnVector(copy, dimension));
        }

        /// <summary>
        /// Computes the sum of two column vectors.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>The sum <paramref name="v1"/> + <paramref name="v2"/>.</returns>
        public static ColumnVector operator + (ColumnVector v1, ColumnVector v2) {
            if (v1 == null) throw new ArgumentNullException(nameof(v1));
            if (v2 == null) throw new ArgumentNullException(nameof(v2));
            if (v1.dimension != v2.dimension) throw new DimensionMismatchException();
            double[] store = VectorAlgorithms.Add(v1.store, v1.offset, v1.stride, v2.store, v2.offset, v2.stride, v1.dimension);
            return (new ColumnVector(store, v1.dimension));
        }

        /// <summary>
        /// Computes the difference of two column vectors.
        /// </summary>
        /// <param name="v1">The first column vector.</param>
        /// <param name="v2">The second column vector.</param>
        /// <returns>The difference <paramref name="v1"/> - <paramref name="v2"/>.</returns>
        public static ColumnVector operator - (ColumnVector v1, ColumnVector v2) {
            if (v1 == null) throw new ArgumentNullException(nameof(v1));
            if (v2 == null) throw new ArgumentNullException(nameof(v2));
            if (v1.dimension != v2.dimension) throw new DimensionMismatchException();
            double[] store = VectorAlgorithms.Subtract(v1.store, v1.offset, v1.stride, v2.store, v2.offset, v2.stride, v1.dimension);
            return (new ColumnVector(store, v1.dimension));
        }

        /// <summary>
        /// Multiplies a column vector by a real, scalar constant.
        /// </summary>
        /// <param name="alpha">The real, scalar constant.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The product &#x3B1;v.</returns>
        public static ColumnVector operator * (double alpha, ColumnVector v) {
            if (v == null) throw new ArgumentNullException(nameof(v));
            double[] store = VectorAlgorithms.Multiply(alpha, v.store, v.offset, v.stride, v.dimension);
            return (new ColumnVector(store, v.dimension));
        }

        /// <summary>
        /// Divides a column vector by a real, scalar constant.
        /// </summary>
        /// <param name="alpha">The real, scalar constant.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The quotient v / &#x3B1;.</returns>
        public static ColumnVector operator / (ColumnVector v, double alpha) {
            if (v == null) throw new ArgumentNullException(nameof(v));
            double[] store = VectorAlgorithms.Multiply(1.0 / alpha, v.store, v.offset, v.stride, v.dimension);
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

    }

}
