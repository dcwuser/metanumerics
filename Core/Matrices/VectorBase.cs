using System;
using System.Collections;
using System.Collections.Generic;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Implements functionality shared between row and column vectors.
    /// </summary>
    /// <seealso cref="RowVector"/>
    /// <seealso cref="ColumnVector"/>
    public abstract class VectorBase : RectangularMatrixBase, IEnumerable, IEnumerable<double>, ICollection<double>, IList<double> {

        internal VectorBase (double[] store, int dimension) {
            if (dimension <= 0) throw new ArgumentOutOfRangeException("dimension");
            this.dimension = dimension;
            this.store = store;
        }

        internal VectorBase (int dimension) : this(new double[dimension], dimension) {
        }

        // this storage is internal so that other routines can access it for fast operations (e.g. multiplication)

        internal int dimension;
        internal double[] store;

        /// <summary>
        /// Gets or sets the specified vector component.
        /// </summary>
        /// <param name="index">The (zero-based) component index.</param>
        /// <returns>The value of the specified vector component.</returns>
        public virtual double this[int index] {
            get {
                if ((index < 0) || (index >= dimension)) throw new ArgumentOutOfRangeException("index");
                return (store[index]);
            }
            set {
                if ((index < 0) || (index >= dimension)) throw new ArgumentOutOfRangeException("index");
                store[index] = value;
            }
        }

        /// <summary>
        /// Gets the dimension of the vector.
        /// </summary>
        public virtual int Dimension {
            get {
                return (dimension);
            }
        }

        // interface implementations

        IEnumerator<double> IEnumerable<double>.GetEnumerator () {
            return ( ((IEnumerable<double>) store).GetEnumerator() );
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (store.GetEnumerator());
        }

        int ICollection<double>.Count {
            get {
                return (dimension);
            }
        }

        bool ICollection<double>.IsReadOnly {
            get {
                // obviously we are not really read-only, but this is consistent with array behavior
                // and, like array, we do not support add and remove operations
                return (true);
            }
        }

        void ICollection<double>.Add (double item) {
            throw new NotSupportedException();
        }

        void ICollection<double>.Clear () {
            throw new NotSupportedException();
        }

        bool ICollection<double>.Contains (double item) {
            return (((ICollection<double>)store).Contains(item));
        }

        void ICollection<double>.CopyTo (double[] array, int arrayIndex) {
            Blas1.dCopy(store, 0, 1, array, arrayIndex, 1, dimension);
        }

        bool ICollection<double>.Remove (double item) {
            throw new NotSupportedException();
        }

        int IList<double>.IndexOf (double item) {
            return (((IList<double>)store).IndexOf(item));
        }

        void IList<double>.Insert (int index, double item) {
            throw new NotSupportedException();
        }

        void IList<double>.RemoveAt (int index) {
            throw new NotSupportedException();
        }

    }

    internal static class VectorAlgorithms {

        public static double[] Copy (double[] store, int dimension) {
            double[] copy = new double[dimension];
            Blas1.dCopy(store, 0, 1, copy, 0, 1, dimension);
            return (copy);
        }

        public static double[] Add (double[] aStore, double[] bStore, int dimension) {
            double[] store = new double[dimension];
            for (int i = 0; i < dimension; i++) {
                store[i] = aStore[i] + bStore[i];
            }
            return(store);
        }

        public static double[] Subtract (double[] aStore, double[] bStore, int dimension) {
            double[] store = new double[dimension];
            for (int i = 0; i < dimension; i++) {
                store[i] = aStore[i] - bStore[i];
            }
            return (store);
        }

        public static double[] Multiply (double alpha, double[] store, int dimension) {
            double[] product = new double[dimension];
            Blas1.dAxpy(alpha, store, 0, 1, product, 0, 1, dimension);
            return (product);
        }


    }

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
        public ColumnVector (int dimension) : base(dimension) {
        }

        public ColumnVector (double[] source) : base(source.Length) {
            Blas1.dCopy(source, 0, 1, store, 0, 1, dimension);
        }

        internal ColumnVector (double[] store, int dimension) : base(store, dimension) {
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
        public ColumnVector Clone () {
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
        public RowVector (int dimension) : base(dimension) {
        }

        public RowVector (double[] source) : base(source.Length) {
            Blas1.dCopy(source, 0, 1, store, 0, 1, dimension);
        }

        internal RowVector (double[] store, int dim) : base(store, dim) {
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
        public RowVector Clone () {
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

        public static RowVector operator * (RowVector v, RectangularMatrixBase A) {

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
