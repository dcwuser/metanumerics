using System;
using System.Collections;
using System.Collections.Generic;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Implements functionality shared between row and column vectors.
    /// </summary>
    /// <seealso cref="RowVector"/>
    /// <seealso cref="ColumnVector"/>
    public abstract class VectorBase : AnyRectangularMatrix, IEnumerable, IEnumerable<double>, ICollection<double>, IList<double> {

        internal VectorBase (double[] store, int dimension) {
            if (dimension <= 0) throw new ArgumentOutOfRangeException("dimension");
            this.dimension = dimension;
            this.store = store;
        }

        internal VectorBase (int dimension) : this(new double[dimension], dimension) {
        }

        internal VectorBase (IList<double> list) {
            if (list == null) throw new ArgumentNullException("list");
            dimension = list.Count;
            store = new double[dimension];
            list.CopyTo(store, 0);
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

        /// <summary>
        /// Computes the magnitude of the vector.
        /// </summary>
        /// <returns>The Euclidean norm of the vector.</returns>
        public virtual double Norm () {
            return(Blas1.dNrm2(store, 0, 1, dimension));
        }

        /// <summary>
        /// Returns the vector elements in an independent array.
        /// </summary>
        /// <returns>An array containing the vector element values.</returns>
        public virtual new double[] ToArray () {
            return(VectorAlgorithms.Copy(store, dimension));
        }

        // interface implementations

        /// <summary>
        /// Gets an enumerator of the vector components.
        /// </summary>
        /// <returns>An enumerator of the vector components.</returns>
        public IEnumerator<double> GetEnumerator () {
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

}
