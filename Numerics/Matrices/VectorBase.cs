using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Implements functionality shared between row and column vectors.
    /// </summary>
    /// <seealso cref="RowVector"/>
    /// <seealso cref="ColumnVector"/>
    public abstract class VectorBase : AnyRectangularMatrix, IEnumerable, IEnumerable<double>, ICollection<double>, IList<double>, IReadOnlyCollection<double>, IReadOnlyList<double> {

        internal VectorBase (double[] store, int offset, int stride, int dimension, bool isReadOnly) : base(isReadOnly) {
            this.store = store;
            this.offset = offset;
            this.stride = stride;
            this.dimension = dimension;
        }

        internal VectorBase (double[] store, int dimension, bool isReadOnly) : this(store, 0, 1, dimension, isReadOnly) { }

        internal VectorBase (double[] store, int dimension) : this(store, dimension, false) { }

        internal VectorBase (int dimension) : this(new double[dimension], dimension) { }

        internal VectorBase (IReadOnlyList<double> list) {
            if (list == null) throw new ArgumentNullException(nameof(list));
            dimension = list.Count;
            store = list.ToArray();
            offset = 0;
            stride = 1;
        }

        // this storage is internal so that other routines can access it for fast operations (e.g. multiplication)

        internal int dimension;
        internal double[] store;
        internal int offset, stride;

        /// <summary>
        /// Gets or sets the specified vector component.
        /// </summary>
        /// <param name="index">The (zero-based) component index.</param>
        /// <returns>The value of the specified vector component.</returns>
        public virtual double this[int index] {
            get {
                if ((index < 0) || (index >= dimension)) throw new ArgumentOutOfRangeException(nameof(index));
                return (store[offset + stride * index]);
            }
            set {
                if ((index < 0) || (index >= dimension)) throw new ArgumentOutOfRangeException(nameof(index));
                if (IsReadOnly) throw new InvalidOperationException();
                store[offset + stride * index] = value;
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
            return(Blas1.dNrm2(store, offset, stride, dimension));
        }

        /// <summary>
        /// Returns the vector elements in an independent array.
        /// </summary>
        /// <returns>An array containing the vector element values.</returns>
        public virtual new double[] ToArray () {
            return(VectorAlgorithms.Copy(store, offset, stride, dimension));
        }

        // interface implementations

        /// <summary>
        /// Gets an enumerator of the vector components.
        /// </summary>
        /// <returns>An enumerator of the vector components.</returns>
        public IEnumerator<double> GetEnumerator () {
            for (int index = 0; index < dimension; index++) {
                yield return (this[index]);
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (store.GetEnumerator());
        }

        int ICollection<double>.Count {
            get {
                return (dimension);
            }
        }

        int IReadOnlyCollection<double>.Count {
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
            for (int index = 0; index < dimension; index++) {
                if (this[index] == item) return (true);
            }
            return (false);
        }

        void ICollection<double>.CopyTo (double[] array, int arrayIndex) {
            Blas1.dCopy(store, offset, stride, array, arrayIndex, 1, dimension);
        }

        bool ICollection<double>.Remove (double item) {
            throw new NotSupportedException();
        }

        int IList<double>.IndexOf (double item) {
            for (int index = 0; index < dimension; index++) {
                if (this[index] == item) return(index);
            }
            return(-1);
        }

        void IList<double>.Insert (int index, double item) {
            throw new NotSupportedException();
        }

        void IList<double>.RemoveAt (int index) {
            throw new NotSupportedException();
        }

    }

}
