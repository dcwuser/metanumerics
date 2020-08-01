using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a column vector of complex numbers.
    /// </summary>
    public sealed class ComplexColumnVector : AnyMatrix<Complex>, IReadOnlyList<Complex>
    {

        /// <summary>
        /// Initializes a new complex column vector with the given dimension.
        /// </summary>
        /// <param name="dimension">The dimension of the vector, which must be positive.</param>
        public ComplexColumnVector(int dimension)
        {
            if (dimension < 1) throw new ArgumentOutOfRangeException(nameof(dimension));
            this.store = new Complex[dimension];
            this.offset = 0;
            this.dimension = dimension;
        }

        /// <summary>
        /// Initializes a new complex column vector with the given components.
        /// </summary>
        /// <param name="components">The components of the vector.</param>
        public ComplexColumnVector(params Complex[] components) {
            if (components == null) throw new ArgumentNullException(nameof(components));
            this.store = components;
            this.offset = 0;
            this.dimension = components.Length;
        }

        internal ComplexColumnVector(Complex[] store, int offset, int dimension, bool isReadOnly) : base(isReadOnly)
        {
            Debug.Assert(store != null);
            Debug.Assert(dimension > 0);
            this.store = store;
            this.offset = offset;
            this.dimension = dimension;
        }

        private readonly Complex[] store;
        private readonly int offset;
        private readonly int dimension;

        /// <summary>
        /// Gets the dimension of the complex column vector.
        /// </summary>
        /// <value>The number of components in the vector.</value>
        public int Dimension
        {
            get
            {
                return (dimension);
            }
        }

        /// <summary>
        /// Gets or sets the specified vector component.
        /// </summary>
        /// <param name="index">The (zero-based) component index.</param>
        /// <value>The value of the specified vector component.</value>
        public Complex this [int index]
        {
            get
            {
                if ((index < 0) || (index >= dimension)) throw new ArgumentOutOfRangeException(nameof(index));
                return (store[offset + index]);
            }
            set
            {
                if ((index < 0) || (index >= dimension)) throw new ArgumentOutOfRangeException(nameof(index));
                if (IsReadOnly) throw new InvalidOperationException();
                store[offset + index] = value;
            }
        }

        /// <inheritdoc />
        public override Complex this[int r, int c]
        {
            get
            {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
                if (c != 0) throw new ArgumentOutOfRangeException(nameof(c));
                return (this[r]);
            }

            set
            {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
                if (c != 0) throw new ArgumentOutOfRangeException(nameof(c));
                this[r] = value;
            }
        }

        /// <inheritdoc />
        public override int ColumnCount
        {
            get
            {
                return (1);
            }
        }

        /// <inheritdoc />
        public override int RowCount
        {
            get
            {
                return (dimension);
            }
        }

        int IReadOnlyCollection<Complex>.Count
        {
            get
            {
                return (dimension);
            }
        }

        /// <summary>
        /// Returns a copy of the complex column vector.
        /// </summary>
        /// <returns>An independent copy of the complex column vector.</returns>
        public ComplexColumnVector Copy ()
        {
            Complex[] copy = new Complex[dimension];
            Blas1.Copy(store, offset, 1, copy, 0, 1, dimension);
            return (new ComplexColumnVector(copy, 0, dimension, false));
        }

        IEnumerator<Complex> IEnumerable<Complex>.GetEnumerator()
        {
            for (int index = 0; index < dimension; index++)
            {
                yield return (this[index]);
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (((IEnumerable<Complex>)this).GetEnumerator());
        }

        /// <summary>
        /// Multiplies a complex column vector by a complex constant.
        /// </summary>
        /// <param name="alpha">The complex constant.</param>
        /// <param name="v">The complex column vector.</param>
        /// <returns>The product &#x3B1;v.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="v"/> is null.</exception>
        public static ComplexColumnVector operator * (Complex alpha, ComplexColumnVector v) {
            if (v == null) throw new ArgumentNullException(nameof(v));
            Complex[] product = new Complex[v.dimension];
            Blas1.zAxpy(alpha, v.store, v.offset, 1, product, 0, 1, v.dimension);
            return (new ComplexColumnVector(product, 0, v.dimension, false));
        }

        /// <summary>
        /// Multiplies a real matrix by a complex column vector.
        /// </summary>
        /// <param name="A">The real matrix.</param>
        /// <param name="v">The complex column vector.</param>
        /// <returns>The product Av.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="v"/> is null.</exception>
        /// <exception cref="DimensionMismatchException">The column count of <paramref name="A"/> does not equal the dimension of <paramref name="v"/>.</exception>
        public static ComplexColumnVector operator * (AnyRectangularMatrix A, ComplexColumnVector v) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (v == null) throw new ArgumentNullException(nameof(v));
            if (A.ColumnCount != v.Dimension) throw new DimensionMismatchException();
            Complex[] product = new Complex[A.RowCount];
            for (int i = 0; i < A.RowCount; i++) {
                Complex s = 0.0;
                for (int j = 0; j < A.ColumnCount; j++) {
                    s += A[i, j] * v.store[j];
                }
                product[i] = s;
            }
            return (new ComplexColumnVector(product, 0, product.Length, false));
        }

    }

}
