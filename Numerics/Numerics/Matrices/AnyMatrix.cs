using System;


namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Describes the form of all matrices.
    /// </summary>
    /// <typeparam name="T">The type of the matrix entries.</typeparam>
    public abstract class AnyMatrix<T> {

        internal AnyMatrix (bool isReadOnly) {
            this.isReadOnly = isReadOnly;
        }

        protected AnyMatrix () : this(false) { }

        /// <summary>
        /// Gets the number of matrix rows.
        /// </summary>
        public abstract int RowCount { get; }

        /// <summary>
        /// Gets the number of matrix columns.
        /// </summary>
        public abstract int ColumnCount { get; }

        /// <summary>
        /// Gets or sets the value of a matrix entry.
        /// </summary>
        /// <param name="r">The (zero-based) row index.</param>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>The value of the <paramref name="r"/>,<paramref name="c"/> matrix entry.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> or <paramref name="c"/> is
        /// outside the valid range.</exception>
        public abstract T this[int r, int c] { get; set; }

        /// <summary>
        /// Sets all matrix entries according to a supplied fill function.
        /// </summary>
        /// <param name="f">The fill function.</param>
        public virtual void Fill (Func<int, int, T> f) {
            if (f == null) throw new ArgumentNullException("f");
            for (int r = 0; r < this.RowCount; r++) {
                for (int c = 0; c < this.ColumnCount; c++) {
                    this[r, c] = f(r, c);
                }
            }
        }

        /// <summary>
        /// Copies the matrix into an array.
        /// </summary>
        /// <returns>A two-dimensional .NET array containing the matrix entries.</returns>
        /// <remarks>
        /// <para>The output array is independent of the matrix. Changes to its elements will not change
        /// the elements of the matrix, and changes to the matrix elements will not be reflected in the array.</para>
        /// </remarks>
        public virtual T[,] ToArray () {
            T[,] result = new T[this.RowCount, this.ColumnCount];
            for (int r = 0; r < this.RowCount; r++) {
                for (int c = 0; c < this.ColumnCount; c++) {
                    result[r, c] = this[r, c];
                }
            }
            return (result);
        }


        /// <summary>
        /// Gets a flag indicating whether the matrix is read-only.
        /// </summary>
        /// <remarks>
        /// <para>Although you can't change the values in a read-only matrix, you can make a writable copy of it.</para>
        /// </remarks>
        public bool IsReadOnly {

            get {
                return (isReadOnly);
            }

            internal set {
                isReadOnly = value;
            }

        }

        private bool isReadOnly = false;

        // to speed access, IsReadOnly shouldn't be virtual

    }

}
