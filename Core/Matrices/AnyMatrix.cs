using System;


namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Describes the form of all matrices.
    /// </summary>
    /// <typeparam name="T">The type of the matrix entries.</typeparam>
    public abstract class AnyMatrix<T> {

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

    }

}
