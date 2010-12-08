using System;


namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Describes all matrices.
    /// </summary>
    /// <typeparam name="T">The type of the matrix entries.</typeparam>
    public abstract class MatrixBase<T> {

        public abstract int RowCount { get; }

        public abstract int ColumnCount { get; }

        public abstract T this[int r, int c] { get; set; }

    }

}
