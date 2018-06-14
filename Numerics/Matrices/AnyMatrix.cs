using System;
using System.Diagnostics;


namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Describes the form of all matrices.
    /// </summary>
    /// <typeparam name="T">The type of the matrix entries.</typeparam>
    public abstract class AnyMatrix<T> : IEquatable<AnyMatrix<T>> where T : IEquatable<T> {

        internal AnyMatrix (bool isReadOnly) {
            this.isReadOnly = isReadOnly;
        }

        /// <summary>
        /// Initializes a new instance of the AnyMatrix class.
        /// </summary>
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
            if (f == null) throw new ArgumentNullException(nameof(f));
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

        // Equality testing
        // It's not worthwhile to re-implement for child classes, because people shouldn't be doing matrix
        // equality comparisons in any code except perhaps test code anyway.

        private static bool InternalEquals (AnyMatrix<T> A, AnyMatrix<T> B) {
            Debug.Assert(!Object.ReferenceEquals(A, null));
            Debug.Assert(!Object.ReferenceEquals(B, null));
            if (A.RowCount != B.RowCount) return (false);
            if (A.ColumnCount != B.ColumnCount) return (false);
            for (int r = 0; r < A.RowCount; r++) {
                for (int c = 0; c < A.ColumnCount; c++) {
                    if (!A[r, c].Equals(B[r, c])) return (false);
                }
            }
            return (true);
        }

        /// <summary>
        /// Determines whether two matrices are equal.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>True if <paramref name="A"/> and <paramref name="B"/> are equal, otherwise false.</returns>
        public static bool operator == (AnyMatrix<T> A, AnyMatrix<T> B) {
            if (Object.ReferenceEquals(A, B)) {
                return (true);
            } else if (Object.ReferenceEquals(A, null) || Object.ReferenceEquals(B, null)) {
                return (false);
            } else {
                return (InternalEquals(A, B));
            }
        }

        /// <summary>
        /// Determines whether two matrices are not equal.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>False if <paramref name="A"/> and <paramref name="B"/> are equal, otherwise true.</returns>
        public static bool operator != (AnyMatrix<T> A, AnyMatrix<T> B) {
            return (!(A == B));
        }

        /// <summary>
        /// Determines whether the given matrix equals the current matrix.
        /// </summary>
        /// <param name="other">The matrix to compare.</param>
        /// <returns>True if the <paramref name="other"/> is equal to the current matrix, otherwise false.</returns>
        public bool Equals (AnyMatrix<T> other) {
            if (Object.ReferenceEquals(other, null)) {
                return (false);
            } else {
                return (InternalEquals(this, other));
            }
        }

        /// <summary>
        /// Determines whether the given object is an equal matrix.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if <paramref name="obj"/> is an equal matrix, otherwise false.</returns>
        public override bool Equals (object obj) {
            return (Equals(obj as AnyMatrix<T>));
        }

        /// <summary>
        /// Not a valid operation.
        /// </summary>
        /// <returns>Throws an <see cref="NotSupportedException"/>.</returns>
        /// <remarks>
        /// <para>The <see cref="Object.GetHashCode"/> method is used to provide a quick equality test when an object
        /// is used as a key in a <see cref="System.Collections.Generic.IDictionary{TKey,TValue}"/> or <see cref="System.Collections.IDictionary"/>.
        /// Since a useful hash code of a matrix would need to involve all its elements, it is not possible to make this a fast operation.
        /// Also, since matrices are not immutable, they should not be used as hash keys. (A matrix might be changed after it
        /// had already been used as a key.) For these reasons, requesting a hash code for a matrix is not supported.
        /// </para>
        /// </remarks>
        /// <exception cref="NotSupportedException">This method always throws this exception.</exception>
        public override int GetHashCode () {
            throw new NotSupportedException();
        }

    }

}
