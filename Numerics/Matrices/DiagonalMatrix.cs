using System;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a diagonal matrix.
    /// </summary>
    /// <remarks>
    /// <para>A diagonal matrix has non-zero entries only on the diagonal.</para>
    /// </remarks>
    public class DiagonalMatrix : AnySquareMatrix {

        /// <summary>
        /// Initializes a new, empty diagonal matrix.
        /// </summary>
        /// <param name="dimension">The dimension of the matrix.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="dimension"/> is non-positive.</exception>
        public DiagonalMatrix (int dimension) {
            if (dimension <= 0) throw new ArgumentOutOfRangeException(nameof(dimension));
            diagonal = new double[dimension];
        }

        /// <summary>
        /// Initializes a new diagonal matrix with the given diagonal entries.
        /// </summary>
        /// <param name="diagonal">The diagonal entries.</param>
        public DiagonalMatrix (params double[] diagonal) {
            if (diagonal == null) throw new ArgumentNullException(nameof(diagonal));
            if (diagonal.Length == 0) throw new ArgumentOutOfRangeException(nameof(diagonal));
            this.diagonal = diagonal;
        }

        internal DiagonalMatrix (double[] diagonal, bool isReadOnly) : base(isReadOnly) {
            Debug.Assert(diagonal != null);
            Debug.Assert(diagonal.Length > 0);
            this.diagonal = diagonal;
        }

        private readonly double[] diagonal;

        /// <inheritdoc/>
        public override int Dimension {
            get {
                return (diagonal.Length);
            }
        }

        /// <inheritdoc/>
        public override double this[int r, int c] {
            get {
                CheckIndexBounds(r, c);
                if (r == c) {
                    return (diagonal[r]);
                } else {
                    return (0.0);
                }
            }
            set {
                CheckIndexBounds(r, c);
                if (r == c) {
                    diagonal[r] = value;
                } else {
                    if (value != 0.0) throw new InvalidOperationException();
                }
            }
        }

        private void CheckIndexBounds (int r, int c) {
            if ((r < 0) || (r >= diagonal.Length)) throw new ArgumentOutOfRangeException(nameof(r));
            if ((c < 0) || (c >= diagonal.Length)) throw new ArgumentOutOfRangeException(nameof(c));
        }

    }


}
