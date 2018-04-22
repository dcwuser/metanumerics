using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a unit matrix.
    /// </summary>
    /// <remarks>
    /// <para>If you need to use the unit matrix in your calculations, this class will do the job with a minimal
    /// memory footprint. Use the static <see cref="OfDimension(int)"/> method to get a unit matrix of the
    /// desired dimension. Since this class can only ever represent a unit matrix, it is read-only,
    /// and any attempts to set elements will result in an <see cref="InvalidOperationException"/>. If you
    /// want to start with a unit matrix and then modify it, you can use the <see cref="ToSquareMatrix"/>,
    /// <see cref="ToSymmetricMatrix"/>, and <see cref="ToDiagonalMatrix"/> methods to produce modify-able
    /// matrices of those sorts with initial unit matrix entries.</para>
    /// </remarks>
    public class UnitMatrix : AnySquareMatrix {

        /// <summary>
        /// Returns a unit matrix of the given dimension.
        /// </summary>
        /// <param name="n">The dimension of the matrix.</param>
        /// <returns>A read-only unit matrix of dimension <paramref name="n"/>.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is non-positive.</exception>
        public static UnitMatrix OfDimension (int n) {
            if (n <= 0) throw new ArgumentOutOfRangeException(nameof(n));
            return (new UnitMatrix(n));
        }

        internal UnitMatrix (int n) : base(true) {
            Debug.Assert(n > 0);
            this.n = n;
        }

        private readonly int n;

        private void CheckIndexBounds (int r, int c) {
            if ((r < 0) || (r >= n)) throw new ArgumentOutOfRangeException(nameof(r));
            if ((c < 0) || (c >= n)) throw new ArgumentOutOfRangeException(nameof(c));
        }

        /// <inheritdoc/>
        public override double this[int r, int c] {
            get {
                CheckIndexBounds(r, c);
                if (r == c) {
                    return (1.0);
                } else {
                    return (0.0);
                }
            }
            set {
                throw new InvalidOperationException();
            }
        }

        /// <inheritdoc/>
        public override int Dimension {
            get {
                return (n);
            }
        }

        /// <summary>
        /// Returns a square matrix with unit matrix entries.
        /// </summary>
        /// <returns>An independent, write-able copy of the unit matrix.</returns>
        public SquareMatrix ToSquareMatrix () {
            double[] storage = SquareMatrixAlgorithms.CreateUnitMatrix(n);
            return (new SquareMatrix(storage, n));
        }

        /// <summary>
        /// Returns a symmetric matrix with unit matrix entries.
        /// </summary>
        /// <returns>An independent, write-able copy of the unit matrix.</returns>
        public SymmetricMatrix ToSymmetricMatrix () {
            SymmetricMatrix S = new SymmetricMatrix(n);
            for (int i = 0; i < n; i++) S[i, i] = 1.0;
            return (S);
        }

        /// <summary>
        /// Returns a diagonal matrix with unit matrix entries.
        /// </summary>
        /// <returns>An independent, write-able copy of the unit matrix.</returns>
        public DiagonalMatrix ToDiagonalMatrix () {
            double[] diagonal = new double[n];
            for (int i = 0; i < diagonal.Length; i++) diagonal[i] = 1.0;
            return (new DiagonalMatrix(diagonal));
        }

        /// <inheritdoc/>
        public override double OneNorm () {
            return (1.0);
        }

        /// <inheritdoc/>
        public override double InfinityNorm () {
            return (1.0);
        }

        /// <inheritdoc/>
        public override double FrobeniusNorm () {
            return (Math.Sqrt(n));
        }

        /// <inheritdoc/>
        public override double MaxNorm () {
            return (1.0);
        }

    }
}
