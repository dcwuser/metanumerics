using Meta.Numerics.Functions;
using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// A permutation matrix.
    /// </summary>
    /// <remarks>
    /// <para>A d-dimensional permutation matrix is a d-dimensional representation of a permutation among d objects. When it right-multiplies a column vector, it
    /// permutes the entries of the vector. When it right-multiples a matrix, it permutes the columns of the matrix. When two permutation matrices are
    /// multiplied, the product is the permutation matrix corresponding to the right-to-left sequential application of the permutations they represent.</para>
    /// <para>The entries of a permutation matrix are all either 0 or 1, and there is exactly one 1 in every row and column.</para>
    /// </remarks>
    public class PermutationMatrix : AnySquareMatrix {

        internal PermutationMatrix(PermutationAsMap permutation) : base(true) {
            this.permutation = permutation;
        }

        private PermutationAsMap permutation;

        /// <inheritdoc/>
        public override int Dimension {
            get {
                return permutation.map.Length;  
            }
        }

        /// <inheritdoc/>
        public override double this[int r, int c] {
            get {
                if ((r < 0) || (r >= permutation.map.Length)) throw new ArgumentOutOfRangeException(nameof(r));
                if ((c < 0) || (c >= permutation.map.Length)) throw new ArgumentOutOfRangeException(nameof(c));
                if (permutation.map[c] == r) {
                    return 1.0;
                } else {
                    return 0.0;
                }
            }
            set {
                throw new InvalidOperationException();
            }
        }

        // Add inverse and transpose

        /// <summary>
        /// Get the inverse of the permutation matrix.
        /// </summary>
        /// <returns>The inverse of the permutation matrix, which is a permutation matrix representing the inverse of the permutation is represents.</returns>
        public PermutationMatrix Inverse() {
            return new PermutationMatrix(permutation.Inverse());
        }

        /// <summary>
        /// The transpose of the permutation matrix.
        /// </summary>
        /// <remarks>
        /// <para>The transpose of a permutation matrix is the same as its inverse.</para>
        /// </remarks>
        public PermutationMatrix Transpose {
            get {
                return new PermutationMatrix(permutation.Inverse());
            }
        }


    }
}
