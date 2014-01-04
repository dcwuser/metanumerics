using System;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a collection of real eigenvalues and eigenvectors.
    /// </summary>
    public sealed class RealEigensystem {

        private readonly int dimension;
        private readonly double[] eigenvalues;
        private readonly double[] eigenvectorStorage;

        internal RealEigensystem (int dimension, double[] eigenvalues, double[] eigenvectorStorage) {
            this.dimension = dimension;
            this.eigenvalues = eigenvalues;
            this.eigenvectorStorage = eigenvectorStorage;
        }

        /// <summary>
        /// Gets the dimension of the eigensystem.
        /// </summary>
        public int Dimension {
            get {
                return (dimension);
            }
        }

        /// <summary>
        /// Gets a specified eigenvalue.
        /// </summary>
        /// <param name="n">The (zero-based) index of the eigenvalue.</param>
        /// <returns>The <paramref name="n"/>th eigenvalue.</returns>
        public double Eigenvalue (int n) {
            if ((n < 0) || (n >= dimension)) throw new ArgumentOutOfRangeException("n");
            return (eigenvalues[n]);
        }

        /// <summary>
        /// Gets a specified eigenvector.
        /// </summary>
        /// <param name="n">The (zero-based) index of the eigenvector.</param>
        /// <returns>The <paramref name="n"/>th eigenvector.</returns>
        /// <remarks>
        /// <para>The returned vector is read-only. If you need to make changes to it, you can call <see cref="ColumnVector.Copy"/> to obtain a writable copy.</para>
        /// </remarks>
        public ColumnVector Eigenvector (int n) {
            if ((n < 0) || (n >= dimension)) throw new ArgumentOutOfRangeException("n");
            return (new ColumnVector(eigenvectorStorage, n * dimension, 1, dimension, true));
        }

        /// <summary>
        /// Gets the transformation matrix that diagonalizes the original matrix.
        /// </summary>
        /// <returns>The orthogonal matrix V such that V<sup>T</sup>AV = D, where A is the orignal matrix and D is diagonal.</returns>
        /// <remarks>
        /// <para>The returned matrix is read-only. If you need to make changes to it, you can call <see cref="SquareMatrix.Copy"/> to obtain a writable copy.</para>
        /// </remarks>
        public SquareMatrix TransformMatrix () {
            return (new SquareMatrix(eigenvectorStorage, dimension, true));
        }

    }

}
