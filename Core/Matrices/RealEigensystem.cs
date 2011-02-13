using System;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a collection of real eigenvalues and eigenvectors.
    /// </summary>
    public sealed class RealEigensystem {

        private int dimension;
        private double[] eigenvalues;
        private double[] eigenvectorStorage;

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
        public ColumnVector Eigenvector (int n) {
            if ((n < 0) || (n >= dimension)) throw new ArgumentOutOfRangeException("n");

            double[] eigenvector = new double[dimension];
            Blas1.dCopy(eigenvectorStorage, n * dimension, 1, eigenvector, 0, 1, dimension);
            return (new ColumnVector(eigenvector, dimension));
        }

        /// <summary>
        /// Gets the transformation matrix that diagonalizes the original matrix.
        /// </summary>
        /// <returns>The orthogonal matrix V such that V<sup>T</sup>AV = D, where D is diagonal.</returns>
        public SquareMatrix Eigentransformation () {
            double[] eigenvectorStorageCopy = MatrixAlgorithms.Copy(eigenvectorStorage, dimension, dimension);
            return (new SquareMatrix(eigenvectorStorageCopy, dimension));
        }

    }

}
