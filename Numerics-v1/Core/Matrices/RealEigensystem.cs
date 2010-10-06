using System;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a collection of real eigenvalues and eigenvectors.
    /// </summary>
    public class RealEigensystem {

        private int dimension;
        private double[] eigenvalues;
        private double[,] eigenvectors;

        internal RealEigensystem (int dimension, double[] eigenvalues, double[,] eigenvectors) {
            this.dimension = dimension;
            this.eigenvalues = eigenvalues;
            this.eigenvectors = eigenvectors;
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
            ColumnVector eigenvector = new ColumnVector(dimension);
            for (int r = 0; r < dimension; r++) {
                eigenvector[r] = eigenvectors[r, n];
            }
            return (eigenvector);
        }

    }

}
