using System;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a collection of complex eigenvalues and eigenvectors.
    /// </summary>
    public sealed class ComplexEigensystem {

        private int dimension;

        private Complex[] eigenvalues;

        private Complex[][] eigenvectors;

        internal ComplexEigensystem (int dimension, Complex[] eigenvalues, Complex[][] eigenvectors) {
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
        /// Gets the specified eigenvalue.
        /// </summary>
        /// <param name="n">The number of the eigenvalue.</param>
        /// <returns>The <paramref name="n"/>th eigenvalue.</returns>
        public Complex Eigenvalue (int n) {
            if ((n < 0) || (n >= dimension)) throw new ArgumentOutOfRangeException("n");
            return (eigenvalues[n]);
        }

        /// <summary>
        /// Gets the specified eigenvector.
        /// </summary>
        /// <param name="n">The number of the eigenvector.</param>
        /// <returns>The <paramref name="n"/>th eigenvector.</returns>
        public Complex[] Eigenvector (int n) {
            if ((n < 0) || (n >= dimension)) throw new ArgumentOutOfRangeException("n");
            Complex[] eigenvector = new Complex[dimension];
            Array.Copy(eigenvectors[n], eigenvector, dimension);
            return (eigenvector);
        }

        /*
        /// <summary>
        /// Gets the matrix transform that diagonalizes the original matrix.
        /// </summary>
        public Complex[,] Transform {
            get {
                return (eigenvectors);
            }
        }
        */

    }

}
