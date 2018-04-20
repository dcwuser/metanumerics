using System;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a collection of real eigenvalues and eigenvectors.
    /// </summary>
    public sealed class RealEigensystem {

        private readonly int dimension;
        private double[] eigenvalues;
        private double[] eigenvectorStorage;

        internal RealEigensystem (int dimension, double[] eigenvalues, double[] eigenvectorStorage) {
            Debug.Assert(dimension > 0);
            Debug.Assert(eigenvalues != null);
            Debug.Assert(eigenvectorStorage != null);
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
            if ((n < 0) || (n >= dimension)) throw new ArgumentOutOfRangeException(nameof(n));
            return (new ColumnVector(eigenvectorStorage, n * dimension, 1, dimension, true));
        }

        /// <summary>
        /// Gets the transformation matrix that diagonalizes the original matrix.
        /// </summary>
        /// <value>The orthogonal matrix V such that V<sup>T</sup>AV = D, where A is the original matrix and D is diagonal.</value>
        /// <remarks>
        /// <para>The returned matrix is read-only. If you need to make changes to it, you can call <see cref="SquareMatrix.Copy"/> to obtain a writable copy.</para>
        /// </remarks>
        public SquareMatrix TransformMatrix {
            get {
                return (new SquareMatrix(eigenvectorStorage, dimension, true));
            }
        }

        /// <summary>
        /// Sort the eigenvalues as specified.
        /// </summary>
        /// <param name="order">The desired ordering.</param>
        public void Sort (OrderBy order) {

            // Create an auxiliary array of indexes to sort
            int[] ranks = new int[dimension];
            for (int i = 0; i < ranks.Length; i++) ranks[i] = i;

            // Sort indexes in the requested order
            Comparison<int> comparison;
            switch (order) {
                case OrderBy.ValueAscending:
                    comparison = (int a, int b) => eigenvalues[a].CompareTo(eigenvalues[b]);
                    break;
                case OrderBy.ValueDescending:
                    comparison = (int a, int b) => eigenvalues[b].CompareTo(eigenvalues[a]);
                    break;
                case OrderBy.MagnitudeAscending:
                    comparison = (int a, int b) => Math.Abs(eigenvalues[a]).CompareTo(Math.Abs(eigenvalues[b]));
                    break;
                case OrderBy.MagnitudeDescending:
                    comparison = (int a, int b) => Math.Abs(eigenvalues[b]).CompareTo(Math.Abs(eigenvalues[a]));
                    break;
                default:
                    throw new NotSupportedException();
            }
            Array.Sort(ranks, comparison);

            // Create new storage in the desired order and discard the old storage
            // This is faster than moving within existing storage since we can move each column just once.
            double[] newEigenvalues = new double[dimension];
            double[] newEigenvectorStorage = new double[dimension * dimension];
            for (int i = 0; i < ranks.Length; i++) {
                newEigenvalues[i] = eigenvalues[ranks[i]];
                Blas1.dCopy(eigenvectorStorage, dimension * ranks[i], 1, newEigenvectorStorage, dimension * i, 1, dimension);
            }
            eigenvalues = newEigenvalues;
            eigenvectorStorage = newEigenvectorStorage;

        }

    }

    /// <summary>
    /// Describes an ordering of eigenvalues.
    /// </summary>
    public enum OrderBy {
        
        /// <summary>
        /// From most negative to most positive.
        /// </summary>
        ValueAscending,
        
        /// <summary>
        /// From most positive to most negative.
        /// </summary>
        ValueDescending,
        
        /// <summary>
        /// From smallest absolute value to largest absolute value.
        /// </summary>
        MagnitudeAscending,
        
        /// <summary>
        /// From largest absolute value to smallest absolute value.
        /// </summary>
        MagnitudeDescending
    
    }

}
