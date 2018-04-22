using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a collection of real eigenvalues and eigenvectors.
    /// </summary>
    public sealed class RealEigendecomposition {

        internal readonly int dimension;
        internal double[] eigenvalues;
        internal double[] eigenvectors;
        private RealEigenpairCollection pairs = null;

        internal RealEigendecomposition (double[] eigenvalues, double[] eigenvectors, int dimension) {
            Debug.Assert(eigenvalues != null);
            Debug.Assert(eigenvectors != null);
            Debug.Assert(dimension > 0);
            this.eigenvalues = eigenvalues;
            this.eigenvectors = eigenvectors;
            this.dimension = dimension;
        }

        /// <summary>
        /// Gets the dimension of the original matrix.
        /// </summary>
        public int Dimension {
            get {
                return (dimension);
            }
        }


        /// <summary>
        /// Gets a list of the eigenvalues and eigenvectors of the original matrix.
        /// </summary>
        public RealEigenpairCollection Eigenpairs
        {
            get
            {
                if (pairs == null) pairs = new RealEigenpairCollection(this);
                return (pairs);
            }
        }

        /// <summary>
        /// Gets the transformation matrix that diagonalizes the original matrix.
        /// </summary>
        /// <value>The read-only, orthogonal matrix V, such that V<sup>T</sup>AV = D, where A is the original matrix and D is the <see cref="DiagonalizedMatrix"/>.</value>
        /// <remarks>
        /// <para>The returned matrix is read-only. If you need to modify it, call <see cref="SquareMatrix.Copy"/> to obtain a writable copy.</para>
        /// </remarks>
        public SquareMatrix TransformMatrix {
            get {
                return (new SquareMatrix(eigenvectors, dimension, true));
            }
        }

        /// <summary>
        /// Gets the diagonalized matrix.
        /// </summary>
        /// <value>A read-only, diagonal matrix D, such that V<sup>T</sup>AV = D, where A is the original matrix and V is the orthogonal <see cref="TransformMatrix"/>.</value>
        public DiagonalMatrix DiagonalizedMatrix {
            get {
                return (new DiagonalMatrix(eigenvalues, true));
            }
        }

        /// <summary>
        /// Computes the determinant of the original matrix.
        /// </summary>
        /// <returns>The determinant of the original matrix.</returns>
        /// <remarks>
        /// <para>The determinant of the original matrix equals the product of the eigenvalues.</para>
        /// </remarks>
        public double Determinant ()
        {
            double det = 1.0;
            foreach (double eigenvalue in eigenvalues) det *= eigenvalue;
            return (det);
        }

    }


    /// <summary>
    /// Contains real-valued eigenvalues and eigenvectors of a matrix.
    /// </summary>
    public sealed class RealEigenpairCollection : IEnumerable, IEnumerable<RealEigenpair>, IReadOnlyCollection<RealEigenpair>, IReadOnlyList<RealEigenpair>
    {

        internal RealEigenpairCollection (RealEigendecomposition system)
        {
            Debug.Assert(system != null);
            this.system = system;
            this.map = GetIntegers(system.dimension);
        }

        private readonly RealEigendecomposition system;
        private readonly int[] map;

        private static int[] GetIntegers (int n)
        {
            int[] integers = new int[n];
            for (int i = 0; i < n; i++) integers[i] = i;
            return (integers);
        }

        /// <summary>
        /// Gets the number of eigenpairs in the collection.
        /// </summary>
        public int Count
        {
            get
            {
                return (map.Length);
            }
        }

        /// <summary>
        /// Gets the eigenpair with the given index.
        /// </summary>
        /// <param name="index">The (zero-based) index.</param>
        /// <returns>The eigenpair with the given index.</returns>
        public RealEigenpair this[int index]
        {
            get
            {
                if ((index < 0) || (index >= map.Length)) throw new ArgumentOutOfRangeException(nameof(index));
                int k = map[index];
                return (new RealEigenpair(system.eigenvalues[k], system.eigenvectors, k * system.dimension, system.dimension));
            }
        }

        IEnumerator<RealEigenpair> IEnumerable<RealEigenpair>.GetEnumerator()
        {
            for (int index = 0; index < map.Length; index++)
            {
                yield return (this[index]);
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (((IEnumerable<RealEigenpair>)this).GetEnumerator());
        }


        /// <summary>
        /// Sort the eigenpairs by eigenvalue as specified.
        /// </summary>
        /// <param name="order">The desired ordering.</param>
        /// <remarks>
        /// <para>Before this method is called, the ordering of the eigenpairs is random and should not be relied upon.</para>
        /// </remarks>
        public void Sort(OrderBy order)
        {
            // Create a comparison to affect the sort
            Comparison<int> comparison;
            switch (order)
            {
                case OrderBy.ValueAscending:
                    comparison = (int a, int b) => system.eigenvalues[a].CompareTo(system.eigenvalues[b]);
                    break;
                case OrderBy.ValueDescending:
                    comparison = (int a, int b) => system.eigenvalues[b].CompareTo(system.eigenvalues[a]);
                    break;
                case OrderBy.MagnitudeAscending:
                    comparison = (int a, int b) => Math.Abs(system.eigenvalues[a]).CompareTo(Math.Abs(system.eigenvalues[b]));
                    break;
                case OrderBy.MagnitudeDescending:
                    comparison = (int a, int b) => Math.Abs(system.eigenvalues[b]).CompareTo(Math.Abs(system.eigenvalues[a]));
                    break;
                default:
                    throw new NotSupportedException();
            }

            // Sort the map
            Array.Sort(map, comparison);
        }

    }

    /// <summary>
    /// Contains a real-valued eigenvector and eigenvalue.
    /// </summary>
    public sealed class RealEigenpair
    {
        internal RealEigenpair(double value, double[] store, int offset, int dimension)
        {
            Debug.Assert(store != null);
            Debug.Assert(offset >= 0);
            Debug.Assert(dimension > 0);
            this.value = value;
            this.store = store;
            this.offset = offset;
            this.dimension = dimension;
        }


        private readonly double value;
        private readonly double[] store;
        private readonly int offset;
        private readonly int dimension;

        /// <summary>
        /// Gets the eigenvalue.
        /// </summary>
        public double Eigenvalue {
            get
            {
                return (value);
            }
        }

        /// <summary>
        /// Gets the eigenvector.
        /// </summary>
        /// <remarks>
        /// <para>The returned vector is read-only. If you need to make changes to it, you can call <see cref="ColumnVector.Copy"/> to obtain a writable copy.</para>
        /// </remarks>
        public ColumnVector Eigenvector {
            get {
                return (new ColumnVector(store, offset, 1, dimension, true));
            }
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
