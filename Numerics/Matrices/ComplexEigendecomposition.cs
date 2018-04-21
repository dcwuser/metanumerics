using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a collection of complex eigenvalues and eigenvectors.
    /// </summary>
    public sealed class ComplexEigendecomposition {

        internal ComplexEigendecomposition (int dimension, Complex[] eigenvalues, Complex[][] eigenvectors) {
            Debug.Assert(dimension > 0);
            Debug.Assert(eigenvalues != null);
            Debug.Assert(eigenvalues.Length == dimension);
            Debug.Assert(eigenvectors != null);
            this.dimension = dimension;
            this.eigenvalues = eigenvalues;
            this.eigenvectors = eigenvectors;
        }

        internal readonly int dimension;
        internal readonly Complex[] eigenvalues;
        internal readonly Complex[][] eigenvectors;

        /// <summary>
        /// Gets the dimension of the eigensystem.
        /// </summary>
        public int Dimension {
            get {
                return (dimension);
            }
        }

        /// <summary>
        /// Gets the eigenvalues and eigenvectors associated with the original matrix.
        /// </summary>
        public ComplexEigenpairCollection Eigenpairs
        {
            get
            {
                return (new ComplexEigenpairCollection(this));
            }
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


    /// <summary>
    /// Contains a collection of complex eigenvalues and eigenvectors.
    /// </summary>
    public sealed class ComplexEigenpairCollection : IEnumerable, IEnumerable<ComplexEigenpair>, IReadOnlyCollection<ComplexEigenpair>, IReadOnlyList<ComplexEigenpair>
    {
        internal ComplexEigenpairCollection (ComplexEigendecomposition system) {
            Debug.Assert(system != null);
            this.system = system;
        }


        private readonly ComplexEigendecomposition system;

        /// <summary>
        /// Gets the number of eigenpairs in the collection.
        /// </summary>
        public int Count
        {
            get
            {
                return (system.dimension);
            }
        }

        /// <summary>
        /// Gets the eigenpair with the given index.
        /// </summary>
        /// <param name="index">The index.</param>
        /// <returns>The eignepair with the given index.</returns>
        public ComplexEigenpair this[int index]
        {
            get
            {
                if ((index < 0) || (index >= system.dimension)) throw new ArgumentOutOfRangeException(nameof(index));
                return (new ComplexEigenpair(system.eigenvalues[index], system.eigenvectors[index], system.dimension));
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (((IEnumerable<ComplexEigenpair>)this).GetEnumerator());
        }

        IEnumerator<ComplexEigenpair> IEnumerable<ComplexEigenpair>.GetEnumerator()
        {
            for (int index = 0; index < this.Count; index++)
            {
                yield return (this[index]);
            }
        }

    }

    /// <summary>
    /// Represents a complex-valued eigenvector and corresponding eigenvalue.
    /// </summary>
    public sealed class ComplexEigenpair
    {

        internal ComplexEigenpair(Complex eigenvalue, Complex[] eigenvector, int dimension)
        {
            Debug.Assert(eigenvector != null);
            Debug.Assert(dimension > 0);
            this.eigenvalue = eigenvalue;
            this.eigenvector = eigenvector;
            this.dimension = dimension;
        }

        private readonly Complex eigenvalue;
        private readonly Complex[] eigenvector;
        private readonly int dimension;


        /// <summary>
        /// Gets the eigenvalue.
        /// </summary>
        public Complex Eigenvalue
        {
            get
            {
                return (eigenvalue);
            }
        }

        /// <summary>
        /// Gets the eigenvector.
        /// </summary>
        /// <remarks>
        /// <para>The returned vector is read-only.
        /// For applications that do not require modifying it, this makes producing it much faster and more memory-effecient.
        /// To obtain a modifyable copy, call <see cref="ComplexColumnVector.Copy"/>.</para>
        /// </remarks>
        public ComplexColumnVector Eigenvector
        {
            get
            {
                return (new ComplexColumnVector(eigenvector, 0, dimension, true));
            }
        }

    }

}
