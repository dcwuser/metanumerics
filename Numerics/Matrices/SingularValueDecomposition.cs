using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Stores the singular value decomposition of a matrix.
    /// </summary>
    /// <remarks>
    /// <para>The singular value decomposition of a matrix represents it as a product of a left orthogonal matrix U, a quasi-diagonal
    /// &#x3A3; matrix, and a right orthogonal matrix V.</para>
    /// <img src="../images/SVDEquation.png" />
    /// <para>Any rectangular matrix has an SVD decomposition. The matrix need not be square. If square, it need not be invertible.
    /// The dimensions of the decomposition factors are illustrated in the following diagram.</para>
    /// <img src="../images/SVDForm.png" />
    /// <para>The elements of the &#x3A3; are called the singular values of the original matrix.</para>
    /// <para>Viewing A as a linear operator, the rows of V<sup>T</sup> (columns of V) form an
    /// orthonormal basis for the domain of the operator, while the columns of U form an orthonormal
    /// basis for the range of the operator. These rows and columns are called, respectively, the right and left singular vectors
    /// of the matrix.</para>
    /// <para>The right singular vectors corresponding to zero singular values span the null-space of A, that is the
    /// set of all x for which Ax = 0. The left singular vectors corresponding to non-zero singular values span the range of A,
    /// that is the space into which all Ax fall; the left singular vectors corresponding to zero singular values span the
    /// complementary space into which no Ax fall. The total number of non-zero singular values is the rank of A.</para>
    /// <para>The SVD can be used to approximate the action of a high-dimensional matrix operator by a lower-rank one. By
    /// keeping only the largest singular values and setting the remaining ones to zero, one obtains a operator that
    /// can applied with fewer operations and approximates the properties of the original operator.</para>
    /// <para>Notice that elements in the last columns of U do not contribute to A since they will be multiplied with elements of zero rows of &#x3A3;.
    /// Many applications use a "thin" or "reduced" form of the SVD, in which the last columns of U and the last rows of &#x3A3; are omitted;
    /// this makes U not square (and therefore not orthogonal) and &#x3A3; square (and diagonal). Since the remaining elements are the same, you
    /// can obtain the thin SVD from this object by simply ignoring the irrelevant elements.</para>
    /// <para>Use the <see cref="RectangularMatrix.SingularValueDecomposition"/> of the <see cref="RectangularMatrix"/> class
    /// to obtain the SVD of an rectangular matrix, or the corresponding <see cref="SquareMatrix.SingularValueDecomposition"/>
    /// method of the <see cref="SquareMatrix"/> class to obtain the SVD of a square matrix.</para>
    /// </remarks>
    /// <see href="http://en.wikipedia.org/wiki/Singular_value_decomposition"/>
    public sealed class SingularValueDecomposition {

        internal SingularValueDecomposition (double[] utStore, double[] wStore, double[] vStore, int rows, int cols) {
            Debug.Assert(utStore != null);
            Debug.Assert(wStore != null);
            Debug.Assert(vStore != null);
            Debug.Assert(rows > 0);
            Debug.Assert(cols > 0);
            this.utStore = utStore;
            this.vStore = vStore;
            this.wStore = wStore;
            this.rows = rows;
            this.cols = cols;
        }

        internal readonly int rows, cols;
        internal readonly double[] utStore;
        internal readonly double[] wStore;
        internal readonly double[] vStore;

        /// <summary>
        /// Gets the number of rows in the original matrix.
        /// </summary>
        public int RowCount {
            get {
                return (rows);
            }
        }

        /// <summary>
        /// Gets the number of columns in the original matrix.
        /// </summary>
        public int ColumnCount {
            get {
                return (cols);
            }
        }

        /// <summary>
        /// Gets the left transform matrix.
        /// </summary>
        /// <value>The matrix U, such that A = U S V<sup>T</sup>.</value>
        /// <remarks>
        /// <para>The returned matrix is read-only. If you need to make changes to it, you can call <see cref="SquareMatrix.Copy"/> to obtain a writable copy.</para>
        /// </remarks>
        public SquareMatrix LeftTransformMatrix {
            get {
                return (new SquareMatrix(utStore, 0, rows, 1, rows, true));
            }
        }

        /// <summary>
        /// Gets the right transform matrix.
        /// </summary>
        /// <value>The matrix V, such that A = U S V<sup>T</sup>.</value>
        /// <remarks>
        /// <para>The returned matrix is read-only. If you need to make changes to it, you can call <see cref="SquareMatrix.Copy"/> to obtain a writable copy.</para>
        /// </remarks>
        public SquareMatrix RightTransformMatrix {
            get {
                return (new SquareMatrix(vStore, cols, true));
            }
        }

        /// <summary>
        /// Gets the number of singular values.
        /// </summary>
        /// <remarks>
        /// <para>For a square matrix, the number of singular values is equal to the dimension of the matrix.
        /// For a rectangular matrix with more rows than columns, the number of singular values is equal to
        /// the number of columns.</para>
        /// </remarks>
        public int Dimension {
            get {
                return (wStore.Length);
            }
        }

        /// <summary>
        /// Gets a collection of all the terms contributing to the singular value decomposition.
        /// </summary>
        public SingularValueContributorCollection Contributors {
            get {
                return (new SingularValueContributorCollection(this));
            }
        }

        /// <summary>
        /// Gets the condition number of the matrix.
        /// </summary>
        /// <remarks>
        /// <para>The condition number is the ratio of the largest singular value to smallest singular value. It is therefore always larger than one.</para>
        /// </remarks>
        public double ConditionNumber {
            get {
                return (wStore[0] / wStore[wStore.Length - 1]);
            }
        }

        private double tolerance = Global.Accuracy;

        /// <summary>
        /// Gets or sets the tolerance with which singular values are distinguished from zero.
        /// </summary>
        /// <remarks>
        /// <para>Some operations offered by singular value decompositions, including rank determination and computation of the pseudo-inverse matrix,
        /// depend on determining whether or not a singular value is zero. Since floating point numbers are only approximate representations of real
        /// numbers, singular values will usually not be exactly zero even for matrices for which they should be, but will instead be very tiny
        /// numbers, on the order of floating point precision (about 10<sup>-16</sup>) as a fraction of the largest singular value. The value
        /// of this property is used to determine how small a singular value must be, as a fraction of the largest singular value, to be considered
        /// zero for these purposes. Usually you will want to maintain its default value.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">The assigned value is outside the range [0,1).</exception>
        public double Tolerance {
            get {
                return (tolerance);
            }
            set {
                if ((value < 0.0) || (value >= 1.0)) throw new ArgumentOutOfRangeException(nameof(value));
                tolerance = value;
            }
        }

        /// <summary>
        /// Computes the rank of the original matrix.
        /// </summary>
        /// <remarks>
        /// <para>The rank of a matrix is the dimension of the space of input vectors which produce non-zero output vectors upon multiplication by the original matrix.</para>
        /// <para>Since this operation depends on identifying zero singular values, the result will depend on the value of the <see cref="Tolerance"/> property.</para>
        /// </remarks>
        public int Rank {
            get {
                double minValue = tolerance * wStore[0];
                int rank = 0;
                for (int i = 0; i < wStore.Length; i++) {
                    if (wStore[i] > minValue) {
                        rank++;
                    } else {
                        break;
                    }
                }
                return (rank);
            }
        }

        /// <summary>
        /// Solves the system of equations Ax = b, where A is the original, square matrix.
        /// </summary>
        /// <param name="rhs">The right-hand-side vector b.</param>
        /// <returns>The solution vector x.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="rhs"/> is null.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="rhs"/> does not have the same dimension as the original matrix.</exception>
        /// <exception cref="InvalidOperationException">The original matrix is not square.</exception>
        /// <remarks>
        /// <para>For singular and nearly-singular matrices, this method operates differently than other solution methods like
        /// <see cref="SquareQRDecomposition.Solve(IReadOnlyList{double})"/> and <see cref="LUDecomposition.Solve(IReadOnlyList{double})"/>.
        /// For singular and nearly-singular matrices, those methods
        /// tend to produce very large or even infinite solution components that delicately cancel to produce the required right-hand-side,
        /// but have little significance to the problem being modeled because they arise from inverting very small singular values
        /// of the original matrix that are dominated by floating point rounding errors. The SVD solution discards those singular
        /// values to obtain a solution vector driven by the dominant, non-singular parts of A. While this solution does not have
        /// the minimum achievable |Ax - b|, it is often more representative of the desired solution to the problem being modeled.
        /// </para>
        /// <para>For original matrices that are not singular or nearly-singular, this method will compute the same solution
        /// as other methods.</para>
        /// <para>This method is only available for square original matrices.</para>
        /// </remarks>
        public ColumnVector Solve (IReadOnlyList<double> rhs) {

            if (rows != cols) throw new InvalidOperationException();
            if (rhs == null) throw new ArgumentNullException(nameof(rhs));
            if (rhs.Count != cols) throw new DimensionMismatchException();

            double[] x = new double[rows];
            double[] y = new double[rows];

            rhs.CopyTo(x, 0);

            Blas2.dGemv(utStore, 0, 1, rows, x, 0, 1, y, 0, 1, rows, rows);

            double minValue = tolerance * wStore[0];
            for (int i = 0; i < wStore.Length; i++) {
                double w = wStore[i];
                if (w > minValue) {
                    y[i] /= w;
                } else {
                    y[i] = 0.0;
                }
            }

            Array.Clear(x, 0, rows);
            Blas2.dGemv(vStore, 0, 1, rows, y, 0, 1, x, 0, 1, rows, rows);

            return (new ColumnVector(x, rows));
        }
    }

    /// <summary>
    /// Contains all the terms contributing to the singular value decomposition of a matrix.
    /// </summary>
    public sealed class SingularValueContributorCollection : IReadOnlyList<SingularValueContributor> {

        internal SingularValueContributorCollection(SingularValueDecomposition system) {
            Debug.Assert(system != null);
            this.system = system;
        }

        private readonly SingularValueDecomposition system;

        /// <summary>
        /// Gets the singular value contributor with the given index.
        /// </summary>
        /// <param name="index">The (zero-based) index.</param>
        /// <returns>The singular value contributor with the given index.</returns>
        /// <remarks>
        /// <para>More significant singular values have lower indexes. The most significant singular
        /// value with have index 0.</para>
        /// </remarks>
        public SingularValueContributor this[int index] {
            get {
                return (new SingularValueContributor(system, index));
            }
        }

        /// <summary>
        /// Gets the number of singular value contributors.
        /// </summary>
        public int Count {
            get {
                return (system.wStore.Length);
            }
        }

        IEnumerator<SingularValueContributor> IEnumerable<SingularValueContributor>.GetEnumerator () {
            for (int i = 0; i < this.Count; i++) {
                yield return (this[i]);
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            throw new NotImplementedException();
        }
    }

    /// <summary>
    /// Contains the value and left and right vectors of one term in a singular value decomposition.
    /// </summary>
    public sealed class SingularValueContributor {

        internal SingularValueContributor (SingularValueDecomposition system, int n) {
            Debug.Assert(system != null);
            Debug.Assert(n >= 0 && n < system.wStore.Length);
            this.system = system;
            this.n = n;
        }

        private readonly SingularValueDecomposition system;
        private readonly int n;

        /// <summary>
        /// Gets the singular value.
        /// </summary>
        public double SingularValue {
            get {
                return (system.wStore[n]);
            }
        }

        /// <summary>
        /// Gets the left singular vector.
        /// </summary>
        /// <remarks>
        /// <para>The returned vector is read-only. If you need to make changes to it, you can call <see cref="ColumnVector.Copy"/> to obtain a writable copy.</para>
        /// </remarks>
        public ColumnVector LeftSingularVector {
            get {
                return (new ColumnVector(system.utStore, n, system.rows, system.rows, true));
            }
        }


        /// <summary>
        /// Gets right singular vector.
        /// </summary>
        /// <remarks>
        /// <para>The returned vector is read-only. If you need to make changes to it, you can call <see cref="ColumnVector.Copy"/> to obtain a writable copy.</para>
        /// </remarks> 
        public ColumnVector RightSingularVector {
            get {
                return (new ColumnVector(system.vStore, n * system.cols, 1, system.cols, true));
            }
        }

    }

}
