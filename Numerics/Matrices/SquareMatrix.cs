using System;
using System.Collections.Generic;
using System.Diagnostics;


namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a square matrix.
    /// </summary>
    [DebuggerDisplay("{RowCount} X {ColumnCount} SquareMatrix")]
    [DebuggerTypeProxy(typeof(AnyRectangularMatrixDebuggerTypeProxy))]
    public sealed class SquareMatrix : AnySquareMatrix {

        private readonly int dimension;
        private readonly double[] store;
        private readonly int offset;
        private readonly int rowStride;
        private readonly int colStride;

        /// <summary>
        /// Initializes a new square matrix.
        /// </summary>
        /// <param name="dimension">The dimension of the matrix, which must be positive.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="dimension"/> &lt; 1.</exception>
        public SquareMatrix (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException(nameof(dimension));
            this.dimension = dimension;
            this.store = MatrixAlgorithms.AllocateStorage(dimension, dimension, ref offset, ref rowStride, ref colStride);
        }

        /// <summary>
        /// Initializes a new square matrix from the given 2D array.
        /// </summary>
        /// <param name="entries">The source 2D array.</param>
        /// <remarks><para>The entries are copied, so changes to <paramref name="entries"/> after the matrix is initialized do not affect matrix entries.</para></remarks>
        /// <exception cref="ArgumentNullException"><paramref name="entries"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The two dimensions of <paramref name="entries"/> are not equal.</exception>
        public SquareMatrix (double[,] entries) {
            if (entries == null) throw new ArgumentNullException(nameof(entries));
            if (entries.GetLength(0) != entries.GetLength(1)) throw new DimensionMismatchException();
            this.dimension = entries.GetLength(0);
            this.store = MatrixAlgorithms.AllocateStorage(dimension, dimension, ref offset, ref rowStride, ref colStride);
            for (int r = 0; r < dimension; r++) {
                for (int c = 0; c < dimension; c++) {
                    this[r, c] = entries[r, c];
                }
            }
        }

        internal SquareMatrix (double[] storage, int offset, int rowStride, int colStride, int dimension, bool isReadOnly) : base(isReadOnly) {
            Debug.Assert(storage != null);
            Debug.Assert(dimension > 0);
            this.store = storage;
            this.dimension = dimension;
            this.offset = offset;
            this.rowStride = rowStride;
            this.colStride = colStride;
        }

        internal SquareMatrix (double[] storage, int dimension, bool isReadOnly) : this(storage, 0, 1, dimension, dimension, isReadOnly) {  }

        internal SquareMatrix (double[] storage, int dimension) : this(storage, dimension, false) { }

        // required methods

        /// <summary>
        /// Gets the dimension of the matrix.
        /// </summary>
        public override int Dimension {
            get {
                return (dimension);
            }
        }

        /// <summary>
        /// Gets or sets an entry of the matrix.
        /// </summary>
        /// <param name="r">The (zero-based) row number.</param>
        /// <param name="c">The (zero-based) column number.</param>
        /// <returns>The value of the specified matrix entry M<sub>r c</sub>.</returns>
        public override double this[int r, int c] {
            get {
                CheckBounds(r, c);
                return (store[MatrixAlgorithms.GetIndex(offset, rowStride, colStride, r, c)]);
            }
            set {
                CheckBounds(r, c);
                if (IsReadOnly) throw new InvalidOperationException();
                store[MatrixAlgorithms.GetIndex(offset, rowStride, colStride, r, c)] = value;
            }
        }

        private void CheckBounds (int r, int c) {
            if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
            if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException(nameof(c));
        }

        /// <inheritdoc />
        public override double OneNorm () {
            return (MatrixAlgorithms.OneNorm(store, offset, rowStride, colStride, dimension, dimension));
        }

        /// <inheritdoc />
        public override double InfinityNorm () {
            return (MatrixAlgorithms.InfinityNorm(store, offset, rowStride, colStride, dimension, dimension));
        }

        /// <summary>
        /// Returns a vector representing a given row of the matrix.
        /// </summary>
        /// <param name="r">The (zero-based) row index to return.</param>
        /// <returns>A vector containing the entries of the specified row.</returns>
        /// <remarks>
        /// <para>The returned vector is independent of the matrix.
        /// If an entry of the returned vector is updated, the corresponding entry of the original matrix will not be updated
        /// as well, Similarly, if an entry in the matrix is updated after this method is called, the corresponding
        /// entry of the vector will not change.</para>
        /// </remarks>
        public override RowVector Row (int r) {
            if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException(nameof(r));
            double[] rStore = new double[dimension];
            Blas1.dCopy(store, offset + rowStride * r, colStride, rStore, 0, 1, dimension);
            return (new RowVector(rStore, dimension));
        }

        /// <summary>
        /// Gets a copy of one column of the the matrix.
        /// </summary>
        /// <param name="c">The (zero-based) column number to return.</param>
        /// <returns>An independent copy of the specified column.</returns>
        /// <remarks>
        /// <para>The returned vector is independent of the matrix.
        /// If an entry of the returned vector is updated, the corresponding entry of the original matrix will not be updated
        /// as well, Similarly, if an entry in the matrix is updated after this method is called, the corresponding
        /// entry of the vector will not change.</para>
        /// </remarks>
        public override ColumnVector Column (int c) {
            if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException(nameof(c));
            double[] cStore = new double[dimension];
            Blas1.dCopy(store, offset + colStride * c, rowStride, cStore, 0, 1, dimension);
            return (new ColumnVector(cStore, dimension));
        }

        /// <summary>
        /// Copies the matrix.
        /// </summary>
        /// <returns>An independent copy of the matrix.</returns>
        public SquareMatrix Copy () {
            double[] copyStore = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, dimension, dimension);
            return (new SquareMatrix(copyStore, dimension));
        }

        /// <summary>
        /// Gets the transpose of the matrix.
        /// </summary>
        /// <value>The matrix transpose M<sup>T</sup>.</value>
        /// <remarks>
        /// <para>The returned transpose matrix is not independent of the original matrix.
        /// Instead, it is a read-only view into the same storage as the original matrix with row an column indices reversed.
        /// This has the advantage that is can be produced with almost no time and memory cost.
        /// It has the disadvantage that any subsequent changes to the original matrix will be reflected in the returned transpose,
        /// which might not be what you expected. If you want an independent, write-able transpose matrix, call <see cref="Copy"/>
        /// on the returned matrix.
        /// </para>
        /// </remarks>
        public SquareMatrix Transpose {
            get {
                return (new SquareMatrix(store, offset, colStride, rowStride, dimension, true));
            }
        }

        /// <summary>
        /// Computes the inverse of the matrix.
        /// </summary>
        /// <returns>The matrix inverse M<sup>-1</sup>.</returns>
        /// <remarks>
        /// <para>The inverse of a matrix M is a matrix M<sup>-1</sup> such that M<sup>-1</sup>M = I, where I is the identity matrix.</para>
        /// <para>If the matrix is singular, inversion is not possible. In that case, this method will fail with a <see cref="DivideByZeroException"/>.</para>
        /// <para>The inversion of a matrix is an O(N<sup>3</sup>) operation.</para>
        /// </remarks>
        /// <exception cref="DivideByZeroException">The matrix is singular.</exception>
        public SquareMatrix Inverse () {
            double[] iStore = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, dimension, dimension);
            SquareMatrixAlgorithms.GaussJordanInvert(iStore, dimension);
            return(new SquareMatrix(iStore, dimension));
        }

        /*
        /// <summary>
        /// Inverts the matrix, in place.
        /// </summary>
        /// <remarks>
        /// <para>This method replaces the elements of M with the elements of M<sup>-1</sup>.</para>
        /// <para>In place matrix inversion saves memory, since separate storage of M and M<sup>-1</sup> is not required.</para></remarks>
        private void InvertInPlace () {
        }
        */
        
        /// <summary>
        /// Computes the LU decomposition of the matrix.
        /// </summary>
        /// <returns>The LU decomposition of the matrix.</returns>
        /// <remarks>
        /// <para>The LU decomposition of a matrix M is a set of matrices L, U, and P such that LU = PM, where L
        /// is lower-left triangular, U is upper-right triangular, and P is a permutation matrix (so that PM is
        /// a row-wise permutation of M).</para>
        /// <para>The LU decomposition of a square matrix is an O(N<sup>3</sup>) operation.</para>
        /// </remarks>
        public LUDecomposition LUDecomposition () {

            // Copy the matrix content
            double[] luStore = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, dimension, dimension);

            // Prepare an initial permutation and parity
            int[] permutation = new int[dimension];
            for (int i = 0; i < permutation.Length; i++) {
                permutation[i] = i;
            }
            int parity = 1;

            // Do the LU decomposition
            SquareMatrixAlgorithms.LUDecompose(luStore, permutation, ref parity, dimension);

            // Package it up and return it
            LUDecomposition LU = new LUDecomposition(luStore, permutation, parity, dimension);
            return (LU);

        }

        /// <summary>
        /// Computes the eigenvalues of the matrix.
        /// </summary>
        /// <returns>The eigenvalues of the matrix.</returns>
        /// <seealso cref="Eigendecomposition"/>
        public Complex[] Eigenvalues () {
            double[] aStore = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, dimension, dimension);
            SquareMatrixAlgorithms.IsolateCheapEigenvalues(aStore, null, dimension);
            SquareMatrixAlgorithms.ReduceToHessenberg(aStore, null, dimension);
            Complex[] eigenvalues = SquareMatrixAlgorithms.ReduceToRealSchurForm(aStore, null, dimension);
            return (eigenvalues);
            
        }

        /// <summary>
        /// Computes the eigenvalues and eigenvectors of the matrix.
        /// </summary>
        /// <returns>A decomposition that makes the eigenvalues and eigenvectors manifest.</returns>
        /// <remarks>
        /// <para>For a generic vector v and matrix M, Mv = u will point in some direction with no particular relationship to v.
        /// The eigenvectors of a matrix M are vectors z that satisfy Mz = &#x3BB;z, i.e. multiplying an eigenvector by a
        /// matrix reproduces the same vector, up to a proportionality constant &#x3BB; called the eigenvalue.</para>
        /// <para>For v to be an eigenvector of M with eigenvalue &#x3BB;, (M - &#x3BB;I)z = 0. But for a matrix to
        /// annihilate any non-zero vector, that matrix's determinant must vanish, so det(M - &#x3BB;I)=0. For a matrix of
        /// order N, this is an equation for the roots of a polynomial of order N. Since an order-N polynomial always has exactly
        /// N roots, an order-N matrix always has exactly N eigenvalues.</para>
        /// <para>Since a polynomial with real coefficients can still have complex roots, a real square matrix can nonetheless
        /// have complex eigenvalues (and corresponding complex eigenvectors). However, again like the complex roots of a real
        /// polynomial, such eigenvalues will always occurs in complex-conjugate pairs.</para>
        /// <para>Although the eigenvalue polynomial ensures that an order-N matrix has N eigenvalues, it can occur that there
        /// are not N corresponding independent eigenvectors. A matrix with fewer eigenvectors than eigenvalues is called
        /// defective. Like singularity, defectiveness represents a delicate balance between the elements of a matrix that can
        /// typically be disturbed by just an infinitesimal perturbation of elements. Because of round-off-error, then, floating-point
        /// algorithms cannot reliably identify defective matrices. Instead, this method will return a full set of eigenvectors,
        /// but some eigenvectors, corresponding to very nearly equal eigenvalues, will be very nearly parallel.</para>
        /// <para>While a generic square matrix can be defective, many subspecies of square matrices are guaranteed not to be.
        /// This includes Markov matrices, orthogonal matrices, and symmetric matrices.</para>
        /// <para>Determining the eigenvalues and eigenvectors of a matrix is an O(N<sup>3</sup>) operation. If you need only the
        /// eigenvalues of a matrix, the <see cref="Eigenvalues"/> method is more efficient.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix"/>
        public ComplexEigendecomposition Eigendecomposition () {

            //double[] aStore = MatrixAlgorithms.Copy(store, dimension, dimension);
            double[] aStore = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, dimension, dimension);

            int[] perm = new int[dimension];
            for (int i = 0; i < perm.Length; i++) perm[i] = i;
            SquareMatrixAlgorithms.IsolateCheapEigenvalues(aStore, perm, dimension);
            double[] qStore = MatrixAlgorithms.AllocateStorage(dimension, dimension);
            for (int i = 0; i < perm.Length; i++) MatrixAlgorithms.SetEntry(qStore, dimension, dimension, perm[i], i, 1.0);
            //double[] qStore = SquareMatrixAlgorithms.CreateUnitMatrix(dimension);

            // Reduce the original matrix to Hessenberg form
            SquareMatrixAlgorithms.ReduceToHessenberg(aStore, qStore, dimension);

            // Reduce the Hessenberg matrix to real Schur form
            SquareMatrixAlgorithms.ReduceToRealSchurForm(aStore, qStore, dimension);

            //SquareMatrix A = new SquareMatrix(aStore, dimension);
            SquareMatrix Q = new SquareMatrix(qStore, dimension);

            Complex[] eigenvalues; Complex[][] eigenvectors;

            // Extract the eigenvalues and eigenvectors of the Schur form matrix
            SquareMatrixAlgorithms.SchurEigensystem(aStore, dimension, out eigenvalues, out eigenvectors);

            // transform eigenvectors of schur form into eigenvectors of original matrix
            // while we are at it, normalize so largest component has value 1
            for (int i = 0; i < dimension; i++) {
                Complex[] v = new Complex[dimension];
                double norm = 0.0;
                for (int j = 0; j < dimension; j++) {
                    for (int k = 0; k < dimension; k++) {
                        v[j] += Q[j, k] * eigenvectors[i][k];
                    }
                    norm = Math.Max(norm, Math.Max(Math.Abs(v[j].Re), Math.Abs(v[j].Im)));
                }
                for (int j = 0; j < dimension; j++) {
                    v[j] = v[j] / norm;
                }
                eigenvectors[i] = v;
            }

            ComplexEigendecomposition eigensystem = new ComplexEigendecomposition(dimension, eigenvalues, eigenvectors);
            return (eigensystem);

        }

        /// <summary>
        /// Computes the singular value decomposition of the square matrix.
        /// </summary>
        /// <returns>The singular value decomposition of the matrix.</returns>
        /// <remarks>
        /// <para>Singular value decomposition is an advanced matrix decomposition technique that can be applied
        /// to all matrices, including non-square and singular square matrices.</para>
        /// </remarks>
        public SingularValueDecomposition SingularValueDecomposition () {
            double[] copy = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, dimension, dimension);
            RectangularMatrix r = new RectangularMatrix(copy, dimension, dimension);
            return (r.SingularValueDecomposition());
        }
       
#if FUTURE

        private static double[] Balance (SquareMatrix A) {

            int d = A.Dimension;

            // a vector to keep track of scale factors
            double[] rhos = new double[d];
            for (int i = 0; i < d; i++) { rhos[i] = 1.0; }

            // iterate over dimensions
            for (int i = 0; i < d; i++) {

                // for each dimension, we will apply a transform Q = diag(1,...,1,rho,1,...,1) where rho is the i'th diagonal element
                // Q A Q^-1 thus multiplies elements of the i'th row by rho and divides elements of the i'th column by rho
                // we choose rho so that (1) the row and column norms are close and (2) rho is an exact power of 2
                double rSum = 0.0;
                double cSum = 0.0;
                for (int j = 0; j < d; j++) {
                    if (i == j) continue;
                    rSum += Math.Abs(A[i, j]);
                    cSum += Math.Abs(A[j, i]);
                }
                if ((rSum == 0.0) || (cSum == 0.0)) continue;
                double rho = ClosestPowerOfTwo(Math.Sqrt(rSum / cSum));
                if (rho != 1.0) {
                    for (int j = 0; j < d; j++) {
                        if (i == j) continue;
                        A[i, j] = A[i, j] * rho;
                        A[j, i] = A[j, i] / rho;
                    }
                    rhos[i] = rhos[i] * rho;
                }
            }

            return (rhos);

        }

        private static double ClosestPowerOfTwo (double x) {
            double y = Math.Round(Math.Log(x) / Math.Log(2.0));
            return (Math.Pow(2.0, y));
        }

#endif

        /// <summary>
        /// Computes a QR decomposition of the matrix.
        /// </summary>
        /// <returns>A QR decomposition of the matrix.</returns>
        public SquareQRDecomposition QRDecomposition () {
            double[] rStore = MatrixAlgorithms.Copy(store, offset, rowStride, colStride, dimension, dimension);
            double[] qtStore = SquareMatrixAlgorithms.CreateUnitMatrix(dimension);
            MatrixAlgorithms.QRDecompose(rStore, qtStore, dimension, dimension);
            return (new SquareQRDecomposition(qtStore, rStore, dimension));
        }

        /// <summary>
        /// Computes the matrix raised to the given power.
        /// </summary>
        /// <param name="n">The power to which to raise the matrix, which must be positive.</param>
        /// <returns>The matrix A<sup>n</sup>.</returns>
        /// <remarks>
        /// <para>This method uses exponentiation-by-squaring, so typically only of order log(n) matrix multiplications are required.
        /// This improves not only the speed with which the matrix power is computed, but also the accuracy.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Exponentiation_by_squaring"/>
        public SquareMatrix Power (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            } else if (n == 0) {
                return (new SquareMatrix(SquareMatrixAlgorithms.CreateUnitMatrix(dimension), dimension));
            } else if (n == 1) {
                return (this.Copy());
            } else {
                double[] pStore = SquareMatrixAlgorithms.Power(store, dimension, n);
                return (new SquareMatrix(pStore, dimension));
            }
        }

        // operators

        /// <summary>
        /// Adds two real, square matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The sum matrix <paramref name="A"/> + <paramref name="B"/>.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="B"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="A"/> is not the same as the dimension of <paramref name="B"/>.</exception>
        public static SquareMatrix operator + (SquareMatrix A, SquareMatrix B) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.dimension != B.dimension) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Add(
                A.store, A.offset, A.rowStride, A.colStride,
                B.store, B.offset, B.rowStride, B.colStride,
                A.dimension, A.dimension
            );
            return (new SquareMatrix(abStore, A.dimension));
        }

        /// <summary>
        /// Computes the difference of two square matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The difference <paramref name="A"/> - <paramref name="B"/>.</returns>
        /// <remarks>
        /// <para>Matrix subtraction is an O(N<sup>2</sup>) process.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="B"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="A"/> is not the same as the dimension of <paramref name="B"/>.</exception>
        public static SquareMatrix operator - (SquareMatrix A, SquareMatrix B) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.dimension != B.dimension) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Subtract(
                A.store, A.offset, A.rowStride, A.colStride,
                B.store, B.offset, B.rowStride, B.colStride,
                A.dimension, A.dimension
            );
            return (new SquareMatrix(abStore, A.dimension));
        }

        /// <summary>
        /// Computes the product of two square matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The product <paramref name="A"/> * <paramref name="B"/>.</returns>
        /// <remarks>
        /// <para>Note that matrix multiplication is not commutative, i.e. M1*M2 is generally not the same as M2*M1.</para>
        /// <para>Matrix multiplication is an O(N<sup>3</sup>) process.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="B"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="A"/> is not the same as the dimension of <paramref name="B"/>.</exception>
        public static SquareMatrix operator * (SquareMatrix A, SquareMatrix B) {
            // this is faster than the base operator, because it knows about the underlying structure
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (B == null) throw new ArgumentNullException(nameof(B));
            if (A.dimension != B.dimension) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Multiply(
                A.store, A.offset, A.rowStride, A.colStride,
                B.store, B.offset, B.rowStride, B.colStride,
                A.dimension, A.dimension, A.dimension
            );
            return (new SquareMatrix(abStore, A.dimension));
        }

        // mixed arithmetic

        /// <summary>
        /// Computes the product of a square matrix and a column vector.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="v">The column vector.</param>
        /// <returns>The column vector Av.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> or <paramref name="v"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="A"/> is not the same as the dimension of <paramref name="v"/>.</exception>
        public static ColumnVector operator * (SquareMatrix A, ColumnVector v) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (v == null) throw new ArgumentNullException(nameof(v));
            if (A.dimension != v.dimension) throw new DimensionMismatchException();
            double[] avStore = new double[A.dimension];
            Blas2.dGemv(A.store, A.offset, A.rowStride, A.colStride, v.store, v.offset, v.stride, avStore, 0, 1, A.dimension, A.dimension);
            return (new ColumnVector(avStore, A.dimension));
        }

        // should above method be in ColumnVector type so that all references to vector fields are in VectorBase, RowVector, and ColumnVector?

        /// <summary>
        /// Multiply a real, square matrix by a real constant.
        /// </summary>
        /// <param name="alpha">The constant.</param>
        /// <param name="A">The matrix.</param>
        /// <returns>The product &#x3B1;A.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> is <see langword="null"/>.</exception>
        public static SquareMatrix operator * (double alpha, SquareMatrix A) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            double[] store = MatrixAlgorithms.Multiply(alpha, A.store, A.offset, A.rowStride, A.colStride, A.dimension, A.dimension);
            return (new SquareMatrix(store, A.dimension));
        }

        /// <summary>
        /// Divides a real, square matrix by a real constant.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="alpha">The constant.</param>
        /// <returns>The quotient A/&#x3B1;.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> is <see langword="null"/>.</exception>
        public static SquareMatrix operator / (SquareMatrix A, double alpha) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            double[] store = MatrixAlgorithms.Multiply(1.0 / alpha, A.store, A.offset, A.rowStride, A.colStride, A.dimension, A.dimension);
            return (new SquareMatrix(store, A.dimension));
        }

        /// <summary>
        /// Negates a real, square matrix.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <returns>The matrix -A.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="A"/> is <see langword="null"/>.</exception>
        public static SquareMatrix operator - (SquareMatrix A) {
            if (A == null) throw new ArgumentNullException(nameof(A));
            double[] store = MatrixAlgorithms.Multiply(-1.0, A.store, A.offset, A.rowStride, A.colStride, A.dimension, A.dimension);
            return (new SquareMatrix(store, A.dimension));
        }

    }



}
