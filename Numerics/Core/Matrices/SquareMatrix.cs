using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

using Meta.Numerics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a square matrix.
    /// </summary>
    public sealed class SquareMatrix : AnySquareMatrix {

        private int dimension;
        private double[] store;

        /// <summary>
        /// Initializes a new square matrix.
        /// </summary>
        /// <param name="dimension">The dimension of the matrix, which must be positive.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="dimension"/> &lt; 1.</exception>
        public SquareMatrix (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException("dimension");
            this.dimension = dimension;
            this.store = MatrixAlgorithms.AllocateStorage(dimension, dimension);
        }

        internal SquareMatrix (double[] storage, int dimension) {
            this.store = storage;
            this.dimension = dimension;
        }

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
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
                return (MatrixAlgorithms.GetEntry(store, dimension, dimension, r, c));
            }
            set {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
                MatrixAlgorithms.SetEntry(store, dimension, dimension, r, c, value);
            }
        }

        /// <inheritdoc />
        public override double OneNorm () {
            return (MatrixAlgorithms.OneNorm(store, dimension, dimension));
        }

        /// <inheritdoc />
        public override double InfinityNorm () {
            return (MatrixAlgorithms.InfinityNorm(store, dimension, dimension));
        }

        // use for fast access
        // this bypasses bounds checks, so make sure the bounds are right!
        /*
        internal double GetEntry (int r, int c) {
            //Console.WriteLine("get from {0}", dimension * c + r);
            return (store[dimension * c + r]);
            //return (values[r * dimension + c]);
        }

        internal void SetEntry (int r, int c, double value) {
            //Console.WriteLine("set {0} at {1}", value, dimension * c + r);
            store[dimension * c + r] = value;
            //values[r * dimension + c] = value;
        }
        */

        /// <summary>
        /// Returns a vector representing a given row of the matrix.
        /// </summary>
        /// <param name="r">The (zero-based) row number to return.</param>
        /// <returns>An independent copy of the specified row.</returns>
        /// <remarks>The returned vector is not linked to the matrix. If an entry in the matrix is updated after this method
        /// is called, the returned object will continue to represent a row of the original, not the updated, matrix. Similiarly,
        /// updates to the elements of the returned vector will not update the original matrix.</remarks>
        public override RowVector Row (int r) {
            if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
            RowVector row = new RowVector(dimension);
            for (int c = 0; c < dimension; c++) {
                row[c] = this[r, c];
            }
            return (row);
        }

        /// <summary>
        /// Gets a copy of one column of the the matrix.
        /// </summary>
        /// <param name="c">The (zero-based) column number to return.</param>
        /// <returns>An independent copy of the specificed column.</returns>
        /// <remarks>The returned vector is not linked to the matrix. If an entry in the matrix is updated after this method
        /// is called, the returned object will continue to represent a row of the original, not the updated, matrix. Similiarly,
        /// updates to the elements of the returned vector will not update the original matrix.</remarks>
        public override ColumnVector Column (int c) {
            if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
            ColumnVector column = new ColumnVector(dimension);
            for (int r = 0; r < dimension; r++) {
                column[r] = this[r, c];
            }
            return (column);
        }

        /// <summary>
        /// Copies the matrix.
        /// </summary>
        /// <returns>An independent copy of the matrix.</returns>
        public SquareMatrix Copy () {
            double[] cStore = MatrixAlgorithms.Copy(store, dimension, dimension);
            return (new SquareMatrix(cStore, dimension));
        }

        /*
        IMatrix IMatrix.Clone () {
            return (Clone());
        }
        */

        /// <summary>
        /// Creates a transpose of the matrix.
        /// </summary>
        /// <returns>The matrix transpose M<sup>T</sup>.</returns>
        public SquareMatrix Transpose () {
            double[] tStore = MatrixAlgorithms.Transpose(store, dimension, dimension);
            return (new SquareMatrix(tStore, dimension));
        }

        /// <summary>
        /// Computes the inverse of the matrix.
        /// </summary>
        /// <returns>The matrix inverse M<sup>-1</sup>.</returns>
        /// <remarks>
        /// <para>The inverse of a matrix M is a matrix M<sup>-1</sup> such that M<sup>-1</sup>M = I, whhere I is the identity matrix.</para>
        /// <para>If the matrix is singular, inversion is not possible. In that case, this method will fail with a <see cref="DivideByZeroException"/>.</para>
        /// <para>The inversion of a matrix is an O(N<sup>3</sup>) operation.</para>
        /// </remarks>
        /// <exception cref="DivideByZeroException">The matrix is singular.</exception>
        public SquareMatrix Inverse () {
            double[] iStore = MatrixAlgorithms.Copy(store, dimension, dimension);
            SquareMatrixAlgorithms.GaussJordanInvert(iStore, dimension);
            return(new SquareMatrix(iStore, dimension));
        }

        /*
        /// <summary>
        /// Inverts the matrix, in place.
        /// </summary>
        /// <remarks>
        /// <para>This method replaces the elements of M with the elements of M<sup>-1</sup>.</para>
        /// <para>In place matrix inversion saves memory, since seperate storage of M and M<sup>-1</sup> is not required.</para></remarks>
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
        public SquareLUDecomposition LUDecomposition () {

            // copy the matrix content
            double[] luStore = MatrixAlgorithms.Copy(store, dimension, dimension);

            // prepare an initial permutation and parity
            int[] permutation = new int[dimension];
            for (int i = 0; i < permutation.Length; i++) {
                permutation[i] = i;
            }
            int parity = 1;

            // do the LU decomposition
            SquareMatrixAlgorithms.LUDecompose(luStore, permutation, ref parity, dimension);

            // package it up and return it
            SquareLUDecomposition LU = new SquareLUDecomposition(luStore, permutation, parity, dimension);

            return (LU);

        }

        /// <summary>
        /// Computes the eigenvalues of the matrix.
        /// </summary>
        /// <returns>The eigenvalues of the matrix.</returns>
        /// <seealso cref="Eigensystem"/>
        public Complex[] Eigenvalues () {

            double[] scratch = MatrixAlgorithms.Copy(store, dimension, dimension);
            SquareMatrixAlgorithms.ReduceToHessenberg(scratch, null, dimension);
            Complex[] lambdas = SquareMatrixAlgorithms.ExtractEigenvalues(scratch, null, dimension);
            return (lambdas);
            
        }

        /// <summary>
        /// Computes the eigenvalues and eigenvectors of the matrix.
        /// </summary>
        /// <returns>A representation of the eigenvalues and eigenvectors of the matrix.</returns>
        /// <remarks>
        /// <para>For a generic vector v and matrix M, Mv = u will point in some direction with no particular relationship to v.
        /// The eigenvectors of a matrix M are vectors z that satisfy Mz = &#x3BB;z, i.e. multiplying an eigenvector by a
        /// matrix reproduces the same vector, up to a prortionality constant &#x3BB; called the eigenvalue.</para>
        /// <para>For v to be an eigenvector of M with eigenvalue &#x3BB;, (M - &#x3BB;I)z = 0. But for a matrix to
        /// anihilate any non-zero vector, that matrix must have determinant, so det(M - &#x3BB;I)=0. For a matrix of
        /// order N, this is an equation for the roots of a polynomial of order N. Since an order-N polynomial always has exactly
        /// N roots, an order-N matrix always has exactly N eigenvalues.</para>
        /// <para>Since a polynomial with real coefficients can still have complex roots, a real square matrix can nonetheless
        /// have complex eigenvalues (and correspondly complex eigenvectors). However, again like the complex roots of a real
        /// polynomial, such eigenvalues will always occurs in complex-conjugate pairs.</para>
        /// <para>Although the eigenvalue polynomial ensures that an order-N matrix has N eigenvalues, it can occur that there
        /// are not N corresponding independent eigenvectors. A matrix with fewer eigenvectors than eigenvalues is called
        /// defective. Like singularity, defectiveness represents a delecate balance between the elements of a matrix that can
        /// typically be disturbed by just an infinitesimal perturbation of elements. Because of round-off-error, then, floating-point
        /// algorithms cannot reliably identify defective matrices. Instead, this method will return a full set of eigenvectors,
        /// but some eigenvectors, corresponding to very nearly equal eigenvalues, will be very nearly parallel.</para>
        /// <para>While a generic square matrix can be defective, many subspecies of square matrices are guaranteed not to be.
        /// This includes Markov matrices, orthogonal matrices, and symmetric matrices.</para>
        /// <para>Determining the eigenvalues and eigenvectors of a matrix is an O(N<sup>3</sup>) operation. If you need only the
        /// eigenvalues of a matrix, the <see cref="Eigenvalues"/> method is more efficient.</para>
        /// </remarks>
        public ComplexEigensystem Eigensystem () {

            double[] aStore = MatrixAlgorithms.Copy(store, dimension, dimension);
            double[] qStore = SquareMatrixAlgorithms.CreateUnitMatrix(dimension);
            SquareMatrixAlgorithms.ReduceToHessenberg(aStore, qStore, dimension);
            Complex[] eigenvalues = SquareMatrixAlgorithms.ExtractEigenvalues(aStore, qStore, dimension);

            SquareMatrix A = new SquareMatrix(aStore, dimension);
            SquareMatrix Q = new SquareMatrix(qStore, dimension);


            //SquareMatrix A = this.Clone();
            //Debug.WriteLine("Start A=");
            //PrintMatrix(A);

            // start with an identity matrix to track the transforms
            //SquareMatrix Q = new SquareMatrix(dimension);
            //for (int i = 0; i < dimension; i++) Q[i, i] = 1.0;

            // balance the matrix
            /*
            Balance(A);
            Debug.WriteLine("Balanced A=");
            PrintMatrix(A);
            Debug.WriteLine("Transform Q=");
            PrintMatrix(Q);
            Debug.WriteLine("Test Q^T M Q = A");
            PrintMatrix(Q.Transpose() * this * Q);
            */

            // reduce to Hessenberg form
            //ReduceToHessenberg(A, Q);
            //Debug.WriteLine("Hessenberg A=");
            //PrintMatrix(A);
            //Debug.WriteLine("Transform Q=");
            //PrintMatrix(Q);
            //Debug.WriteLine("Test Q^T M Q = A");
            //PrintMatrix(Q.Transpose() * this * Q);
            
            // reduce to Schur form, extracting eigenvalues as we go
            //Complex[] eigenvalues = ExtractEigenvalues(A, Q);
            //Debug.WriteLine("Schur A=");
            //PrintMatrix(A);
            //Debug.WriteLine("Transform Q=");
            //PrintMatrix(Q);
            //Debug.WriteLine("Test Q^T M Q = A");
            //PrintMatrix(Q.Transpose() * this * Q);

            // get eigenvectors
            Complex[,] eigenvectors = ExtractEigenvectors(A, Q, eigenvalues);
            NormalizeEigenvectors(eigenvectors);

            ComplexEigensystem eigensystem = new ComplexEigensystem(dimension, eigenvalues, eigenvectors);
            return (eigensystem);

        }

        // given A in real Schur form (i.e. upper triangular except for 2 X 2 blocks along the diagonal), and Q
        // that got us there, extract eigenvectors of A and apply Q to transform them to the eigenvectors of the
        // original matrix

        private static Complex[,] ExtractEigenvectors (SquareMatrix A, SquareMatrix Q, Complex[] e) {

            // a store which will be used to store each eigenvalue of A
            int dim = A.Dimension;

            // a store for the eigenvectors of the the original matrix,
            // which are th
            Complex[,] X = new Complex[dim, dim];

            // get eigenvectors of A
            for (int k = 0; k < dim; k++) {

                // find the kth eigenvector of the (nearly) upper triangular matrix A
                Complex[] b;
                int imax;
                //Console.WriteLine("Computing eigenvector {0}", k);

                if ((k > 0) && (A[k, k - 1] != 0.0)) {
                    // an extra element to the left
                    b = new Complex[k + 1];
                    b[k] = 1.0;
                    b[k - 1] = A[k - 1, k] / (e[k] - A[k - 1, k - 1]);
                    //Console.WriteLine("b[{0}] = {1}", k, b[k]);
                    //Console.WriteLine("b[{0}] = {1}", k - 1, b[k - 1]);
                    imax = k - 2;
                } else if (((k + 1) < dim) && (A[k + 1, k] != 0.0)) {
                    // an extra element below
                    b = new Complex[k + 2];
                    b[k + 1] = 1.0;
                    b[k] = A[k, k + 1] / (e[k] - A[k, k]);
                    //b[k + 1] = A[k+1,k] / (e[k] - A[k+1,k+1]);
                    //Console.WriteLine("b[{0}] = {1}", k+1, b[k+1]);
                    //Console.WriteLine("b[{0}] = {1}", k, b[k]);
                    imax = k - 1;
                } else {
                    // the pure upper triangular case
                    b = new Complex[k + 1];
                    b[k] = 1.0;
                    imax = k - 1;
                }

                for (int i = imax; i >= 0; i--) {
                    //Console.WriteLine("Component {0}", i);
                    if ((i == 0) || (A[i, i - 1] == 0.0)) {
                        // system is pure tridiagonal, so solution is straightforward
                        Complex s = 0.0;
                        //Console.WriteLine("Start from A[{0},{1}] = {2}", i, k, s);
                        for (int j = i + 1; j < b.Length; j++) {
                            s += A[i, j] * b[j];
                            //Console.WriteLine("Add A[{0},{1}] * B[{2}] = {3} * {4}", i, j, j, A[i, j], b[j]);
                        }
                        Complex t = e[k] - A[i, i];
                        if (s == 0.0) {
                            // deal with trivial Shur form; this arises e.g. for decomposition of unit matrix
                            // without this, we get zero divided by zero, which is NaN
                            b[i] = 0.0;
                        } else {
                            b[i] = s / t;
                        }
                    } else {
                        // system has a sub-diagonal element, so solution is a little more complex
                        Complex s1 = 0.0;
                        Complex s2 = 0.0;
                        for (int j = i + 1; j < b.Length; j++) {
                            s1 += A[i - 1, j] * b[j];
                            s2 += A[i, j] * b[j];
                        }
                        Complex t1 = e[k] - A[i - 1, i - 1];
                        Complex t2 = e[k] - A[i, i];
                        b[i] = (s2 + A[i,i-1] * s1 / t1) / (t2 - A[i,i-1] * A[i-1,i]/ t1);
                    }
                    //Console.WriteLine("b[{0}] = {1}", i, b[i]);
                }

                // transform it to the original basis
                for (int i = 0; i < dim; i++) {
                    Complex x = 0.0;
                    for (int j = 0; j < b.Length; j++) {
                        x += Q[i, j] * b[j];
                    }
                    X[i,k] = x;
                    //Console.Write("{0}  ", x);
                }
                //Console.WriteLine();
                //Console.ReadLine();

            }

            //

            return (X);

        }

        // renormalize eigenvectors so that their 2-norm is unity

        private static void NormalizeEigenvectors (Complex[,] Z) {

            int d = Z.GetLength(0);

            // loop over eigenvectors
            for (int n = 0; n < d; n++) {

                // find the normalization factor
                double x = 0.0;
                for (int i = 0; i < d; i++) {
                    Complex z = Z[i,n];
                    x += z.Re * z.Re + z.Im * z.Im;
                }
                x = Math.Sqrt(x);

                // divide by it
                for (int i = 0; i < d; i++) {
                    Z[i, n] = Z[i, n] / x;
                }

            }

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
            double[] copy = MatrixAlgorithms.Copy(store, dimension, dimension);
            RectangularMatrix r = new RectangularMatrix(copy, dimension, dimension);
            return (r.SingularValueDecomposition());
        }

#if PAST

        private static void WilkersonOneStep (double[] store, int dimension, int a, int n) {

            //Console.WriteLine("Before:");
            //Write(store, dimension, dimension);

            int x = n - 2;
            int y = n - 1;
            int z = n;
            double Axx = MatrixAlgorithms.GetEntry(store, dimension, dimension, x, x);
            double Axy = MatrixAlgorithms.GetEntry(store, dimension, dimension, x, y);
            double Axz = MatrixAlgorithms.GetEntry(store, dimension, dimension, x, z);
            double Ayx = MatrixAlgorithms.GetEntry(store, dimension, dimension, y, x);
            double Ayy = MatrixAlgorithms.GetEntry(store, dimension, dimension, y, y);
            double Ayz = MatrixAlgorithms.GetEntry(store, dimension, dimension, y, z);
            double Azy = MatrixAlgorithms.GetEntry(store, dimension, dimension, z, y);
            double Azz = MatrixAlgorithms.GetEntry(store, dimension, dimension, z, z);

            
            double[] roots = CubicRealRoots(
                -(Axx + Ayy + Azz),
                Axx * Ayy + Axx * Azz + Ayy * Azz - Axy * Ayx - Ayz * Ayz,
                Axy * Ayx * Azz + Axx * Ayz * Azy - Axx * Ayy * Azz - Axz * Azy * Ayx
            );

            double mu = roots[0];
            for (int i = 1; i < roots.Length; i++) {
                if (Math.Abs(roots[i] - Azz) < Math.Abs(mu - Azz)) mu = roots[i];
            }
            //Console.WriteLine("mu = {0}", mu);
            
            ShiftedQRStep(store, dimension, a, n, mu);

            //Console.WriteLine("After:");
            //Write(store, dimension, dimension);

        }

        private static double[] CubicRealRoots (double a, double b, double c) {
            double Q = (a * a - 3.0 * b) / 9.0;
            double R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

            double Q3 = Q * Q * Q;
            double R2 = R * R;

            //double S = R * R / (Q * Q * Q);
            //Console.WriteLine("Q={0}, R={1} S={2}", Q, R, S);

            if (Q3 > R2) {
                double t = Math.Acos(R/Math.Sqrt(Q3));
                double x0 = -a / 3.0;
                double x1 = -2.0 * Math.Sqrt(Q);
                return (new double[] {
                    x0 + x1 * Math.Cos(t / 3.0),
                    x0 + x1 * Math.Cos((t + 2.0 * Math.PI) / 3.0),
                    x0 + x1 * Math.Cos((t - 2.0 * Math.PI) / 3.0)
                });
            } else {
                double A = Math.Pow(Math.Abs(R) + Math.Sqrt(R2 - Q3), 1.0 / 3.0);
                if (R >= 0.0) A = -A;

                double B;
                if (A != 0.0) {
                    B = Q / A;
                } else {
                    B = 0.0;
                }
                return (new double[] { A + B - a / 3.0 });

                throw new NotImplementedException();
            }
        }

        // a Householder reflection matrix is a rank-1 update to the identity.
        //   P = I - b v v^T
        // Unitarity requires b = 2 / |v|^2. To anihilate all but the first components of vector x,
        //   P x = a e_1
        // we must have
        //   v = x +/- |x| e_1
        // that is, all elements the same except the first, from which we have either added or subtracted |x|.
        // (We choose the sign that avoids canelation.) This makes
        //   a = -/+ |x|
        // 

        public static void GenerateHouseholderReflection (double[] store, int offset, int stride, int count, out double a) {
            double xm = Blas1.dNrm2(store, offset, stride, count);
            if (xm == 0.0) {
                a = 0.0;
            } else {
                double x0 = store[offset];
                double u0;
                if (x0 < 0.0) {
                    // subtract |x| from (negative) x_0, avoiding cancelation (and making it more negative)
                    u0 = x0 - xm;
                    a = xm;
                } else {
                    // add |x| to (positive) x_0, avoiding cancelation (and making it more positive)
                    u0 = x0 + xm;
                    a = -xm;
                }
                store[offset] = u0;
                // rescale u to make |u|^2 = 2 and avoid having to multiply later
                // this will not result in division by zero, because xm > 0 and u0 > xm
                double um = Math.Sqrt(xm * Math.Abs(u0));
                Blas1.dScal(1.0 / um, store, offset, stride, count);
            }
        }

        public static double[] GenerateHouseholderReflection (double[] store, int offset, int stride, int count) {

            // copy the vector, so we have seperate storage for the original column and the Householder vector
            double[] v = new double[count];
            Blas1.dCopy(store, offset, stride, v, 0, 1, count);

            // replace v with the Householder vector, and the first element of the original column by the value produced
            GenerateHouseholderReflection(v, 0, 1, count, out store[offset]);

            // zero the other column elements
            for (int i = 1; i < count; i++) {
                offset += stride;
                store[offset] = 0.0;
            }

            return (v);
        }

        public static void ApplyHouseholderReflection (
            double[] uStore, int uOffset, int uStride,
            double[] yStore, int yOffset, int yStride,
            int count
        ) {
            // (1 - b v v^T) y = y - (b v^T y) v
            double s = Blas1.dDot(uStore, uOffset, uStride, yStore, yOffset, yStride, count);
            Blas1.dAxpy(-s, uStore, uOffset, uStride, yStore, yOffset, yStride, count);
        }


        // XXXXX    XXXXX
        // XXXXX    XXXXX
        // 00XXX -> X0XXX ->
        // 000XX    X00XX
        // 000XX    000XX

        public static void Write (double[] store, int rows, int cols) {
            Console.WriteLine("M");
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < cols; c++) {
                    Console.Write("  {0}", store[rows * c + r]);
                }
                Console.WriteLine();
            }
        }

        public static void Write (SquareMatrix A) {
            Console.WriteLine("A");
            for (int r = 0; r < A.Dimension; r++) {
                for (int c = 0; c < A.Dimension; c++) {
                    Console.Write("  {0}", A[r,c]);
                }
                Console.WriteLine();
            }
        }
#endif
       
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
            double[] rStore = MatrixAlgorithms.Copy(store, dimension, dimension);
            double[] qtStore = SquareMatrixAlgorithms.CreateUnitMatrix(dimension);
            MatrixAlgorithms.QRDecompose(rStore, qtStore, dimension, dimension);
            return (new SquareQRDecomposition(qtStore, rStore, dimension));
        }

        // operators

        /// <summary>
        /// Adds two real, square matrices.
        /// </summary>
        /// <param name="A">The first matrix.</param>
        /// <param name="B">The second matrix.</param>
        /// <returns>The sum matrix <paramref name="A"/> + <paramref name="B"/>.</returns>
        public static SquareMatrix operator + (SquareMatrix A, SquareMatrix B) {
            if (A == null) throw new ArgumentNullException("A");
            if (B == null) throw new ArgumentNullException("B");
            if (A.dimension != B.dimension) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Add(A.store, B.store, A.dimension, A.dimension);
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
        public static SquareMatrix operator - (SquareMatrix A, SquareMatrix B) {
            if (A == null) throw new ArgumentNullException("A");
            if (B == null) throw new ArgumentNullException("B");
            if (A.dimension != B.dimension) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Subtract(A.store, B.store, A.dimension, A.dimension);
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
        public static SquareMatrix operator * (SquareMatrix A, SquareMatrix B) {
            // this is faster than the base operator, because it knows about the underlying structure
            if (A == null) throw new ArgumentNullException("A");
            if (B == null) throw new ArgumentNullException("B");
            if (A.dimension != B.dimension) throw new DimensionMismatchException();
            double[] abStore = MatrixAlgorithms.Multiply(A.store, A.dimension, A.dimension, B.store, B.dimension, B.dimension);
            return (new SquareMatrix(abStore, A.dimension));
        }

        // mixed arithmetic

        /// <summary>
        /// Multiply a real, square matrix by a real constant.
        /// </summary>
        /// <param name="alpha">The constant.</param>
        /// <param name="A">The matrix.</param>
        /// <returns>The product aA.</returns>
        public static SquareMatrix operator * (double alpha, SquareMatrix A) {
            if (A == null) throw new ArgumentNullException("A");
            double[] store = MatrixAlgorithms.Multiply(alpha, A.store, A.dimension, A.dimension);
            return (new SquareMatrix(store, A.dimension));
        }

        /// <summary>
        /// Negates a real, square matrix.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <returns>The matrix -A.</returns>
        public static SquareMatrix operator - (SquareMatrix A) {
            if (A == null) throw new ArgumentNullException("A");
            double[] store = MatrixAlgorithms.Multiply(-1.0, A.store, A.dimension, A.dimension);
            return (new SquareMatrix(store, A.dimension));
        }

        /*

        /// <summary>
        /// Computes the the quotient of a square matrix and a real number.
        /// </summary>
        /// <param name="M">The matrix.</param>
        /// <param name="x">The real number.</param>
        /// <returns>The quotient <paramref name="M"/>/<paramref name="x"/>.</returns>
        public static SquareMatrix operator / (SquareMatrix M, double x) {
            return (Multiply(1.0 / x, M));
        }
        */
    }



}
