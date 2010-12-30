using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

using Meta.Numerics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a square matrix.
    /// </summary>
    public sealed class SquareMatrix : SquareMatrixBase {

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
            
            //SquareMatrix A = new SquareMatrix(scratch, dimension);
            //Complex[] eigenvalues = ExtractEigenvalues(A, null);
            
            /*
            //SquareMatrix A = this.Clone();
            //ReduceToHessenberg(A, null);

            s.Stop();
            Console.WriteLine(s.ElapsedMilliseconds);

            //A.Write(Console.Out);
            Complex[] eigenvalues = ExtractEigenvalues(A, null);
            return (eigenvalues);
            */

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

#if PAST

        private static void ReduceToHessenberg (SquareMatrix A, SquareMatrix Q) {

            // note the dimension of the problem
            int d = A.Dimension;

            // iterate over columns
            for (int k=0; k<(d-2); k++) {

                // determine the Householder transform P that will zero the elements of the column below the sub-diagonal
                double[] u = new double[d-1-k];
                for (int rp=0; rp < u.Length; rp++) {
                    int r = k + 1 + rp;
                    u[rp] = A[r,k]; 
                }
                double x = ComputeHouseholderVector(ref u);

                // compute P * A
                // use P = ( 1 - u u^T) so P * A = A - u (u^T * A) = A - u v^T where v^T = u^T A
                for (int c = k + 1; c < d; c++) {
                    double v = 0.0;
                    for (int rp = 0; rp < u.Length; rp++) {
                        int r = k + 1 + rp;
                        v += u[rp] * A[r, c];
                    }
                    for (int rp = 0; rp < u.Length; rp++) {
                        int r = k + 1 + rp;
                        A[r, c] -= u[rp] * v;
                    }
                }

                // note that in the outer loop above c would naively go from 0 to dimension, but for c < k, we know that
                // A[r,c] = 0 for r > k + 1, because those columns are already in Hessenberg form, so we will get v = 0 for them.

                // for c = k, we know that the transform will make the subdiagonal element x and the elements below it zero;
                // we can save computational effort and avoid round-off error by just setting them to that result.
                A[k + 1, k] = x;
                for (int r = k + 2; r < d; r++) {
                    A[r, k] = 0.0;
                }

                // compute A * P
                // use P = ( 1 - u u^T) so A * P = A - (A u) u^T = A - v u^T where v = A u
                for (int r = 0; r < d; r++) {
                    double v = 0.0;
                    for (int cp = 0; cp < u.Length; cp++) {
                        int c = k + 1 + cp;
                        v += A[r, c] * u[cp];
                        //if (A[r, c] == 0.0) Console.WriteLine("{0},{1} => 0", r, c);
                    }
                    for (int cp = 0; cp < u.Length; cp++) {
                        int c = k + 1 + cp;
                        A[r, c] -= v * u[cp];
                    }
                }

                //PrintMatrix(A);
                //Console.ReadLine();

                // if we we are keeping track of transformations, record the transform that got us here
                if (Q != null) {
                    // compute Q * P
                    for (int r = 0; r < d; r++) {
                        double v = 0.0;
                        for (int cp = 0; cp < u.Length; cp++) {
                            int c = k + 1 + cp;
                            v += Q[r, c] * u[cp];
                        }
                        for (int cp = 0; cp < u.Length; cp++) {
                            int c = k + 1 + cp;
                            Q[r, c] -= v * u[cp];
                        }
                    }
                }

                // done with that column

            }

            // done with all columns

        }


        // updates a transform matrix Q with a (right) Householder transform
        // this is used to keep track of accumulated transforms

        /*
        private void HouseholderTransform (ref double[,] Q, int offset, double[] u) {
            // compute Q * P
            for (int r = 0; r < dimension; r++) {
                double v = 0.0;
                for (int cp = 0; cp < u.Length; cp++) {
                    int c = offset + cp;
                    v += Q[r, c] * u[cp];
                }
                for (int cp = 0; cp < u.Length; cp++) {
                    int c = offset + cp;
                    Q[r, c] -= v * u[cp];
                }
            }
        }
        */

        // on return, upper Hessenberg part of the matrix (i.e. the upper right triangle and the first subdiagonal)
        // contains H and the returned matrix contains a unitary transform Q such that Q^(-1) A Q = H.
        /*
        public double[,] HessenbergReduce () {

            // start with an identity matrix to track the transforms
            double[,] Q = new double[dimension, dimension];
            for (int i = 0; i < dimension; i++) Q[i, i] = 1.0;

            // iterate over columns,
            // doing a Householder transform for each column to eliminate the entires in rows below the sub-diagonal
            for (int c = 0; c < (dimension - 1); c++) {

                // determine the sub-diagonal norm
                double x = 0.0;
                for (int r = c + 1; r < dimension; r++) {
                    double v = values[r, c];
                    x += v * v;
                }
                x = Math.Sqrt(x);

                // if the norm is zero, the column in already Hessenberg
                if (x == 0.0) continue;

                // choose the first component of the Householder transform vector to minimize roundoff error
                double s, d;
                if (values[c + 1, c] < 0) {
                    values[c + 1, c] -= x;
                    d = x;
                    s = -x * values[c + 1, c];
                } else {
                    values[c + 1, c] += x;
                    d = -x;
                    s = x * values[c + 1, c];
                }

                // if s=0, we have a problem

                // set A=AH
                for (int i = 0; i < dimension; i++) {
                    double sum = 0.0;
                    for (int k = c + 1; k < dimension; k++) {
                        sum += values[i, k] * values[k, c];
                    }
                    sum = sum / s;
                    for (int j = c + 1; j < dimension; j++) {
                        values[i, j] -= sum * values[j, c];
                    }
                }
                // set A=HA (= HAH)
                for (int j = c + 1; j < dimension; j++) {
                    double sum = 0.0;
                    for (int k = c + 1; k < dimension; k++) {
                        sum += values[k, c] * values[k, j];
                    }
                    sum = sum / s;
                    for (int i = c + 1; i < dimension; i++) {
                        values[i, j] -= values[i, c] * sum;
                    }
                }

                // update Q=QH to reflect the transform
                for (int i = 0; i < dimension; i++) {
                    double sum = 0.0;
                    for (int k = c + 1; k < dimension; k++) {
                        sum += Q[i, k] * values[k, c];
                    }
                    sum = sum / s;
                    for (int j = c + 1; j < dimension; j++) {
                        Q[i, j] -= sum * values[j, c];
                    }
                }

                // update subdiagonal elements in column c
                values[c + 1, c] = d;
                for (int i = c + 2; i < dimension; i++) {
                    values[i, c] = 0.0;
                }
            }

            return (Q);

        }
        */

        private static Complex[] ExtractEigenvalues (SquareMatrix A, SquareMatrix Q) {

            // note the dimension
            int dim = A.Dimension;

            // keep track of eigenvalues found
            Complex[] values = new Complex[dim]; // storage for the eigenvalues
            int c = 0; // count of eigenvalues found

            // keep track of the active area
            int a = 0; // the upper-left index of the active area of the  matrix
            int n = dim - 1; // the lower-right index of the active area of the matrix
            // the effective dimension de = n - a

            // maximum number of iterations
            int countMax = dim * 30;

			// keep track of iterations
			int count = 0;

            while (true) {

                // check whether we have all the eigenvalues
                if (c >= dim) break;

                //Debug.WriteLine(String.Format("iteration count = {0}", count));
                //Debug.WriteLine("A = ");
                //PrintMatrix(A);
                
                //Console.WriteLine("count = {0}", count);
                //A.Write(Console.Out);

                // zero negligible sub-diagonal elements
                for (int r = a+1; r <= n; r++) {
                    double x = Math.Abs(A[r, r]) + Math.Abs(A[r - 1, r - 1]);
                    if ((x + A[r, r - 1]) == x) {
                        A[r, r - 1] = 0.0;
                        if (a == (r - 1)) {
                            // we have uncovered an eigenvalue at the top
                            values[a] = A[a,a]; // record the eigenvalue
                            c++; // one more eigenvalue
                            //Debug.WriteLine(String.Format("Got eigenvalue {0} from top: {1}", c, values[a]));
                            a++; // active area shrinks by one from the top
                            //count = 0; // reset the iteration count
                        } else if (a == (r - 2)) {
                            // we have isolated a 2x2 matrix at the top
                            // compute its eigenvalues (move this to a subroutine)
                            double tr = A[a, a] + A[a + 1, a + 1];
                            double det = A[a, a] * A[a + 1, a + 1] - A[a, a + 1] * A[a + 1, a];
                            ExtractTwoByTwoEigenvalues(tr, det, out values[a], out values[a + 1]);
                            c += 2; // eigenvalue count increases by two
                            //Debug.WriteLine(String.Format("Got eigenvalues up to {0} from the top: {1}, {2}", c, values[a], values[a+1]));
                            a += 2; // active area shrinks by two from the top
                            //count = 0; // reset the iteration count
                        }
                    }
                }
                //Debug.WriteLine(String.Format("a={0}", a));

                // check again
                if (c >= dim) break;

                if ((n == 0) || (A[n, n - 1] == 0.0)) {
                    // we have isolated a single eigenvalue in the lower-right corner
                    values[n] = A[n, n]; // record the eigenvalue
                    c++; // one more eigenvalue
                    //Debug.WriteLine(String.Format("Got eigenvalue {0} from bottom: {1}", c, values[n]));
                    n--; // active area decreases by one from bottom
                    //count = 0; // reset the iteration count
                } else {
                    // look at the lower 2X2 matrix
                    int m = n - 1;
                    double tr = A[m, m] + A[n, n];
                    double det = A[m, m] * A[n, n] - A[m, n] * A[n, m];
                    // construct the eigenvalues of the 2 X 2 lower-right matrix
                    if ((m == 0) || (A[m, m - 1] == 0.0)) {
                        // we have isolated a 2 X 2 matrix

                        // compute its eigenvalues
                        ExtractTwoByTwoEigenvalues(tr, det, out values[m], out values[n]);
                        c += 2; // two more eigenvalues
                        //Debug.WriteLine(String.Format("Got eigenvalues up to {0} from bottom: {1}, {2}", c, values[m], values[n]));
                        n -= 2; // the active area decreases by two from the bottom
                        //count = 0; // reset the iteration count

                    } else {

                        // an ad hoc shift if we are not converging
                        // i see no evidence that this acomplishes anything!
                        /*
                        if ((count % 8) == 7) {
                            //Console.WriteLine("was tr={0} det={1}", tr, det);
                            tr = A[m, m-1] + A[n, n-1];
                            det = 0.0;
                            //Console.WriteLine("now tr={0} det={1}", tr, det);
                        }
                        */

                        // do a Francis implicit QR step to reduce sub-diagonal elements
                        FrancisTwoStep(A, Q, a, n, tr, det);

                        // up the iteration count
                        count++;
                        //if (count > (countMax-5)) {
                        //    PrintMatrix(A);
                        //}
                        if (count > countMax) {
                            throw new NonconvergenceException();
                        }

                    }

                }
            }

            // we have reduced the dimension to zero, so we have all the eigenvalues

            return (values);
        }

        // get the eigenvalues of a real 2x2 matrix

        private static void ExtractTwoByTwoEigenvalues (double tr, double det, out Complex e1, out Complex e2) {
            double p = tr / 2.0;
            double q2 = p * p - det;
            if (q2 < 0.0) {
                // if the descriminant is negative, the eigenvalues are complex
                double q = Math.Sqrt(-q2);
                e1 = new Complex(p, q);
                e2 = new Complex(p, -q);
            } else {
                // otherwise they are real
                double q = Math.Sqrt(q2);
                e1 = p + q;
                e2 = p - q;
            }
        }

        // do a Francis implicit QR step on A
        // the active area is between indices a and n
        // the sum of the eigenvalues is tr, the product is det
        // keep track of the transformations in Q, which may be null

        public static void FrancisTwoStep (SquareMatrix A, SquareMatrix Q, int a, int n, double tr, double det) {

            //Write(A);

            int dim = A.Dimension;
            //a = 0; // temporary
            int m = n - 1;

            // compute the first column of A' = A^2 - 2 Re(\lambda) A + |\lambda|^2 I
            // note 2 Re(\lambda) = tr and |\lambda|^2 = det
            // because A is Hessenberg, only the first three elements of this column are non-zero

            double[] u = new double[3];
            u[0] = A[a, a] * (A[a, a] - tr) + A[a, a+1] * A[a+1, a] + det;
            u[1] = A[a+1, a] * (A[a, a] + A[a+1, a+1] - tr);
            u[2] = A[a+1, a] * A[a+2, a+1];
            //Debug.WriteLine("u = ");
            //PrintVector(u);

            // the implicit Q theorem says that the Q required to bring A' back to Hessenberg form
            // is determined by its first column, so we can go ahead and apply

            // compute the corresponding Householder transform P that would zero that column
            ComputeHouseholderVector(ref u);

            // compute P * A
            // (PA)_rc = A_rc - u_r (u_k A_kc)
            for (int c = a; c < dim; c++) {
                double v = 0.0;
                for (int rp = 0; rp < 3; rp++) {
                    int r = a + rp;
                    v += u[rp] * A[r, c];
                }
                for (int rp = 0; rp < 3; rp++) {
                    int r = a + rp;
                    A[r, c] -= u[rp] * v;
                }
            }
            // compute A * P
            // (AP)_rc = A_rc - (A_rk u_k) u_c
            for (int r = 0; r < dim; r++) {
                double v = 0.0;
                for (int cp = 0; cp < 3; cp++) {
                    int c = a + cp;
                    v += A[r, c] * u[cp];
                }
                for (int cp = 0; cp < 3; cp++) {
                    int c = a + cp;
                    A[r, c] -= v * u[cp];
                }
            }
            // actually, there is a maximum r we needn't go beyond; replace r < dim by this limit
            // i think we go to at most row a+3, i.e. r < a + 4
            // this same logic applies to loop below, too

            // if we are keeping track of transformations, compute Q * P
            if (Q != null) {
                for (int r = 0; r < dim; r++) {
                    double v = 0.0;
                    for (int cp = 0; cp < u.Length; cp++) {
                        int c = a + cp;
                        v += Q[r, c] * u[cp];
                    }
                    for (int cp = 0; cp < u.Length; cp++) {
                        int c = a + cp;
                        Q[r, c] -= v * u[cp];
                    }
                }
            }

            //Write(A);

            // A is now no longer Hessenberg; it has two extra elements in the first column and one extra
            // element in the second. use a series of Householder transforms to fix this

            // the following logic is effectively the same as calling HessenbergReduce(A, Q)
            // however, we reproduce that logic here so we can take advantage of our knowledge that
            // the matrix is already Hessenberg with the exception of three specific elements

            for (int k = a; k < (n - 2); k++) {

                // compute the householder to zero the elements of this column below the sub-diagonal
                u[0] = A[k + 1, k];
                u[1] = A[k + 2, k];
                u[2] = A[k + 3, k];
                //PrintVector(u);
                double x = ComputeHouseholderVector(ref u);

                // compute P * A
                for (int c = k + 1; c < dim; c++) {
                    double v = 0.0;
                    for (int rp = 0; rp < 3; rp++) {
                        int r = k + 1 + rp;
                        v += u[rp] * A[r, c];
                    }
                    for (int rp = 0; rp < 3; rp++) {
                        int r = k + 1 + rp;
                        A[r, c] -= u[rp] * v;
                    }
                }

                // set zeroed elements
                A[k + 1, k] = x;
                for (int r = k + 2; r < dim; r++) {
                    A[r, k] = 0.0;
                }

                // compute A * P
                for (int r = 0; r < dim; r++) {
                    double v = 0.0;
                    for (int cp = 0; cp < 3; cp++) {
                        int c = k + 1 + cp;
                        v += A[r, c] * u[cp];
                        //if (A[r, c] == 0.0) Debug.WriteLine(String.Format("{0},{1} => 0", r, c));
                    }
                    for (int cp = 0; cp < 3; cp++) {
                        int c = k + 1 + cp;
                        A[r, c] -= v * u[cp];
                    }
                }

                // if we are keeping track of transforms, compute Q * P
                if (Q != null) {
                    for (int r = 0; r < dim; r++) {
                        double v = 0.0;
                        for (int cp = 0; cp < u.Length; cp++) {
                            int c = k + 1 + cp;
                            v += Q[r, c] * u[cp];
                        }
                        for (int cp = 0; cp < u.Length; cp++) {
                            int c = k + 1 + cp;
                            Q[r, c] -= v * u[cp];
                        }
                    }
                }

            }

            //Write(A);

            // now only one non-Hessenberg element remains, in the lower-right corner
            // this can be zerod with a 2D Hessenberg transform or a 2D rotation

            double[] up = new double[2];
            up[0] = A[m, n - 2];
            up[1] = A[n, n - 2];
            double xp = ComputeHouseholderVector(ref up);

            // compute P * A
            for (int c = n - 1; c < dim; c++) {
                double v = 0.0;
                for (int rp = 0; rp < 2; rp++) {
                    int r = n - 1 + rp;
                    v += up[rp] * A[r, c];
                }
                for (int rp = 0; rp < 2; rp++) {
                    int r = n - 1 + rp;
                    A[r, c] -= up[rp] * v;
                }
            }

            // set zeroed elements
            A[m, n - 2] = xp;
            A[n, n - 2] = 0.0;

            // compute A * P
            for (int r = 0; r < dim; r++) {
                double v = 0.0;
                for (int cp = 0; cp < 2; cp++) {
                    int c = n - 1 + cp;
                    v += A[r, c] * up[cp];
                }
                for (int cp = 0; cp < 2; cp++) {
                    int c = n - 1 + cp;
                    A[r, c] -= v * up[cp];
                }
            }

            // if we are keeping track of transforms, compute Q*P
            if (Q != null) {
                for (int r = 0; r < dim; r++) {
                    double v = 0.0;
                    for (int cp = 0; cp < up.Length; cp++) {
                        int c = n - 1 + cp;
                        v += Q[r, c] * up[cp];
                    }
                    for (int cp = 0; cp < up.Length; cp++) {
                        int c = n - 1 + cp;
                        Q[r, c] -= v * up[cp];
                    }
                }
            }

            // finished Francis step; A is now once again Hessenberg,
            // hopefully with smaller sub-diagonal elements

            //Write(A);

        }

#endif

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
        /*
        public static Complex[] ExtractEigenvalues (double[] aStore, double[] qStore, int dimension) {

            int count = 0;

            // keep track of extracted eigenvalues
            Complex[] lambdas = new Complex[dimension];

            double sum_old = Double.PositiveInfinity;

            int n = dimension - 1;
            while (n >= 0) {

                //Write(aStore, dimension, dimension);

                // find the lowest decoupled, unreduced block
                int a = n;
                double sum = 0.0;
                while (a > 0) {
                    double f = Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a)) +
                        Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a - 1, a - 1));
                    double e = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a - 1);
                    if ((f + e) == f) {
                        //Console.WriteLine("zero at {0}", a);
                        MatrixAlgorithms.SetEntry(aStore, dimension, dimension, a, a - 1, 0.0);
                        break;
                    }
                    sum += Math.Abs(e);
                    a--;
                }

                // check for maximum numbers of iterations without finding an eigenvalue
                if (count > 32) {
                    Console.WriteLine("max iter");
                    throw new NonconvergenceException();
                }

                // we are working between a and n, and our block is at least 3 wide
                Console.WriteLine("a={0} n={1} sum={2}", a, n, sum);

                // reduce if possible, otherwise do a iteration step

                if (a == n) {
                    // 1 X 1 block isolated
                    Console.WriteLine("one block");
                    lambdas[a] = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a);
                    n -= 1;
                    count = 0;
                    sum_old = Double.PositiveInfinity;
                } else if (a == (n - 1)) {
                    // 2 X 2 block isolated
                    Console.WriteLine("two blocks");
                    double Aaa = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a);
                    double Aba = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a + 1, a);
                    double Aab = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a + 1);
                    double Abb = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a + 1, a + 1);
                    // eigenvalues are given by the quadratic formula
                    double c = (Aaa + Abb) / 2.0;
                    double d = Aaa * Abb - Aba * Aab;
                    double q2 = c * c - d;
                    if (q2 >= 0.0) {
                        double q = Math.Sqrt(q2);
                        if (c >= 0.0) {
                            lambdas[a] = c + q;
                        } else {
                            lambdas[a] = c - q;
                        }
                        // i tried to do this as c + Math.Sign(c) * q, but when c = 0, Math.Sign(c) = 0, not 1 
                        lambdas[a + 1] = d / lambdas[a];
                    } else {
                        double q = Math.Sqrt(-q2);
                        lambdas[a] = new Complex(c, q);
                        lambdas[a + 1] = new Complex(c, -q);
                    }

                    n -= 2;
                    count = 0;
                    sum_old = Double.PositiveInfinity;
                } else {
                    Console.WriteLine("iterate");
                    //if (1.25 * sum < sum_old) {
                        Console.WriteLine("francis count={0}", count);

                        // use the lower-left 2 X 2 matrix to generate an approximate eigenvalue pair

                        int m = n - 1;
                        double Amm = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, m);
                        double Amn = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, n);
                        double Anm = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, m);
                        double Ann = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, n);

                        // the trace of this 2 X 2 matrix is the sum of its eigenvalues, and its determinate is their product

                        double tr = Amm + Ann;
                        double det = Amm * Ann - Amn * Anm;

                        if ((count == 8) || (count == 16) || (count == 24)) {
                            double w = Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, n - 1)) +
                                Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, m - 1));
                            Console.WriteLine("ad hoc w={0}", w);
                            tr = 2.0 * w;
                            det = w * w;
                        }
                        Console.WriteLine("tr={0} det={1}", tr, det);

                        FrancisTwoStep(aStore, qStore, dimension, a, n, tr, det);
                        sum_old = sum;
                    //} else {
                    //    Console.WriteLine("wilkerson");
                    //    WilkersonOneStep(aStore, dimension, a, n);
                    //    sum_old = Double.PositiveInfinity;
                    //}
                    count++;
                }

            }

            return (lambdas);


        }

        public static void FrancisTwoStep (double[] aStore, double[] qStore, int dimension, int a, int n, double sum, double product) {

            //Console.WriteLine("Start {0} {1}", a, n);
            //Write(aStore, dimension, dimension);

            // to apply the implicit Q theorem, get the first column of what (A - mu I)(A - mu* I) would be
            double[] v = GenerateFirstColumn(aStore, dimension, a, n, sum, product);

            int d = n - a + 1;
            Debug.Assert(d >= 3);

            // generate a householder reflection for this column, and apply it to the matrix
            double e;
            GenerateHouseholderReflection(v, 0, 1, 3, out e);

            // determine P * A
            for (int c = a; c < dimension; c++) {
                ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + a, 1, 3);
            }
            // determine A * P
            for (int r = 0; r <= Math.Min(a + 3, n); r++) {
                ApplyHouseholderReflection(v, 0, 1, aStore, dimension * a + r, dimension, 3);
            }
            // determine Q * P
            if (qStore != null) {
                for (int r = 0; r < dimension; r++) {
                    ApplyHouseholderReflection(v, 0, 1, qStore, dimension * a + r, dimension, 3);
                }
            }

            //Write(aStore, dimension, dimension);

            // The matrix now looks like this:
            // XXXXXX    ABCDEF    AAAXXX
            // XXXXXX    ABCDEF    BBBXXX
            // 0XXXXX -> ABCDEF -> CCCXXX
            // 00XXXX    00XXXX    DDDXXX
            // 000XXX    000XXX    000XXX
            // 0000XX    0000XX    0000XX
            // Note there are three non-Hessenberg elements. These constitute a "bulge". Now "chase the bulge" down the matrix, using
            // Householder relfections to zero the non-Hessenberg elements of each column in turn.

            for (int k = a; k <= (n - 3); k++) {

                // determine the required Householder reflection
                v[0] = aStore[dimension * k + (k + 1)];
                v[1] = aStore[dimension * k + (k + 2)];
                v[2] = aStore[dimension * k + (k + 3)];
                GenerateHouseholderReflection(v, 0, 1, 3, out e);

                //v = GenerateHouseholderReflection(aStore, dimension * k + k + 1, 1, 3);

                // determine P * A and A * P
                aStore[dimension * k + (k + 1)] = e;
                aStore[dimension * k + (k + 2)] = 0.0;
                aStore[dimension * k + (k + 3)] = 0.0;
                for (int c = k + 1; c < dimension; c++) {
                    ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + (k + 1), 1, 3);
                }
                for (int r = 0; r < Math.Max(k + 4, dimension); r++) {
                    ApplyHouseholderReflection(v, 0, 1, aStore, dimension * (k + 1) + r, dimension, 3);
                }
                // if tracking eigenvalues, determine Q * P
                if (qStore != null) {
                    for (int r = 0; r < dimension; r++) {
                        ApplyHouseholderReflection(v, 0, 1, qStore, dimension * (k + 1) + r, dimension, 3);
                    }
                }

            }

            //Write(aStore, dimension, dimension);

            // restoring Householder form in the last column requires only a length-2 Householder reflection
            // a Givens rotation would be another possibility here
            int l = n - 2;
            v[0] = aStore[dimension * l + (l + 1)];
            v[1] = aStore[dimension * l + (l + 2)];
            GenerateHouseholderReflection(v, 0, 1, 2, out e);
            aStore[dimension * l + (l + 1)] = e;
            aStore[dimension * l + (l + 2)] = 0.0;
            for (int c = l + 1; c < dimension; c++) {
                ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + (l + 1), 1, 2);
            }
            for (int r = 0; r <= n; r++) {
                ApplyHouseholderReflection(v, 0, 1, aStore, dimension * (l + 1) + r, dimension, 2);
            }
            if (qStore != null) {
                for (int r = 0; r < dimension; r++) {
                    ApplyHouseholderReflection(v, 0, 1, qStore, dimension * (l + 1) + r, dimension, 2);
                }
            }
            //double cos, sin;
            //GenerateGivensRotation(ref aStore[dimension * l + l + 1], ref aStore[dimension * l + l + 2], out cos, out sin);

            //Write(aStore, dimension, dimension);
            //Console.WriteLine("Finish");
        }

        public static double[] GenerateFirstColumn (double[] aStore, int dimension, int a, int n, double sum, double product) {

            // compute the first column of M = (A - \lambda I) (A - \lambda^* I) = A^2 - 2 Re(\lambda) + |\lambda|^2 I
            // where 2 Re(\lambda) = sum of approximate eigenvalues, |\lambda|^2 = product of approximate eigenvalues
            // since A is Hessenberg, this has only three entries

            // first get the relevant A entries, which for the first column are from the leading Hesenberg entries

            int b = a + 1;
            int c = a + 2;
            double Aaa = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a);
            double Aab = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, b);
            double Aba = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, b, a);
            double Abb = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, b, b);
            double Acb = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, c, b);

            // Aba should not be zero, or the problem would have been reduced
            Debug.Assert(Aba != 0.0);

            // use these entries and the sum and product of eigenvalues of get first column

            double[] u = new double[3];
            u[0] = Aaa * (Aaa - sum) + Aba * Aab + product;
            u[1] = Aba * (Aaa + Abb - sum);
            u[2] = Aba * Acb;

            return (u);

        }

        public static double[] GenerateFirstColumn (double[] aStore, int dimension, int a, int n) {

            // use the lower-left 2 X 2 matrix to generate an approximate eigenvalue pair

            int m = n - 1;
            double Amm = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, m);
            double Amn = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, n);
            double Anm = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, m);
            double Ann = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, n);

            // the trace of this 2 X 2 matrix is the sum of its eigenvalues, and its determinate is their product

            double tr = Amm + Ann;
            double det = Amm * Ann - Amn * Anm;
            //Console.WriteLine("tr = {0}, det = {1}", tr, det);

            return (GenerateFirstColumn(aStore, dimension, a, n, tr, det));

        }
         */

#if PAST

        // generate a Givens rotation matrix, such that
        //   (  c  s ) ( a ) = ( r )
        //   ( -s  c ) ( b )   ( 0 )
        // such rotations are used to zero subdiagonal elements

        // given a and b, this method computes c and s, then sets a to r and b to 0

        private static void GenerateGivensRotation (ref double a, ref double b, out double c, out double s) {

            if (b == 0.0) {
                c = 1.0;
                s = 0.0;
            } else {
                double r = MoreMath.Hypot(a, b);
                c = a / r;
                s = b / r;
                a = r;
                b = 0.0;
            }

        }

        private static void ApplyGivensRotation (
            double c, double s,
            double[] store1, int offset1, int stride1,
            double[] store2, int offset2, int stride2,
            int count)
        {
            int p1 = offset1;
            int p2 = offset2;
            for (int i = 0; i < count; i++) {
                double t = c * store1[p1] + s * store2[p2];
                store2[p2] = -s * store1[p1] + c * store2[p2];
                store1[p1] = t;
                p1 += stride1;
                p2 += stride2;
            }


        }

        // 1 0 0
        // 1 1 1
        // 0 1 1
        // is a good example, with eigenvalues 0, 1, and 2 and a simple QR decomposition

        // Start with Hessenberg matrix. Apply G1, G2, G3, etc. from right in order to zero sub-diagonal elements.
        // Then apply G1, G2, G3 from left to complete transform. We don't want to apply G1 from left immediately after
        // applying it from right, because it would create non-Hessenberg elements (they would become zero eventually,
        // but we would have to track them in the meantime because they would interact with other elements). That means
        // we have to keep track of the rotation angles on a stack as we apply them from the right so we have them when
        // we are ready to apply from the left.

        // XXXX    ABCD    XXXX    XXXX    AAXX    XAAX    XXAA
        // XXXX    0BCD    0ABC    0XXX    BBXX    XBBX    XXBB
        // 0XXX -> 0XXX -> 00BC -> 00AB -> 00XX -> 0CCX -> XXCC
        // 00XX    00XX    00XX    000B    000X    000X    XXDD

        // a is the starting index and n is the final index

        public static void UnshiftedQRStep (double[] store, int dimension, int a, int n) {

            // the count of rotations, which is one more than the dimension of the submatrix
            int c = n - a;
            int d = c + 1;

            // create storage for angles
            double[] cos = new double[c];
            double[] sin = new double[c];

            // determine Givens rotations and apply them from the left (i.e. on columns)
            for (int i = 0; i < c; i++) {
                int k = a + i;
                GenerateGivensRotation(ref store[dimension * k + k], ref store[dimension * k + k + 1], out cos[i], out sin[i]);
                //Console.WriteLine("{0} {1}", i, cos[i]);
                ApplyGivensRotation(cos[i], sin[i], store, dimension * (k + 1) + k, dimension, store, dimension * (k + 1) + (k + 1), dimension, dimension - (k + 1));
            }

            // now apply them from the right (i.e. on rows)
            for (int i = 0; i < c; i++) {
                int k = a + i;
                ApplyGivensRotation(cos[i], sin[i], store, dimension * k, 1, store, dimension * (k + 1), 1, k + 2);
            }
        }

        public static void UnshiftedQRStep (double[,] H, int n) {

            double[] c = new double[n];
            double[] s = new double[n];
            for (int k = 0; k < n - 1; k++) {
                GenerateGivensRotation(ref H[k, k], ref H[k + 1, k], out c[k], out s[k]);
                Console.WriteLine("{0} {1}", k, c[k]);
                for (int i = k + 1; i < n; i++) {
                    double t = c[k] * H[k, i] + s[k] * H[k + 1, i];
                    H[k + 1, i] = -s[k] * H[k, i] + c[k] * H[k + 1, i];
                    H[k, i] = t;
                }
            }
            for (int k = 0; k < n - 1; k++) {
                for (int i = 0; i <= k + 1; i++) {
                    double t = c[k] * H[i, k] + s[k] * H[i, k + 1];
                    H[i, k + 1] = -s[k] * H[i, k] + c[k] * H[i, k + 1];
                    H[i, k] = t;
                }
            }
        }

        public static void ShiftedQRStep (double[,] H, double mu, int n) {
            for (int k = 0; k < n; k++) {
                H[k, k] -= mu;
            }
            UnshiftedQRStep(H, n);
            for (int k = 0; k < n; k++) {
                H[k, k] += mu;
            }
        }

        public static void ShiftedQRStep (double[] store, int dimension, int a, int n, double mu) {
            for (int k = a; k <= n; k++) {
                store[dimension * k + k] -= mu;
            }
            UnshiftedQRStep(store, dimension, a, n);
            for (int k = a; k <= n; k++) {
                store[dimension * k + k] += mu;
            }
        }


        public static void GolubKahanStep (double[] a, double [] b, int dimension, int f, int n) {

            // t is a temporary storage register for copying, c contains the current non-diagonal element
            double t, c;
            
            // compute the shift
            double mu = 0.0;

            // compute the Givens rotation that would zero the first column of B^T B, which has just two entries
            double p = a[f] * a[f] - mu;
            double q = a[f] * b[f];
            double cos, sin;
            GenerateGivensRotation(ref p, ref q, out cos, out sin);

            // apply it to B (we don't call ApplyGivensRotation because 
            t = cos * a[f] - sin * b[f]; b[f] = sin * a[f] + cos * b[f]; a[f] = t;
            c = -sin * a[f + 1];
            a[f + 1] = cos * a[f + 1];

            // now "chase the bulge" down the bidiagonal; it alternates sides



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
#if PAST

        // given a matrix column, returns the properly normalized vector defining the Householder transform that zeros it
        // if v is the returned vector, T = I - v v^T is the transform, applied to A as T A T
        // the computation requires 2 N flops
        private static double ComputeHouseholderVector (ref double[] v) {

            // compute column norm
            double x = 0.0;
            for (int i = 0; i < v.Length; i++) {
                x += v[i] * v[i];
            }
            x = Math.Sqrt(x);

            // create the un-normalized Householder vector and compute the normalization factor
            double norm;
            if (v[0] < 0.0) {
                norm = Math.Sqrt(x * (x - v[0]));
                v[0] -= x;
            } else {
                norm = Math.Sqrt(x * (x + v[0]));
                v[0] += x;
                x = -x;
            }

            // normalize the transform
            if (norm != 0.0) {
                for (int i = 0; i < v.Length; i++) {
                    v[i] = v[i] / norm;
                }
            }

            // return the value of the first component of P v that remains after the transform 
            return (x);
        }

        // Start of the reduction to Hessenberg form: Given the matrix A partitioned
        //   ( a_11 a_12 )
        //   ( a_21 A_22 )
        // where a_11 is a number, a_12 is a row vector, a_21 is a column vector, and A_22 is a matrix
        // and a Householder transform P a_21 = x e1 then PAP is
        //   ( a_11 a_12^T P )
        //   ( x e1 P A_22 P )
        // Continue the reduction to Hessenberg form: Given a part-Hessenberg matrix A partitioned
        //   ( A_11 a_12 A_13 )
        //   (  0   a_22 a_23 )
        //   (  0   a_32 A_33 )
        // where A_11 is a Hessenberg matrix, a_12 is a column vector, A_13 is a matrix, a_22 is a number,
        // a_23 is a row vector, a_32 is a column vector, and A_33 is a matrix
        // and a Householder transform P a_32 = x e1 than PAP is
        //   ( A_11 a_12 A_13 P   )
        //   (  0   a_22 a_23 P   )
        //   (  0   x e1 P A_33 P )
        // Doing Householder transform: Given P = 1 - u u^T,
        //   A P = A - (A u) u^T = A - v u^T where v = A u
        // This is a O(N^2) operation.

        private double ComputeHouseHolderTransform2 (double[] v) {

            //PrintVector(v);

            // compute column norm
            double y = 0.0;
            for (int i = 1; i < v.Length; i++) {
                y += MoreMath.Pow2(v[i]);
            }
            double x = Math.Sqrt(MoreMath.Pow2(v[0]) + y);
            Debug.WriteLine(String.Format("x={0}", x));

            // create the un-normalized Householder vector and compute the normalization factor
            if (v[0] < 0.0) x = -x;
            v[0] -= x;
            double norm = Math.Sqrt(MoreMath.Pow2(v[0]) + y);
            //Debug.WriteLine(String.Format("x={0}", x));

            //PrintVector(v);
            //Debug.WriteLine(String.Format("norm={0}", norm));

            // normalize the transform
            if (norm != 0.0) {
                for (int i = 0; i < v.Length; i++) {
                    v[i] = v[i] / norm;
                }
            }
            //PrintVector(v);

            return (x);
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
        internal static SquareMatrix Multiply (double x, SquareMatrix M) {
            SquareMatrix N = new SquareMatrix(M.Dimension);
            for (int r = 0; r < N.Dimension; r++) {
                for (int c = 0; c < N.Dimension; c++) {
                    N[r, c] = x * M[r, c];
                }
            }
            return (N);
        }

        /// <summary>
        /// Computes the product of a real number and a square matrix.
        /// </summary>
        /// <param name="x">The real number.</param>
        /// <param name="M">The matrix.</param>
        /// <returns>The product of <paramref name="x"/> and <paramref name="M"/>.</returns>
        public static SquareMatrix operator * (double x, SquareMatrix M) {
            return (Multiply(x, M));
        }

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


    internal static class SquareMatrixAlgorithms {

        public static double[] CreateUnitMatrix (int dimension) {
            double[] store = MatrixAlgorithms.AllocateStorage(dimension, dimension);
            for (int i = 0; i < dimension; i++) {
                store[dimension * i + i] = 1.0;
            }
            return (store);
        }

        // on input:
        // store contains the matrix in column-major order (must have the length dimension^2)
        // permutation contains the row permutation (typically 0, 1, 2, ..., dimension - 1; must have the length dimension^2)
        // parity contains the parity of the row permutation (typically 1; must be 1 or -1)
        // dimension contains the dimension of the matrix (must be non-negative)
        // on output:
        // store is replaced by the L and U matrices, L in the lower-left triangle (with 1s along the diagonal), U in the upper-right triangle
        // permutation is replaced by the row permutation, and parity be the parity of that permutation

        // A = PLU

        public static void LUDecompose (double[] store, int[] permutation, ref int parity, int dimension) {

            for (int d = 0; d < dimension; d++) {

                int pivotRow = -1;
                double pivotValue = 0.0;

                for (int r = d; r < dimension; r++) {
                    int a0 = dimension * d + r;
                    double t = store[a0] - Blas1.dDot(store, r, dimension, store, dimension * d, 1, d);
                    store[a0] = t;
                    if (Math.Abs(t) > Math.Abs(pivotValue)) {
                        pivotRow = r;
                        pivotValue = t;
                    }
                }
                //Console.WriteLine("pivot row = {0}, pivot value = {1}", pivotRow, pivotValue);
                //pivotValue = aStore[dimension * d + d];

                if (pivotValue == 0.0) throw new DivideByZeroException();

                if (pivotRow != d) {
                    // switch rows
                    Blas1.dSwap(store, d, dimension, store, pivotRow, dimension, dimension);
                    int t = permutation[pivotRow];
                    permutation[pivotRow] = permutation[d];
                    permutation[d] = t;
                    parity = -parity;
                }

                Blas1.dScal(1.0 / pivotValue, store, dimension * d + d + 1, 1, dimension - d - 1);

                for (int c = d + 1; c < dimension; c++) {
                    double t = Blas1.dDot(store, d, dimension, store, dimension * c, 1, d);
                    store[dimension * c + d] -= t;
                }

            }

        }

        // inverts the matrix in place
        // the in-place-ness makes this a bit confusing

        public static void GaussJordanInvert (double[] store, int dimension) {

            // keep track of row exchanges
            int[] ps = new int[dimension];

            // iterate over dimensions
            for (int k = 0; k < dimension; k++) {

                // look for a pivot in the kth column on any lower row
                int p = k;
                double q = MatrixAlgorithms.GetEntry(store, dimension, dimension, k, k);
                
                for (int r = k + 1; r < dimension; r++) {
                    double s = MatrixAlgorithms.GetEntry(store, dimension, dimension, r, k);
                    if (Math.Abs(s) > Math.Abs(q)) {
                        p = r;
                        q = s;
                    }
                }
                ps[k] = p;

                // if no non-zero pivot is found, the matrix is singular and cannot be inverted
                if (q == 0.0) throw new DivideByZeroException();

                // if the best pivot was on a lower row, swap it into the kth row
                if (p != k) {
                    Blas1.dSwap(store, k, dimension, store, p, dimension, dimension);
                }

                // divide the pivot row by the pivot element, so the diagonal element becomes unity
                MatrixAlgorithms.SetEntry(store, dimension, dimension, k, k, 1.0);
                Blas1.dScal(1.0 / q, store, k, dimension, dimension);

                // add factors to the pivot row to zero all off-diagonal elements in the kth column
                for (int r = 0; r < dimension; r++) {
                    if (r == k) continue;
                    double a = MatrixAlgorithms.GetEntry(store, dimension, dimension, r, k);
                    MatrixAlgorithms.SetEntry(store, dimension, dimension, r, k, 0.0);
                    Blas1.dAxpy(-a, store, k, dimension, store, r, dimension, dimension);
                }

            }

            // unscramble exchanges
            for (int k = dimension - 1; k >= 0; k--) {
                int p = ps[k];
                if (p != k) Blas1.dSwap(store, dimension * p, 1, store, dimension * k, 1, dimension);
            }

        }

        public static void SolveLowerLeftTriangular (double[] tStore, double[] x, int xOffset, int dimension) {
            for (int r = 0; r < dimension; r++) {
                double t = Blas1.dDot(tStore, r, dimension, x, xOffset, 1, r);
                x[xOffset + r] -= t;
            }
        }

        public static void SolveUpperRightTriangular (double[] tStore, double[] x, int xOffset, int dimension) {
            for (int r = dimension - 1; r >= 0; r--) {
                int i0 = dimension * r + r;
                double t = Blas1.dDot(tStore, i0 + dimension, dimension, x, xOffset + r + 1, 1, dimension - r - 1);
                x[xOffset + r] = (x[xOffset + r] - t) / tStore[i0];
            }
        }

        // The reduction to Hessenberg form via a similiarity transform proceeds as follows. A Householder reflection can
        // zero the non-Householder elements in the first column. Since this is a similiarity transform, we have to apply
        // it from the other side as well. The letters indicate which elements are mixed by each transform.

        // XXXXX    XXXXX    XFFFF
        // XXXXX    ABCDE    AGGGG    
        // XXXXX -> 0BCDE -> 0HHHH
        // XXXXX    0BCDE    0IIII
        // XXXXX    0BCDE    0JJJJ

        // Note that if we had tried to zero all the subdiagonal elements, rather than just the (sub+1)-diagonal elements, then the
        // applicaion of the Householder reflection from the other side would have messed with our zeros. But because we only tried
        // for Hessenberg, not upper-right-triangular, form, our zeros are safe. We then do it again for the next column.

        // XXXXX    XXXXX    XXEEE
        // XXXXX    XXXXX    XXFFF
        // 0XXXX -> 0ABCD -> 0AGGG
        // 0XXXX    00BCD    00HHH
        // 0XXXX    00BCD    00III

        // And so on until we are done.

        public static void ReduceToHessenberg (double[] aStore, double[] vStore, int dim) {

            int dm2 = dim - 2;
            for (int k = 0; k < dm2; k++) {

                // determine a Householder reflection to zero the (sub+1)-diagonal elements of the current column
                int offset = dim * k + (k + 1);
                int length = dim - (k + 1);
                double a;
                GenerateHouseholderReflection(aStore, offset, 1, length, out a);

                // determine P * A
                for (int c = k + 1; c < dim; c++) {
                    ApplyHouseholderReflection(aStore, offset, 1, aStore, dim * c + (k + 1), 1, length);
                }
                // determine A * P
                for (int r = 0; r < dim; r++) {
                    ApplyHouseholderReflection(aStore, offset, 1, aStore, dim * (k + 1) + r, dim, length);

                }

                // if we are keeping track of the transformation, determine V * P
                if (vStore != null) {
                    for (int r = 0; r < dim; r++) {
                        ApplyHouseholderReflection(aStore, offset, 1, vStore, dim * (k + 1) + r, dim, length);
                    }
                }

                // We are done with our Householder vector, so we can overwrite the space we were using for it with the values produced by
                // the reflection.
                aStore[offset] = a;
                for (int i = 1; i < length; i++) {
                    aStore[offset + i] = 0.0;
                }
                //aStore[dim * k + (k + 1)] = a;
                //for (int i = k + 2; i < dim; i++) {
                //    MatrixAlgorithms.SetEntry(aStore, dim, dim, i, k, 0.0);
                //}
                // Not having done this earlier is not a problem because the values in that column do not come into play when computing
                // the effects of the transform on the other entries, i.e. it doesn't "mix" with the other entries

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

        public static void ApplyHouseholderReflection (
            double[] uStore, int uOffset, int uStride,
            double[] yStore, int yOffset, int yStride,
            int count
        ) {
            // (1 - b v v^T) y = y - (b v^T y) v
            double s = Blas1.dDot(uStore, uOffset, uStride, yStore, yOffset, yStride, count);
            Blas1.dAxpy(-s, uStore, uOffset, uStride, yStore, yOffset, yStride, count);
        }

        // EIGENVALUE ALGORITHMS 

        public static Complex[] ExtractEigenvalues (double[] aStore, double[] qStore, int dimension) {

            int count = 0;

            // keep track of extracted eigenvalues
            Complex[] lambdas = new Complex[dimension];

            double sum_old = Double.PositiveInfinity;

            int n = dimension - 1;
            while (n >= 0) {

                //Write(aStore, dimension, dimension);

                // find the lowest decoupled, unreduced block
                int a = n;
                double sum = 0.0;
                while (a > 0) {
                    double f = Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a)) +
                        Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a - 1, a - 1));
                    double e = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a - 1);
                    if ((f + e) == f) {
                        //Console.WriteLine("zero at {0}", a);
                        MatrixAlgorithms.SetEntry(aStore, dimension, dimension, a, a - 1, 0.0);
                        break;
                    }
                    sum += Math.Abs(e);
                    a--;
                }

                // check for maximum numbers of iterations without finding an eigenvalue
                if (count > 32) {
                    Console.WriteLine("max iter");
                    throw new NonconvergenceException();
                }

                // we are working between a and n, and our block is at least 3 wide
                Console.WriteLine("a={0} n={1} sum={2}", a, n, sum);

                // reduce if possible, otherwise do a iteration step

                if (a == n) {
                    // 1 X 1 block isolated
                    Console.WriteLine("one block");
                    lambdas[a] = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a);
                    n -= 1;
                    count = 0;
                    sum_old = Double.PositiveInfinity;
                } else if (a == (n - 1)) {
                    // 2 X 2 block isolated
                    Console.WriteLine("two blocks");
                    double Aaa = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a);
                    double Aba = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a + 1, a);
                    double Aab = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a + 1);
                    double Abb = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a + 1, a + 1);
                    // eigenvalues are given by the quadratic formula
                    double c = (Aaa + Abb) / 2.0;
                    double d = Aaa * Abb - Aba * Aab;
                    double q2 = c * c - d;
                    if (q2 >= 0.0) {
                        double q = Math.Sqrt(q2);
                        if (c >= 0.0) {
                            lambdas[a] = c + q;
                        } else {
                            lambdas[a] = c - q;
                        }
                        // i tried to do this as c + Math.Sign(c) * q, but when c = 0, Math.Sign(c) = 0, not 1 
                        lambdas[a + 1] = d / lambdas[a];
                    } else {
                        double q = Math.Sqrt(-q2);
                        lambdas[a] = new Complex(c, q);
                        lambdas[a + 1] = new Complex(c, -q);
                    }

                    n -= 2;
                    count = 0;
                    sum_old = Double.PositiveInfinity;
                } else {
                    Console.WriteLine("iterate");
                    //if (1.25 * sum < sum_old) {
                    Console.WriteLine("francis count={0}", count);

                    // use the lower-left 2 X 2 matrix to generate an approximate eigenvalue pair

                    int m = n - 1;
                    double Amm = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, m);
                    double Amn = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, n);
                    double Anm = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, m);
                    double Ann = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, n);

                    // the trace of this 2 X 2 matrix is the sum of its eigenvalues, and its determinate is their product

                    double tr = Amm + Ann;
                    double det = Amm * Ann - Amn * Anm;

                    if ((count == 8) || (count == 16) || (count == 24)) {
                        double w = Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, n - 1)) +
                            Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, m - 1));
                        Console.WriteLine("ad hoc w={0}", w);
                        tr = 2.0 * w;
                        det = w * w;
                    }
                    Console.WriteLine("tr={0} det={1}", tr, det);

                    FrancisTwoStep(aStore, qStore, dimension, a, n, tr, det);
                    sum_old = sum;
                    //} else {
                    //    Console.WriteLine("wilkerson");
                    //    WilkersonOneStep(aStore, dimension, a, n);
                    //    sum_old = Double.PositiveInfinity;
                    //}
                    count++;
                }

            }

            return (lambdas);


        }

        private static void FrancisTwoStep (double[] aStore, double[] qStore, int dimension, int a, int n, double sum, double product) {

            //Console.WriteLine("Start {0} {1}", a, n);
            //Write(aStore, dimension, dimension);

            // to apply the implicit Q theorem, get the first column of what (A - mu I)(A - mu* I) would be
            double[] v = GenerateFirstColumn(aStore, dimension, a, n, sum, product);

            int d = n - a + 1;
            Debug.Assert(d >= 3);

            // generate a householder reflection for this column, and apply it to the matrix
            double e;
            GenerateHouseholderReflection(v, 0, 1, 3, out e);

            // determine P * A
            for (int c = a; c < dimension; c++) {
                ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + a, 1, 3);
            }
            // determine A * P
            for (int r = 0; r <= Math.Min(a + 3, n); r++) {
                ApplyHouseholderReflection(v, 0, 1, aStore, dimension * a + r, dimension, 3);
            }
            // determine Q * P
            if (qStore != null) {
                for (int r = 0; r < dimension; r++) {
                    ApplyHouseholderReflection(v, 0, 1, qStore, dimension * a + r, dimension, 3);
                }
            }

            //Write(aStore, dimension, dimension);

            // The matrix now looks like this:
            // XXXXXX    ABCDEF    AAAXXX
            // XXXXXX    ABCDEF    BBBXXX
            // 0XXXXX -> ABCDEF -> CCCXXX
            // 00XXXX    00XXXX    DDDXXX
            // 000XXX    000XXX    000XXX
            // 0000XX    0000XX    0000XX
            // Note there are three non-Hessenberg elements. These constitute a "bulge". Now "chase the bulge" down the matrix, using
            // Householder relfections to zero the non-Hessenberg elements of each column in turn.

            for (int k = a; k <= (n - 3); k++) {

                // determine the required Householder reflection
                v[0] = aStore[dimension * k + (k + 1)];
                v[1] = aStore[dimension * k + (k + 2)];
                v[2] = aStore[dimension * k + (k + 3)];
                GenerateHouseholderReflection(v, 0, 1, 3, out e);

                //v = GenerateHouseholderReflection(aStore, dimension * k + k + 1, 1, 3);

                // determine P * A and A * P
                aStore[dimension * k + (k + 1)] = e;
                aStore[dimension * k + (k + 2)] = 0.0;
                aStore[dimension * k + (k + 3)] = 0.0;
                for (int c = k + 1; c < dimension; c++) {
                    ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + (k + 1), 1, 3);
                }
                for (int r = 0; r < Math.Max(k + 4, dimension); r++) {
                    ApplyHouseholderReflection(v, 0, 1, aStore, dimension * (k + 1) + r, dimension, 3);
                }
                // if tracking eigenvalues, determine Q * P
                if (qStore != null) {
                    for (int r = 0; r < dimension; r++) {
                        ApplyHouseholderReflection(v, 0, 1, qStore, dimension * (k + 1) + r, dimension, 3);
                    }
                }

            }

            //Write(aStore, dimension, dimension);

            // restoring Householder form in the last column requires only a length-2 Householder reflection
            // a Givens rotation would be another possibility here
            int l = n - 2;
            v[0] = aStore[dimension * l + (l + 1)];
            v[1] = aStore[dimension * l + (l + 2)];
            GenerateHouseholderReflection(v, 0, 1, 2, out e);
            aStore[dimension * l + (l + 1)] = e;
            aStore[dimension * l + (l + 2)] = 0.0;
            for (int c = l + 1; c < dimension; c++) {
                ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + (l + 1), 1, 2);
            }
            for (int r = 0; r <= n; r++) {
                ApplyHouseholderReflection(v, 0, 1, aStore, dimension * (l + 1) + r, dimension, 2);
            }
            if (qStore != null) {
                for (int r = 0; r < dimension; r++) {
                    ApplyHouseholderReflection(v, 0, 1, qStore, dimension * (l + 1) + r, dimension, 2);
                }
            }
            //double cos, sin;
            //GenerateGivensRotation(ref aStore[dimension * l + l + 1], ref aStore[dimension * l + l + 2], out cos, out sin);

            //Write(aStore, dimension, dimension);
            //Console.WriteLine("Finish");
        }

        private static double[] GenerateFirstColumn (double[] aStore, int dimension, int a, int n, double sum, double product) {

            // compute the first column of M = (A - \lambda I) (A - \lambda^* I) = A^2 - 2 Re(\lambda) + |\lambda|^2 I
            // where 2 Re(\lambda) = sum of approximate eigenvalues, |\lambda|^2 = product of approximate eigenvalues
            // since A is Hessenberg, this has only three entries

            // first get the relevant A entries, which for the first column are from the leading Hesenberg entries

            int b = a + 1;
            int c = a + 2;
            double Aaa = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a);
            double Aab = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, b);
            double Aba = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, b, a);
            double Abb = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, b, b);
            double Acb = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, c, b);

            // Aba should not be zero, or the problem would have been reduced
            Debug.Assert(Aba != 0.0);

            // use these entries and the sum and product of eigenvalues of get first column

            double[] u = new double[3];
            u[0] = Aaa * (Aaa - sum) + Aba * Aab + product;
            u[1] = Aba * (Aaa + Abb - sum);
            u[2] = Aba * Acb;

            return (u);

        }

        // SVD

        // Bidiagonalization uses seperate left and right Householder reflections to bring a matrix into a form of bandwidth two.

        // XXXXX    ABCDE    XA000    XX000    XX000
        // XXXXX    0BCDE    0BBBB    0ABCD    0XA00
        // XXXXX    0BCDE    0CCCC    00BCD    00BBB
        // XXXXX -> 0BCDE -> 0DDDD -> 00BCD -> 00CCC -> etc
        // XXXXX    0BCDE    0EEEE    00BCD    00DDD
        // XXXXX    0BCDE    0FFFF    00BCD    00EEE
        // XXXXX    0BCDE    0GGGG    00BCD    00FFF


        // The end result is B = U A V, where U and V are different (so this isn't a similiarity transform) orthogonal matrices.

        // Note that we can't get fully diagonal with this approach, because if the right transform tried to zero the first superdiagonal
        // element, it would interfere with the zeros already created.

        // This method gradually replaces the matrix elements with the elements of the Householder reflections used.
        // The resulting diagonal and first superdiagonal elements are returned in d and e.

        public static void Bidiagonlize (double[] store, int rows, int cols) {

            Debug.Assert(rows >= cols);

            for (int k = 0; k < cols; k++) {

                // generate a householder reflection which, when applied from the left, will zero sub-diagonal elements in the kth column
                // use those elements to store the reflection vector; store the resulting diagonal element seperately
                double a;
                GenerateHouseholderReflection(store, rows * k + k, 1, rows - k, out a);

                // apply that reflection to all the columns; only subsequent ones are affected since preceeding ones
                // contain only zeros in the subdiagonal rows
                for (int c = k + 1; c < cols; c++) {
                    ApplyHouseholderReflection(store, rows * k + k, 1, store, rows * c + k, 1, rows - k);
                }
                store[rows * k + k] = a;

                if ((k + 1) < cols) {
                    // generate a Householder reflection which, when applied from the right, will zero (super+1)-diagonal elements in the kth row
                    // again, store the elements of the reflection vector in the zeroed matrix elements and store the resulting superdiagonal element seperately
                    double b;
                    GenerateHouseholderReflection(store, rows * (k + 1) + k, rows, cols - (k + 1), out b);
                    // apply the reflection to all the rows; only subsequent ones are affected since the preceeding ones contain only zeros in the
                    // affected columns; the already-zeroed column elements are not distrubed because those columns are not affected
                    // this restriction is why we cannot fully diagonalize using this transform; if we tried to zero all the super-diagonal
                    for (int r = k + 1; r < rows; r++) {
                        ApplyHouseholderReflection(store, rows * (k + 1) + k, rows, store, rows * (k + 1) + r, rows, cols - (k + 1));
                    }
                    store[rows * (k + 1) + k] = b;
                }


            }



        }

    }

    /*
     * Some test Hessenberg and eigenvalue problems:
     * 
     * (  1  -4   0   3 )
     * ( -1   3   0  -2 )
     * (  3  -7  -2   6 )
     * (  0   4   0  -2 )
     * 
     * has eigenvalues (2,-2,1,-1)
     */

}
