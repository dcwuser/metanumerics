using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

    /// <summary>
    /// Represents a symmetric matrix.
    /// </summary>
    public sealed class SymmetricMatrix : ISymmetricMatrix {

        private int dimension;
        private double[][] values;

        // storage is in lower triangular form, i.e. values[r][c] with c <= r

        /// <summary>
        /// Instantiates a new symmetric matrix.
        /// </summary>
        /// <param name="dimension">The dimension of the matrix, which must be positive.</param>
        public SymmetricMatrix (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException("dimension");
            this.dimension = dimension;
            values = new double[dimension][];
            for (int r = 0; r < dimension; r++) {
                values[r] = new double[r + 1];
            }
        }

        /// <summary>
        /// Gets the dimension of the matrix.
        /// </summary>
        public int Dimension {
            get {
                return (dimension);
            }
        }

        int IMatrix.RowCount {
            get {
                return (dimension);
            }
        }

        int IMatrix.ColumnCount {
            get {
                return (dimension);
            }
        }

        /// <summary>
        /// Gets or sets an element of the matrix.
        /// </summary>
        /// <param name="r">The (zero-based) row number.</param>
        /// <param name="c">The (zero-based) column number.</param>
        /// <returns>The value of the specified matrix entry M<sub>r c</sub>.</returns>
        /// <remarks>
        /// <para>The set operation preserves the symmetry of the matrix; when entry M<sub>r c</sub> is changed, entry
        /// M<sub>c r</sub> is updated automatically.</para>
        /// </remarks>
        public double this[int r, int c] {
            get {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
                if (c < r) {
                    return (values[r][c]);
                } else {
                    return (values[c][r]);
                }
            }
            set {
                if ((r < 0) || (r >= dimension)) throw new ArgumentOutOfRangeException("r");
                if ((c < 0) || (c >= dimension)) throw new ArgumentOutOfRangeException("c");
                if (c < r) {
                    values[r][c] = value;
                } else {
                    values[c][r] = value;
                }
            }
        }

        private double GetEntry (int r, int c) {
            Debug.Assert(c <= r);
            return (values[r][c]);
        }

        private void SetEntry (int r, int c, double value) {
            Debug.Assert(c <= r);
            values[r][c] = value;
        }

        /// <summary>
        /// Returns an independent copy of the matrix.
        /// </summary>
        /// <returns></returns>
        public SymmetricMatrix Clone () {
            SymmetricMatrix clone = new SymmetricMatrix(Dimension);
            for (int r = 0; r < Dimension; r++) {
                for (int c = 0; c <= r; c++) {
                    clone[r,c] = this[r,c];
                }
            }
            return (clone);
        }

        IMatrix IMatrix.Clone () {
            return (Clone());
        }

        IMatrix IMatrix.Transpose () {
            return (Clone());
        }

        /// <summary>
        /// Computes the trace of the matrix.
        /// </summary>
        /// <returns>The trace of the matrix tr(M).</returns>
        public double Trace () {
            double tr = 0.0;
            for (int i = 0; i < Dimension; i++) {
                tr += this[i,i];
            }
            return (tr);
        }

        /// <summary>
        /// Computes the inverse of the matrix.
        /// </summary>
        /// <returns>The matrix inverse M<sup>-1</sup>.</returns>
        public SymmetricMatrix Inverse () {

            // this implementation re-uses the generic square matrix algorithm
            // replace it with one based on LDL decomposition, which exploits the symmetry

            SquareMatrix M = new SquareMatrix(dimension);
            for (int r = 0; r < dimension; r++) {
                for (int c = 0; c < dimension; c++) {
                    M[r, c] = this[r, c];
                }
            }

            SquareMatrix MI = M.Inverse();

            SymmetricMatrix I = new SymmetricMatrix(dimension);
            for (int r = 0; r < dimension; r++) {
                for (int c = 0; c <= r; c++) {
                    I[r, c] = MI[r, c];
                }
            }

            return (I);
        }

        /// <summary>
        /// Returns the Cholesky decomposition of the matrix.
        /// </summary>
        /// <returns>The Cholesky decomposition of the matrix, or null is the matrix is not positive definite.</returns>
        /// <remarks>
        /// <para>A Cholesky decomposition is a special decomposition that is possible only for positive definite matrices.
        /// (A positive definite matrix has x<sup>T</sup>Mx > 0 for any vector x. Equivilently, M is positive definite if
        /// all its eigenvalues are positive.)</para>
        /// <para>The Cholesky decomposition represents M = C C<sup>T</sup>, where C is lower-left triangular (and thus C<sup>T</sup>
        /// is upper-right triangular. It is basically an LU decomposition where the L and U factors are related by transposition.
        /// Since the M is produced by multiplying C "by itself", the matrix C is sometimes call the "square root" of M.</para>
        /// <para>Cholesky decomposition is an O(N<sup>3</sup>) operation. It is about a factor of two faster than LU decomposition,
        /// so it is a faster way to obtain inverses, determinates, etc. if you know that M is positive definite.</para>
        /// <para>The fastest way to test whether your matrix is positive definite is attempt a Cholesky decomposition. If this
        /// method returns null, M is not positive definite.</para>
        /// </remarks>
        public CholeskyDecomposition CholeskyDecomposition () {
            SymmetricMatrix LU = new SymmetricMatrix(dimension);
            for (int i = 0; i < dimension; i++) {
                double p = this[i, i];
                for (int j = 0; j < i; j++) {
                    p -= Math.Pow(LU[i, j], 2);
                }

                // If a pivot is negative, the matrix is not positive definite
                if (p <= 0) return (null);
                LU[i, i] = Math.Sqrt(p);

                for (int j = i + 1; j < dimension; j++) {
                    LU[j, i] = this[i, j];
                    for (int k = 0; k < i; k++) {
                        LU[j, i] -= LU[i, k] * LU[j, k];
                    }
                    LU[j, i] = LU[j, i] / LU[i, i];
                }
            }

            CholeskyDecomposition CD = new CholeskyDecomposition(LU);
            return (CD);
        }


#if FUTURE

        public static AasenDecomposition LTLDecompose3 (SymmetricMatrix M) {

            // Aasen's method, now with pivoting

            int n = M.Dimension;

            double[] a = new double[n];
            double[] b = new double[n - 1];
            SquareMatrix L = new SquareMatrix(n);
            for (int i = 0; i < n; i++) L[i, i] = 1.0;

            // working space for d'th column of H = T L^T
            double[] h = new double[n];

            // first row
            a[0] = M[0, 0];
            if (n > 1) b[0] = M[1, 0];
            for (int i = 2; i < n; i++) L[i, 1] = M[i, 0] / b[0];

            PrintLTLMatrices(h, a, b, L);

            // second row
            if (n > 1) {
                a[1] = M[1, 1];
                if (n > 2) b[1] = M[2, 1] - L[2, 1] * a[1];
                for (int i = 3; i < n; i++) L[i, 2] = (M[i, 1] - L[i, 1] * a[1]) / b[1];
            }

            PrintLTLMatrices(h, a, b, L);

            for (int d = 0; d < n; d++) {

                Console.WriteLine("d = {0}", d);

                // compute h (d'th row of T L^T)
                if (d == 0) {
                    h[0] = M[0, 0];
                } else if (d == 1) {
                    h[0] = b[0];
                    h[1] = M[1, 1];
                } else {
                    h[0] = b[0] * L[d, 1];
                    h[1] = b[1] * L[d, 2] + a[1] * L[d, 1];
                    h[d] = M[d, d] - L[d, 1] * h[1];
                    for (int i = 2; i < d; i++) {
                        h[i] = b[i] * L[d, i + 1] + a[i] * L[d, i] + b[i - 1] * L[d, i - 1];
                        h[d] -= L[d, i] * h[i];
                    }
                }

                // compute alpha (d'th diagonal element of T)
                if ((d == 0) || (d == 1)) {
                    a[d] = h[d];
                } else {
                    a[d] = h[d] - b[d - 1] * L[d, d - 1];
                }

                Console.WriteLine("before pivot");
                PrintMatrix(M);
                PrintMatrix(L);

                // find the pivot
                if (d < (n - 1)) {
                    int p = d + 1;
                    double q = M[p, d];
                    for (int i = d + 2; i < n; i++) {
                        if (Math.Abs(M[i, d]) > Math.Abs(q)) {
                            p = i;
                            q = M[i, d];
                        }
                    }

                    Console.WriteLine("pivot = {0}", p);

                    // symmetricly permute the pivot element to M[d+1,d]
                    if (p != d + 1) {

                        // symmetricly permute the pivot element to M[d+1, d]
                        // we have to be a bit careful here, because some permutations will be done
                        // automatically by our SymmetricMatrix class due to symmetry
                        for (int i = 0; i < n; i++) {
                            if ((i == p) || (i == d + 1)) continue;
                            double t = M[d + 1, i];
                            M[d + 1, i] = M[p, i];
                            M[p, i] = t;
                        }
                        double tt = M[d + 1, d + 1];
                        M[d + 1, d + 1] = M[p, p];
                        M[p, p] = tt;

                        // also reorder the affected previously computed elements of L 
                        for (int i = 1; i <= d; i++) {
                            double t = L[d + 1, i];
                            L[d + 1, i] = L[p, i];
                            L[p, i] = t;
                        }

                        Console.WriteLine("after pivot");
                        PrintMatrix(M);
                        PrintMatrix(L);

                    }

                }

                // compute beta (d'th subdiagonal element of T)
                if (d < (n - 1)) {
                    b[d] = M[d + 1, d];
                    for (int i = 0; i <= d; i++) {
                        Console.WriteLine("n={0} d={1} i={2}", n, d, i);
                        b[d] -= L[d + 1, i] * h[i];
                    }
                }

                // compute (d+1)'th column of L
                for (int i = d + 2; i < n; i++) {
                    L[i, d + 1] = M[i, d];
                    for (int j = 0; j <= d; j++) L[i, d + 1] -= L[i, j] * h[j];
                    L[i, d + 1] = L[i, d + 1] / b[d];
                }

                PrintLTLMatrices(h, a, b, L);

            }


            Console.WriteLine("Reconstruct");
            SymmetricMatrix T = new SymmetricMatrix(n);
            for (int i = 0; i < n; i++) {
                T[i, i] = a[i];
            }
            for (int i = 0; i < (n - 1); i++) {
                T[i + 1, i] = b[i];
            }
            SquareMatrix A = L * T * L.Transpose();
            PrintMatrix(A);


            SymmetricMatrix D = new SymmetricMatrix(n);
            for (int i = 0; i < n; i++) {
                D[i, i] = a[i];
            }
            for (int i = 0; i < (n - 1); i++) {
                D[i + 1, i] = b[i];
            }
            for (int c = 1; c < (n - 1); c++) {
                for (int r = c + 1; r < n; r++) {
                    D[r, c - 1] = L[r, c];
                }
            }
            AasenDecomposition LTL = new AasenDecomposition(D);
            return (LTL);

        }

        private static void PrintMatrix (IMatrix M) {
            for (int i = 0; i < M.RowCount; i++) {
                for (int j = 0; j < M.ColumnCount; j++) {
                    Console.Write(" {0}", M[i, j]);
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        private static void PrintLTLMatrices (double[] h, double[] a, double[] b, IMatrix L) {

            Console.WriteLine("h=");
            for (int i = 0; i < h.Length; i++) {
                Console.Write("  {0}", h[i]);
            }
            Console.WriteLine();

            Console.WriteLine("alpha=");
            for (int i = 0; i < a.Length; i++) {
                Console.Write("  {0}", a[i]);
            }
            Console.WriteLine();

            Console.WriteLine("beta=");
            for (int i = 0; i < b.Length; i++) {
                Console.Write("  {0}", b[i]);
            }
            Console.WriteLine();

            Console.WriteLine("L=");
            for (int i = 0; i < L.RowCount; i++) {
                for (int j = 0 ; j <L.ColumnCount; j++) {
                    Console.Write(" {0}", L[i, j]);
                }
                Console.WriteLine();
            }
            Console.WriteLine();
        }

#endif


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
        /// <para>An alternative way of expressing the same relationship is to say that the eigenvalues of a matrix are its
        /// diagonal elements when the matrix is expressed in a basis that diagonalizes it. That is, given Z such that Z<sup>-1</sup>MZ = D,
        /// where D is diagonal, the columns of Z are the eigenvectors of M and the diagonal elements of D are the eigenvalues.</para>
        /// <para>Note that the eigenvectors of a matrix are not entirely unique. Given an eigenvector z, any scaled vector &#x3B1;z
        /// is an eigenvector with the same eigenvalue, so eigenvectors are at most unique up to a rescaling. If an eigenvalue
        /// is degenerate, i.e. there are two or more linearly independent eigenvectors with the same eigenvalue, then any linear
        /// combination of the eigenvectors is also an eigenvector with that eigenvalue, and in fact any set of vectors that span the
        /// same subspace could be taken as the eigenvector set corresponding to that eigenvalue.</para>
        /// <para>The eigenvectors of a symmetric matrix are always orthogonal and the eigenvalues are always real. The transformation
        /// matrix Z is thus orthogonal (Z<sup>-1</sup> = Z<sup>T</sup>).</para>
        /// <para>Finding the eigenvalues and eigenvectors of a symmetric matrix is an O(N<sup>3</sup>) operation.</para>
        /// </remarks>
        public RealEigensystem Eigensystem () {
            SymmetricMatrix clone = this.Clone();

            double[,] v = new double[dimension,dimension];
            for (int i = 0; i < dimension; i++) {
                v[i, i] = 1.0;
            }

            clone.Jacobi(v);

            double[] e = new double[dimension];
            for (int i = 0; i < dimension; i++) {
                e[i] = clone.values[i][i];
            }

            return (new RealEigensystem(dimension, e, v));

        }

        /// <summary>
        /// Computes the eigenvalues of the matrix. 
        /// </summary>
        /// <returns></returns>
        public double[] Eigenvalues () {
            SymmetricMatrix clone = this.Clone();
            clone.Jacobi(null);

            double[] e = new double[dimension];
            for (int i = 0; i < dimension; i++) {
                e[i] = clone.values[i][i];
            }

            return (e);
        }

        private void Jacobi (double[,] V) {

        	int count = 0;
			double sum;
				
			do {

				// check for non-convergence
				count++;
				if (count > 50) throw new NonconvergenceException();

                //Console.WriteLine("count = {0}", count);
                //PrintMatrix(this);
                //if (V!= null) PrintMatrix(V);

				// sweep over off-diagonal elements p,q where q < p
				for (int p=0; p<dimension; p++) {
					for (int q=0; q<p; q++) {

                        //Console.WriteLine("n={0}, p={1}, q={2}", count, p, q);

						// skip if already zero
						double M_pq = values[p][q];
						if (M_pq == 0.0) continue;

						// compute the rotation
                        // -pi/4 < phi < pi/4 is the rotation 
                        // s, c, and t are its sine, cosine, and tangent
                        // theta = cot(2 phi) and t2 = tan(phi/2)
						double theta = 0.5 * ( values[q][q] - values[p][p] ) / M_pq;
						double t;
						if (theta > 0) {
							t = 1.0 / ( Math.Sqrt(theta*theta + 1.0) + theta );
						} else {
							t = -1.0 / ( Math.Sqrt(theta*theta + 1.0) - theta );
						}
						double c = 1.0 / Math.Sqrt(1.0 + t*t);
						double s = t * c;
						double t2 = s / ( 1.0 + c);
                        //Console.WriteLine("s={0}, c={1}, t2={2}", s, c, t2);

						// do the rotation
						values[p][p] += -t * M_pq;
						values[q][q] += t * M_pq;
						for (int r=0; r<q; r++) {
							double M_pr = values[p][r];
                            //double M_pr = GetEntry(p, r);
                            double M_qr = values[q][r];
                            //double M_qr = GetEntry(q, r);
                            values[p][r] = M_pr - s * (M_qr + t2 * M_pr);
                            //SetEntry(p, r, M_pr - s * (M_qr + t2 * M_pr));
							values[q][r] = M_qr + s * ( M_pr - t2 * M_qr );
                            //SetEntry(q, r, M_qr + s * (M_pr - t2 * M_qr));
						}
						for (int r=q+1; r<p; r++) {
							double M_pr = values[p][r];
							double M_rq = values[r][q];
							values[p][r] = M_pr - s * ( M_rq + t2 * M_pr );
							values[r][q] = M_rq + s * ( M_pr - t2 * M_rq );
						}
						for (int r=p+1; r<dimension; r++) {
							double M_rp = values[r][p];
							double M_rq = values[r][q];
							values[r][p] = M_rp - s * ( M_rq + t2 * M_rp );
							values[r][q] = M_rq + s * ( M_rp - t2 * M_rq );
						}
						values[p][q] = 0.0;

                        //PrintMatrix(this);

                        // accumulate the rotations, if we are keeping track of them
                        if (V != null) {
                            for (int r = 0; r < dimension; r++) {
                                double V_rp = V[r, p];
                                double V_rq = V[r, q];
                                V[r, p] = c * V_rp - s * V_rq;
                                V[r, q] = s * V_rp + c * V_rq;
                            }
                        }

					}
				}

				// sum off-diagonal elements
				sum = 0;
				for (int p=0; p<dimension; p++) {
					for (int q=0; q<p; q++) {
						double M_pq = values[p][q];
						if (M_pq > 0) {
							sum += M_pq;
						} else {
							sum -= M_pq;
						}
					}
				}

			} while (sum > 0.0);

        }

        /*
        internal static void PrintMatrix (IMatrix M) {
            for (int r = 0; r < M.RowCount; r++) {
                for (int c = 0; c < M.ColumnCount; c++) {
                    Console.Write("{0} ", M[r, c]);
                }
                Console.WriteLine();
            }
            if (M.RowCount == 3) {
                double det =
                    M[0, 0] * M[1, 1] * M[2, 2] + M[0, 1] * M[1, 2] * M[2, 0] + M[0, 2] * M[1, 0] * M[2, 1] -
                    M[0, 2] * M[1, 1] * M[2, 0] - M[0, 0] * M[1, 2] * M[2, 1] - M[0, 1] * M[1, 0] * M[2, 2];
                Console.WriteLine("det = {0}", det);
            }
            Console.WriteLine("--");
        }
        */

        /*
        internal static void PrintMatrix (double[,] M) {
            for (int r = 0; r < M.GetLength(0); r++) {
                for (int c = 0; c < M.GetLength(1); c++) {
                    Console.Write("{0} ", M[r, c]);
                }
                Console.WriteLine();
            }
            Console.WriteLine("--");
        }
        */

        // operators

        internal static bool Equals (SymmetricMatrix M1, SymmetricMatrix M2) {
            if (Object.ReferenceEquals(M1, null)) {
                if (Object.ReferenceEquals(M2, null)) {
                    // both are null
                    return(true);
                } else {
                    // only M1 is null
                    return (false);
                }
            } else {
                if (Object.ReferenceEquals(M2, null)) {
                    // only M2 is null
                    return(false);
                } else {
                    // both are non-null
                    if (M1.Dimension != M2.Dimension) throw new DimensionMismatchException();
                    for (int r = 0; r < M1.Dimension; r++) {
                        for (int c = 0; c <= r; c++) {
                            if (M1[r, c] != M2[r, c]) return (false);
                        }
                    }
                    return (true);
                }
            }
        }

        /// <summary>
        /// Determines whether two symmetric matrices are equal.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>True if <paramref name="M1"/> and <paramref name="M2"/> are equal, otherwise false.</returns>
        public static bool operator == (SymmetricMatrix M1, SymmetricMatrix M2) {
            return (Equals(M1, M2));
        }

        /// <summary>
        /// Determines whether two symmetric matrices are not equal.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>False if <paramref name="M1"/> and <paramref name="M2"/> are equal, otherwise true.</returns>
        public static bool operator != (SymmetricMatrix M1, SymmetricMatrix M2) {
            return (!Equals(M1, M2));
        }

        /// <summary>
        /// Determines whether the given object is an equal matrix.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if <paramref name="obj"/> is an equal matrix, otherwise false.</returns>
        public override bool Equals (object obj) {
            return (Matrix.Equals(this, obj as IMatrix));
        }

        // arithmetic operators

        internal static SymmetricMatrix Add (SymmetricMatrix M1, SymmetricMatrix M2) {
            if (M1.Dimension != M2.Dimension) throw new DimensionMismatchException();
            SymmetricMatrix N = new SymmetricMatrix(M1.Dimension);
            for (int r = 0; r < N.Dimension; r++) {
                for (int c = 0; c <= r; c++) {
                    N[r, c] = M1[r, c] + M2[r, c];
                }
            }
            return (N);
        }

        /// <summary>
        /// Adds two symmetric matrices.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>The sum <paramref name="M1"/> + <paramref name="M2"/>.</returns>
        public static SymmetricMatrix operator + (SymmetricMatrix M1, SymmetricMatrix M2) {
            return (Add(M1, M2));
        }

        internal static SymmetricMatrix Subtract (SymmetricMatrix M1, SymmetricMatrix M2) {
            if (M1.Dimension != M2.Dimension) throw new DimensionMismatchException();
            SymmetricMatrix N = new SymmetricMatrix(M1.Dimension);
            for (int r = 0; r < N.Dimension; r++) {
                for (int c = 0; c <= r; c++) {
                    N[r, c] = M1[r, c] - M2[r, c];
                }
            }
            return (N);
        }

        /// <summary>
        /// Subtracts two symmetric matrices.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>The difference <paramref name="M1"/> - <paramref name="M2"/>.</returns>
        public static SymmetricMatrix operator - (SymmetricMatrix M1, SymmetricMatrix M2) {
            return (Subtract(M1, M2));
        }

        /// <summary>
        /// Multiplies two symmetric matrices.
        /// </summary>
        /// <param name="M1">The first matrix.</param>
        /// <param name="M2">The second matrix.</param>
        /// <returns>The matrix product <paramref name="M1"/> <paramref name="M2"/>.</returns>
        /// <remarks><para>Note that the product of two symmetric matrix is, in general, not itself symmetric.</para></remarks>
        public static SquareMatrix operator * (SymmetricMatrix M1, SymmetricMatrix M2) {
            return (SquareMatrix.Multiply(M1, M2));
        }

        // arithmatic operators with doubles

        private static SymmetricMatrix Multiply (double x, SymmetricMatrix M) {
            int d = M.Dimension;
            SymmetricMatrix N = new SymmetricMatrix(d);
            for (int r = 0; r < d; r++) {
                for (int c = 0; c <= r; c++) {
                    N[r, c] = x * M[r, c];
                }
            }
            return (N);
        }

        /// <summary>
        /// Multiplies a symmetric matrix by a real factor.
        /// </summary>
        /// <param name="f">The factor.</param>
        /// <param name="M">The matrix.</param>
        /// <returns>The product of the matrix and the factor.</returns>
        public static SymmetricMatrix operator * (double f, SymmetricMatrix M) {
            return (Multiply(f, M));
        }

        /// <summary>
        /// Computes the the quotient of a symmetric matrix and a real number.
        /// </summary>
        /// <param name="M">The matrix.</param>
        /// <param name="x">The real number.</param>
        /// <returns>The quotient <paramref name="M"/>/<paramref name="x"/>.</returns>
        public static SymmetricMatrix operator / (SymmetricMatrix M, double x) {
            return (Multiply(1.0 / x, M));
        }

    }

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

    /// <summary>
    /// Represents the Cholesky Decomposition of a symmetric, positive definite matrix. 
    /// </summary>
    /// <seealso cref="SymmetricMatrix.CholeskyDecomposition"/>
    public class CholeskyDecomposition : ISquareDecomposition {

        internal SymmetricMatrix sqrtM;

        /// <summary>
        /// Gets the dimension of the system.
        /// </summary>
        public virtual int Dimension {
            get {
                return (sqrtM.Dimension);
            }
        }

        /*
        public virtual LowerTriangularMatrix LeftFactor {
            get {
                LowerTriangularMatrix L = new LowerTriangularMatrix(Dimension);
                for (int r = 0; r < Dimension; r++) {
                    for (int c = 0; c <= r; c++) {
                        L[r, c] = sqrtM[r, c];
                    }
                }
                return (L);
            }
        }

        public virtual UpperTriangularMatrix RightFactor {
            get {
                UpperTriangularMatrix R = new UpperTriangularMatrix(Dimension);
                for (int c = 0; c < Dimension; c++) {
                    for (int r = 0; r <= c; r++) {
                        R[r, c] = sqrtM[r, c];
                    }
                }
                return (R);
            }
        }
        */
 
        /// <summary>
        /// Computes the solution vector that, when multiplied by the original matrix, produces the given left-hand side vector.
        /// </summary>
        /// <param name="rhs">The right-hand-side vector.</param>
        /// <returns>The left-hand-side (solution) vector.</returns>
        public virtual ColumnVector Solve (IList<double> rhs) {
            if (rhs == null) throw new ArgumentNullException("rhs");
            if (rhs.Count != Dimension) throw new InvalidOperationException();

            // Determine Ly = x
            double[] y = new double[Dimension];
            for (int i = 0; i < Dimension; i++) {
                y[i] = rhs[i];
                for (int j = 0; j < i; j++) {
                    y[i] -= sqrtM[j, i] * y[j];
                }
                y[i] = y[i] / sqrtM[i, i];
            }

            // Determine L^{T} z = y
            double[] z = new double[Dimension];

            // Re-use y for z, since we don't need a value after it's set
            for (int i = (Dimension - 1); i >= 0; i--) {
                z[i] = y[i];
                for (int j = (Dimension - 1); j > i; j--) {
                    z[i] -= sqrtM[i, j] * z[j];
                }
                z[i] = z[i] / sqrtM[i, i];
            }
            return (new ColumnVector(z));
        }

        /// <summary>
        /// Computes the inverse of the original matrix.
        /// </summary>
        /// <returns>M<sup>-1</sup></returns>
        public virtual SymmetricMatrix Inverse() {
                SymmetricMatrix MI = new SymmetricMatrix(Dimension);

                // do each column as a RHS
                for (int c = 0; c < Dimension; c++) {
                    MI[c, c] = 1.0 / sqrtM[c, c];
                    for (int r = c + 1; r < Dimension; r++) {
                        MI[r, c] = 0.0;
                        for (int i = c; i < r; i++) {
                            MI[r, c] -= sqrtM[r, i] * MI[i, c];
                        }
                        MI[r, c] = MI[r, c] / sqrtM[r, r];
                    }
                }

                // unnecessary?
                for (int c = 0; c < Dimension; c++) {
                    for (int r = (Dimension - 1); r >= c; r--) {
                        for (int i = r + 1; i < Dimension; i++) {
                            MI[r, c] -= sqrtM[r, i] * MI[i, c];
                        }
                        MI[r, c] = MI[r, c] / sqrtM[r, r];
                    }
                }
                // end unnecessary?

                return (MI);
        }

        ISquareMatrix ISquareDecomposition.Inverse () {
            return (Inverse());
        }

        /// <summary>
        /// Computes the determinant of the original matrix.
        /// </summary>
        /// <returns>det M</returns>
        public virtual double Determinant () {
                double lnDet = 0.0;
                for (int i = 0; i < Dimension; i++) {
                    lnDet += 2.0 * Math.Log(sqrtM[i, i]);
                }
                return (Math.Exp(lnDet));
        }

        internal CholeskyDecomposition (SymmetricMatrix sqrtM) {
            this.sqrtM = sqrtM;
        }

    }

#if FUTURE

    public class AasenDecomposition : ISquareDecomposition {

        // D is a packed representation of L and T. Since L and T have the form
        //     ( 1       )          ( A B     )
        //     ( 0 1     )          ( B A B   )
        // L = ( 0 X 1   )      T = (   B A B )
        //     ( 0 X X 1 )          (     B A )
        // we can fit all entries in
        //     ( A       )
        //     ( B A     )
        // D = ( X B A   )
        //     ( X X B A )

        private SymmetricMatrix D;

        internal AasenDecomposition (SymmetricMatrix D) {
            this.D = D;
        }

        /// <summary>
        /// Get the dimension of the original matrix.
        /// </summary>
        public int Dimension {
            get {
                return (D.Dimension);
            }
        }

        /// <summary>
        /// Computes the determinant of the original matrix.
        /// </summary>
        /// <returns>det M</returns>
        public double Determinant () {

            // the determinant is just the determinant of T, which can be computed by the recursion
            // det A_{n,n} = a_{n, n} det A_{n-1,n-1} - a_{n, n-1} a_{n-1, n} det A_{n-2,n-2}

            int n = Dimension;

            double a0 = D[0, 0];

            if (n == 1) return (a0);

            double a1 = D[0, 0] * D[1, 1] - D[0, 1] * D[1, 0];

            for (int i = 2; i < n; i++) {
                double a2 = D[i, i] * a1 - D[i - 1, i] * D[i, i - 1] * a0;
                a0 = a1;
                a1 = a2;
            }

            return (a1);
        }

        public ColumnVector Solve (IList<double> rhs) {


            int n = Dimension;

            double[] DD = new double[n];
            for (int i = 0; i < n; i++) {
                DD[i] = D[i, i];
            }

            double[] DL = new double[n - 1];
            double[] DU = new double[n - 1];
            for (int i = 0; i < (n - 1); i++) {
                DL[i] = D[i + 1, i];
                DU[i] = D[i + 1, i];
            }

            double[] UU = new double[n - 2];

            int[] permutation = new int[n];
            for (int i = 0; i < n; i++) {
                permutation[i] = i;
            }
            /*
            for (int i = 0; i < (n-1); i++) {

                // the diagonal is the pivot candidate
                double q = D[i];

                if (q == 0.0) throw new DivideByZeroException();

                // superdiagonal stays the same
                // subdiagonal gets divided by the pivot
                DL[i] = DL[i] / q;
                // product gets subtracted from the diagonal
                DD[i + 1] = DD[i + 1] - DU[i] * DL[i];

                // decide whether to pivot
                if (Math.Abs(DD[i]) >= Math.Abs(DL[i])) {

                    // diagonal is larger than subdiagonal; do not pivot

                    if (D[i] == 0.0) throw new DivideByZeroException();

                    double t = DL[i] / D[i];
                    DL[i] = t;
                    DD[i+1] = DD[i+1] - t * DU[i];

                } else {
                    // diagonal is smaller than subdiagonal; pivot

                    // switch row d and row d+1
                    double t1 = DD[i]; DD[i] = DL[i]; DL[i] = t1;
                    double t2 = DU[i]; DU[i] = DD[i + 1]; DD[i + 1] = t2;
                    UU[i] = DU[i + 1]; DU[i + 1] = 0.0;

                    double q = D[i];

                    DL[i] = DL[i] / q;
                    DU[i] = DU[i] - UU[i - 1] * DL[i - 1];
                    DD[i


                }

            }
            */

            throw new NotImplementedException();
        }

        /// <summary>
        /// Compute the inverse of the original matrix.
        /// </summary>
        /// <returns>M<sup>-1</sup></returns>
        public SymmetricMatrix Inverse () {
            throw new NotImplementedException();
        }

        ISquareMatrix ISquareDecomposition.Inverse () {
            return(Inverse());
        }

    }

#endif

}
