using System;
using System.Diagnostics;

namespace Meta.Numerics.Matrices {

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

        public static void IsolateCheapEigenvalues (double[] aStore, int[] perm, int dimension) {

            int n = dimension;
            while (n > 0) {
                int r = FindZeroRow(aStore, dimension, n);
                if (r < 0) {
                    break;
                } else {
                    n--;
                    SwapIndexes(aStore, perm, dimension, r, n);
                }
            }

            int m = 0;
            while (m < n) {
                int c = FindZeroColumn(aStore, dimension, m);
                if (c < 0) {
                    break;
                } else {
                    SwapIndexes(aStore, perm, dimension, c, m);
                    m++;
                }
            }

        }

        private static int FindZeroRow (double[] aStore, int dimension, int n) {
            for (int r = n - 1; r >= 0; r--) {
                if (IsZeroRow(aStore, dimension, r, n)) return (r);
            }
            return (-1);
        }

        private static bool IsZeroRow (double[] aStore, int dimension, int r, int n) {
            for (int c = 0; c < n; c++) {
                if ((r != c) && (aStore[r + dimension * c] != 0.0)) return (false);
            }
            return (true);
        }

        private static int FindZeroColumn (double[] aStore, int dimension, int m) {
            for (int c = m; c < dimension; c++) {
                if (IsZeroColumn(aStore, dimension, c, m)) return(c);
            }
            return (-1);
        }

        private static bool IsZeroColumn (double[] aStore, int dimension, int c, int m) {
            for (int r = m; r < dimension; r++) {
                if ((r != c) && (aStore[r + dimension * c] != 0.0)) return (false);
            }
            return (true);
        }

        private static void SwapIndexes (double[] aStore, int[] perm, int dimension, int p, int q) {
            if (p == q) return;
            Blas1.dSwap(aStore, p, dimension, aStore, q, dimension, dimension);
            Blas1.dSwap(aStore, dimension * p, 1, aStore, dimension * q, 1, dimension);
            if (perm != null) Global.Swap<int>(ref perm[p], ref perm[q]);
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
                VectorAlgorithms.GenerateHouseholderReflection(aStore, offset, 1, length, out a);

                // determine P * A
                for (int c = k + 1; c < dim; c++) {
                    VectorAlgorithms.ApplyHouseholderReflection(aStore, offset, 1, aStore, dim * c + (k + 1), 1, length);
                }
                // determine A * P
                for (int r = 0; r < dim; r++) {
                    VectorAlgorithms.ApplyHouseholderReflection(aStore, offset, 1, aStore, dim * (k + 1) + r, dim, length);

                }

                // if we are keeping track of the transformation, determine V * P
                if (vStore != null) {
                    for (int r = 0; r < dim; r++) {
                        VectorAlgorithms.ApplyHouseholderReflection(aStore, offset, 1, vStore, dim * (k + 1) + r, dim, length);
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

        // EIGENVALUE ALGORITHMS
 
        // 

        public static Complex[] ReduceToRealSchurForm (double[] aStore, double[] qStore, int dimension) {

            // keep track of eigenvalues
            Complex[] eigenvalues = new Complex[dimension];

            // keep track of interations
            int count = 0;

            // isolate the upper and lower boundaries of the curent subproblem
            // p is the upper index, n the lower index

            int n = dimension - 1;
            while (n >= 0) {

                // move up the matrix from endpoint n, looking for subdiagonal elements negligible compared to the neighboring diagonal elements
                // if we find one, that reduces the problem to the submatrix between p and n
                int p = n;
                while (p > 0) {

                    double d = Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, p, p)) + Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, p - 1, p - 1));
                    double e = Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, p, p - 1));

                    if (d + e == d) {

                      //  double f = d * Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, p, p) - MatrixAlgorithms.GetEntry(aStore, dimension, dimension, p - 1, p - 1));
                      //  double g = e * Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, p - 1, p));

                      //  if (f + g == f) {
                            MatrixAlgorithms.SetEntry(aStore, dimension, dimension, p, p - 1, 0.0);
                            break;
                      //  }
                    }

                    p--;
                }

                /*
                if (p == n) {
                    // one eigenvalue
                    // reduce n and re-set count
                } else {
                    // compute m, matrix elements
                    // reduce n and re-set count
                    if (p == m) {
                        // two eigenvalues
                    } else {
                        // do step
                    }
                }
                */

                //Print(a, dimension, p, n);

                // get the (indexes for) the entries in the trailing 2x2 matrix
                int m = n - 1;
                int ammi = MatrixAlgorithms.GetIndex(dimension, dimension, m, m);
                int amni = MatrixAlgorithms.GetIndex(dimension, dimension, m, n);
                int anmi = MatrixAlgorithms.GetIndex(dimension, dimension, n, m);
                int anni = MatrixAlgorithms.GetIndex(dimension, dimension, n, n);

                if (n - p > 1) {

                    count++;
                    if (count > 32) throw new NonconvergenceException();

                    double tr = aStore[ammi] + aStore[anni];
                    double det = aStore[ammi] * aStore[anni] - aStore[amni] * aStore[anmi];

                    // ad hoc shift
                    if ((count == 8) || (count == 16) || (count == 24)) {
                        double w = Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, n - 1)) +
                            Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n - 1, n - 2));
                        tr = 2.0 * w;
                        det = w * w;
                    }

                    FrancisTwoStep(aStore, qStore, dimension, p, n, tr, det);
                    
                } else {

                    if (p == n) {

                        eigenvalues[n] = aStore[anni];
 
                    } else if (p == m) {
                        
                        double sn, cn;
                        TwoByTwoRealSchur(ref aStore[ammi], ref aStore[amni], ref aStore[anmi], ref aStore[anni], out sn, out cn, out eigenvalues[m], out eigenvalues[n]);

                        // Multiply A from left by the rotation matrix R
                        for (int cc = p + 2; cc < dimension; cc++) {
                            int i = MatrixAlgorithms.GetIndex(dimension, dimension, p, cc);
                            int j = MatrixAlgorithms.GetIndex(dimension, dimension, p + 1, cc);
                            double t = aStore[i];
                            aStore[i] = cn * t + sn * aStore[j];
                            aStore[j] = cn * aStore[j] - sn * t;
                        }
                        // Multiply A from the right by R^T
                        for (int rr = 0; rr < p; rr++) {
                            int i = MatrixAlgorithms.GetIndex(dimension, dimension, rr, p);
                            int j = MatrixAlgorithms.GetIndex(dimension, dimension, rr, p + 1);
                            double t = aStore[i];
                            aStore[i] = cn * t + sn * aStore[j];
                            aStore[j] = cn * aStore[j] - sn * t;
                        }
                        // Multiply Q^T from the left by R
                        if (qStore != null) {
                            for (int rr = 0; rr < dimension; rr++) {
                                int i = MatrixAlgorithms.GetIndex(dimension, dimension, rr, p);
                                int j = MatrixAlgorithms.GetIndex(dimension, dimension, rr, p + 1);
                                double t = qStore[i];
                                qStore[i] = cn * t + sn * qStore[j];
                                qStore[j] = cn * qStore[j] - sn * t;
                            }
                        }
                        
                    }
                    
                    n = p - 1;
                    count = 0;
                }
            }

            return (eigenvalues);

        }

#if PAST

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
                        MatrixAlgorithms.SetEntry(aStore, dimension, dimension, a, a - 1, 0.0);
                        break;
                    }
                    sum += Math.Abs(e);
                    a--;
                }

                /*
                Console.WriteLine("count = {0}", count);
                double[] qt = MatrixAlgorithms.Transpose(qStore, dimension, dimension);
                double[] qa = MatrixAlgorithms.Multiply(qStore, dimension, dimension, aStore, dimension, dimension);
                double[] qaqt = MatrixAlgorithms.Multiply(qa, dimension, dimension, qt, dimension, dimension);
                MatrixAlgorithms.PrintMatrix(qaqt, dimension, dimension);
                */

                // check for maximum numbers of iterations without finding an eigenvalue
                if (count > 32) {
                    throw new NonconvergenceException();
                }

                // we are working between a and n, and our block is at least 3 wide

                // reduce if possible, otherwise do a iteration step

                if (a == n) {
                    // 1 X 1 block isolated
                    lambdas[a] = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a);
                    n -= 1;
                    count = 0;
                    sum_old = Double.PositiveInfinity;
                } else if (a == (n - 1)) {
                    // 2 X 2 block isolated
                    double Aaa = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a);
                    double Aba = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a + 1, a);
                    double Aab = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a, a + 1);
                    double Abb = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, a + 1, a + 1);
                    // eigenvalues are given by the quadratic formula
                    //   e = \frac{a_11 + a22}{2} \pm q
                    // descriminant q^2 = ( \frac{a_11 - a_22}{2} )^2 + a_21 a_1
                    double c = (Aaa + Abb) / 2.0;
                    double d = Aaa * Abb - Aba * Aab;
                    double q2 = c * c - d;
                    if (q2 >= 0.0) {
                        // eigenvalues are real
                        
                        double q = Math.Sqrt(q2);
                        if (c >= 0.0) {
                            lambdas[a] = c + q;
                        } else {
                            lambdas[a] = c - q;
                        }
                        // i tried to do this as c + Math.Sign(c) * q, but when c = 0, Math.Sign(c) = 0, not 1 
                        lambdas[a + 1] = d / lambdas[a];
                        
                        /*
                        double sn, cn;
                        TwoByTwoSchur(ref Aaa, ref Aab, ref Aba, ref Abb, out sn, out cn);
                        MatrixAlgorithms.SetEntry(aStore, dimension, dimension, a, a, Aaa);
                        MatrixAlgorithms.SetEntry(aStore, dimension, dimension, a, a + 1, Aab);
                        MatrixAlgorithms.SetEntry(aStore, dimension, dimension, a + 1, a, Aba);
                        MatrixAlgorithms.SetEntry(aStore, dimension, dimension, a + 1, a + 1, Abb);

                        // Multiply A from left by the rotation matrix R
                        for (int cc = a + 2; cc < dimension; cc++) {
                            int i = MatrixAlgorithms.GetIndex(dimension, dimension, a, cc);
                            int j = MatrixAlgorithms.GetIndex(dimension, dimension, a + 1, cc);
                            double t = aStore[i];
                            aStore[i] = cn * t + sn * aStore[j];
                            aStore[j] = cn * aStore[j] - sn * t;
                        }
                        // Multiply A from the right by R^T
                        for (int rr = 0; rr < a; rr++) {
                            int i = MatrixAlgorithms.GetIndex(dimension, dimension, rr, a);
                            int j = MatrixAlgorithms.GetIndex(dimension, dimension, rr, a + 1);
                            double t = aStore[i];
                            aStore[i] = cn * t + sn * aStore[j];
                            aStore[j] = cn * aStore[j] - sn * t;
                        }
                        // Multiply Q^T from the left by R
                        if (qStore != null) {
                            for (int cc = 0; cc < dimension; cc++) {
                                int i = MatrixAlgorithms.GetIndex(dimension, dimension, a, cc);
                                int j = MatrixAlgorithms.GetIndex(dimension, dimension, a + 1, cc);
                                double t = qStore[i];
                                qStore[i] = cn * t + sn * qStore[j];
                                qStore[j] = cn * qStore[j] - sn * t;
                            }
                        }

                        lambdas[a] = Aaa;
                        lambdas[a + 1] = Abb;
                        */
                    } else {
                        double q = Math.Sqrt(-q2);
                        lambdas[a] = new Complex(c, q);
                        lambdas[a + 1] = new Complex(c, -q);
                    }

                    n -= 2;
                    count = 0;
                    sum_old = Double.PositiveInfinity;
                } else {

                    // use the lower-left 2 X 2 matrix to generate an approximate eigenvalue pair

                    int m = n - 1;
                    double Amm = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, m);
                    double Amn = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, n);
                    double Anm = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, m);
                    double Ann = MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, n);

                    // the trace of this 2 X 2 matrix is the sum of its eigenvalues, and its determinate is their product

                    double tr = Amm + Ann;
                    double det = Amm * Ann - Amn * Anm;

                    // ad hoc shift
                    if ((count == 8) || (count == 16) || (count == 24)) {
                        double w = Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, n, n - 1)) +
                            Math.Abs(MatrixAlgorithms.GetEntry(aStore, dimension, dimension, m, m - 1));
                        tr = 2.0 * w;
                        det = w * w;
                    }

                    FrancisTwoStep(aStore, qStore, dimension, a, n, tr, det);
                    sum_old = sum;
                    count++;
                }

            }

            return (lambdas);


        }

#endif

        private static void FrancisTwoStep (double[] aStore, double[] qStore, int dimension, int a, int n, double sum, double product) {

            // to apply the implicit Q theorem, get the first column of what (A - mu I)(A - mu* I) would be
            double[] v = GenerateFirstColumn(aStore, dimension, a, n, sum, product);

            int d = n - a + 1;
            Debug.Assert(d >= 3);

            // generate a householder reflection for this column, and apply it to the matrix
            double e;
            VectorAlgorithms.GenerateHouseholderReflection(v, 0, 1, 3, out e);

            // determine P * A
            for (int c = a; c < dimension; c++) {
                VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + a, 1, 3);
            }
            // determine A * P
            for (int r = 0; r <= Math.Min(a + 3, n); r++) {
                VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, aStore, dimension * a + r, dimension, 3);
            }
            // determine Q * P
            if (qStore != null) {
                for (int r = 0; r < dimension; r++) {
                    VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, qStore, dimension * a + r, dimension, 3);
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
                VectorAlgorithms.GenerateHouseholderReflection(v, 0, 1, 3, out e);

                //v = GenerateHouseholderReflection(aStore, dimension * k + k + 1, 1, 3);

                // determine P * A and A * P
                aStore[dimension * k + (k + 1)] = e;
                aStore[dimension * k + (k + 2)] = 0.0;
                aStore[dimension * k + (k + 3)] = 0.0;
                for (int c = k + 1; c < dimension; c++) {
                    VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + (k + 1), 1, 3);
                }
                for (int r = 0; r < Math.Max(k + 4, dimension); r++) {
                    VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, aStore, dimension * (k + 1) + r, dimension, 3);
                }
                // if tracking eigenvalues, determine Q * P
                if (qStore != null) {
                    for (int r = 0; r < dimension; r++) {
                        VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, qStore, dimension * (k + 1) + r, dimension, 3);
                    }
                }

            }

            //Write(aStore, dimension, dimension);

            // restoring Householder form in the last column requires only a length-2 Householder reflection
            // a Givens rotation would be another possibility here
            int l = n - 2;
            v[0] = aStore[dimension * l + (l + 1)];
            v[1] = aStore[dimension * l + (l + 2)];
            VectorAlgorithms.GenerateHouseholderReflection(v, 0, 1, 2, out e);
            aStore[dimension * l + (l + 1)] = e;
            aStore[dimension * l + (l + 2)] = 0.0;
            for (int c = l + 1; c < dimension; c++) {
                VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, aStore, dimension * c + (l + 1), 1, 2);
            }
            for (int r = 0; r <= n; r++) {
                VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, aStore, dimension * (l + 1) + r, dimension, 2);
            }
            if (qStore != null) {
                for (int r = 0; r < dimension; r++) {
                    VectorAlgorithms.ApplyHouseholderReflection(v, 0, 1, qStore, dimension * (l + 1) + r, dimension, 2);
                }
            }
            //double cos, sin;
            //GenerateGivensRotation(ref aStore[dimension * l + l + 1], ref aStore[dimension * l + l + 2], out cos, out sin);

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

        // We want to find the rotation that brings a 2 X 2 matrix into triangular form.
        //   (  c  s ) ( a11  a12 ) ( c  -s )
        //   ( -s  c ) ( a21  a22 ) ( s   c )
        // Multiplying out gives
        //   ( c^2 a11 + cs a21 + cs a12 + s^2 a22  c^2 a12 - cs a11 + cs a22 - s^2 a21 )
        //   ( c^2 a21 - cs a11 + cs a22 - s^2 a12  c^2 a22 - cs a21 - cs a12 + s^2 a11 )

        // We want a21' = 0. Since the rotation must vanish when a21 = 0 and grow as a21 grows, make the Ansatz that s ~ a21, i.e.
        //   s = \frac{a_{21}}{\sqrt{a_{21}^2 + b^2}}  c = \frac{b}{\sqrt{a_{21}^2 + b^2}}
        // It is then straightforward to derive
        //   b = \frac{a_{11} - a_{22} \pm q}{2}
        // where q^2 = ( a_{11} - a_{22} )^2 + 4 a12 a21 is the same descriminant as appears in the eigenvalue problem, so this rotation
        // exists iff the eigenvalues are real.

        // A useful test case is
        //   ( 1 4 )
        //   ( 2 3 )
        // with eigenvalues 5, -1 and s = -\sqrt{1/5}.

        private static void TwoByTwoRealSchur (ref double a11, ref double a12, ref double a21, ref double a22, out double s, out double c, out Complex e1, out Complex e2) {

            // compute some quantities we will use 
            double u = a11 + a22;
            double v = a11 - a22;
            double w = a11 * a22 - a12 * a21;
            double q2 = v * v + 4.0 * a12 * a21;
            // note u is the trace and w is the determinant

            if (q2 >= 0.0) {

                // the descriminant is positive so the eigenvalues are real
                double q = Math.Sqrt(q2);

                // find the rotation sets a21' = 0
                // in the equation for b, choose the sign so as to minimize cancelation
                double b = (v >= 0.0) ? (v + q) / 2.0 : (v - q) / 2.0;
                double rho = MoreMath.Hypot(a21, b);
                s = a21 / rho;
                c = b / rho;

                // Note that a12' - a21' = a12 - a21, and since a21' = 0, a12' = a12 - a21
                a12 = a12 - a21;
                a21 = 0.0;

                // the eigenvalues are (u \pm q) / 2
                // we avoid cancelation by computing the non-canceling one first and
                // computing the other using the fact that their product equals the determinant
                if (u >= 0.0) {
                    a11 = (u + q) / 2.0;
                    a22 = w / a11;
                } else {
                    a22 = (u - q) / 2.0;
                    a11 = w / a22;
                }
                // we have placed the u + q eigenvalue in the a11 slot. if v >= 0, that is where the rotation puts it
                // if v < 0, the u - q eigenvalue belongs these, so we need to swap them
                if (v < 0) { double t = a11; a11 = a22; a22 = t; }

                e1 = a11;
                e2 = a22;

            } else {

                // In the q2 < 0 case, we can't zero a21. But we can rotate so a11' = a22'
                /*
                double r = -(a12 + a21) / v;
                double t = Math.Sign(r) / (Math.Abs(r) + MoreMath.Hypot(1.0, r));
                c = 1.0 / MoreMath.Hypot(1.0, t);
                s = t * c;

                // Since rotations preserve the trace, the equal diagonal elements must equal the average of the previous elements 
                a11 = u / 2.0;
                a22 = a11;
                */
                double q = Math.Sqrt(-q2);
                e1 = new Complex(u / 2.0, q / 2.0);
                e2 = e1.Conjugate;

                s = 0.0;
                c = 1.0;

            }
        }

        public static void SchurEigensystem (double[] aStore, int n, out Complex[] eigenvalues, out Complex[][] eigenvectors) {

            eigenvalues = new Complex[n];
            eigenvectors = new Complex[n][];

            int i = n - 1;
            while (i >= 0) {
                // find the eigenvector corresponding to the ith eigenvalue

                if ((i == 0) || (MatrixAlgorithms.GetEntry(aStore, n, n, i, i-1) == 0.0)) {

                    // We have no sub-diagonal element on the ith row. Since
                    //  (XXXXX) ( X )   ( X )
                    //  (XXXXX) ( X )   ( X )
                    //  (  aXX) ( 1 ) = ( a )
                    //  (   XX) ( 0 )   ( 0 )
                    //  (   XX) ( 0 )   ( 0 )
                    // a vector with v[i] = 1, v[j > i] = 0, and v[j < i] properly chosen will be an eigenvector.

                    // the eigenvalue is the ith diagonal element
                    double e = MatrixAlgorithms.GetEntry(aStore, n, n, i, i);

                    // we take the ith component of the eigenvector to be one; higher components must be zero
                    Complex[] v = new Complex[n];
                    v[i] = 1.0;

                    CompleteEigenvector(aStore, n, v, e, i, i - 1);
                    //WriteTransformedEigenvector(i, e, v, n, Q);
                    eigenvectors[i] = v;
                    eigenvalues[i] = e;

                    i--;

                } else {

                    // We have a 2x2 sub-block
                    //   ( X X X )
                    //   (  [X X])
                    //   (  [X X])
                    // Find its eigenvalues

                    Complex e1, e2;
                    TwoByTwoEigenvalues(
                        MatrixAlgorithms.GetEntry(aStore, n, n, i - 1, i - 1),
                        MatrixAlgorithms.GetEntry(aStore, n, n, i - 1, i),
                        MatrixAlgorithms.GetEntry(aStore, n, n, i, i - 1),
                        MatrixAlgorithms.GetEntry(aStore, n, n, i, i),
                        out e1, out e2
                    );
                    //TwoByTwoEigenvalues(a[i - 1, i - 1], a[i - 1, i], a[i, i - 1], a[i, i], out e1, out e2);

                    // For each eigenvalue, compute the trailing eigenvector components and then the remaining components.
                    // For the trailing components, since 
                    //   (a11 a12) ( v1 ) = \lambda ( v1 )
                    //   (a21 a22) ( v2 )           ( v2 )
                    // It follows that
                    //   v1 = (\lambda - a22) / a21 * v2
                    // and we normalize by choosing v2 = 1. Since we are only here if a21 != 0, there is no danger of
                    // dividing by zero, but there may well be a danger that (\lambda - a22) looses significant
                    // precision due to cancelation.

                    Complex[] v = new Complex[n];
                    v[i] = 1.0;
                    //v[i - 1] = (e1 - a[i, i]) / a[i, i - 1];
                    v[i - 1] = (e1 - MatrixAlgorithms.GetEntry(aStore, n, n, i, i)) / MatrixAlgorithms.GetEntry(aStore, n, n, i, i - 1);
                    CompleteEigenvector(aStore, n, v, e1, i, i - 2);
                    //WriteTransformedEigenvector(i, e1, v, n, Q);
                    eigenvectors[i] = v;
                    eigenvalues[i] = e1;

                    v = new Complex[n];
                    v[i] = 1.0;
                    //v[i - 1] = (e2 - a[i, i]) / a[i, i - 1];
                    v[i - 1] = (e2 - MatrixAlgorithms.GetEntry(aStore, n, n, i, i)) / MatrixAlgorithms.GetEntry(aStore, n, n, i, i - 1);
                    CompleteEigenvector(aStore, n, v, e2, i, i - 2);
                    //WriteTransformedEigenvector(i - 1, e2, v, n, Q);
                    eigenvectors[i - 1] = v;
                    eigenvalues[i - 1] = e2;

                    i -= 2;

                }

            }

        }

        // v is assumed to 

        private static void CompleteEigenvector (double[] aStore, int n, Complex[] v, Complex e, int i, int j) {

            while (j >= 0) {

                if ((j == 0) || (MatrixAlgorithms.GetEntry(aStore, n, n, j, j - 1) == 0.0)) {
                    Complex t = 0.0;
                    for (int k = j + 1; k <= i; k++) {
                        t += MatrixAlgorithms.GetEntry(aStore, n, n, j, k) * v[k];
                        //t += a[j, k] * v[k];
                    }

                    if (t == 0.0) {
                        // this branch exists to avoid undefined 0.0 / 0.0
                        v[j] = 0.0;
                    } else {
                        //v[j] = t / (e - a[j, j]);
                        v[j] = t / (e - MatrixAlgorithms.GetEntry(aStore, n, n, j, j));
                    }

                    j--;
                } else {
                    Complex t1 = 0.0;
                    Complex t2 = 0.0;
                    for (int k = j + 1; k <= i; k++) {
                        t1 -= MatrixAlgorithms.GetEntry(aStore, n, n, j - 1, k) * v[k];
                        t2 -= MatrixAlgorithms.GetEntry(aStore, n, n, j, k) * v[k];
                        //t1 -= a[j - 1, k] * v[k];
                        //t2 -= a[j, k] * v[k];
                    }
                    TwoByTwoSystem(
                        MatrixAlgorithms.GetEntry(aStore, n, n, j - 1, j - 1) - e,
                        MatrixAlgorithms.GetEntry(aStore, n, n, j - 1, j),
                        MatrixAlgorithms.GetEntry(aStore, n, n, j, j - 1),
                        MatrixAlgorithms.GetEntry(aStore, n, n, j, j) - e,
                        ref t1, ref t2
                    );
                    //TwoByTwoSystem(a[j - 1, j - 1] - e, a[j - 1, j], a[j, j - 1], a[j, j] - e, ref t1, ref t2);
                    v[j - 1] = t1;
                    v[j] = t2;

                    j -= 2;
                }
            }

        }


        // The eigenvalues of a 2x2 matrix
        //   (a11 a12)
        //   (a21 a22)
        // Are given by the quadratic equation \det(A - \lambda 1) = 0. To avoid cancelation, we do not use the +/-
        // form of the quadratic formula, but instead take the root that does not involve cancelation, then use
        // \lambda_1 \lambda_2 = \det A to determine the other root.

        private static void TwoByTwoEigenvalues (double a11, double a12, double a21, double a22, out Complex e1, out Complex e2) {

            double u = a11 + a22;
            double v = a11 - a22;
            double x = a12 * a21;
            double q2 = v * v + 4.0 * x;
            //double q2 = u * u - 4.0 * (a11 * a22 - a12 * a21);
            if (q2 >= 0.0) {
                double q = Math.Sqrt(q2);
                double det = a11 * a22 - x;
                e1 = (u >= 0) ? (u + q) / 2.0 : (u - q) / 2.0;
                e2 = det / e1;
                // because we compute e1 without cancelation, e1 can only be zero if u = 0 and q = 0
            } else {
                double q = Math.Sqrt(-q2);
                e1 = new Complex(u, q) / 2.0;
                e2 = new Complex(u, -q) / 2.0;
            }
        }

        // Solve the 2 x 2 system
        //   ( a11 a12 ) ( x1 ) = ( b1 )
        //   ( a21 a22 ) ( x2 )   ( b2 )
        // 

        private static void TwoByTwoSystem (Complex a11, double a12, double a21, Complex a22, ref Complex x1, ref Complex x2) {

            Complex det = a11 * a22 - a12 * a21;

            Complex t = x1;
            x1 = (a22 * t - a12 * x2) / det;
            x2 = (a11 * x2 - a21 * t) / det;

        }

        public static double[] Power (double[] aStore, int dimension, int power) {

            // We take powers via the exponentiation-by-squaring algorithm.
            // This is not strictly optimal, but it is very simple, is optimal in most cases (e.g. all n<15),
            // and is nearly optimal (e.g. 6 multiplies instead of 5 for n=15) even when it is not perfectly optimal.

            if (power == 1) {
                return (aStore);
                // Returning the given storage for A^1 instead of copying it saves space and time and works fine when this point is reached via recursion,
                // since that storage will just be read for a multiply, not returned to the calling user. But it would be bad to do so when the user requests
                // A^1, since the returned matrix would not be independent of the original. So don't call this method to compute A^1 for the user.
            } else {
                if (power % 2 == 0) {
                    // return (A * A)^{n/2}
                    double[] bStore = MatrixAlgorithms.Multiply(aStore, dimension, dimension, aStore, dimension, dimension);
                    double[] cStore = Power(bStore, dimension, power / 2);
                    return (cStore);
                } else {
                    // return A * (A * A)^{(n-1)/2}
                    double[] bStore = MatrixAlgorithms.Multiply(aStore, dimension, dimension, aStore, dimension, dimension);
                    double[] cStore = Power(bStore, dimension, (power - 1) / 2);
                    double[] dStore = MatrixAlgorithms.Multiply(aStore, dimension, dimension, cStore, dimension, dimension);
                    return (dStore);
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
