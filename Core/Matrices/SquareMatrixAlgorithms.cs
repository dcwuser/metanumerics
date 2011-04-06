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

        // A Householder reflection matrix is a rank-1 update to the identity.
        //   P = I - b v v^T
        // Unitarity requires b = 2 / |v|^2. To anihilate all but the first components of vector x,
        //   P x = a e_1
        // we must have
        //   v = x +/- |x| e_1
        // that is, all elements the same except the first, from which we have either added or subtracted |x|. This makes
        //   a = -/+ |x|
        // There are two way to handle the sign. One is to choose the sign that avoids cancelation when calculating v_1,
        // i.e. + for positive x and - for negative x. This works fine, but makes a negative for positive x, which is
        // weird (even 1 0 0 gets turned into -1 0 0). An alternative is to write

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
                        MatrixAlgorithms.SetEntry(aStore, dimension, dimension, a, a - 1, 0.0);
                        break;
                    }
                    sum += Math.Abs(e);
                    a--;
                }

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

        private static void FrancisTwoStep (double[] aStore, double[] qStore, int dimension, int a, int n, double sum, double product) {

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
