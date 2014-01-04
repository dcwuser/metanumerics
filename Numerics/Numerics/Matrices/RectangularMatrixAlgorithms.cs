using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace Meta.Numerics.Matrices {

    // The new design idea is as follows:
    // 1. The Matrix class is public. It is basically a pure shim over the storage format. To actually do all the operations it makes public,
    //    it calls into the MatrixAlgorithm class.
    // 2. The MatrixAlgorithm class is internal. It contains all the algorithmic logic for doing matrix operations on the storage format. It
    //    calls into BLAS to perform optomized low-level operations.
    // 3. The Blas classes have optomized low-level operations.
    // (1) and (2) should be repeated for each kind of matrix: square, symmetric, tridiagonal, etc.

    internal static class MatrixAlgorithms {

        /*
        public static void PrintMatrix (double[] store, int nRows, int nCols) {
            for (int r = 0; r < nRows; r++) {
                for (int c = 0; c < nCols; c++) {
                    Console.Write("{0}  ", GetEntry(store, nRows, nCols, r, c));
                }
                Console.WriteLine();
            }
        }
        */

        public static double[] AllocateStorage (int nRows, int nCols) {
            return (new double[nRows * nCols]);
        }

        // basic access: O(1) operations

        public static int GetIndex (int nRows, int nCols, int r, int c) {
            return (nRows * c + r);
        }

        public static void SetEntry (double[] store, int nRows, int nCols, int r, int c, double value) {
            store[nRows * c + r] = value;
        }

        public static double GetEntry (double[] store, int nRows, int nCols, int r, int c) {
            return (store[nRows * c + r]);
        }

        // transformations and simple arithmetic: O(N^2) operations

        public static double[] Copy (double[] store, int nRows, int nColumns) {
            double[] cloneStore = new double[store.Length];
            store.CopyTo(cloneStore, 0);
            return (cloneStore);
        }

        public static double[] Transpose (double[] store, int nRows, int nColumns) {
            double[] tStore = new double[store.Length];
            for (int r = 0; r < nRows; r++) {
                for (int c = 0; c < nColumns; c++) {
                    tStore[nColumns * r + c] = store[nRows * c + r];
                }
            }
            return (tStore);
        }

        public static double OneNorm (double[] store, int nRows, int nColumns) {
            double norm = 0.0;
            for (int c = 0; c < nColumns; c++) {
                double csum = Blas1.dNrm1(store, c * nRows, 1, nRows);
                if (csum > norm) norm = csum;
            }
            return (norm);
        }

        public static double InfinityNorm (double[] store, int nRows, int nColumns) {
            double norm = 0.0;
            for (int r = 0; r < nRows; r++) {
                double rsum = Blas1.dNrm1(store, r, nRows, nColumns);
                if (rsum > norm) norm = rsum;
            }
            return (norm);
        }

        public static double[] Add (double[] xStore, double[] yStore, int nRows, int nCols) {
            double[] store = new double[nRows * nCols];
            for (int i = 0; i < store.Length; i++) {
                store[i] = xStore[i] + yStore[i];
            }
            return (store);
        }

        public static double[] Subtract (double[] xStore, double[] yStore, int nRows, int nCols) {
            double[] store = new double[nRows * nCols];
            for (int i = 0; i < store.Length; i++) {
                store[i] = xStore[i] - yStore[i];
            }
            return (store);
        }

        public static double[] Multiply (double alpha, double[] store, int nRows, int nCols) {
            double[] result = new double[nRows * nCols];
            for (int i = 0; i < store.Length; i++) {
                result[i] = alpha * store[i];
            }
            return (result);
        }

        public static double[] Multiply (double[] aStore, double[] xStore, int nRows, int nCols) {
            double[] yStore = new double[nCols];
            Blas2.dGemv(aStore, 0, 1, nRows, xStore, 0, 1, yStore, 0, 1, nRows, nCols);
            return (yStore);
        }

        // multiplication: an O(N^3) algorithm

        public static double[] Multiply (double[] aStore, int aRows, int aCols, double[] bStore, int bRows, int bCols) {

            double[] abStore = new double[aRows * bCols];
            for (int j = 0; j < bCols; j++) {
                int aRows_j = aRows * j;
                int bRows_j = bRows * j;
                for (int k = 0; k < aCols; k++) {
                    double t = bStore[bRows_j + k];
                    if (t != 0.0) {
                        int aRows_k = aRows * k;
                        for (int i = 0; i < aRows; i++) {
                            abStore[aRows_j + i] += t * aStore[aRows_k + i];
                        }
                    }
                }
            }

            return (abStore);

        }

        public static void SolveUpperRightTriangular (double[] tStore, int tRows, int tCols, double[] xStore, int xOffset) {
            Debug.Assert(tRows >= tCols);
            // pre-calculate some useful quantities
            int trp1 = tRows + 1;
            int tcm1 = tCols - 1;
            // initialize the x index (to its highest value) the the diagonal t index (to its highest value)
            int xIndex = xOffset + tcm1;
            int tIndex = trp1 * tcm1;
            // do the first, trivial back-substitution
            xStore[xIndex] /= tStore[tIndex];
            // do the remaining back-substitutions, each of which requires subtracting n previously computed x-values
            for (int n = 1; n <= tcm1; n++) {
                xIndex -= 1;
                tIndex -= trp1;
                double t = Blas1.dDot(tStore, tIndex + tRows, tRows, xStore, xIndex + 1, 1, n);
                xStore[xIndex] = (xStore[xIndex] - t) / tStore[tIndex];
            }
        }

        // more complicated operations: also O(N^3)

        // QR decomposition write A = Q R, where Q is orthogonal (Q^T = Q^{-1}) and R is upper-right triangular.

        // Given a QR decomposition, it is easy to solve Ax=b by writing (QR)x=b, Q(Rx)=b, Rx=y and Qy=b.
        // Then y = Q^T b (an N^2 multiplication) and R x = y is a triangular system (also N^2 to solve).

        // To form a QR decomposition, we apply householder reflections to zero each column in turn.

        // XXXXX    ABCDE    XXXXX    XXXXX    XXXXX    XXXXX
        // XXXXX    0BCDE    0ABCD    0XXXX    0XXXX    0XXXX
        // XXXXX    0BCDE    00BCD    00ABC    00XXX    00XXX
        // XXXXX -> 0BCDE -> 00BCD -> 000BC -> 000AB -> 000XX
        // XXXXX    0BCDE    00BCD    000BC    0000B    0000A    
        // XXXXX    0BCDE    00BCD    000BC    0000B    00000
        // XXXXX    0BCDE    00BCD    000BC    0000B    00000

        // This gives us (Q_N \cdots Q_2 Q_1) A = R, so A = (Q_N \cdots Q_2 Q_1)^T R = (Q_1^T Q_2^T \cdots Q_N^T) R

        public static void QRDecompose (double[] store, double[] qtStore, int rows, int cols) {

            double[] vStore = new double[rows];

            // loop over columns
            for (int k = 0; k < cols; k++) {

                int offset = rows * k + k;
                int length = rows - k;
                double a;
                SquareMatrixAlgorithms.GenerateHouseholderReflection(store, offset, 1, length, out a);

                // apply P to the other columns of A
                for (int c = k + 1; c < cols; c++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, offset, 1, store, rows * c + k, 1, length);
                }

                // apply P to Q to accumulate the transform
                // since Q is rows X rows, its column count is just rows
                for (int c = 0; c < rows; c++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, offset, 1, qtStore, rows * c + k, 1, length);
                }

                // we are done with the Householder vector now, so we can zero the first column
                store[offset] = a;
                for (int i = 1; i < length; i++) {
                    store[offset + i] = 0.0;
                }

            }

        }

        // ** SVD **

        public static void ExtractSingularValues (double[] a, double[] b, double[] uStore, double[] vStore, int rows, int cols) {

            // n is the upper limit index
            int n = cols - 1;

            int count = 0;
            while (count < 32) {

                // move the lower boundary up as far as we can
                while ((n > 0) && (Math.Abs(b[n - 1]) <= (2.0E-16) * Math.Abs(a[n - 1] + a[n]))) {
                    b[n - 1] = 0.0;
                    if (a[n] < 0.0) {
                        a[n] = -a[n];
                        if (vStore != null) Blas1.dScal(-1.0, vStore, n * cols, 1, cols);
                    }
                    n--;
                    count = 0;
                }

                if (n == 0) return;

                // p is the upper limit index
                // find any zeros on the diagonal and reduce the problem
                int p = n;
                while (p > 0) {
                    if (Math.Abs(b[p - 1]) <= (2.0E-16) * (Math.Abs(a[p]) + Math.Abs(a[p - 1]))) break;
                    p--;
                    // do we need a better zero test?
                    if (a[p] == 0.0) {
                        // 
                        ChaseBidiagonalZero(a, b, p, n, uStore, rows);
                        p++;
                        break;
                    }
                }
                // diagonal zeros will cause the Golulb-Kahn step to reproduce the same matrix,
                // and eventually we will fail with non-convergence

                GolubKahn(a, b, p, n, uStore, vStore, rows, cols);
                count++;

            }

            throw new NonconvergenceException();

        }

        // A Golub-Kahn step in the SVD is like a Francis step in a QR eigendecomposition. After each step, the superdiagonal
        // elements of the bidiagonal matrix become smaller, making it converge to diagonal form.

        // Finding the singular values of A is equivilent to finding the eigenvalues of A^T A. For a bidiagonal matrix,
        // B^T B is tridiagonal. A QR step on B^T B would find the eigenvalues of its lower right 2 X 2 part, determine
        // the transformation that would bring (B^T B - e I)

        // |XX   |    |AA   |    |ABC  |    |XA0  |    |XX   |    |XX   |    |XX   |    |XX   |    |XX   |
        // | XX  |    |BBX  |    |0BC  |    | BB  |    | ABC |    | XA0 |    | XX  |    | XX  |    | XX  |
        // |  XX | -> |  XX | -> |  XX | -> | CCX | -> | 0BC | -> |  BB | -> |  ABC| -> |  XA0| -> |  XX |
        // |   XX|    |   XX|    |   XX|    |   XX|    |   XX|    |  CCX|    |  0BC|    |   BB|    |   AB|
        // |    X|    |    X|    |    X|    |    X|    |    X|    |    X|    |    X|    |   CC|    |   0B|

        private static void GolubKahn (double[] a, double[] b, int p, int n, double[] uStore, double[] vStore, int rows, int cols) {

            //WriteBidiagonal(a, b, 0.0, n);
            if (p == n) return;

            int m = n - 1;

            // find the eigenvalues of the bottom 2 X 2 of the tridiagonal matrix B^T * B, which has the form
            //   | a_m^2 + b_{m-1}^2      a_m b_m     |  =  | T_mm T_mn |
            //   |     a_m b_m          a_n^2 + b_m^2 |     | T_mn T_nn |
            // and pick the eigenvalue closer to the bottom value
            double T_mm = MoreMath.Pow2(a[m]);
            if (m > p) T_mm += MoreMath.Pow2(b[m - 1]);
            double T_nn = MoreMath.Pow2(a[n]) + MoreMath.Pow2(b[m]);
            double T_mn = a[m] * b[m];
            double dd = MoreMath.Hypot(T_mm - T_nn, 2.0 * T_mn);
            //double dd = Math.Sqrt(Math.Pow(T_mm - T_nn, 2) + 4.0 * Math.Pow(T_mn, 2));
            double e1 = (T_mm + T_nn + dd) / 2.0;
            double e2 = (T_mm + T_nn - dd) / 2.0;
            // shouldn't we compute e2 = T_mn / e1 instead?
            double e;
            if (Math.Abs(e1 - T_nn) < Math.Abs(e2 - T_nn)) {
                e = e1;
            } else {
                e = e2;
            }

            // find the givens rotation that would zero the first column of (B^T B - e I)

            // variables for cosine and sine of rotation angles
            double c, s;

            double u = MoreMath.Pow2(a[p]) - e;
            double v = a[p] * b[p];
            GenerateGivensRotation(ref u, ref v, out c, out s);

            // apply that rotation directly as B * G

            // a temporary variable to store the one non-biadiagonal element
            double t = 0.0;

            ApplyGivensRotation(c, s, ref a[p], ref b[p]);
            ApplyGivensRotation(c, s, ref t, ref a[p + 1]);

            if (vStore != null) {
                for (int j = 0; j < cols; j++) {
                    ApplyGivensRotation(c, s, ref vStore[p * cols + j], ref vStore[(p + 1) * cols + j]);
                }
            }

            //WriteBidiagonal(a, b, t, n);

            // "chase the bulge" down B
            for (int i = p; i < m; i++) {

                // G * B to zero sub-diagonal t
                GenerateGivensRotation(ref a[i], ref t, out c, out s);
                ApplyGivensRotation(c, s, ref b[i], ref a[i + 1]);
                ApplyGivensRotation(c, s, ref t, ref b[i + 1]);

                if (uStore != null) {
                    for (int j = 0; j < rows; j++) {
                        ApplyGivensRotation(c, s, ref uStore[j * rows + i], ref uStore[j * rows + i + 1]);
                    }
                }

                // B * G to zero (super+1)-diagonal t
                GenerateGivensRotation(ref b[i], ref t, out c, out s);
                ApplyGivensRotation(c, s, ref a[i + 1], ref b[i + 1]);
                ApplyGivensRotation(c, s, ref t, ref a[i + 2]);

                //WriteBidiagonal(a, b, t, n);

                if (vStore != null) {
                    for (int j = 0; j < cols; j++) {
                        ApplyGivensRotation(c, s, ref vStore[(i+1) * cols + j], ref vStore[(i + 2) * cols + j]);
                    }
                }

            }

            // G * B to remove one last sub-diagonal t
            GenerateGivensRotation(ref a[m], ref t, out c, out s);
            ApplyGivensRotation(c, s, ref b[m], ref a[n]);

            if (uStore != null) {
                for (int j = 0; j < rows; j++) {
                    ApplyGivensRotation(c, s, ref uStore[j * rows + m], ref uStore[j * rows + n]);
                }
            }

            //WriteBidiagonal(a, b, t, n);

        }

        // If there is a diagonal zero, then a Golub-Kahn step will leave the matrix unchanged. To deal with this, we "chase" the zero away
        // with a series of Givens rotations. Each rotation zeros a superdiagonal element on row with the diagonal zero by mixing with a lower row,
        // moving the non-zero element over by one column.

        // |XX   |    [XX   ]    [XX   ]    [XX   ]
        // | 0X  |    [ 00B ]    [ 000B]    [ 0000]
        // |  XX | -> [  AB ] -> [  XX ] -> [  XX ]
        // |   XX|    [   XX]    [   AB]    [   XX]
        // |    X|    [    X]    [    X]    [    A]

        private static void ChaseBidiagonalZero (double[] a, double[] b, int p, int n, double[] uStore, int rows) {

            Debug.Assert(a[p] == 0.0);

            double s, c;
            /*
            GenerateGivensRotation(ref a[q], ref b[p], out c, out s);

            // continue to chase the zero to the right by mixing with lower rows
            while (q < n) {
                double t = 0.0;
                ApplyGivensRotation(c, s, ref b[q], ref t);
                q++;
                GenerateGivensRotation(ref a[q], ref t, out c, out s);
            }
            */

            // a variable to hold the off-biadiagonal element
            // it starts off as b[p], which gets zeroed by the first operation
            double t = b[p]; b[p] = 0.0;

            for (int q = p + 1; q <= n; q++) {
                // t is non-zero
                GenerateGivensRotation(ref a[q], ref t, out c, out s);
                // t is zero
                if (q < n) ApplyGivensRotation(c, s, ref b[q], ref t);
                // t is non-zero

                // keep track of rotations
                if (uStore != null) {
                    for (int j = 0; j < rows; j++) {
                        ApplyGivensRotation(c, s, ref uStore[j * rows + q], ref uStore[j * rows + p]);
                    }
                }
            }
            

        }

        // [  c  s ] [ u ] = [ u' ]
        // [ -s  c ] [ v ] = [ v' ]

        private static void ApplyGivensRotation (double c, double s, ref double u, ref double v) {
            double t = u;
            u = c * t + s * v;
            v = c * v - s * t;
        }

        // A Givens rotation is a simple 2D rotation to make a vector horizontal, i.e. to zero one component.

        // [  c  s ] [ a ] = [ r ]
        // [ -s  c ] [ b ] = [ 0 ]

        // This routine replaces a and b with r and 0, and returns c and s.

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

        // The reslting diagonal and first superdiagonal elements are returned in a and b. The matrix itself is gradually
        // overwritten with the elements of the Householder reflections used. When we are done
        // the matrix is replaced by Householder vector components as follows:

        // APPPP
        // ABQQQ
        // ABCRR
        // ABCDS
        // ABCDE
        // ABCDE

        // where A, B, C, D, E are the Householder reflections applied from the left to zero columns, and P, Q, R, S are
        // the Householder reflections applied from the right to zero rows.

        public static void Bidiagonalize (double[] store, int rows, int cols, out double[] a, out double[] b) {

            //Debug.Assert(rows >= cols);
            a = new double[cols];
            b = new double[cols - 1];


            for (int k = 0; k < cols; k++) {

                // generate a householder reflection which, when applied from the left, will zero sub-diagonal elements in the kth column
                // use those elements to store the reflection vector; store the resulting diagonal element seperately
                SquareMatrixAlgorithms.GenerateHouseholderReflection(store, rows * k + k, 1, rows - k, out a[k]);

                // apply that reflection to all the columns; only subsequent ones are affected since preceeding ones
                // contain only zeros in the subdiagonal rows
                for (int c = k + 1; c < cols; c++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, rows * k + k, 1, store, rows * c + k, 1, rows - k);
                }

                if ((k + 1) < cols) {
                    // generate a Householder reflection which, when applied from the right, will zero (super+1)-diagonal elements in the kth row
                    // again, store the elements of the reflection vector in the zeroed matrix elements and store the resulting superdiagonal element seperately
                    SquareMatrixAlgorithms.GenerateHouseholderReflection(store, rows * (k + 1) + k, rows, cols - (k + 1), out b[k]);

                    // apply the reflection to all the rows; only subsequent ones are affected since the preceeding ones contain only zeros in the
                    // affected columns; the already-zeroed column elements are not distrubed because those columns are not affected
                    // this restriction is why we cannot fully diagonalize using this transform
                    for (int r = k + 1; r < rows; r++) {
                        SquareMatrixAlgorithms.ApplyHouseholderReflection(store, rows * (k + 1) + k, rows, store, rows * (k + 1) + r, rows, cols - (k + 1));
                    }
                }


            }

        }

        public static double[] AccumulateBidiagonalV (double[] store, int rows, int cols) {
            //  Q_1 * ... (Q_{n-1} * (Q_n * 1))
            double[] result = SquareMatrixAlgorithms.CreateUnitMatrix(cols);
            for (int k = cols - 2; k >= 0; k--) {
                // apply Householder reflection to each column from the left
                for (int j = k + 1; j < cols; j++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, (k + 1) * rows + k, rows, result, j * cols + (k + 1), 1, cols - k - 1);
                }
            }
            return (result);

        }

        public static double[] AccumulateBidiagonalU (double[] store, int rows, int cols) {
            // ((1 *  Q_n) * Q_{n-1}) ... * Q_1
            double[] result = SquareMatrixAlgorithms.CreateUnitMatrix(rows);

            // iterate over Householder reflections
            for (int k = cols - 1; k >= 0; k--) {
                // apply Householder reflection to each row from the right
                for (int j = k; j < rows; j++) {
                    SquareMatrixAlgorithms.ApplyHouseholderReflection(store, k * rows + k, 1, result, k * rows + j, rows, rows - k);
                }
            }

            return (result);
        }

        public static void SortValues (double[] values, double[] utStore, double[] vStore, int rows, int cols) {

            // this is a selection sort, an O(N^2) sort which requires fewer swaps than an insertion sort or an O(N ln N) sort

            // loop over ranks
            for (int i = 0; i < cols; i++) {

                // find the next largest value
                int j = i;
                double t = values[i];
                for (int k = i + 1; k < cols; k++) {
                    if (values[k] > t) {
                        j = k;
                        t = values[k];
                    }
                }

                // if necessary swap it with the current element
                if (j != i) {
                    Global.Swap(ref values[i], ref values[j]);
                    if (utStore != null) Blas1.dSwap(utStore, i, rows, utStore, j, rows, rows);
                    if (vStore != null) Blas1.dSwap(vStore, i * cols, 1, vStore, j * cols, 1, cols);
                }

            }

        }

    }

}
