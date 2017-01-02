using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Matrices {

    // BLAS (Basic Linear Algebra Subsysem) operations
    
    // I tried making these IList<double> instead of double[], but performance degraded significantly. Arrays appers to be intrinsically faster.

    // BLAS Level 1 functions support O(N) operations with O(1) auxiluary storage requirements, basically scalar/vector operations
    // BLAS Level 2 functions support O(N^2) operations with O(N) auxiluary storage requirements, basically matrix/vector operations
    // BLAS Level 3 functions support O(N^3) operations with O(N^2) auxiluary storage requirements, basically matrix/matrix operations

    internal static class Blas1 {

        // y <- x
        public static void dCopy (double[] xStore, int xOffset, int xStride, double[] yStore, int yOffset, int yStride, int count) {
            if ((xStride == 1) && (yStride == 1)) {
                Array.Copy(xStore, xOffset, yStore, yOffset, count);
            } else {
                int n = 0;
                int x = xOffset;
                int y = yOffset;
                while (n < count) {
                    yStore[y] = xStore[x];
                    n++;
                    x += xStride;
                    y += yStride;
                }
            }
        }

        // x <-> y
        public static void dSwap (double[] xStore, int xOffset, int xStride, double[] yStore, int yOffset, int yStride, int count) {
            int n = 0;
            int x = xOffset;
            int y = yOffset;
            while (n < count) {
                double t = xStore[x];
                xStore[x] = yStore[y];
                yStore[y] = t;
                n++;
                x += xStride;
                y += yStride;
            }
        }

        // |x| = \sqrt( \sum_i x_i^2 )
        public static double dNrm2 (double[] store, int offset, int stride, int count) {

            // This is essentially Blue's algorithm. Divide the double range into numbers whose square would
            // overflow (large), numbers whose square would underflow (small), and numbers in betweeen (medium).
            // Sum each seperately, scaling large and small to avoid to avoid over- and under-flow before squaring.
            // Then combine sums when finished. In "usual" circumstances (medium numbers only), this requires
            // no more flops than a straight sum-of-squares, just a few extra comparisons.

            double small = 0.0;
            double medium = 0.0;
            double large = 0.0;

            int i = offset;
            for (int n = 0; n < count; n++) {
                double x = Math.Abs(store[i]);
                if (x != 0.0) {
                    if (x < smallLimit) {
                        small += MoreMath.Sqr(x / smallLimit);
                    } else if (x < largeLimit) {
                        medium += x * x;
                    } else {
                        large += MoreMath.Sqr(x / largeLimit);
                    }
                }
                i += stride;
            }
            if (large > 0.0) {
                return (MoreMath.Hypot(Math.Sqrt(large) * largeLimit, medium));
            } else if (small > 0.0) {
                return (MoreMath.Hypot(Math.Sqrt(small) * smallLimit, medium));
            } else {
                return (Math.Sqrt(medium));
            }
        }

        private static readonly double largeLimit = MoreMath.Pow(2.0, 508);

        private static readonly double smallLimit = 1.0 / largeLimit;

        // sum_i |x_i|
        public static double dNrm1 (double[] store, int offset, int stride, int count) {
            double m = 0.0;
            int n = 0;
            int i = offset;
            while (n < count) {
                m += Math.Abs(store[i]);
                n++;
                i += stride;
            }
            return (m);
        }

        // a^T b
        public static double dDot (double[] aStore, int aOffset, int aStride, double[] bStore, int bOffset, int bStride, int count) {
            double m = 0.0;
            int n = 0;
            int a = aOffset;
            int b = bOffset;
            while (n < count) {
                m += aStore[a] * bStore[b];
                n++;
                a += aStride;
                b += bStride;
            }
            return (m);
        }

        // x <- a x
        public static void dScal (double alpha, double[] store, int offset, int stride, int count) {
            int n = 0;
            int i = offset;
            while (n < count) {
                store[i] *= alpha;
                n++;
                i += stride;
            }
        }

        // y <- a x + y
        public static void dAxpy (double alpha, double[] xStore, int xOffset, int xStride, double[] yStore, int yOffset, int yStride, int count) {
            int n = 0;
            int x = xOffset;
            int y = yOffset;
            while (n < count) {
                yStore[y] += alpha * xStore[x];
                n++;
                x += xStride;
                y += yStride;
            }
        }

    }

    internal static class Blas2 {

        // y <- A x + y

        public static void dGemv (
            double[] aStore, int aOffset, int aRowStride, int aColStride,
            double[] xStore, int xOffset, int xStride,
            double[] yStore, int yOffset, int yStride,
            int rows, int cols
        ) {

            int aIndex = aOffset;
            int yIndex = yOffset;
            for (int n = 0; n < rows; n++) {
                yStore[yIndex] += Blas1.dDot(aStore, aIndex, aColStride, xStore, xOffset, xStride, cols);
                aIndex += aRowStride;
                yIndex += yStride;
            }

        }

        // Solve a triangular system
        //   aIsUpper indicates whether A is upper/right or lower/left
        //   aIsUnit indicates whether A's diagonal elements are all 1 or not
        // The relevent elements of xStore are overwritten by the solution vector

        public static void dTrsv (
            bool aIsUpper, bool aIsUnit, double[] aStore, int aOffset, int aRowStride, int aColStride,
            double[] xStore, int xOffset, int xStride,
            int count
        ) {

            int aDiagonalStride = aColStride + aRowStride;

            // If A is upper triangular, start at the bottom/right and reverse the direction of progression.
            // This does not affect the passed in variables because integers are pass-by-value.
            if (aIsUpper) {
                aOffset = aOffset + (count - 1) * aDiagonalStride;
                aRowStride = -aRowStride;
                aColStride = -aColStride;
                aDiagonalStride = -aDiagonalStride;
                xOffset = xOffset + (count - 1) * xStride;
                xStride = -xStride;
            }

            // The index of the first row element to be subtracted in the numerator.
            int aRowIndex = aOffset;
            // The index of the x element to be solved for.
            int xIndex = xOffset;

            // We hoist the aIsUnit test outside the loop and pay the cost of a little repeated code in
            // order to avoid an unnecessary per-loop test. It's possible the compiler could do this for us.
            if (aIsUnit) {
                for (int n = 0; n < count; n++) {
                    xStore[xIndex] -= Blas1.dDot(aStore, aRowIndex, aColStride, xStore, xOffset, xStride, n);
                    aRowIndex += aRowStride;
                    xIndex += xStride;
                }
            } else {
                int aDiagonalIndex = aOffset;
                for (int n = 0; n < count; n++) {
                    xStore[xIndex] -= Blas1.dDot(aStore, aRowIndex, aColStride, xStore, xOffset, xStride, n);
                    xStore[xIndex] /= aStore[aDiagonalIndex];
                    aRowIndex += aRowStride;
                    aDiagonalIndex += aDiagonalStride;
                    xIndex += xStride;
                }
            }

        }
    }

}
