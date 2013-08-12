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

        // |x| = sqrt( sum_i x_i^2 )
        public static double dNrm2 (double[] store, int offset, int stride, int count) {
            double m = 0.0;
            int n = 0;
            int i = offset;
            while (n < count) {
                double x = store[i];
                m += x * x;
                n++;
                i += stride;
            }
            return (Math.Sqrt(m));
        }

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

    }

}
