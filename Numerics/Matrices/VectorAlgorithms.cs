using System;

namespace Meta.Numerics.Matrices {

    internal static class VectorAlgorithms {

        public static double[] Copy (double[] store, int offset, int stride, int dimension) {
            double[] copy = new double[dimension];
            Blas1.dCopy(store, offset, stride, copy, 0, 1, dimension);
            return (copy);
        }

        public static double[] Add (double[] aStore, int aOffset, int aStride, double[] bStore, int bOffset, int bStride, int dimension) {
            double[] store = new double[dimension];
            Blas1.dCopy(aStore, aOffset, aStride, store, 0, 1, dimension);
            Blas1.dAxpy(1.0, bStore, bOffset, bStride, store, 0, 1, dimension);
            return (store);
        }

        public static double[] Subtract (double[] aStore, int aOffset, int aStride, double[] bStore, int bOffset, int bStride, int dimension) {
            double[] store = new double[dimension];
            Blas1.dCopy(aStore, aOffset, aStride, store, 0, 1, dimension);
            Blas1.dAxpy(-1.0, bStore, bOffset, bStride, store, 0, 1, dimension);
            return (store);
        }

        public static double[] Multiply (double alpha, double[] store, int offset, int stride, int dimension) {
            double[] product = new double[dimension];
            Blas1.dAxpy(alpha, store, offset, stride, product, 0, 1, dimension);
            return (product);
        }


        // A Householder reflection matrix is a rank-1 update to the identity.
        //   P = I - b v v^T
        // Unitarity requires b = 2 / |v|^2. To anihilate all but the first components of vector x,
        //   P x = a e_1
        // we must have
        //   v = x +/- |x| e_1
        // that is, all elements the same except the first, from which we have either added or subtracted |x|. This makes
        //   a = -/+ |x|
        // There are two way to handle the sign. One is to simply choose the sign that avoids cancelation when calculating v_1,
        // i.e. + for positive x_1 and - for negative x_1. This works fine, but makes a negative for positive x_1, which is
        // weird-looking (1 0 0 gets turned into -1 0 0). An alternative is to choose the negative sign even for positive x_1,
        // but to avoid cancelation write
        //   v_1 = x_1 - |x| = ( x_1 - |x| ) ( x_1 + |x|) / ( x_1 + |x|) = ( x_1^2 - |x|^2 ) / ( x_1 + |x| )
        //       = - ( x_2^2 + \cdots + x_n^2 ) / ( x_1 + |x| )
        // We have now moved to the second method. Note that v is the same as x except for the first element.

        public static void GenerateHouseholderReflection (double[] store, int offset, int stride, int count, out double a) {

            double x0 = store[offset];

            // Compute |x| and u_0.
            double xm, u0;
            if (x0 < 0.0) {
                xm = Blas1.dNrm2(store, offset, stride, count);
                u0 = x0 - xm;
            } else {
                // This method of computing ym and xm does incur and extra square root compared to naively computing x_2 + \cdots + x_n^2,
                // but doing it this way allows us to use dNrm2's over/under-flow prevention when we have large/small elements.
                double ym = Blas1.dNrm2(store, offset + stride, stride, count - 1);
                xm = MoreMath.Hypot(x0, ym);
                // Writing ym / (x0 + xm) * ym instead of ym * ym / (x0 + xm) prevents over/under-flow for large/small ym. Note this will
                // be 0 / 0 = NaN when xm = 0.
                u0 = -ym / (x0 + xm) * ym;
            }

            // Set result element.
            a = xm;

            // If |x| = 0 there is nothing to do; we could have done this a little earlier but the compiler requires us to set a before returning.
            if (xm == 0.0) return;

            // Set the new value of u_0
            store[offset] = u0;

            // Rescale to make b = 1.
            double um = Math.Sqrt(xm * Math.Abs(u0));
            if (um > 0.0) Blas1.dScal(1.0 / um, store, offset, stride, count);

        }


        // Apply a Householder transfrom defined by v to the vector x (which may be
        // a column of a matrix, if H is applied from the left, or a row of a matrix
        // if H is applied from the right). On exit v is unchanged, x is changed.
        // We have H = I - \beta v v^T, so H x = x - (\beta v^T x) v.

        public static void ApplyHouseholderReflection (
            double[] uStore, int uOffset, int uStride,
            double[] yStore, int yOffset, int yStride,
            int count
        ) {
            double s = Blas1.dDot(uStore, uOffset, uStride, yStore, yOffset, yStride, count);
            Blas1.dAxpy(-s, uStore, uOffset, uStride, yStore, yOffset, yStride, count);
        }

        public static void Zero (double[] store, int offset, int stride, int count) {
            if (stride == 1 && count > 64) {
                // Information on internet indicates that Array.Clear is faster than setting
                // individual elements to zero for size larger than ~75
                Array.Clear(store, offset, count);
            } else {
                int n = 0;
                int i = offset;
                while (n < count) {
                    store[i] = 0.0;
                    n++;
                    i += stride;
                }
            }
        }

    }

}
