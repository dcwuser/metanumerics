using System;

namespace Meta.Numerics.Extended {

    /// <summary>
    /// Contains method for computed advanced functions to double double accuracy.
    /// </summary>
    public static class AdvancedDoubleDoubleMath {

        /// <summary>
        /// Computes the error function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of erf(x).</returns>
        public static DoubleDouble Erf (DoubleDouble x) {
            if (x < 0.0) {
                return (-Erf(-x));
            } else if (x < 4.0) {
                return (Erf_Series(x));
            } else {
                return (1.0 - Erfc_ContinuedFraction(x));
            }
        }

        /// <summary>
        /// Computes the complementatry error function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of erfc(x).</returns>
        public static DoubleDouble Erfc (DoubleDouble x) {
            if (x < 4.0) {
                return (1.0 - Erf_Series(x));
            } else {
                return (Erfc_ContinuedFraction(x));
            }
        }

        private static DoubleDouble Erf_Series (DoubleDouble x) {
            DoubleDouble mx2 = -(x * x);
            DoubleDouble t = DoubleDouble.One;
            DoubleDouble f = t;
            for (int k = 1; k < 100; k++) {
                DoubleDouble f_old = f;
                t *= mx2 / k;
                f += t / (2 * k + 1);
                if (f == f_old) {
                    return (2.0 / DoubleDouble.Sqrt(DoubleDouble.Pi)) * x * f;
                }
            }
            throw new InvalidOperationException();
        }

        private static DoubleDouble Erfc_ContinuedFraction (DoubleDouble x) {
            DoubleDouble x2 = x * x;
            double aa = 1.0;
            DoubleDouble bb = x2 + 0.5;
            DoubleDouble D = 1.0 / bb;
            DoubleDouble Df = aa / bb;
            DoubleDouble f = Df;
            for (int k = 1; k < 100; k++) {
                DoubleDouble f_old = f;
                aa = -k * (k - 0.5);
                bb += 2.0;
                D = 1.0 / (bb + aa * D);
                Df = (bb * D - 1.0) * Df;
                f += Df;
                if (f == f_old) {
                    DoubleDouble e = DoubleDouble.Exp(-x2);
                    DoubleDouble g = x / DoubleDouble.Sqrt(DoubleDouble.Pi);
                    return (e * f * g);
                }
            }
            throw new InvalidOperationException();
        }
    }

}
