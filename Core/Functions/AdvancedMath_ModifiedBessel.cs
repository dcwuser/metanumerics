using System;
using System.Diagnostics;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    // mostly direct and straightforward adaptations of the Bessel routines

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns></returns>
        public static double BesselI (int n, double x) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (x < 0) throw new ArgumentOutOfRangeException("x");

            if (x < 4.0 + 2.0 * Math.Sqrt(n)) {
                return (BesselI_Series(n, x));
            } else if (x > 20.0 + n * n / 4.0) {
                return (Math.Exp(x) * BesselI_Reduced_Asymptotic(n, x));
            } else {
                throw new NotImplementedException();
            }
        }

        private static double BesselI_Series (double nu, double x) {
            double dJ = Math.Pow(0.5 * x, nu) / AdvancedMath.Gamma(nu + 1.0);
            double J = dJ;
            double xx = x * x / 4.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double J_old = J;
                dJ = dJ * xx / (nu + k) / k;
                J += dJ;
                if (J == J_old) {
                    return (J);
                }
            }
            throw new NonconvergenceException();
        }

        private static double BesselI_Reduced_Asymptotic (double nu, double x) {
            throw new NotImplementedException();
        }

    }

}
