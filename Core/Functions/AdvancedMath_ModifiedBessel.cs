using System;
using System.Diagnostics;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    // mostly direct and straightforward adaptations of the Bessel routines

    public static partial class AdvancedMath {

#if FUTURE

        /// <summary>
        /// Computes the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of I<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The modified Bessel function are essentially the Bessel functions with an imaginary argument.
        /// In particular, I<sub>n</sub>(x) = (-1)<sup>n</sup> J<sub>n</sub>(i x).</para>
        /// </remarks>
        /// <seealso cref="BesselJ(int,double)"/>
        public static double BesselI (int n, double x) {

            if (n < 0) {
                return (BesselI(-n, x));
            }

            if (x < 0) {
                double I = BesselI(n, -x);
                if (n % 2 != 0) I = -I;
                return (I);
            }

            if (x < 4.0 + 2.0 * Math.Sqrt(n)) {
                return (BesselI_Series(n, x));
            } else if (x > 20.0 + n * n / 4.0) {
                return (Math.Exp(x) * BesselI_Reduced_Asymptotic(n, x));
            } else {
                return (Math.Exp(x) * BesselI_Reduced_Miller(n, x));
            }
        }

        private static double BesselI_Series (double nu, double x) {

            if (x == 0.0) {
                if (nu == 0.0) {
                    return (1.0);
                } else {
                    return (0.0);
                }
            } else {
                double dI = Math.Pow(x / 2.0, nu) / AdvancedMath.Gamma(nu + 1.0);
                double I = dI;
                double xx = x * x / 4.0;
                for (int k = 1; k < Global.SeriesMax; k++) {
                    double I_old = I;
                    dI = dI * xx / (nu + k) / k;
                    I += dI;
                    if (I == I_old) {
                        return (I);
                    }
                }
                throw new NonconvergenceException();
            }
        }

        private static double BesselI_Reduced_Miller (int n, double x) {

            // use I_{k-1} = (2k/x) I_k + I_{k+1} to iterate down from a high k
            // use I_0 + 2 I_1 + 2 I_2 + 2 I_3 + 2 I_4 = e^x to normalize

            // this should be stable for all x, since I_k increases as k decreases

            // assume zero values at a high order
            double IP = 0.0;
            double I = 1.0;

            // do the renormalization sum as we go
            double sum = 0.0;

            // how much higher in order we start depends on how many significant figures we need
            // 30 is emperically enough for double precision
            // iterate down to I_n
            for (int k = n + 40; k > n; k--) {
                sum += I;
                double IM = (2 * k / x) * I + IP;
                IP = I;
                I = IM;
            }

            // remember I_n; we will renormalize it to get the answer
            double In = I;

            // iterate down to I_0
            for (int k = n; k > 0; k--) {
                sum += I;
                double IM = (2 * k / x) * I + IP;
                IP = I;
                I = IM;
            }

            // add the I_0 term to sum
            sum = 2.0 * sum + I;

            // since we are dividing by e^x, sum should be 1, so divide all terms by the actual sum

            return (In / sum);

        }

        private static double BesselI_Reduced_Asymptotic (double nu, double x) {

            // asmyptotic expansion for modified Bessel I
            // development is ~ (nu^2 / x)

            double mu = 4.0 * nu * nu;
            double xx = 8.0 * x;

            double dI = 1.0;
            double I = dI;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double I_old = I;
                int kk = 2 * k - 1;
                dI = - dI * (mu - kk * kk) / xx / k;
                I += dI;
                if (I == I_old) {
                    return (I / Math.Sqrt(2.0 * Math.PI * x));
                }
            }

            // can rewrite as dI = - dI * a * b / xx / k, where a+=2 and b -= 2 for each iteration

            throw new NonconvergenceException();
        }

        // Neuman series no good for computing K0; it relies on a near-perfect cancelation between the I0 term
        // and the higher terms to achieve an exponentially supressed small value

#endif

    }

}
