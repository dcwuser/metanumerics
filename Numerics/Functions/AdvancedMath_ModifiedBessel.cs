using System;
using System.Diagnostics;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    // mostly direct and straightforward adaptations of the Bessel routines

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes modified cylindrical Bessel functions.
        /// </summary>
        /// <param name="nu">The order, which must be non-negative.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The values of I, I', K, and K' for the given order and argument.</returns>
        /// <remarks>
        /// <para>The modified Bessel functions fulfill a differential equation similar to the Bessel differential equation.</para>
        /// <img src="../images/ModifiedBesselODE.png" />
        /// </remarks>
        public static SolutionPair ModifiedBessel (double nu, double x) {
            if (nu < 0.0) throw new ArgumentOutOfRangeException(nameof(nu));
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

            if (x == 0.0) {
                return (Bessel_Zero(nu));
            } else if (x < 2.0) {
                
                // use series to determine I and I'
                double I, IP;
                ModifiedBesselI_Series(nu, x, out I, out IP);

                // use series to determine K and K' at -1/2 <= mu <= 1/2 that is an integer offset from nu
                int n = (int) Math.Round(nu);
                double mu = nu - n;

                double K, K1;
                ModifiedBesselK_Series(mu, x, out K, out K1);
                double KP = (mu / x) * K - K1;

                // recurr K, K' upward to order nu
                ModifiedBesselK_RecurrUpward(mu, x, ref K, ref KP, n);

                // return the result
                return (new SolutionPair(I, IP, K, KP));

            } else if (x > 32.0 + nu * nu / 2.0) {

                double sI, sIP, sK, sKP;
                ModifiedBessel_Asymptotic(nu, x, out sI, out sIP, out sK, out sKP);
                double e = Math.Exp(x);
                return (new SolutionPair(e * sI, e * sIP, sK / e, sKP / e));

            } else {

                // find 0 <= mu < 1 with same fractional part as nu
                // this is necessary because CF2 does not produce K with good accuracy except at very low orders
                int n = (int) Math.Floor(nu);
                double mu = nu - n;

                // compute K, K' at this point (which is beyond the turning point because mu is small) using CF2
                double sK, g;
                ModifiedBessel_CF_K(mu, x, out sK, out g);
                double sKP = g * sK;

                // recurr upward to order nu
                ModifiedBesselK_RecurrUpward(mu, x, ref sK, ref sKP, n);

                // determine I'/I at the desired point
                double f = ModifiedBessel_CF1(nu, x);

                // Use the Wronskian relationship K I' - I K' = 1/x to determine I and I' separately
                double sI = 1.0 / (f * sK - sKP) / x;
                double sIP = f * sI;

                double e = Math.Exp(x);
                return (new SolutionPair(e * sI, e * sIP, sK / e, sKP / e));

            }

        }


        // Unlike the analogous recurrence for J, Y, this recurrence, as written, is only good for K, not I
        // I and e^{i\pi\nu} K share the same recurrence, but we have written it here so as to avoid the e^{i\pi\nu} factor

        private static void ModifiedBesselK_RecurrUpward (double mu, double x, ref double K, ref double KP, int n) {

            for (int i = 0; i < n; i++) {
                double t = K;
                K = (mu / x) * K - KP;
                mu += 1.0;
                KP = -(t + (mu / x) * K);
            }

        }

        /// <summary>
        /// Computes the regular modified cylindrical Bessel function.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of I<sub>&#x3BD;</sub>(x).</returns>
        /// <remarks>
        /// <para>The modified Bessel functions appear as the solutions of hyperbolic differential equations with
        /// cylindrical or circular symmetry, for example the conduction of heat through a cylindrical pipe.</para>
        /// <para>The regular modified Bessel functions are related to the Bessel functions with pure imaginary arguments.</para>
        /// <img src="../images/BesselIBesselJRelation.png" />
        /// <para>The regular modified Bessel functions increase monotonically and exponentially from the origin.</para>
        /// </remarks>
        /// <seealso cref="ModifiedBesselK"/>
        public static double ModifiedBesselI (double nu, double x) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

            // Reflect negative nu to positive nu.
            if (nu < 0.0) {
                if (Math.Round(nu) == nu) {
                    return (ModifiedBesselI(-nu, x));
                } else {
                    return (ModifiedBesselI(-nu, x) + 2.0 / Math.PI * MoreMath.SinPi(-nu) * ModifiedBesselK(-nu, x));
                }
            }

            if (x == 0.0) {
                if (nu == 0.0) {
                    return (1.0);
                } else {
                    return (0.0);
                }
            } else if (x <= 2.0 * Math.Sqrt(nu + 1.0)) {
                // Close to the origin, use the power series.
                return (ModifiedBesselI_Series(nu, x));
            } else if (x >= 32.0 + 0.5 * nu * nu) {
                // Far from the origin, use the asymptotic expansion.
                double sI, sIP, sK, sKP;
                ModifiedBessel_Asymptotic(nu, x, out sI, out sIP, out sK, out sKP);
                return (Math.Exp(x) * sI);
            } else {
                // In the intermediate region, we will use CF1 + CF2 in a slightly different way than for the Bessel functions

                // find 0 <= mu < 1 with same fractional part as nu
                // this is necessary because CF2 does not produce K with good accuracy except at very low orders
                int n = (int) Math.Floor(nu);
                double mu = nu - n;

                // compute K, K' at this point (which is beyond the turning point because mu is small) using CF2
                double sK, g;
                ModifiedBessel_CF_K(mu, x, out sK, out g);
                double sKP = g * sK;

                ModifiedBesselK_RecurrUpward(mu, x, ref sK, ref sKP, n);

                // determine I'/I at the desired point
                double f = ModifiedBessel_CF1(nu, x);

                // use the Wronskian to determine I from (I'/I), (K'/K), and K
                double sI = 1.0 / (f * sK - sKP) / x;

                return (Math.Exp(x) * sI);
            }

        }

        /// <summary>
        /// Computes the irregular modified cylindrical Bessel function.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of K<sub>&#x3BD;</sub>(x).</returns>
        /// <remarks>
        /// <para>The modified Bessel functions are related to the Bessel functions with pure imaginary arguments.</para>
        /// <para>The irregular modified Bessel function decreases monotonically and exponentially from the origin.</para>
        /// </remarks>
        public static double ModifiedBesselK (double nu, double x) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

            if (nu < 0.0) return (ModifiedBesselK(-nu, x));

            if (x == 0.0) {
                return (Double.NegativeInfinity);
            } else if (x < 2.0) {
                // For small x, determine the value at small mu using the series and recurr up.
                int n = (int) Math.Round(nu);
                double mu = nu - n;

                double K, K1;
                ModifiedBesselK_Series(mu, x, out K, out K1);
                if (n == 0) {
                    return (K);
                } else if (n == 1) {
                    return (K1);
                } else {
                    double KP = (mu / x) * K - K1;
                    ModifiedBesselK_RecurrUpward(mu, x, ref K, ref KP, n);
                    return (K);
                }
            } else if (x > 32.0 + nu * nu / 2.0) {
                // For large x, use the asymptotic series.
                double sI, sIP, sK, sKP;
                ModifiedBessel_Asymptotic(nu, x, out sI, out sIP, out sK, out sKP);
                return(Math.Exp(-x) * sK);
            } else {
                // In the intermediate region, determine the value at small mu using continued fraction and recurr up.
                int n = (int) Math.Round(nu);
                double mu = nu - n;

                double sK, g;
                ModifiedBessel_CF_K(mu, x, out sK, out g);
                double sKP = g * sK;
                ModifiedBesselK_RecurrUpward(mu, x, ref sK, ref sKP, n);
                return (Math.Exp(-x) * sK);
            }

        }

        // series near the origin; this is entirely analogous to the Bessel series near the origin
        // it has a corresponding radius of rapid convergence, x < 4 + 2 Sqrt(nu)

        // This is exactly the same as BesselJ_Series with xx -> -xx.
        // We could even factor this out into a common method with an additional parameter.

        private static void ModifiedBesselI_Series (double nu, double x, out double I, out double IP) {

            Debug.Assert(x > 0.0);

            double x2 = 0.5 * x;
            double xx = x2 * x2;
            double dI = AdvancedMath.PowerOverFactorial(x2, nu);
            I = dI;
            IP = nu * dI;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double I_old = I;
                double IP_old = IP;
                dI *= xx / (k * (nu + k));
                I += dI;
                IP += (2 * k + nu) * dI;
                if ((I == I_old) && (IP == IP_old)) {
                    IP = IP / x;
                    return;
                }
            }

            throw new NonconvergenceException();
        }

        private static double ModifiedBesselI_Series (double nu, double x) {

            double x2 = x / 2.0;
            double dI = Math.Pow(x2, nu) / AdvancedMath.Gamma(nu + 1.0);
            double I = dI;
            double xx = x2 * x2;
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

        // good for x > 32 + nu^2 / 2
        // this returns scaled values; I = e^(x) sI, K = e^(-x) sK

        // Asymptotic expansions
        //   I_{\nu}(x) = \frac{e^x}{\sqrt{2 \pi x}}  \left[ 1 - \frac{\mu-1}{8x} + \frac{(\mu-1)(\mu-9)}{2! (8x)^2} + \cdots \right]
        //   K_{\nu}(x) = \sqrt{\frac{\pi}{2x}}e^{-x} \left[ 1 + \frac{\mu-1}{8x} + \frac{(\mu-1)(\mu-9)}{2! (8x)^2} + \cdots \right]
        // where \mu = 4 \nu^2. Derivatives
        //   I_{\nu}'(x) = \frac{e^x}{\sqrt{2 \pi x}}   \left[ 1 - \frac{\mu+3}{8x} + \frac{(\mu-1)(\mu+15)}{2! (8x)^2} + \cdots \right]
        //   K_{\nu}'(x) = -\sqrt{\frac{\pi}{2x}}e^{-x} \left[ 1 + \frac{\mu+3}{8x} + \frac{(\mu-1)(\mu+15)}{2! (8x)^2} + \cdots \right]
        // Note series differ only by alternating vs. same sign of terms.

        private static void ModifiedBessel_Asymptotic (double nu, double x, out double sI, out double sIP, out double sK, out double sKP) {

            // precompute some values we will use
            double mu = 4.0 * nu * nu;
            double xx = 8.0 * x;

            // initialize the series
            double I = 1.0; double IP = 1.0;
            double K = 1.0; double KP = 1.0;

            // intialize the term value
            double t = 1.0;

            for (int k = 1; k < Global.SeriesMax; k++) {

                double I_old = I; double K_old = K;
                double IP_old = IP; double KP_old = KP;

                // determine next term values
                int k2 = 2 * k;
                t /= k * xx;
                double tp = (mu + (k2 * k2 - 1)) * t;
                //t *= (mu - (k2-1) * (k2-1));
                t *= (2.0 * (nu - k) + 1.0) * (2.0 * (nu + k) - 1.0);

                // add them, with alternating-sign for I series and same-sign for K series
                if (k % 2 == 0) {
                    I += t;
                    IP += tp;
                } else {
                    I -= t;
                    IP -= tp;
                }
                K += t;
                KP += tp;

                // check for convergence
                if ((I == I_old) && (K == K_old) && (IP == IP_old) && (KP == KP_old)) {
                    double fI = Math.Sqrt(Global.TwoPI * x);
                    double fK = Math.Sqrt(Global.HalfPI / x);
                    sI = I / fI;
                    sIP = IP / fI;
                    sK = K * fK;
                    sKP = -KP * fK;
                    return;
                }

            }

            throw new NonconvergenceException();

        }

        // compute I'/I via continued fraction
        // this is analogous to CF1 for normal Bessel functions
        // it converges quickly for x < nu and slowly for x > nu

        private static double ModifiedBessel_CF1 (double nu, double x) {
            double q = 2.0 * (nu + 1.0);
            double a = 1.0;
            double b = (q / x);
            double D = 1.0 / b;
            double Df = a / b;
            double f = nu / x + Df;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double f_old = f;

                q += 2.0;
                a = 1.0;
                b = q / x;
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;

                if (f == f_old) {
                    return (f);
                }
            }

            throw new NonconvergenceException();
        }

        // the continued fraction here is 

        private static void ModifiedBessel_CF_K (double nu, double x, out double K, out double g) {

            double a = 1.0;
            double b = 2.0 * (1.0 + x);

            double D = 1.0 / b;
            double Df = a / b;
            double f = Df;

            // Thompson-Barnett do normalizing sum using these auxiliary variables (see Numerical Recipes)

            double C = (0.5 - nu) * (0.5 + nu);
            double q0 = 0.0;
            double q1 = 1.0;
            double Q = C * q1;
            double S = 1.0 + Q * Df;

            for (int k = 2; k < Global.SeriesMax; k++) {

                double f_old = f;

                // Lenz method for z1/z0

                a = -(k - nu - 0.5) * (k + nu - 0.5);
                b += 2.0;

                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;

                // Thompson-Barnett sum trick

                double S_old = S;
                
                C = -a / k * C;
                double q2 = (q0 - (b - 2.0) * q1) / a;
                Q += C * q2;
                S += Q * Df;

                if ((S == S_old) && (f == f_old)) {
                    g = -((x + 0.5) + (nu + 0.5) * (nu - 0.5) * f) / x;
                    K = Math.Sqrt(Math.PI / 2.0 / x) / S;
                    return;
                }

                q0 = q1; q1 = q2;

            }

            throw new NonconvergenceException();
        }

        // use Wronskian I' K - I K' = 1/x, plus values of I'/I, K'/K, and K, to compute I, I', K, and K'

        // Temme's series for K0 and K1 for small nu and small x
        // applies only for -1/2 <= nu <= 1/2
        // accurate to all but last couple digits for x < 2; converges further out but accuracy deteriorates rapidly
        // initial constant determination shared with BesselY_Series; we should factor out the common parts

        private static void ModifiedBesselK_Series (double nu, double x, out double K0, out double K1) {

            Debug.Assert((-0.5 <= nu) && (nu <= 0.5));

            if (x == 0.0) {
                K0 = Double.PositiveInfinity;
                K1 = Double.PositiveInfinity;
                return;
            }

            // polynomial variable and coefficient

            double xx = x * x / 4.0;
            double cx = 1.0;

            // determine initial p, q

            double GP = Gamma(1.0 + nu);
            double GM = Gamma(1.0 - nu);

            double y = x / 2.0;
            double z = Math.Pow(y, nu);

            double p = 0.5 / z * GP;
            double q = 0.5 * z * GM;

            // determine initial f; this is rather complicated

            double ln = -Math.Log(y);

            double s = nu * ln;
            double C1 = Math.Cosh(s);
            double C2 = (s == 0.0) ? 1.0 : Math.Sinh(s) / s;

            double G1 = (Math.Abs(nu) < 0.25) ? NewG(1.0 - nu, 2.0 * nu) : (1.0 / GM - 1.0 / GP) / (2.0 * nu);
            double G2 = (1.0 / GM + 1.0 / GP) / 2.0;
            double F = (nu == 0.0) ? 1.0 : (nu * Math.PI) / MoreMath.SinPi(nu);

            double f = F * (C1 * G1 + C2 * G2 * ln);

            // run the series

            K0 = cx * f;
            K1 = cx * p;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double K0_old = K0;
                double K1_old = K1;

                cx *= xx / k;

                double kpnu = k + nu;
                double kmnu = k - nu;
                f = (k * f + p + q) / (kpnu * kmnu);
                p = p / kmnu;
                q = q / kpnu;
                double h = p - k * f;

                K0 += cx * f;
                K1 += cx * h;

                if ((K0 == K0_old) && (K1 == K1_old)) {
                    K1 = 2.0 / x * K1;
                    return;
                }
            }
            throw new NonconvergenceException();

        }

        /// <summary>
        /// Computes the modified Bessel function of the first kind.
        /// </summary>
        /// <param name="n">The order.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of I<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The modified Bessel function are essentially the Bessel functions with an imaginary argument.
        /// In particular, I<sub>n</sub>(x) = (-1)<sup>n</sup> J<sub>n</sub>(i x).</para>
        /// </remarks>
        /// <seealso cref="BesselJ(int,double)"/>
        public static double ScaledBesselI (int n, double x) {

            if (n < 0) {
                return (ScaledBesselI(-n, x));
            }

            if (x < 0) {
                double I = ScaledBesselI(n, -x);
                if (n % 2 != 0) I = -I;
                return (I);
            }

            if (x <= 2.0 * Math.Sqrt(n + 1)) {
                return (Math.Exp(-x) * BesselI_Series(n, x));
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
                double xh = 0.5 * x;
                double dI = 1.0;
                double I = dI;
                double xx = xh * xh;
                for (int k = 1; k < Global.SeriesMax; k++) {
                    double I_old = I;
                    dI = dI * xx / (k * (nu + k));
                    I += dI;
                    if (I == I_old) {
                        return (PowerOverFactorial(xh, nu) * I);
                    }
                }
                throw new NonconvergenceException();
            }
        }

        private static double BesselI_Reduced_Miller (int n, double x) {

            // use I_{k-1} = (2k/x) I_k + I_{k+1} to iterate down from a high k
            // use I_0 + 2 I_1 + 2 I_2 + 2 I_3 + 2 I_4 = e^x to normalize

            // this should be stable for all x, since I_k increases as k decreases

            int m = Bessel_Miller_Limit(n, x);

            // assume zero values at a high order
            double IP = 0.0;
            double I = 1.0;

            // do the renormalization sum as we go
            double sum = 0.0;

            // iterate down to 0
            double In = Double.NaN;
            for (int k = m; k > 0; k--) {
                if (k == n) {
                    In = I;
                }
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

    }

}
