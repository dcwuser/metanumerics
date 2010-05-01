using System;
using System.Diagnostics;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    // mostly direct and straightforward adaptations of the Bessel routines

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the regular modified cynlindrical Bessel function.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of I<sub>&#x3BD;</sub>(x).</returns>
        /// <remarks>
        /// <para>The modified Bessel functions appear as the solutions of hyperbolic differential equations with
        /// cylindrical or circular symmetry, for example the conduction of heat through a cylindrical pipe.</para>
        /// <para>The regular modified Bessel functions are related to the Bessel fuctions with pure imaginary arguments.</para>
        /// <img src="../images/BesselIBesselJRelation.png" />
        /// <para>The regular modified Bessel functions increase monotonically and exponentially from the origin.</para>
        /// </remarks>
        /// <seealso cref="ModifiedBesselK"/>
        public static double ModifiedBesselI (double nu, double x) {

            // reflect negative nu to positive nu
            if (nu < 0.0) {
                if (Math.Round(nu) == nu) {
                    return (ModifiedBesselI(-nu, x));
                } else {
                    return (ModifiedBesselI(-nu, x) + 2.0 / Math.PI * Sin(0.0, -nu) * ModifiedBesselK(-nu, x));
                }
            }

            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            if (x < 4.0 + 2.0 * Math.Sqrt(nu)) {
                // close to the origin, use series
                return (ModifiedBesselI_Series(nu, x));
            } else if (x > 32.0 + nu * nu / 2.0) {
                // far from the origin, use asymptotic expansion
                double sI, sK;
                ModifiedBessel_Asymptotic(nu, x, out sI, out sK);
                return (Math.Exp(x) * sI);
            } else {
                // beyond the turning point, we will use CF1 + CF2 in a slightly different way than for the Bessel functions

                // find 0 <= mu < 1 with same fractional part as nu
                // this is necessary because CF2 does not produce K with good accuracy except at very low orders
                int n = (int) Math.Floor(nu);
                double mu = nu - n;

                // compute K, K' at this point (which is beyond the turning point because mu is small) using CF2
                //double K, KP, g;
                double K, g;
                ModifiedBessel_CF_K(mu, x, out K, out g);
                //KP = g * K;

                // recurse K, K' upward to nu
                for (int k = 0; k < n; k++) {
                    //double Kp1 = -KP + (mu / x) * K;
                    double Kp1 = (mu / x - g) * K;
                    mu += 1.0;
                    //KP = -K - (mu / x) * Kp1;
                    g = -(mu / x + K / Kp1);
                    K = Kp1;
                }

                // determine I'/I at the desired point
                double f = ModifiedBessel_CF1(nu, x);

                // use the Wronskian to determine I from (I'/I), (K'/K), and K
                //double I = 1.0 / (f * K - KP) / x;
                double I = 1.0 / (f - g) / K / x;

                return (I);

            }

        }

        /// <summary>
        /// Computes the irregular modified cynlindrical Bessel function.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of K<sub>&#x3BD;</sub>(x).</returns>
        /// <remarks>
        /// <para>The modified Bessel functions are related to the Bessel fuctions with pure imaginary arguments.</para>
        /// <para>The irregular modified Bessel function decreases monotonically and exponentially from the origin.</para>
        /// </remarks>
        public static double ModifiedBesselK (double nu, double x) {

            if (nu < 0.0) return (ModifiedBesselI(-nu, x));

            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            if (x > 32.0 + nu * nu / 2.0) {
                // for large x, use asymptotic series
                double sI, sK;
                ModifiedBessel_Asymptotic(nu, x, out sI, out sK);
                return(Math.Exp(-x) * sK);
            } else {
                // otherwise, reduce to a problem with -1/2 <= mu <= 1/2 and recurse up

                // determine mu
                int n = (int)Math.Floor(nu + 0.5);
                double mu = nu - n;

                // determine K and K' for mu
                double K, g;
                if (x < 2.0) {
                    double Kp1;
                    ModifiedBesselK_Series(mu, x, out K, out Kp1);
                    g = -Kp1 / K + mu / x;
                } else {
                    ModifiedBessel_CF_K(mu, x, out K, out g);
                }

                // recurse upward to nu
                for (int k = 0; k < n; k++) {
                    double Kp1 = (mu / x - g) * K;
                    mu += 1.0;
                    g = -(mu / x + K / Kp1);
                    K = Kp1;
                }

                return (K);

            }

        }

        // series near the origin; this is entirely analogous to the Bessel series near the origin
        // it has a corresponding radius of rapid convergence, x < 4 + 2 Sqrt(nu)

        private static double ModifiedBesselI_Series (double nu, double x) {

            if (x == 0.0) {
                // handle x = 0 specially
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

        // good for x > 32 + nu^2 / 2
        // this returns scaled values; I = e^(x) sI, K = e^(-x) sK

        private static void ModifiedBessel_Asymptotic (double nu, double x, out double sI, out double sK) {

            double fI = 1.0;
            double fK = 1.0;
            double df = 1.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double fI_old = fI; double fK_old = fK;

                df = df * (2.0 * (nu - k) + 1.0) * (2.0 * (nu + k) - 1.0) / k / (8.0 * x);

                if (k % 2 == 0) {
                    fI += df;
                } else {
                    fI -= df;
                }

                fK += df;

                if ((fI == fI_old) && (fK == fK_old)) {
                    sI = fI / Math.Sqrt(2.0 * Math.PI * x);
                    sK = fK / Math.Sqrt(2.0 * x / Math.PI);
                    //lnI = x - Math.Log(2.0 * Math.PI * x) / 2.0 + Math.Log(fI);
                    //lnK = -x - Math.Log(2.0 * x / Math.PI) / 2.0 + Math.Log(fK);
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

            // Thompson-Barnett do normalizing sum using these auxiluary variables (see Numerical Recipies)

            double C = (0.5 - nu) * (0.5 + nu);
            double q0 = 0.0;
            double q1 = 1.0;
            double Q = C * q1;
            double S = 1.0 + Q * Df;
            //decimal S = (decimal) (1.0 + Q * Df);

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
                //S += (decimal) (Q * Df);
                S += Q * Df;

                if ((S == S_old) && (f == f_old)) {
                    g = -((x + 0.5) + (nu + 0.5) * (nu - 0.5) * f) / x;
                    K = Math.Sqrt(Math.PI / 2.0 / x) * Math.Exp(-x) / ((double) S);
                    return;
                }

                q0 = q1; q1 = q2;

            }

            throw new NonconvergenceException();
        }

        // use Wronskian I' K - I K' = 1/x, plus values of I'/I, K'/K, and K, to compute I, I', K, and K'

        /*
        private static SolutionPair ModifiedBessel_Steed (double nu, double x) {

            double f = ModifiedBessel_CF1(nu, x);

            double K, g;
            ModifiedBessel_CF_K(nu, x, out K, out g);

            double KP = g * K;

            double I = 1.0 / K / (f - g) / x;

            double IP = f * I;

            return (new SolutionPair(nu, x, I, IP, K, KP));

        }
        */

        // Temme's series for K0 and K1 for small nu and small x
        // applies only for -1/2 <= nu <= 1/2
        // accurate to all but last couple digits for x < 2; converges further out but accuracy deteriorates rapidly
        // initial constant determination shared with BesselY_Series; we should factor out the common parts

        private static void ModifiedBesselK_Series (double nu, double x, out double K0, out double K1) {

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

            double p = Math.Pow(x / 2.0, -nu) * GP / 2.0;
            double q = Math.Pow(x / 2.0, +nu) * GM / 2.0;

            // determine initial f; this is rather complicated

            double s = nu * Math.Log(2.0 / x);
            double C1, C2;
            if (Math.Abs(s) < 1.0e-3) {
                C1 = 1.0 + s * s / 2.0 + s * s * s * s / 24;
                C2 = 1.0 + s * s / 6.0 + s * s * s * s / 120.0;
            } else {
                C1 = Math.Cosh(s);
                C2 = Math.Sinh(s) / s;
            }

            double G1, G2, F;
            if (Math.Abs(nu) < 1.0e-5) {
                G1 = -EulerGamma;
                G2 = 1.0;
                F = 1.0;
            } else {
                G1 = (1.0 / GM - 1.0 / GP) / (2.0 * nu);
                G2 = (1.0 / GM + 1.0 / GP) / 2.0;
                F = (nu * Math.PI) / Math.Sin(nu * Math.PI);
            }

            double f = F * (C1 * G1 + C2 * G2 * Math.Log(2.0 / x));

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

                // Console.WriteLine(String.Format("k = {0}, Y0={1}", k, -Y0));

                if ((K0 == K0_old) && (K1 == K1_old)) {
                    K1 = 2.0 / x * K1;
                    return;
                }
            }
            throw new NonconvergenceException();

        }

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

        /// <summary>
        /// Computes the Airy function of the first kind.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value Ai(x).</returns>
        /// <remarks>
        /// <para>The Airy functions appear in quantum mechanics in the semiclassical WKB solution to the wave functions in a potential.</para>
        /// <para>For negative arguments, Ai(x) is oscilatory. For positive arguments, it decreases exponentially with increasing x.</para>
        /// </remarks>
        /// <seealso cref="AiryBi"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Airy_functions" />
        public static double AiryAi (double x) {

            if (Math.Abs(x) < 2.0) {
                // close to the origin, use a power series
                // the definitions in terms of bessel functions becomes 0 X infinity at x=0, so we can't use them
                return (AiryAi_Series(x));
            } else if (x > 0.0) {
                double y = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
                return (Math.Sqrt(x / 3.0) / Math.PI * ModifiedBesselK(1.0 / 3.0, y));
            } else {
                // change in future to call a function which returns J and Y together
                double y = 2.0 / 3.0 * Math.Pow(-x, 3.0 / 2.0);
                double J = BesselJ(1.0 / 3.0, y);
                double Y = BesselY(1.0 / 3.0, y);
                return (Math.Sqrt(-x) / 2.0 * (J - Y / Math.Sqrt(3.0)));
            }
        }

        // problematic for positive x once exponential decay sets in; don't use for x > 2

        private static double AiryAi_Series (double x) {

            double p = 1.0 / (Math.Pow(3.0, 2.0 / 3.0) * AdvancedMath.Gamma(2.0 / 3.0));
            double q = x / (Math.Pow(3.0, 1.0 / 3.0) * AdvancedMath.Gamma(1.0 / 3.0));
            double f = p - q;

            double x3 = x * x * x;
            for (int k = 0; k < Global.SeriesMax; k+=3) {
                double f_old = f;
                p *= x3 / ((k + 2)*(k + 3));
                q *= x3 / ((k + 3)*(k + 4));
                f += p - q;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();

        }

        /// <summary>
        /// Computes the Airy function of the second kind.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value Bi(x).</returns>
        /// <remarks>
        /// <para>The Airy functions appear in quantum mechanics in the semiclassical WKB solution to the wave functions in a potential.</para>
        /// <para>For negative arguments, Bi(x) is oscilatory. For positive arguments, it increases exponentially with increasing x.</para>
        /// <para>While the notation Bi(x) was chosen simply as a natural complement to Ai(x), it has influenced the common nomenclature for
        /// this function, which is now often called the "Bairy function".</para>
        /// </remarks>        
        /// <seealso cref="AiryAi"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Airy_functions" />
        public static double AiryBi (double x) {
            if (Math.Abs(x) < 2.0) {
                return (AiryBi_Series(x));
            } else if (x > 0.0) {
                // change to use a function that returns I and K together
                double y = 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0);
                double I = ModifiedBesselI(1.0 / 3.0, y);
                double K = ModifiedBesselK(1.0 / 3.0, y);
                return (Math.Sqrt(x) * (2.0 / Math.Sqrt(3.0) * I + K / Math.PI));
            } else {
                // change to use a function that returns J and Y together
                double y = 2.0 / 3.0 * Math.Pow(-x, 3.0 / 2.0);
                double J = BesselJ(1.0 / 3.0, y);
                double Y = BesselY(1.0 / 3.0, y);
                return (-Math.Sqrt(-x) / 2.0 * (J / Math.Sqrt(3.0) + Y));
            }
        }

        private static double AiryBi_Series (double x) {

            double p = 1.0 / (Math.Pow(3.0, 1.0 / 6.0) * AdvancedMath.Gamma(2.0 / 3.0));
            double q = x * (Math.Pow(3.0, 1.0 / 6.0) / AdvancedMath.Gamma(1.0 / 3.0));
            double f = p + q;

            double x3 = x * x * x;
            for (int k = 0; k < Global.SeriesMax; k += 3) {
                double f_old = f;
                p *= x3 / ((k + 2) * (k + 3));
                q *= x3 / ((k + 3) * (k + 4));
                f += p + q;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();

        }

    }

}
