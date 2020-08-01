using System;
using System.Diagnostics;

using Meta.Numerics.Analysis;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        // **** Integer order Bessel functions ****

        /// <summary>
        /// Computes the regular Bessel function for integer orders.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of J<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>For information on the cylindrical Bessel functions, see <see cref="AdvancedMath.Bessel"/>.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
        public static double BesselJ (int n, double x) {

            // Relate negative n to positive n.
            if (n < 0) {
                if ((n % 2) == 0) {
                    return (BesselJ(-n, x));
                } else {
                    return (-BesselJ(-n, x));
                }
            }

            // Relate negative x to positive x.
            if (x < 0) {
                if ((n % 2) == 0) {
                    return (BesselJ(n, -x));
                } else {
                    return (-BesselJ(n, -x));
                }
            }

            if (x <= Math.Sqrt(2 * (n + 1))) {
                // If x is small enough, use the series.
                // The transition point is chosen so that 2nd term cannot overwhelm 1st.
                return (BesselJ_Series(n, x));
            } else if (x >= 32.0 + 0.5 * n * n) {
                // If x is large enough, use the asymptotic expansion.
                return (Bessel_Asymptotic(n, x).FirstSolutionValue);
            } else if ((x > n) && (n > 32)) {
                // We are in the transition region \sqrt{n} < x < n^2. If x > n we are still allowed to recurr upward,
                // and as long as we can find an n below us for which x lies in the asymptotic region, we have a value
                // to recur upward from. In general we expect this to be better than Miller's algorithm, because
                // it will require fewer recurrance steps.

                int m = (int) Math.Floor(Math.Sqrt(2.0 * (x - 32.0)));
                Debug.Assert(m <= n);

                SolutionPair s = Bessel_Asymptotic(m, x);
                double J = s.FirstSolutionValue;
                double JP = s.FirstSolutionDerivative;
                Bessel_RecurrUpward(m, x, ref J, ref JP, n - m);
                return (J);
            } else {
                // We are in the transition region; x is too large for the series but still less than n. In this region
                // upward recurrance is unstable, since even a tiny mixture of the Y solution will come to dominate as n increases.

                // Instead we will use Miller's algorithm: recurr downward from far above and use a sum rule to normalize the result
                // The sum rule we will use is:
                //   \sum_{k=-\infty}^{\infty} J_{2k}(x) = J_0(x) + 2 J_2(x) + 2 J_4(x) + \cdots = 1
                return (Bessel_Miller(n, -1, x).FirstSolutionValue);
                //return (BesselJ_Miller(n, x));
            }
        }

        /// <summary>
        /// Computes the irregular Bessel function for integer orders.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of Y<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>For information on the cylindrical Bessel functions, see <see cref="AdvancedMath.Bessel"/>.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
        public static double BesselY (int n, double x) {

            // Relate negative n to positive n.
            if (n < 0) {
                if ((n % 2) == 0) {
                    return (BesselY(-n, x));
                } else {
                    return (-BesselY(-n, x));
                }
            }

            // Relate negative x to positive x.
            if (x < 0) {
                if ((n % 2) == 0) {
                    return (BesselY(n, -x));
                } else {
                    return (-BesselY(n, -x));
                }
            }


            if (x == 0.0) {
                // special case at origin
                return (Double.NegativeInfinity);
            } else if (x < 4.0) {
                // for small enough x, use the power series for Y0 and Y1, then recurr upward

                // get Y0 and Y1
                double Y0, Y1;
                BesselY_Series(x, out Y0, out Y1);

                if (n == 0) {
                    return (Y0);
                } else if (n == 1) {
                    return (Y1);
                } else {
                    double Ym1 = Y0;
                    double Y = Y1;
                    for (int k = 1; k < n; k++) {
                        double Yp1 = (2 * k) / x * Y - Ym1;
                        Ym1 = Y;
                        Y = Yp1;
                        if (Double.IsInfinity(Y)) break;
                    }
                    return (Y);
                }

            } else if (x >= (32.0 + 0.5 * n * n)) {
                // far enough out, use the asymptotic series
                return (Bessel_Asymptotic(n, x).SecondSolutionValue);
            } else {
                // use Miller's algorithm to compute a tower of J's
                // use the tower of J's to compute Y0
                // use the Wronskian to compute Y1
                // recurr Y upward to desired order
                return (Bessel_Miller(-1, n, x).SecondSolutionValue);
                //return (BesselY_Miller(n, x));
            }

        }

        // This implementation of the Miller algorithm starts from a high enough order to accurately
        // determine J_{nJ}, iterates down to J_0, and uses a sum rule to normalize correctly.
        // It then uses another sum rule to determine Y_0 and iterates up to Y_{nY}.
        // Coding this way allows one routine to cover the case of compute just one J, just one Y,
        // or both Y and J and deratives, at the cost of just a few extra flops and branch tests
        // that would not be necessary in invidually coded routines.

        private static SolutionPair Bessel_Miller (int nJ, int nY, double x) {

            // We can compute J alone or Y alone, but if we are computing both, it must be for the same order.
            Debug.Assert((nJ < 0) || (nY < 0) || (nJ == nY));

            // Start at a high enough order that values are neglible relative to the order sought.
            int m = Bessel_Miller_Limit(Math.Max(1, nJ), x);
            Debug.Assert(m > nJ);

            // Recurr J down to 0, keeping track of some sums we will use for re-scaling.
            // Starting values are arbitrary, but should be small since they will grow,
            // and we choose consistent with J_{m+1} = 0 to reflect assumption of rapid decrease.
            double J = Bessel_Tower_Start;
            double JP = m / x * Bessel_Tower_Start;
            double s = 0.0;
            double Y = 0.0;
            double Jn = Double.NaN;
            double JPn = Double.NaN;
            while (m > 0) {
                if (m % 2 == 0) {
                    s += J;
                    int mh = m / 2;
                    double dy = J / mh;
                    if (mh % 2 == 0) dy = -dy;
                    Y += dy;
                }
                double t = J;
                J = m / x * J + JP;
                m--;
                JP = m / x * J - t;
                // If there is a J we are trying to compute, remember it on our way to zero.
                if (m == nJ) {
                    Jn = J;
                    JPn = JP;
                }
            }

            // Use 1 = J_0 + 2 J_2 + 2 J_4 + 2 J_6 + 2 J_8 + \cdots to re-scale Js.
            s = J + 2.0 * s;

            // Re-scale the Js we want.
            Jn /= s;
            JPn /= s;

            // If we don't want any Ys, we are done.
            if (nY < 0) return new SolutionPair(Jn, JPn, Double.NaN, Double.NaN);

            // Find J_0 and J_0'.
            J /= s;
            JP /= s;

            // Use Y_0 = (2 / \pi) \left[ ( \ln(x/2) + \gamma ) J_0 + J_2 - 1/2 J_4 + 1/3 J_6 - 1/4 J_8 + \cdots \right]
            // to find Y_0.
            Y /= s;
            Y = 2.0 / Math.PI * ((Math.Log(0.5 * x) + AdvancedMath.EulerGamma) * J + 2.0 * Y);

            // Use Wronskian J_0 Y_0 - Y_0 J_0' = \frac{2}{\pi x} to find Y_0'.
            double YP = (2.0 / Math.PI / x + Y * JP) / J;

            // Recurr Y up to the desired n.
            while (m < nY) {
                double t = Y;
                Y = m / x * Y - YP;
                m++;
                YP = t - m / x * Y;
                if (Double.IsNegativeInfinity(Y)) break;
            }

            return new SolutionPair(JP, JPn, Y, YP);

        }

        private static double Bessel_Tower_Start = Math.Pow(2.0, -512);

        // this is just an integer version of the series we implement below for doubles;
        // having an integer-specific version is slightly faster

        private static double BesselJ_Series (int n, double x) {
            Debug.Assert(n >= 0);
            Debug.Assert(x >= 0.0);
            double z = 0.5 * x;
            double dJ = AdvancedMath.PowerOverFactorial(z, n);
            double J = dJ;
            double zz = -z * z;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double J_old = J;
                dJ *= zz / ((n + k) * k);
                J += dJ;
                if (J == J_old) {
                    return (J);
                }
            }
            throw new NonconvergenceException();
        }

        // Miller's algorithm requires that we start with arbitrary values at an L far above the one sought, then recurr downward
        // to the sought L and use a relation to normalize the result. Since J decreases with increasing L in the region of interest,
        // the initial values don't matter as long as L is far enough above that their contribution to the normaliztion is supressed
        // below the target accuracy. How far above is sufficient?

        // There doesn't appear to be a universally accepted answer to this question. NR uses a formula of the norm m = n + c \sqrt{n},
        // which they claim can be justified via some hueristic manipulation of limiting expressions that I cannot reproduce.
        // Other literature, and my own experience, implies life isn't so simple. Olver gives the most convincing answer:
        // run the recurrence up until it has produced a factor larger than the error supression you require, then start back
        // down from there. This is the approach we adopt.

        private static int Bessel_Miller_Limit (int n, double x) {
            double Fm = 0.0;
            double F = 1.0;
            while (Math.Abs(F) < 1.0E20) {
                double Fp = (2 * n) / x * F - Fm;
                Fm = F;
                F = Fp;
                n++;
            }
            return (n);
        }

        // a series expansion in x, fastest for small x, but works pretty well quite far out
        // 10 terms at x~1.0, 30 terms at x~10.0, but accuracy in the high digits suffers that far out
        // this organization of the series is due to Temme

        private static void BesselY_Series (double x, out double Y0, out double Y1) {

            if (x == 0.0) {
                Y0 = Double.NegativeInfinity;
                Y1 = Double.NegativeInfinity;
                return;
            }

            double xx = -x * x / 4.0;
            double cx = 1.0;

            double f = (2.0 / Math.PI) * (Math.Log(2.0 / x) - EulerGamma);
            double p = 1.0 / Math.PI;
            double h = p;

            Y0 = cx * f;
            Y1 = cx * h;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double Y0_old = Y0;
                double Y1_old = Y1;

                cx = cx * xx / k;

                f = f / k + 2 * p / (k * k);
                p = p / k;
                h = p - k * f;

                Y0 += cx * f;
                Y1 += cx * h;

                if ((Y0 == Y0_old) && (Y1 == Y1_old)) {
                    Y0 = -Y0;
                    Y1 = -2.0 * Y1 / x;
                    return;
                }
            }
            throw new NonconvergenceException();

        }

        // this is the Temme series for non-integer -0.5 < nu < 0.5

        private static void BesselY_Series (double nu, double x, out double Y0, out double Y1) {

            Debug.Assert(Math.Abs(nu) <= 0.5);
            Debug.Assert(x > 0.0);

            // polynomial variable and coefficient

            double xx = -x * x / 4.0;
            double cx = 1.0;

            // determine initial p, q

            double GP = Gamma(1.0 + nu);
            double GM = Gamma(1.0 - nu);

            double p = Math.Pow(x / 2.0, -nu) * GP / Math.PI;
            double q = Math.Pow(x / 2.0, +nu) * GM / Math.PI;

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

            double G1, G2, F, cq;
            if (Math.Abs(nu) < 1.0e-5) {
                G1 = -EulerGamma;
                G2 = 1.0;
                F = 2.0 / Math.PI;
                cq = 0.0;
            } else {
                G1 = (1.0 / GM - 1.0 / GP) / (2.0 * nu);
                G2 = (1.0 / GM + 1.0 / GP) / 2.0;
                F = 2.0 * nu / Math.Sin(nu * Math.PI);
                cq = (2.0 / nu) * Math.Pow(Math.Sin(nu * Math.PI / 2.0), 2.0);
            }

            double f = F * (C1 * G1 + C2 * G2 * Math.Log(2.0 / x));

            // run the series

            Y0 = cx * (f + cq * q);
            Y1 = cx * p;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double Y0_old = Y0;
                double Y1_old = Y1;

                cx = cx * xx / k;

                double kpnu = k + nu;
                double kmnu = k - nu;
                f = (k * f + p + q) / (kpnu * kmnu);
                p = p / kmnu;
                q = q / kpnu;

                double g = f + cq * q;
                double h = p - k * g;

                Y0 += cx * g;
                Y1 += cx * h;

                if ((Y0 == Y0_old) && (Y1 == Y1_old)) {
                    Y0 = -Y0;
                    Y1 = -2.0 * Y1 / x;
                    return;
                }
            }

            throw new NonconvergenceException();

        }


        // **** Real order Bessel functions ****

        // Series development of Bessel J
        //   J_{\nu}(x) = \sum_{k=0}^{\infty} \frac{( )^{\nu + 2k}}{\Gamma(\nu + 2k + 1)}
        // As can be seen by comparing leading term and first correction, this is good for z <~ \max(1 , \sqrt{\nu})

        // For nu=0, it requires 10 terms at x~1.0, 25 terms at x~10.0, but accuracy in the last few digits suffers that far out
        // Gets better for higher nu (for nu=10, only 20 terms are required at x~10.0 and all digits are good), but only with sqrt(nu)

        private static double BesselJ_Series (double nu, double x) {
            double z = 0.5 * x;
            double dJ = AdvancedMath.PowerOverFactorial(z, nu);
            double J = dJ;
            double zz = -z * z;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double J_old = J;
                dJ *= zz / (k * (nu + k));
                J += dJ;
                if (J == J_old) {
                    return (J);
                }
            }
            throw new NonconvergenceException();
        }

        // Regular Bessel function series about origin
        //   J_{\nu}(x) = (x/2)^{\nu} \sum{k=0}^{\infty} \frac{(-x^2/4)^k}{k! \Gamma(\nu + k + 1)}
        // This can be turned into a series for J' by term-by-term differentiation
        // For nu=0, x~1 it requires ~10 terms, x~4 requires ~17 terms, x~8 requires ~24 terms
        // Since nu appears in the denominator, higher nu converges faster
        // Since first correction term is ~\frac{x^2}{4\nu}, the radius of fast convergence should grow ~2\sqrt{\nu}

        private static void BesselJ_Series (double nu, double x, out double J, out double JP) {

            // We divide by x to deal with derivative, so handle x = 0 separately.
            Debug.Assert(x > 0.0);

            // Evaluate the series of J and J', knowing x != 0
            // Note that the J' series is just with J series with teach term having one less power,
            // and each term multipled by the power it has in the J series.
            double x2 = 0.5 * x;
            double x22 = -x2 * x2;
            double dJ = AdvancedMath.PowerOverFactorial(x2, nu);
            J = dJ;
            JP = nu * dJ;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double J_old = J;
                double JP_old = JP;
                dJ *= x22 / ( k * (nu + k));
                J += dJ;
                JP += (nu + 2 * k) * dJ;
                if ((J == J_old) && (JP == JP_old)) {
                    JP = JP / x;
                    return;
                }
            }
            // The J series alone requires 5 ops per term; J' evaulation adds 2 ops per term.
            throw new NonconvergenceException();

        }

        // computes the ratio f=J_{nu+1}/J_{nu} via a continued fraction
        // converges rapidly for x < Sqrt(nu(nu+1)) ~ nu, and slowly for larger x
        // problematic for very small x, because we divide by x, so very close to the origin use series instead
        // note that J'/J = nu/x - f so this can also be used to compute the ratio of J to its derivative
        private static double Bessel_CF1 (double nu, double x, out int sign) {

            // seed with the first fraction
            sign = 1;
            double Am2 = 0.0;
            double Am1 = 1.0;
            double Bm2 = 1.0;
            double Bm1 = 2.0 * (nu + 1.0) / x;
            double f = Am1 / Bm1;

            // for x > nu, we don't expect the series to even start converging for about (x - nu) steps
            int kmax = Global.SeriesMax;
            if (x > nu) kmax += (int) Math.Ceiling(x - nu);

            // move to higher fractions
            for (int k = 2; k < 2 * kmax; k++) {
                double a = -1.0;
                double b = 2.0 * (nu + k) / x;
                double A = b * Am1 + a * Am2;
                double B = b * Bm1 + a * Bm2;

                //Console.WriteLine("k={0}, A={1}, B={2}, f={3}, sign={4}", k, A, B, f, sign);

                // renormalize and check for convergence
                if (B != 0.0) {
                    if (B < 0.0) sign = -sign;

                    Am1 = Am1 / B;
                    A = A / B;
                    Bm1 = Bm1 / B;
                    B = 1.0;

                    if (A == f) {
                        return (f);
                    } else {
                        f = A;
                    }
                }

                // prepare for next cycle
                Am2 = Am1;
                Am1 = A;
                Bm2 = Bm1;
                Bm1 = B;
            }

            /*
            double a = 1.0;
            double b = 2.0*(nu+1.0)/x;
            double D = 1.0 / b;
            double Df = a / b;
            double f = Df;
            sign = 1;
            for (int k = 2; k < seriesMax; k++) {
                double f_old = f;
                a = -1.0;
                b = 2.0*(nu+k)/x;
                D = b + a * D;
                if (D == 0.0) D = 1.0e-30;
                D = 1.0 / D;
                Df = (b * D - 1.0) * Df;
                f += Df;
                if (D < 0.0) sign = -sign;
                Console.WriteLine("k={0}, D={1}, Df = {2}", k, D, Df);
                if (f == f_old) return (k);
            }
            */
            throw new NonconvergenceException();
        }

        // computes the ratio (J' + i Y')/(J + i Y) via a continued fraction
        // converges rapidly for x > Sqrt(nu(nu+1)), but does not converge for lower x
        private static Complex Bessel_CF2 (double nu, double x) {
            double nu2 = nu * nu;
            double a = 0.25 - nu2;
            Complex b = 2.0 * new Complex(x, 1.0);
            Complex D = 1.0 / b;
            Complex Df = a / b;
            Complex f = Df;
            for (int k = 3; k < Global.SeriesMax; k += 2) {
                Complex f_old = f;
                a = k * k / 4.0 - nu2;
                //b.Im = b.Im + 2.0;
                b = b + new Complex(0.0, 2.0);
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;
                if (f == f_old) return (-0.5 / x + ComplexMath.I + ComplexMath.I * f / x);
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes the regular Bessel function for real orders.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of J<sub>&#x3BD;</sub>(x).</returns>
        /// <remarks>
        /// <para>For information on the cylindrical Bessel functions, see <see cref="AdvancedMath.Bessel"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
        public static double BesselJ (double nu, double x) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

            // Use reflection to turn negative orders into positive orders.
            if (nu < 0.0) {
                double mu = -nu;
                return (MoreMath.CosPi(mu) * BesselJ(mu, x) - MoreMath.SinPi(mu) * BesselY(mu, x));
            }

            if (x <= Math.Sqrt(2.0 * (nu + 1.0))) {
                // We are close enough to origin to use the series.
                return (BesselJ_Series(nu, x));
            } else if (x >= 32.0 + 0.5 * nu * nu) {
                // We are far enough from origin to use the asymptotic expansion.
                SolutionPair result = Bessel_Asymptotic(nu, x);
                return (result.FirstSolutionValue);
            } else if (x >= nu) {
                // We are far enough from origin to evaluate CF2, so use Steed's method.
                SolutionPair result = Bessel_Steed(nu, x);
                return (result.FirstSolutionValue);
            } else {
                // We have x < nu, but x is still not small enough to use the series; this only occurs for nu >~ 6.
                // To handle this case, compute J_{nu+1}/J_{nu}, recurse down to mu where mu ~ x < nu; use
                // Steed's method to evaluate J_{mu} and re-normalize J_{nu}.

                // for example, for mu = 16, x = 12, we can't evaluate CF2 because x < 16; so assume J_16=1, compute J_17 / J_16,
                // and recurse down to J_11; since 12 > 11, we can compute CF2 and get J_11 and Y_11; comparing this
                // with our value of J_11 gives the proper renormalization factor for J_16

                int sign;
                double r = nu / x - Bessel_CF1(nu, x, out sign);

                double J = sign * Bessel_Tower_Start;
                //double J = 1.0 * sign;
                double JP = r * J;
                double mu = nu;
                while (mu > x) {
                    double t = J;
                    J = mu / x * J + JP;
                    mu -= 1.0;
                    JP = mu / x * J - t;
                    // Not sure if this is really okay in all cases.
                    if (Double.IsInfinity(J)) return (0.0);
                }

                Complex z = Bessel_CF2(mu, x);
                SolutionPair result = Bessel_Steed(JP / J, z, 2.0 / Math.PI / x, Math.Sign(J));
                //SolutionPair result = Bessel_Steed(mu / x - Jp1 / J, z, 2.0 / Math.PI / x, Math.Sign(J));

                return (result.FirstSolutionValue / J * Bessel_Tower_Start);

            }
        }


        /// <summary>
        /// Computes the irregual Bessel function for real orders.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of Y<sub>&#x3BD;</sub>(x).</returns>
        /// <remarks>
        /// <para>For information on the cylindrical Bessel functions, see <see cref="AdvancedMath.Bessel"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="Bessel"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
        public static double BesselY (double nu, double x) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

            if (nu < 0.0) {
                return (MoreMath.CosPi(nu) * BesselY(-nu, x) - MoreMath.SinPi(nu) * BesselJ(-nu, x));
            }

            if (x == 0.0) {
                return (Double.NegativeInfinity);
            } else if (x >= 32.0 + nu * nu / 2.0) {

                // far from the origin, use the asymptotic expansion

                SolutionPair result = Bessel_Asymptotic(nu, x);
                return (result.SecondSolutionValue);
            } else if (x < 4.0) {

                // close to the origin, use the Temme series

                double Y0, Y1;
                if (nu <= 0.5) {
                    BesselY_Series(nu, x, out Y0, out Y1);
                    return (Y0);
                } else if (nu <= 1.5) {
                    BesselY_Series(nu - 1.0, x, out Y0, out Y1);
                    return (Y1);
                } else {

                    // compute for mu -0.5 <= mu <= 0.5
                    double mu = nu - Math.Round(nu);
                    BesselY_Series(mu, x, out Y0, out Y1);

                    // recurr upward to nu
                    int n = ((int) Math.Round(nu - mu)) - 1;
                    for (int i = 0; i < n; i++) {
                        // use recurrence on Y, Y':
                        //   Y_{mu+1} = (2 mu/x) Y_{mu} - Y_{mu}
                        mu += 1.0;
                        double t = Y0;
                        Y0 = Y1;
                        Y1 = (2.0 * mu / x) * Y1 - t;
                        if (Double.IsNegativeInfinity(Y1)) break;
                    }
                    return (Y1);

                }

            } else if (x > nu) {
                SolutionPair result = Bessel_Steed(nu, x);
                return (result.SecondSolutionValue);
            } else {

                // we have 4 < x < nu; evaluate at mu~x and recurr upward to nu

                // figure out mu
                int n = (int) Math.Ceiling(nu - x);
                double mu = nu - n;

                // evaluate at mu
                SolutionPair result = Bessel_Steed(mu, x);
                double Y = result.SecondSolutionValue;
                double YP = result.SecondSolutionDerivative;


                // recurr upward to nu
                Bessel_RecurrUpward(mu, x, ref Y, ref YP, n);

                return (Y);

            }

        }

        /// <summary>
        /// Computes both solutions of the Bessel differential equation.
        /// </summary>
        /// <param name="nu">The order, which must be non-negative.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The values of J<sub>nu</sub>(x), J'<sub>nu</sub>(x), Y<sub>nu</sub>(x), and Y'<sub>nu</sub>(x).</returns>
        /// <remarks>
        /// <para>Bessel functions often occur in the analysis of physical phenomena with cylindrical symmetry.
        /// They satisfy the differential equation:</para>
        /// <img src="../images/BesselODE.png" />
        /// <para>Since this is a second order linear equation, it has two linearly independent solutions. The regular Bessel function
        /// J<sub>&#x3BD;</sub>(x), which is finite at the origin, and the irregular Bessel function Y<sub>&#x3BD;</sub>(x), which diverges
        /// at the origin. Far from the origin, both functions are oscilatory.</para>
        /// <para>This method simultaneously computes both Bessel functions and their derivatives. If you need both J and Y, it is faster to call this method once than to call
        /// <see cref="BesselJ(double,double)"/> and <see cref="BesselY(double,double)"/> seperately. If on, the other hand, you need only J or only Y, it is faster to
        /// call the appropriate method to compute the one you need.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="nu"/> or <paramref name="x"/> is negative.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
        public static SolutionPair Bessel (double nu, double x) {

            if (nu < 0.0) throw new ArgumentOutOfRangeException(nameof(nu));
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

            if (x == 0.0) {
                return Bessel_Zero(nu, +1);
            } else if (x < 4.0) {

                // Close to the origin, use the series for Y and J.
                // The range does not grow with increasing order because our series for Y is only for small orders.

                // for J, we can use the series directly
                double J, JP;
                BesselJ_Series(nu, x, out J, out JP);

                // for Y, we only have a series for -0.5 < mu < 0.5, so pick an order in that range offset from our
                // desired order by an integer number of recurrence steps. Evaluate at the mu order and recurr upward to nu.
                double mu = nu - Math.Round(nu);
                int n = (int) Math.Round(nu - mu);
                double Y, YP1;
                BesselY_Series(mu, x, out Y, out YP1);
                double YP = (mu / x) * Y - YP1;
                Bessel_RecurrUpward(mu, x, ref Y, ref YP, n);

                // return J, J', Y, Y' in a solution pair
                return (new SolutionPair(J, JP, Y, YP));

            } else if (x > 32.0 + nu * nu / 2.0) {
                // Far enough out, use the asymptotic series.
                return (Bessel_Asymptotic(nu, x));
            } else if (x > nu) {
                // In the intermediate region, use Steed's method directly if x > nu.
                return (Bessel_Steed(nu, x));
            } else {

                // If x < nu, we can't use Steed's method directly because CF2 does not converge.
                // So pick a mu ~ x < nu, offset from nu by an integer; use Steed's method there,
                // then recurr up to nu.

                // this only occurs for 4 < x < nu, e.g. nu, x = 16, 8

                int n = (int) Math.Ceiling(nu - x);
                double mu = nu - n;
                SolutionPair result = Bessel_Steed(mu, x);
                double Y = result.SecondSolutionValue;
                double YP = result.SecondSolutionDerivative;
                Bessel_RecurrUpward(mu, x, ref Y, ref YP, n);

                // being in a region with x < nu, we can evaluate CF1 to get J'/J
                // this ratio, together with Y and Y' and the Wronskian, gives us J' and J seperately
                int sign;
                double r = nu / x - Bessel_CF1(nu, x, out sign);
                double J = (2.0 / Math.PI / x) / (YP - r * Y);
                double JP = r * J;

                // return the result
                return (new SolutionPair(J, JP, Y, YP));

            }

        }

        // We have the recurrance
        //   F_{\nu+1}(x) = \frac{\nu}{x} F_{\nu}(x) - F_{\nu}'(x)
        // and the relationship
        //   F_{\nu}'(x) = F_{\nu-1}(x) - \frac{\nu}{x} F_{\nu}(x)
        // We can use this to recurr upward arbitrarily.

        // While mathematically valid for F being either J or Y, it is unstable for J in the upward direction in the region x < nu
        // This is because Y >> J in that region, so a very small admixture of Y into F, introduced by roundoff error, will
        // quickly dwarf the desired J solution.

        // This method recurrs upward n times, starting at the given nu value.

        private static void Bessel_RecurrUpward (double nu, double x, ref double F, ref double FP, int n) {

            for (int i = 0; i < n; i++) {
                double t = F;
                F = (nu / x) * F - FP;
                nu += 1.0;
                FP = t - (nu / x) * F;
                if (Double.IsInfinity(F)) break;
            }

        }

        // CF1, CF2, and Wronskian to determine Bessel function values
        private static SolutionPair Bessel_Steed (double nu, double x) {

            // CF2 only converges for x > nu
            Debug.Assert(x + 1.0 > nu);

            // Wronskian
            double W = 2.0 / Math.PI / x;

            // CF1
            int sign;
            double r = nu / x - Bessel_CF1(nu, x, out sign);

            // CF2
            Complex z = Bessel_CF2(nu, x);

            // use these results to compute J, JP, Y, and YP
            // what about zeros?

            double g = (z.Re - r) / z.Im;

            double J = Math.Sqrt(W / (z.Im + g * (z.Re - r)));
            if (sign < 0) J = -J;

            double JP = r * J;

            double Y = g * J;

            double YP = Y * (z.Re + z.Im / g);

            return new SolutionPair(J, JP, Y, YP);

        }

        // Evaluate J, J', Y, Y' via Steed's method
        // r = J'/J, z = ()/(), W = J Y' - Y J', and sign is the sign of J
        private static SolutionPair Bessel_Steed (double r, Complex z, double W, int sign) {

            double g = (z.Re - r) / z.Im;

            double J = Math.Sqrt(W / (z.Im + g * (z.Re - r)));
            if (sign < 0) J = -J;

            double JP = r * J;

            double Y = g * J;

            double YP = Y * (z.Re + z.Im / g);

            //return new SolutionPair(0, 0.0, J, JP, Y, YP);
            return new SolutionPair(J, JP, Y, YP);

        }

        // Hankel's asymptotic expansions for Bessel function (A&S 9.2)
        //   J = \sqrt{\frac{2}{\pi x}} \left[ P \cos\omega - Q \sin\omega \right]
        //   Y = \sqrt{\frac{2}{\pi x}} \left[ P \sin\omega + Q \cos\omega \right]
        // where \omega = x - \left( \nu / 2 + 1 / 4 \right) \pi and \mu = 4 \nu^2 and
        //   P = 1 - \frac{(\mu-1)(\mu-9)}{2! (8x)^2} + \frac{(\mu-1)(\mu-9)(\mu-25)(\mu-49}{4! (8x)^4} + \cdots
        //   Q = \frac{(\mu-1)}{8x} - \frac{(\mu-1)(\mu-9)(\mu-25)}{3! (8x)^3} + \cdots
        // Derivatives have similiar expressions
        //   J' = - \sqrt{\frac{2}{\pi x}} \left[ R \sin\omega + S \cos\omega \right]
        //   Y' = \sqrt{\frac{2}{\pi x}} \left[ R \cos\omega - S \sin\omega \right]
        // where
        //   R = 1 - \frac{(\mu-1)(\mu+15)}{2! (8x)^2} + \cdots
        //   S = \frac{(\mu+3)}{8x} - \frac{(\mu-1)(\mu - 9)(\mu+35)}{3! (8x)^3} + \cdots

        // For nu=0, this series converges to full precision in about 10 terms at x~100, and in about 25 terms even as low as x~25.
        // It fails to converge to full precision for x <~ 25.
        // Since the first correction term is ~ (4 \nu^2)^2 / (8 x)^2 ~ (\nu^2 / 2 x)^2, the minimum x should grow like \nu^2 / 2

        private static SolutionPair Bessel_Asymptotic (double nu, double x) {

            Debug.Assert(nu >= 0.0);
            Debug.Assert(x > 0.0);

            // pre-compute factors of nu and x as they appear in the series
            double mu = 4.0 * nu * nu;
            double xx = 8.0 * x;

            // initialize P and Q
            double P = 1.0; double R = 1.0;
            double Q = 0.0; double S = 0.0;

            // k is the current term number, k2 is (2k - 1), and t is the value of the current term
            int k = 0;
            int k2 = -1;
            double t = 1.0;

            while (true) {

                double Q_old = Q; double P_old = P;
                double R_old = R; double S_old = S;

                k++; k2 += 2;
                t /= k * xx;
                S += (mu + k2 * (k2 + 2)) * t;
                t *= (mu - k2 * k2);
                Q += t;

                k++; k2 += 2;
                t /= -k * xx;
                R += (mu + k2 * (k2 + 2)) * t;
                t *= (mu - k2 * k2);
                P += t;

                if ((P == P_old) && (Q == Q_old) && (R == R_old) && (S == S_old)) {
                    break;
                }

                if (k > Global.SeriesMax) throw new NonconvergenceException();

            }

            // The computation of \sin\omega and \cos\omega is fairly delicate. Since x is large,
            // it's important to use MoreMath.Sin and MoreMath.Cos to ensure correctness for
            // large arguments. Furthermore, since the shift by 1/2 ( \nu + 1/2) \pi is a
            // multiple of \pi and could be dwarfed by or have significant cancellation with x,
            // it's better to use the trig addition formulas to handle it seperately using
            // SinPi and CosPi functions.

            // This works well, but it took a long time to get to this.

            double a = 0.5 * (nu + 0.5);
            double sa = MoreMath.SinPi(a);
            double ca = MoreMath.CosPi(a);

            double sx = MoreMath.Sin(x);
            double cx = MoreMath.Cos(x);

            // It would be great to have a SinAndCos function so as not to do the range reduction twice.

            double s1 = sx * ca - cx * sa;
            double c1 = cx * ca + sx * sa;
            
            // Assemble the solution
            double N = Math.Sqrt(2.0 / Math.PI / x);
            SolutionPair result = new SolutionPair(
                N * (c1 * P - s1 * Q), -N * (R * s1 + S * c1),
                N * (s1 * P + c1 * Q), N * (R * c1 - S * s1)
            );

            return (result);

        }

        // Behavior at zero argument is slightly complicated, depending on \nu.
        // The limits for J and I are the same. The limits for Y and K differ
        // only in sign. Therefore we use the same method for both.

        private static SolutionPair Bessel_Zero (double nu, int s) {

            Debug.Assert(nu >= 0.0);

            double J, JP;
            if (nu == 0.0) {
                J = 1.0;
                JP = 0.0;
            } else {
                J = 0.0;
                if (nu < 1.0) {
                    JP = Double.PositiveInfinity;
                } else if (nu == 1.0) {
                    JP = 1.0 / 2.0;
                } else {
                    JP = 0.0;
                }
            }

            double Y = s * Double.NegativeInfinity;
            double YP = s * Double.PositiveInfinity;

            return (new SolutionPair(J, JP, Y, YP));

        }

        // **** Spherical Bessel functions ****

        

        // Functions needed for uniform asymptotic expansions

        /// <summary>
        /// Computes the requested zero of the regular Bessel function.
        /// </summary>
        /// <param name="nu">The order, which must be non-negative.</param>
        /// <param name="k">The index of the zero, which must be positive.</param>
        /// <returns>The <paramref name="k"/>th value of x for which J<sub>nu</sub>(x) = 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="nu"/> is negative or <paramref name="k"/> is non-positive.</exception>
        public static double BesselJZero (double nu, int k) {
            if (nu < 0.0) throw new ArgumentOutOfRangeException(nameof(nu));
            if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

            double jMin, jMax;
            if (k == 1) {
                // Rayleigh inequalities
                // Cited in Ismail and Muldoon, Bounds for the small real and purely imaginary zeros of
                // Bessel and related functions, Methods and Applications of Analysis 2 (1995) p. 1-21
                // http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.140.288&rep=rep1&type=pdf
                // They cite Watson p. 502, but I don't find them there.
                //jMin = 2.0 * Math.Sqrt((nu + 1.0) * Math.Sqrt(nu + 2.0));
                // Lee Lorch, "Some Inequalities for the First Positive Zeros of the Bessel Functions",
                // SIAM J. Math. Anal. 24 (1993) 814
                jMin = Math.Sqrt((nu + 1.0) * (nu + 5.0));
                jMax = ApproximateBesselJZero(nu + 1, k);
                double tt = Math.Sqrt(2.0 * (nu + 1.0) * (nu + 3.0));
                if (tt < jMax) jMax = tt;
            } else {
                // Use the interlacing property j_{\nu + 1, k - 1} < j_{\nu, k} < j_{\nu + 1, k}
                // and the approximation function to get a bound. The interlacing property is exact,
                // but our approximation functions are not, so it is concievable this could go wrong.
                jMin = ApproximateBesselJZero(nu + 1, k - 1);
                jMax = ApproximateBesselJZero(nu + 1, k);
            }
            Debug.Assert(jMin < jMax);
            // We should use internal method that takes best guess. Also, add derivative and use Newton's method.
            double j = FunctionMath.FindZero(x => BesselJ(nu, x), Interval.FromEndpoints(jMin, jMax));
            return (j);
        }

        private static double ApproximateBesselJZero (double nu, int k) {
            if (nu <= 2 * k) {
                // Since k >= 1, this covers all nu < 2.
                return (ApproximateBesselJZero_LargeIndex(nu, k));
            } else if (k == 1) {
                // Although it is better for large nu, I have verified that this gives
                // multi-digit accuracy for nu = 2, the lowest order for which it will be used.
                double cbrt_nu = Math.Pow(nu, 1.0 / 3.0);
                return (nu + 1.8557571 * cbrt_nu + 1.033150 / cbrt_nu);
            } else {
                // We use the complicated uniform asymptotic expansion only when we must.
                return (ApproximateBesselJZero_LargeOrder(nu, k));
            }
        }

        private static double ApproximateBesselJZero_LargeIndex (double nu, int k) {
            // I verified emperically that, over 0 < \nu < 256, I get 4-digit accuracy
            // for k ~ \nu and 2-digit accuracy for k ~ 1/2 \nu. 
            double mu = MoreMath.Sqr(2.0 * nu);
            double t = (k + 0.5 * nu - 0.25) * Math.PI;
            double t8 = 8.0 * t;
            double t82 = t8 * t8;
            return (t - (mu - 1.0) / t8 * (1.0 + 4.0 / 3.0 * (7.0 * mu - 31.0) / t82));
        }

        private static double ApproximateBesselJZero_LargeOrder (double nu, int k) {
            double a = AiryAiZero(k);
            double zeta = a / Math.Pow(nu, 2.0 / 3.0);
            double z = ZFromZeta(zeta, out double c1);
            return (z * (nu + c1 / nu));
        }

        private static double ZetaFromZ (double z) {
            Debug.Assert(z > 0.0);
            if (z < 0.75) {
                double s = Math.Sqrt((1.0 - z) * (1.0 + z));
                return (Math.Pow(3.0 / 2.0 * (Math.Log((1.0 + s) / z) - s), 2.0 / 3.0));
            } else if (z < 1.25) {
                double y = 1.0 - z;
                double c = Math.Pow(2.0, 1.0 / 3.0);
                // Need more terms
                return (c * y * (1.0 + 3.0 / 10.0 * y + 32.0 / 175.0 * y * y + 1037.0 / 7875 * y * y * y));
            } else {
                double s = Math.Sqrt((z - 1.0) * (z + 1.0));
                return (-Math.Pow(3.0 / 2.0 * (s - Math.Acos(1.0 / z)), 2.0 / 3.0));
            }
        }
        
        // This is an inversion of the zeta-from-z function.
        // Since it is used only in the approximate root expressions, we
        // don't try to achieve full accuracy. We achieve 6-8 digit
        // accuracy over most of the range.

        private static double ZFromZeta (double zeta, out double c1) {
            if (zeta < -0.375) {

                // Make an approximation based on the 1st two terms
                // of the small z (corresponding to large negative zeta)
                // development.
                double y = 2.0 / 3.0 * Math.Pow(-zeta, 3.0 / 2.0);
                double g = y + Math.PI / 2.0;
                double z = g - 1.0 / (2.0 * g);
                Debug.Assert(z > 1.0);

                // Do one or two Newton cycles to improve the value.
                for (int k = 0; k < 1; k++) {
                    double q = Math.Sqrt(z * z - 1.0);
                    double f = q - Math.Acos(1.0 / z);
                    double fPrime = (z - 1.0 / z) / q;
                    double dz = (y - f) / fPrime;
                    z += dz;
                }
                Debug.Assert(z > 1.0);

                // Compute c1
                double s = z * z - 1.0;
                double h2 = Math.Sqrt(-4.0 * zeta / s);
                double b0 = -5.0 / 48.0 / (zeta * zeta) + 1.0 / Math.Sqrt(-zeta) * (
                    5.0 / 24.0 / Math.Pow(s, 3.0 / 2.0) + 1.0 / 8.0 / Math.Sqrt(s)
                );
                c1 = h2 * b0 / 2.0;
                return (z);
            } else if (zeta < 0.375) {
                // This is an inversion of the series for \zeta(1+e).
                double cbrt2 = Math.Pow(2.0, 1.0 / 3.0);
                double s = zeta / cbrt2;
                double z = 1.0 + s * (-1.0 + s * (3.0 / 10.0 + s * (1.0 / 350.0 + s * (-479.0 / 63000.0 + s * (-20231.0 / 8085000.0)))));
                // Iteration to improve near z ~ 1 could actually reduce accuracy,
                // because of cancellation in terms of f near z ~ 1.
                c1 = 1.0 / 70.0 + s * (23.0 / 1575.0 + s * 838.0 / 121275.0);
                return (z);
            } else {
                throw new NotImplementedException();
            }
        }

    }

}
