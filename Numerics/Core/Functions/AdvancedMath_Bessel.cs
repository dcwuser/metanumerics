using System;
using System.Diagnostics;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        // integer Bessel functions
        // appear in cases of cylindrical symmetry

        // **** Integer order Bessel functions ****

        /// <summary>
        /// Computes the regular Bessel function for integer orders.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of J<sub><paramref name="n"/></sub>(<paramref name="x"/>).</returns>
        /// <remarks>
        /// <para>The Bessel functions of integral order occur in solutions to the wave equations with cylindrical symmetry. The
        /// regular Bessel functions are finite at the origin, and thus occur in situations where the wave equation is satisfied
        /// at the origin.</para></remarks>
        public static double BesselJ (int n, double x) {

            // relate negative n to positive n
            if (n < 0) {
                if ((n % 2) == 0) {
                    return (BesselJ(-n, x));
                } else {
                    return (-BesselJ(-n, x));
                }
            }

            // relate negative x to positive x
            if (x < 0) {
                if ((n % 2) == 0) {
                    return (BesselJ(n, -x));
                } else {
                    return (-BesselJ(n, -x));
                }
            }

            if (x < (3.0 + 2.0 * Math.Sqrt(n))) {
                // if x is small enough, use series
                return (BesselJ_Series(n, x));
            } else if (x > (30.0 + n * n / 2.0)) {
                // if x is large enough, use asymptotic expansion
                return (Bessel_Asymptotic(n, x).Regular);
            } else {
                // in the transition region, use Miller's algorithm:
                // recurr downward and use the sum rule to normalize the result

                // the height of the tower required depends on how many digits accuracy we need
                int nmax = n + 40 * ((int)Math.Ceiling(Math.Sqrt(n + 1)));

                double Jp1 = 0.0;
                double J = 1.0 / Math.Sqrt(x); // correct order of magnitude
                double sum = 0.0;
                for (int k = nmax; k > n; k--) {
                    if ((k % 2) == 0) sum += J;
                    double Jm1 = (2 * k) / x * J - Jp1;
                    Jp1 = J;
                    J = Jm1;
                }
                double Jn = J;
                for (int k = n; k > 0; k--) {
                    if ((k % 2) == 0) sum += J;
                    double Jm1 = (2 * k) / x * J - Jp1;
                    Jp1 = J;
                    J = Jm1;
                }
                sum = 2.0 * sum + J;
                return (Jn / sum);
            }
        }

        /// <summary>
        /// Computes the irregular Bessel function for integer orders.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of Y<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The Bessel functions of integral order occur in solutions to the wave equations with cylindrical symmetry. The
        /// irregular Bessel functions diverge at the origin, and thus occur in situations where the region in which the wave
        /// equation is satisfied not not include the origin.</para></remarks>
        public static double BesselY (int n, double x) {

            // relate negative n to positive n
            if (n < 0) {
                if ((n % 2) == 0) {
                    return (BesselY(-n, x));
                } else {
                    return (-BesselY(-n, x));
                }
            }

            // relate negative x to positive x
            if (x < 0) {
                if ((n % 2) == 0) {
                    return (BesselY(n, -x));
                } else {
                    return (-BesselY(n, -x));
                }
            }


            if (x == 0.0) {
                return (Double.NegativeInfinity);
            } else if (x < 5.0) {

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
                    }
                    return (Y);
                }

            } else if (x > (30.0 + 0.5 * n * n)) {
                return (Bessel_Asymptotic(n, x).Irregular);
            } else {

                // use Miller's algorithm to compute a tower of J's
                // use the tower of J's to compute Y0
                // use the Wronskian to compute Y1

                // this doesn't seem to be working; use the real routine for now
                return (BesselY((double) n, x));
                
                /*
                double sum1 = 0.0;
                double sum2 = 0.0;

                double Jp1 = 0.0;
                double J = 1.0 / Math.Sqrt(x);
                for (int k = 40; k > 0; k--) {
                    if ((k % 2) == 0) {
                        sum1 += J;
                        if ((k % 4) == 0) {
                            sum2 += J / k;
                        } else {
                            sum2 -= J / k;
                        }
                    }
                    double Jm1 = (2 * k) / x * J - Jp1;
                    Jp1 = J;
                    J = Jm1;
                }
                sum1 = 2.0 * sum1 + J;
                double J0 = J / sum1;
                double J1 = Jp1 / sum1;
                sum2 = sum2 / sum1;
                double Y0 = (2.0 * J0 / Math.PI) * (Math.Log(x / 2.0) + EulerGamma) - (8.0 / Math.PI) * sum2;
                double Y1 = (J1 * Y0 - (2.0/Math.PI/x))/J0;

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
                    }
                    return (Y);
                }
                */

            }
            
        }

        private static double BesselJ_Series (int n, double x) {
            double dJ = Math.Pow(x / 2.0, n) / AdvancedIntegerMath.Factorial(n);
            double J = dJ;
            double zz = - x * x / 4.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double J_old = J;
                dJ = dJ * zz / ((n + k) * k);
                J += dJ;
                if (J == J_old) {
                    return (J);
                }
            }
            throw new NonconvergenceException();
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

            if (x == 0.0) {
                Y0 = Double.NegativeInfinity;
                Y1 = Double.NegativeInfinity;
                return;
            }

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

                // Console.WriteLine(String.Format("k = {0}, Y0={1}", k, -Y0));

                if ((Y0 == Y0_old) && (Y1 == Y1_old)) {
                    Y0 = -Y0;
                    Y1 = -2.0 * Y1 / x;
                    return;
                }
            }
            throw new NonconvergenceException();

        }


        // **** Real order Bessel functions ****

        // works for z << Sqrt(nu), say z < Sqrt(nu)+1
        // a series expansion in x, but works pretty well quite far out
        // for nu=0, requires 10 terms at x~1.0, 25 terms at x~10.0, but accuracy in the last few digits suffers that far out
        // gets better for higher nu (for nu=10, only 20 terms are required at x~10.0 and all digits are good), but only with Sqrt(nu)

        private static double BesselJ_Series (double nu, double x) {
            double dJ = Math.Pow(0.5 * x, nu) / AdvancedMath.Gamma(nu + 1.0);
            double J = dJ;
            double zz = -0.25 * x * x;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double J_old = J;
                dJ = dJ * zz / (nu + k) / k;
                J += dJ;
                if (J == J_old) {
                    return (J);
                }
            }
            throw new NonconvergenceException();
        }

        // computes the ratio f=J_{nu+1}/J_{nu} via a continued fraction
        // converges rapidly for x < Sqrt(nu(nu+1)), and slowly for larger x
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

            // move to higher fractions
            for (int k = 2; k < Global.SeriesMax; k++) {
                double a = -1.0;
                double b = 2.0 * (nu +k ) / x;
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
            for (int k = 3; k < Global.SeriesMax; k+=2) {
                Complex f_old = f;
                a = k * k / 4.0 - nu2;
                b.Im = b.Im + 2.0;
                D = 1.0 / (b + a * D);
                Df = (b * D - 1.0) * Df;
                f += Df;
                if (f == f_old) return (-0.5/x + ComplexMath.I + ComplexMath.I * f / x);
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes the regular Bessel function for real orders.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of J<sub>&#x3BD;</sub>(x).</returns>
        public static double BesselJ (double nu, double x) {
            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            if (nu < 0.0) {
                double mu = -nu;
                return (Math.Cos(Math.PI * mu) * BesselJ(mu, x) - Math.Sin(Math.PI * mu) * BesselY(mu, x));
            }

            if (x < (4.0 + Math.Sqrt(nu))) {
                // we are close enough to origin to use series
                return(BesselJ_Series(nu, x));
            } else if (x > 30.0 + 0.5 * nu * nu) {
                // we are far enough from origin to use the asymptotic expansion
                SolutionPair result = Bessel_Asymptotic(nu, x);
                return (result.Regular);
            } else {
                if (x > nu) {
                    // we are far enough from origin to evaluate CF2, so use Steed's method
                    SolutionPair result = Bessel_Steed(nu, x);
                    return (result.Regular);
                } else {
                    // we have x < nu, but x is still not small enough to use series; this only occurs for nu >~ 6
                    // to handle this case, compute J_{nu+1}/J_{nu}, recurse down to mu where mu ~ x < nu; use
                    // Steed's method to evaluate J_{mu} and re-normalize J_{nu}.

                    // for example, for mu = 16, x = 12, we can't evaluate CF2 because x < 16; so assume J_16=1, compute J_17 / J_16,
                    // and recurse down to J_11; since 12 > 11, we can compute CF2 and get J_11 and Y_11; comparing this
                    // with our value of J_11 gives the proper renormalization factor for J_16

                    int sign;
                    double r = Bessel_CF1(nu, x, out sign);
                    //Console.WriteLine("r={0}", r);

                    double J = 1.0 * sign;
                    double Jp1 = r * J;
                    double mu = nu;
                    while (mu > x) {
                        double Jm1 = (2 * mu) / x * J - Jp1;
                        Jp1 = J;
                        J = Jm1;
                        mu = mu - 1.0;
                    }
                    //Console.WriteLine("derived J = {0} at mu={1}", J, mu);

                    Complex z = Bessel_CF2(mu, x);
                    SolutionPair result = Bessel_Steed(mu / x - Jp1 / J, z, 2.0 / Math.PI / x , Math.Sign(J));
                    //Console.WriteLine("actual J = {0} at mu={1}", result.Regular, mu);
                    //Console.WriteLine("scale factor {0}", result.Regular / J);

                    return (result.Regular / J);

                }
            }
        }


        /// <summary>
        /// Computes the irregual Bessel function for real orders.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of Y<sub>&#x3BD;</sub>(x).</returns>
        public static double BesselY (double nu, double x) {
            if (nu < 0.0) throw new ArgumentOutOfRangeException("nu");
            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            if (nu < 0.0) {
                double np = nu * Math.PI;
                return (Math.Cos(np) * BesselY(-nu, x) - Math.Sin(np) * BesselJ(-nu, x));
            }

            if (x == 0.0) {
                return (Double.NegativeInfinity);
            } else if (x > 30.0 + nu * nu) {

                // far from the origin, use the asymptotic expansion

                SolutionPair result = Bessel_Asymptotic(nu, x);
                return (result.Irregular);
            } else if (x < 3.0) {

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
                    }
                    return (Y1);

                }

            } else if (x > nu) {
                SolutionPair result = Bessel_Steed(nu, x);
                return (result.Irregular);
            } else {

                // we have 3 < x < nu; evaluate at mu~x and recurr upward to nu

                // figure out mu
                int n = (int) Math.Ceiling(nu - x);
                double mu = nu - n;

                // evaluate at mu
                SolutionPair result = Bessel_Steed(mu, x);
                double Y = result.Irregular;
                double YP = result.IrregularPrime;

                // recurr upward to nu
                for (int i = 0; i < n; i++) {
                    // use recurrence on Y, Y':
                    //   Y_{mu+1} = (mu/x) Y_{mu} - Y_{mu}'
                    //   Y_{mu}' = Y_{mu-1} - (mu/x) Y_{mu}
                    double t = Y;
                    Y = (mu / x) * t - YP;
                    mu += 1.0;
                    YP = t - (mu / x) * Y;
                }

                return (Y);

            }

        }

        // CF1, CF2, and Wronskian to determine Bessel function values
        private static SolutionPair Bessel_Steed (double nu, double x) {

            // CF2 only converges for x > nu
            Debug.Assert(x > nu);

            // Wronskian
            double W = 2.0 / Math.PI / x;

            // CF1
            int sign;
            double r = nu / x - Bessel_CF1(nu, x, out sign);

            // CF2
            Complex z = Bessel_CF2(nu, x);

            // use these results to compute J, JP, Y, and YP
            // what about zeros?

            double g = (z.Re - r)/z.Im;

            double J = Math.Sqrt(W / (z.Im + g * (z.Re - r)));
            if (sign < 0) J = -J;

            double JP = r * J;

            double Y = g * J;

            double YP = Y * (z.Re + z.Im / g);

            return new SolutionPair(nu, x, J, JP, Y, YP);

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

            return new SolutionPair(0, 0.0, J, JP, Y, YP);

        }

        // Hankel's asymptotic expansions
        // for nu=0, converges in about 10 terms at x~100, and in about 25 terms even as low as x~25 at full precision,
        // but fails to converge at all for lower x; minimum x should grow about like 1/2 * nu^2

        private static SolutionPair Bessel_Asymptotic (double nu, double x) {

            double n2 = 4.0 * nu * nu;
            double xx = 8.0 * x;

            double a = 1.0;
            double b = 0.0;

            int k = 0;
            int kn = -1;
            double d = 1.0;

            while (true) {

                k++;
                kn += 2;
                d = d * (n2 - kn * kn) / k / xx;
                double b_old = b;
                b = b_old + d;

                k++;
                kn += 2;
                d = - d * (n2 - kn * kn) / k / xx;
                double a_old = a;
                a = a_old + d;

                if ((b == b_old) && (a == a_old)) {
                    break;
                }

                if (k > Global.SeriesMax) throw new NonconvergenceException();

            }

            //double t = x - ( 0.5 * nu + 0.25 ) * Math.PI;
            //double s = Math.Sin(t);
            //double c = Math.Cos(t);
            double s = Sin(x, -(nu + 0.5)/4.0);
            double c = Cos(x, -(nu + 0.5)/4.0);
            //Console.WriteLine("t={0}  s={1}  c={2}", t, s, c);
            //Console.WriteLine("t={0} sa={1} ca={2}", t, sa, ca);

            double N = Math.Sqrt(2.0 / Math.PI / x);

            //Debug.WriteLine(String.Format("x={0} nu={1} ca={2} sb={3}", x, nu, c * a, s * b));
            SolutionPair result = new SolutionPair(nu, x, N * (c * a - s * b), 0.0, N * (s * a + c * b), 0.0);
            //BesselResult result = new BesselResult(nu, x, N * Math.Sin(t - p), 0.0, N * Math.Sin(t + p), 0.0);

            return (result);

        }

        // **** Spherical Bessel functions ****

        /// <summary>
        /// Computes the regular spherical Bessel function of integer order.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns> The value of j<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The spherical Bessel functions occur in solutions to the wave equations with spherical symmetry. The
        /// regular sperhical Bessel functions are finite at the origin, and thus occur in situations where the wave equation is satisfied
        /// at the origin.</para>
        /// <para>The regular spherical Bessel functions are related to the regular Bessel functions of half-integer order by
        /// j<sub>n</sub>(x) = Sqrt(&#x3C0;/2x) J<sub>n+1/2</sub>(x).</para></remarks>
        /// <seealso cref="BesselJ(double,double)"/>
        public static double SphericalBesselJ (int n, double x) {

            if (n < 0) {
                if ((n % 2) == 0) {
                    return(SphericalBesselY(-n-1, x));
                } else {
                    return(-SphericalBesselY(-n-1, x));
                }
            } else if (n == 0) {
                return (SphericalBesselJ_Zero(x));
            } else if (n == 1) {
                return (SphericalBesselJ_One(x));
            } else {

                if (x < 2.0 + 2.0 * Math.Sqrt(n)) {
                    // close enough to the origin, use the power series
                    return (SphericalBesselJ_Series(n, x));
                } else if (x > (30.0 + n * n / 2.0)) {
                    // far enough from the origin, use the asymptotic expansion
                    return (Math.Sqrt(Math.PI/2.0/x) * Bessel_Asymptotic(n+0.5, x).Regular);
                } else {
                    // in the transition region, use Miller's algorithm
                    return (SphericalBesselJ_Miller(n, x));
                }
            }
        }

        /// <summary>
        /// Computes the irregular spherical Bessel function of integer order.
        /// </summary>
        /// <param name="n">The order parameter.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of y<sub>n</sub>(x).</returns>
        public static double SphericalBesselY (int n, double x) {

            if (n < 0) {
                if ((n % 2) == 0) {
                    return(-SphericalBesselJ(-n-1, x));
                } else {
                    return(SphericalBesselJ(-n-1, x));
                }
            } else if (n == 0) {
                return (SphericalBesselY_Zero(x));
            } else if (n == 1) {
                return (SphericalBesselY_One(x));
            } else {
                if (x < (2.0 + Math.Sqrt(n))) {
                    return (SphericalBesselY_Series(n, x));
                } else if (x > (30.0 + 0.5 * n * n)) {
                    // if x is large enough, use asymptotic expansion
                    return (Math.Sqrt(Math.PI / 2.0 / x) * Bessel_Asymptotic(n + 0.5, x).Irregular);
                } else {
                    // move up using the recursion relation
                    double ym1 = SphericalBesselY_Zero(x);
                    double y = SphericalBesselY_One(x);
                    for (int k = 1; k < n; k++) {
                        double yp1 = (2 * k + 1) / x * y - ym1;
                        ym1 = y;
                        y = yp1;
                    }
                    return (y);
                }
            }

        }

        private static double SphericalBesselJ_SeriesZero (double x) {
            double xx = x * x / 2.0;
            double dj = 1.0;
            double j = dj;
            for (int i = 1; i < Global.SeriesMax; i++) {
                double j_old = j;
                dj = - dj * xx / i / (2 * i + 1);
                j = j_old + dj;
                if (j == j_old) return (j);
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselJ_Zero (double x) {
            if (Math.Abs(x) < 0.25) {
                return (SphericalBesselJ_SeriesZero(x));
            } else {
                return (Math.Sin(x) / x);
            }
        }

        private static double SphericalBesselJ_SeriesOne (double x) {
            double xx = x * x / 2.0;
            double dj = x / 3.0;
            double j = dj;
            for (int i = 1; i < Global.SeriesMax; i++) {
                double j_old = j;
                dj = -dj * xx / i / (2 * i + 3);
                j = j_old + dj;
                if (j == j_old) return (j);
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselJ_One (double x) {
            if (Math.Abs(x) < 1.0) {
                return(SphericalBesselJ_SeriesOne(x));
            } else if (Math.Abs(x) > 100.0) {
                return (Math.Sqrt(Math.PI / 2.0 / x) * Bessel_Asymptotic(1.5, x).Regular);
            } else {
                return((Math.Sin(x)/x - Math.Cos(x))/x);
            }
        }

        private static double SphericalBesselJ_Series (int n, double x) {
            //Console.WriteLine("Series n={0} x={1}", n, x);
            double xx = x*x/2.0;
            double df = Math.Exp(n * Math.Log(x) - AdvancedIntegerMath.LogDoubleFactorial(2 * n + 1));
            //double df = Math.Pow(x, n) / AdvancedIntegerMath.DoubleFactorial(2*n+1);
            double f = df;
            for (int i = 1; i < Global.SeriesMax; i++) {
                double f_old = f;
                df = -df * xx / i / (2 * (n + i) + 1);
                //Console.WriteLine("f={0}, df={1}", f, df);
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselY_Series (int n, double x) {
            double xx = x * x / 2.0;
            double df = - AdvancedIntegerMath.DoubleFactorial(2 * n - 1) / Math.Pow(x, n + 1);
            double f = df;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                df = -df * xx / k / (2 * (k - n) - 1);
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        /*
        private static double DoubleFactorial (int n) {
            if (n < 20) {
                double f = 1.0;
                for (int i = n; i > 1; i -= 2) {
                    f = f * i;
                }
                return (f);
            } else {
                int n2 = n / 2;
                return (Math.Exp(LogGamma(n + 1) - n2 * Math.Log(2.0) - LogGamma(n2 + 1)));
            }
        }
        */

        private static double SphericalBesselY_SeriesZero (double x) {
            if (x == 0) return (Double.NegativeInfinity);
            double xx = x * x / 2.0;
            double dy = - 1.0 / x;
            double y = dy;
            for (int i = 1; i < Global.SeriesMax; i++) {
                double y_old = y;
                dy = - dy * xx / i / (2 * i - 1);
                y = y_old + dy;
                if (y == y_old) return (y);
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselY_Zero (double x) {
            if (Math.Abs(x) < 0.25) {
                return (SphericalBesselY_SeriesZero(x));
            } else {
                return (-Math.Cos(x) / x);
            }
        }

        private static double SphericalBesselY_SeriesOne (double x) {
            if (x == 0) return (Double.NegativeInfinity);
            double xx = x * x / 2.0;
            double dy = - 1.0 / (x * x);
            double y = dy;
            for (int i = 1; i < Global.SeriesMax; i++) {
                double y_old = y;
                dy = - dy * xx / i / (2 * i - 3);
                y = y_old + dy;
                if (y == y_old) return (y);
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselY_One (double x) {
            if (Math.Abs(x) < 1.0) {
                return (SphericalBesselY_SeriesOne(x));
            } else if (Math.Abs(x) > 100.0) {
                return (Math.Sqrt(Math.PI / 2.0 / x) * Bessel_Asymptotic(1.5, x).Irregular);
            } else {
                return(-(Math.Cos(x)/x + Math.Sin(x))/x);
            }
        }

        // Miller's method assumes a value at some high N and recurrs downward
        // the result is then normalized using a sum relation or a known value

        private static double SphericalBesselJ_Miller (int n, double x) {

            // pick starting value for the downward recursion that
            // takes us well into the x < nu regime, where J_nu decreases with nu
            int kmax = n;
            if (x > n) kmax = (int) Math.Ceiling(x);
            kmax += 50; // since J_(nu+1)/J_(nu) ~ 1/2 for x~v, taking N steps supresses by 2^(N) = 10^(16) at N ~ 50

            double jp1 = 0.0;
            double j = 1.0 / ((double) kmax); // look for a better guess

            // recur downward to order zero
            // the recurrence j_{k-1} = (2k+1)/x * j_k - j_{k+1} is stable in this direction
            for (int k = kmax; k > n; k--) {
                double jm1 = (2 * k + 1) / x * j - jp1;
                jp1 = j;
                j = jm1;
                //Console.WriteLine("{0} {1}", k, j);
            }
            double jn = j;
            for (int k = n; k > 0; k--) {
                double jm1 = (2 * k + 1) / x * j - jp1;
                jp1 = j;
                j = jm1;
                //Console.WriteLine("{0} {1}", k, j);
            }

            // compute the value we should have got and use it to normalize our result
            double j0 = SphericalBesselJ_Zero(x);
            //Console.WriteLine("j0={0}", j0);
            return ((j0 / j) * jn);

        }

        /*
        private static BesselResult SphericalBessel_Zero (double x) {
            if (Math.Abs(x) < 0.25) {
                return (new BesselResult(0, x, SphericalBesselJ_SeriesZero(x), 0.0, 0.0, 0.0));
            } else {
                double s = Math.Sin(x) / x;
                double c = Math.Cos(x) / x;
                return (new BesselResult(0, x, s, s / x - c, -c, -c / x - s));
            }
        }
        */

        /*
        // an asymptotic series good for the phase angle when x >> nu^2, say x > 20 + nu^2
        private static double Bessel_AsymptoticPhase (double nu, double x) {
            double mu = 4.0 * nu * nu;
            double x2 = 16.0 * x * x;
            double theta = x - (0.5 * nu + 0.25) * Math.PI;
            theta += (mu - 1.0) / 2.0 / (4.0 * x);
            theta += (mu - 1.0) * (mu - 25.0) / 6.0 / Math.Pow(4.0 * x, 3.0);
            theta += (mu - 1.0) * (mu * mu - 114.0 * mu + 1073.0) / 5.0 / Math.Pow(4.0 * x, 5.0);
            return (theta);
        }
        */

    }


    internal struct SolutionPair {

        private double nu, x, j, jPrime, y, yPrime;

        public double Order {
            get {
                return (nu);
            }
        }

        public double Argument {
            get {
                return (x);
            }
        }

        public double Regular {
            get {
                return (j);
            }
            set {
                j = value;
            }
        }

        public double RegularPrime {
            get {
                return (jPrime);
            }
            set {
                jPrime = value;
            }
        }

        public double Irregular {
            get {
                return (y);
            }
            set {
                y = value;
            }
        }

        public double IrregularPrime {
            get {
                return (yPrime);
            }
            set {
                yPrime = value;
            }
        }

        /*
        public double Wronskian {
            get {
                return (j * yPrime - y * jPrime);
            }
        }
        */

        internal SolutionPair (double nu, double x, double j, double jPrime, double y, double yPrime) {
            this.nu = nu;
            this.x = x;
            this.j = j;
            this.jPrime = jPrime;
            this.y = y;
            this.yPrime = yPrime;
        }

    }

}
