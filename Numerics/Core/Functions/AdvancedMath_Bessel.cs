using System;
using System.Diagnostics;

using Meta.Numerics;

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
        /// <para>The Bessel functions of integral order occur in solutions to the wave equations with cylindrical symmetry. The
        /// regular Bessel functions are finite at the origin, and thus occur in situations where the wave equation is satisfied
        /// at the origin.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
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

            if (x < 4.0 + 2.0 * Math.Sqrt(n)) {
                // if x is small enough, use series
                return (BesselJ_Series(n, x));
            } else if (x > 32.0 + n * n / 2.0) {
                // if x is large enough, use asymptotic expansion
                return (Bessel_Asymptotic(n, x).FirstSolutionValue);
            } else if ((x > n) && (n > 32)) {
                // in this region, we can recurr upward from a smaller order m for which x is in the asymptotic region
                // recurrence is stable for x > n (and therefore also x > m) because J and Y are of the same magnitude
                // this is better than Miller's algorithm, because Miller's algorithm would require us to start from m >> x

                int m = (int) Math.Floor(Math.Sqrt(2.0 * (x - 32.0)));
                Debug.Assert(m <= n);

                SolutionPair s = Bessel_Asymptotic(m, x);
                double J = s.FirstSolutionValue; double JP = s.FirstSolutionDerivative;
                Bessel_RecurrUpward(m, x, ref J, ref JP, n - m);
                return (J);
            } else {
                // we are in the transition region; x is too large for the series but still less than n, so
                // we can't use upward recurrence

                // in the transition region, use Miller's algorithm:
                // recurr downward and use the sum rule to normalize the result

                // the height of the tower required depends on how many digits accuracy we need
                int nmax = n + 32 + 32 * (int) Math.Ceiling(Math.Sqrt(n+1));

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
        /// equation is satisfied not not include the origin.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
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
                    }
                    return (Y);
                }

            } else if (x > (32.0 + n * n / 2.0)) {
                // far enough out, use the asymptotic series
                return (Bessel_Asymptotic(n, x).SecondSolutionValue);
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

        // this is just an integer version of the series we implement below for doubles;
        // having an integer-specific version is slightly faster

        private static double BesselJ_Series (int n, double x) {
            double dJ = MoreMath.Pow(x / 2.0, n) / AdvancedIntegerMath.Factorial(n);
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

            Debug.Assert(Math.Abs(nu) <= 0.5);

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
            double z = x / 2.0;
            double dJ = Math.Pow(z, nu) / AdvancedMath.Gamma(nu + 1.0);
            double J = dJ;
            double zz = -z * z;
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

        // Regular Bessel function series about origin
        //   J_{\nu}(x) = (x/2)^{\nu} \sum{k=0}^{\infty} \frac{(-x^2/4)^k}{k! \Gamma(\nu + k + 1)}
        // This can be turned into a series for J' by term-by-term differentiation
        // For nu=0, x~1 it requires ~10 terms, x~4 requires ~17 terms, x~8 requires ~24 terms
        // Since nu appears in the denominator, higher nu converges faster
        // Since first correction term is ~\frac{x^2}{4\nu}, the radius of fast convergence should grow ~2\sqrt{\nu}

        private static void BesselJ_Series (double nu, double x, out double J, out double JP) {

            if (x == 0.0) {
                // treat the x = 0 case specially, becase various multipliciations and divisions by zero would give NaNs
                if (nu == 0.0) {
                    J = 1.0;
                    JP = 0.0;
                } else if (nu < 1.0) {
                    J = 0.0;
                    JP = Double.PositiveInfinity;
                } else if (nu == 1.0) {
                    J = 0.0;
                    JP = 0.5;
                } else {
                    J = 0.0;
                    JP = 0.0;
                }
            } else {
                // evaluate the series of J and J', knowing x != 0
                // note that the J' series is just with J series with teach term having one less power, and each term multipled
                // by the power it has in the J series
                double x2 = x / 2.0;
                double x22 = - x2 * x2;
                double dJ;
                if (nu < 128.0) {
                    dJ = Math.Pow(x2, nu) / AdvancedMath.Gamma(nu + 1.0);
                } else {
                    // if nu is very big, use log gamma
                    dJ = Math.Exp(nu * Math.Log(x2) - AdvancedMath.LogGamma(nu + 1.0));
                }
                J = dJ; JP = nu * dJ; 
                for (int k = 1; k < Global.SeriesMax; k++) {
                    double J_old = J; double JP_old = JP;
                    dJ *= x22 / k / (nu + k);
                    J += dJ; JP += (nu + 2 * k) * dJ;
                    if ((J == J_old) && (JP == JP_old)) {
                        JP = JP / x;
                        return;
                    }
                }
                // the J series alone requires 5 ops per term; J' evaulation adds 2 ops per term
                throw new NonconvergenceException();
            }
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
            for (int k = 2; k < kmax; k++) {
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
            for (int k = 3; k < Global.SeriesMax; k+=2) {
                Complex f_old = f;
                a = k * k / 4.0 - nu2;
                //b.Im = b.Im + 2.0;
                b = b + new Complex(0.0, 2.0);
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
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of J<sub>&#x3BD;</sub>(x).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
        public static double BesselJ (double nu, double x) {

            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            // use reflection to turn negative orders into positive orders
            if (nu < 0.0) {
                double mu = -nu;
                return (Math.Cos(Math.PI * mu) * BesselJ(mu, x) - Math.Sin(Math.PI * mu) * BesselY(mu, x));
            }

            if (x < 4.0 + Math.Sqrt(nu)) {
                // we are close enough to origin to use series
                return(BesselJ_Series(nu, x));
            } else if (x > 32.0 + nu * nu / 2.0) {
                // we are far enough from origin to use the asymptotic expansion
                SolutionPair result = Bessel_Asymptotic(nu, x);
                return (result.FirstSolutionValue);
            } else if (x > nu) {
                // we are far enough from origin to evaluate CF2, so use Steed's method
                SolutionPair result = Bessel_Steed(nu, x);
                return (result.FirstSolutionValue);
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

                return (result.FirstSolutionValue / J);

            }
        }


        /// <summary>
        /// Computes the irregual Bessel function for real orders.
        /// </summary>
        /// <param name="nu">The order parameter.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of Y<sub>&#x3BD;</sub>(x).</returns>
        /// <remarks>
        /// <para>For information on Bessel functions, see <see cref="Bessel"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="Bessel"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
        public static double BesselY (double nu, double x) {

            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            if (nu < 0.0) {
                double np = nu * Math.PI;
                return (Math.Cos(np) * BesselY(-nu, x) - Math.Sin(np) * BesselJ(-nu, x));
            }

            if (x == 0.0) {
                return (Double.NegativeInfinity);
            } else if (x > 32.0 + nu * nu / 2.0) {

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
                /*
                for (int i = 0; i < n; i++) {
                    // use recurrence on Y, Y':
                    //   Y_{mu+1} = (mu/x) Y_{mu} - Y_{mu}'
                    //   Y_{mu}' = Y_{mu-1} - (mu/x) Y_{mu}
                    double t = Y;
                    Y = (mu / x) * t - YP;
                    mu += 1.0;
                    YP = t - (mu / x) * Y;
                }
                */
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
        /// <para>Bessel functions often occur in physical phenomenon with cylindrical symmetry. They satisfy the differential equation</para>
        /// <img src="../images/BesselODE.png" />
        /// <para>Since this is second order linear equation, is has two linearly independent solutions. The regular Bessel function
        /// J<sub>nu</sub>(x), which is regular at the origin, and the irregular Bessel function Y<sub>nu</sub>(x), which diverges
        /// at the origin.</para>
        /// <para>This method simultaneously computes both Bessel functions and their derivatives. If you need both J and Y, it is faster to call this method once than to call
        /// <see cref="BesselJ(double,double)"/> and <see cref="BesselY(double,double)"/> seperately. If on, the other hand, you need only J or only Y, it is faster to
        /// call the appropirate method to compute only the one you need.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="nu"/> or <paramref name="x"/> is negative.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Bessel_function"/>
        public static SolutionPair Bessel (double nu, double x) {

            if (nu < 0.0) throw new ArgumentOutOfRangeException("nu");
            if (x < 0.0) throw new ArgumentOutOfRangeException("x");

            if (x < 4.0) {

                // close to the origin, use the series for Y and J
                // the range does not grow with increasing order because our series for Y is only for small orders

                // for J, we can use the series directly
                double J, JP; BesselJ_Series(nu, x, out J, out JP);

                // for Y, we only have a series for -0.5 < mu < 0.5, so pick an order in that range offset from our
                // desired order by an integer number of recurrence steps. Evaluate at the mu order and recurr upward to nu.
                double mu = nu - Math.Round(nu);
                int n = (int) Math.Round(nu - mu);
                double Y, YP1;
                BesselY_Series(mu, x, out Y, out YP1);
                double YP = (mu / x) * Y - YP1;
                Bessel_RecurrUpward(mu, x, ref Y, ref YP, n);

                // return J, J', Y, Y' in a solution pair
                return(new SolutionPair(J, JP, Y, YP));

            } else if (x > 32.0 + nu * nu / 2.0) {
                // far enough out, use the asymptotic series
                return (Bessel_Asymptotic(nu, x));
            } else if (x > nu) {
                // in the intermediate region, use Steed's method directly if x > nu
                return (Bessel_Steed(nu, x));
            } else {

                // if x < nu, we can't use Steed's method directly because CF2 does not converge
                // so pick a mu ~ x < nu, offset from nu by an integer; use Steed's method there, then recurse up to nu

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

            double g = (z.Re - r)/z.Im;

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
        //   J = \sqrt{\frac{2}{\pi x}} \left[ P \cos\phi - Q \sin\phi \right]
        //   Y = \sqrt{\frac{2}{\pi x}} \left[ P \sin\phi + Q \cos\phi \right]
        // where \phi = x - \left( \nu / 2 + 1 / 4 \right) \pi and \mu = 4 \nu^2 and
        //   P = 1 - \frac{(\mu-1)(\mu-9)}{2! (8x)^2} + \frac{(\mu-1)(\mu-9)(\mu-25)(\mu-49}{4! (8x)^4} + \cdots
        //   Q = \frac{(\mu-1)}{8x} - \frac{(\mu-1)(\mu-9)(\mu-25)}{3! (8x)^3} + \cdots
        // Derivatives have similiar expressions
        //   J' = - \sqrt{\frac{2}{\pi x}} \left[ R \sin\phi + S \cos\phi \right]
        //   Y' = \sqrt{\frac{2}{\pi x}} \left[ R \cos\phi - S \sin\phi \right]
        // where
        //   R = 1 - \frac{(\mu-1)(\mu+15)}{2! (8x)^2} + \cdots
        //   S = \frac{(\mu+3)}{8x} - \frac{(\mu-1)(\mu - 9)(\mu+35)}{3! (8x)^3} + \cdots

        // For nu=0, this series converges to full precision in about 10 terms at x~100, and in about 25 terms even as low as x~25
        // It fails to converge at all for lower x <~ 25
        // Since the first correction term is ~ (4 \nu^2)^2 / (8 x)^2 ~ (\nu^2 / 2 x)^2, the minimum x should grow like \nu^2 / 2

        private static SolutionPair Bessel_Asymptotic (double nu, double x) {

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

                if ((Q == Q_old) && (P == P_old) && (R == R_old) && (S == S_old)) {
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

            SolutionPair result = new SolutionPair(N * (c * P - s * Q), -N * (R * s + S * c), N * (s * P + c * Q), N * (R * c - S * s));

            return (result);

        }

        // large order expansion using Airy functions
        /*
        public static double BesselJ_Airy (double nu, double x) {

            double z = Bessel_Zeta(x / nu);
            double nu3 = Math.Pow(nu, 1.0 / 3.0);
            double nu23 = nu3 * nu3;
            double x = nu23 * z;

            double ai, aip; // airy functions

            throw new NotImplementedException();

        }
        */

        /*
        public static double Bessel_Zeta (double z) {
            if (z < 1.0) {
                double sz = Math.Sqrt((1.0 + z) * (1.0 - z));
                return (Math.Log((1.0 + sz) / z) - sz);
            } else {
                double sz = Math.Sqrt((z + 1.0) * (z - 1.0));
                return (sz - Math.Acos(1.0 / z));
            }
        }
        */

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
        /// <seealso cref="SphericalBesselY" />
        /// <seealso cref="BesselJ(double,double)"/>
        /// <seealso href="http://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html" />
        public static double SphericalBesselJ (int n, double x) {

            if (x < 0.0) {
                if (n % 2 == 0) {
                    return (SphericalBesselJ(n, -x));
                } else {
                    return (-SphericalBesselJ(n, -x));
                }
            }

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
                } else if (x > (32.0 + n * n / 2.0)) {
                    // far enough from the origin, use the asymptotic expansion
                    return (Math.Sqrt(Global.HalfPI / x) * Bessel_Asymptotic(n+0.5, x).FirstSolutionValue);
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
        /// <seealso cref="SphericalBesselJ"/>
        /// <seealso href="http://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html" />
        public static double SphericalBesselY (int n, double x) {

            if (x < 0.0) {
                if (n % 2 == 0) {
                    return (-SphericalBesselY(n, -x));
                } else {
                    return (SphericalBesselY(n, -x));
                }
            }

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
                    return (Math.Sqrt(Global.HalfPI / x) * Bessel_Asymptotic(n + 0.5, x).SecondSolutionValue);
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
            if (Math.Abs(x) < 0.5) {
                return(SphericalBesselJ_SeriesOne(x));
            } else if (Math.Abs(x) > 100.0) {
                return (Math.Sqrt(Global.HalfPI / x) * Bessel_Asymptotic(1.5, x).FirstSolutionValue);
            } else {
                return((Math.Sin(x) / x - Math.Cos(x)) / x);
            }
        }

        private static double SphericalBesselJ_Series (int n, double x) {
            double xx = x * x / 2.0;
            double df = Math.Exp(n * Math.Log(x) - AdvancedIntegerMath.LogDoubleFactorial(2 * n + 1));
            double f = df;
            for (int i = 1; i < Global.SeriesMax; i++) {
                double f_old = f;
                df = -df * xx / i / (2 * (n + i) + 1);
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselY_Series (int n, double x) {
            double xx = x * x / 2.0;
            double df = - AdvancedIntegerMath.DoubleFactorial(2 * n - 1) / MoreMath.Pow(x, n + 1);
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
                return (Math.Sqrt(Global.HalfPI / x) * Bessel_Asymptotic(1.5, x).SecondSolutionValue);
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
            }
            double jn = j;
            for (int k = n; k > 0; k--) {
                double jm1 = (2 * k + 1) / x * j - jp1;
                jp1 = j;
                j = jm1;
            }

            // compute the value we should have got and use it to normalize our result
            double j0 = SphericalBesselJ_Zero(x);
            return ((j0 / j) * jn);

        }

    }


    /// <summary>
    /// Contains a pair of solutions to a differential equation.
    /// </summary>
    /// <remarks>
    /// <para>Any linear second order differential equation has two independent solutions. For example,
    /// the Bessel differential equation (<see cref="AdvancedMath.Bessel"/>) has solutions J and Y,
    /// the Coulomb wave equation has solutions F and G,
    /// and the Airy differential equation has solutions Ai and Bi.</para>
    /// <para>A solution pair structure contains values for both solutions and for their derivatives. It is often useful to
    /// have all this information together when fitting boundary conditions.</para>
    /// <para>Which solution is considered the first and which is considered the second is
    /// a matter of convention. When one solution is regular (finite) at the origin and the other is not, we take the regular solution
    /// to be the first.</para>
    /// </remarks>
    public struct SolutionPair {

        private double j, jPrime, y, yPrime;

        /// <summary>
        /// Gets the value of the first solution.
        /// </summary>
        public double FirstSolutionValue {
            get {
                return (j);
            }
            internal set {
                j = value;
            }
        }

        /// <summary>
        /// Gets the derivative of the first solution.
        /// </summary>
        public double FirstSolutionDerivative {
            get {
                return (jPrime);
            }
            internal set {
                jPrime = value;
            }
        }

        /// <summary>
        /// Gets the value of the second solution.
        /// </summary>
        public double SecondSolutionValue {
            get {
                return (y);
            }
            internal set {
                y = value;
            }
        }

        /// <summary>
        /// Gets the derivative of the second solution.
        /// </summary>
        public double SecondSolutionDerivative {
            get {
                return (yPrime);
            }
            internal set {
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
        
       internal SolutionPair (double j, double jPrime, double y, double yPrime) {
            this.j = j;
            this.jPrime = jPrime;
            this.y = y;
            this.yPrime = yPrime;
        }

    }

}
