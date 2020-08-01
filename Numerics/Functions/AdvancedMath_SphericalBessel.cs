using System;
using System.Diagnostics;

namespace Meta.Numerics.Functions {


    public static partial class AdvancedMath {

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
                    return SphericalBesselJ(n, -x);
                } else {
                    return -SphericalBesselJ(n, -x);
                }
            }

            if (n < 0) {
                if ((n % 2) == 0) {
                    return SphericalBesselY(-n - 1, x);
                } else {
                    return -SphericalBesselY(-n - 1, x);
                }
            } else if (n == 0) {
                return SphericalBesselJ_Zero(x);
            } else if (n == 1) {
                return SphericalBesselJ_One(x);
            } else {
                if (x <= Math.Sqrt(2 * n + 3)) {
                    // Close enough to the origin, use the power series.
                    return SphericalBesselJ_Series(n, x);
                } else if (x < (32.0 + 0.5 * n * n)) {
                    // In the transition region, use Miller's algorithm.
                    return SphericalBesselJ_Miller(n, x);
                } else {
                    // Far enough from the origin, use the asymptotic expansion.
                    return Math.Sqrt(Math.PI / 2.0 / x) * Bessel_Asymptotic(n + 0.5, x).FirstSolutionValue;
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
                    return -SphericalBesselY(n, -x);
                } else {
                    return SphericalBesselY(n, -x);
                }
            }

            if (n < 0) {
                if ((n % 2) == 0) {
                    return -SphericalBesselJ(-n - 1, x);
                } else {
                    return SphericalBesselJ(-n - 1, x);
                }
            } else if (n == 0) {
                return SphericalBesselY_Zero(x);
            } else if (n == 1) {
                SphericalBesselY_ZeroAndOne(x, out double _, out double y1);
                return y1;
            } else {
                if (x < (2.0 + Math.Sqrt(n))) {
                    return SphericalBesselY_Series(n, x);
                } else if (x > (30.0 + 0.5 * n * n)) {
                    // if x is large enough, use asymptotic expansion
                    return (Math.Sqrt(Math.PI / 2.0 / x) * Bessel_Asymptotic(n + 0.5, x).SecondSolutionValue);
                } else {
                    // Recurse upward from y0 and y1.
                    // This is okay, even in transition region, because y grows with order.
                    SphericalBesselY_ZeroAndOne(x, out double ym1, out double y);
                    for (int k = 1; k < n; k++) {
                        double yp1 = (2 * k + 1) / x * y - ym1;
                        ym1 = y;
                        y = yp1;
                        // We have to break immediately when we hit infinity because otherwise the next iteration
                        // will produce infinity - infinity = NaN.
                        if (Double.IsInfinity(y)) break;
                    }
                    return y;
                }
            }

        }

        private static double SphericalBesselJ_Zero (double x) {
            return MoreMath.Sinc(x);
        }

        private static double SphericalBesselJ_SeriesOne (double x) {
            double xx = -0.5 * x * x;
            double dj = x / 3.0;
            double j = dj;
            for (int i = 1; i < Global.SeriesMax; i++) {
                double j_old = j;
                dj *= xx / (i * (2 * i + 3));
                j += dj;
                if (j == j_old) return j;
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselJ_One (double x) {
            if (Math.Abs(x) < 0.5) {
                // Close to the origin, sinc and cosine terms cancel, so use the series.
                return SphericalBesselJ_SeriesOne(x);
            } else if (Math.Abs(x) < 100.0) {
                // There is no need to use sinc because we know we're not using at x=0.
                return (MoreMath.Sin(x) / x - MoreMath.Cos(x)) / x;
            } else {
               // Do I need to fall back to this? Isn't intermediate expression good enough since using MoreMath?
                return Math.Sqrt(Math.PI / 2.0 / x) * Bessel_Asymptotic(1.5, x).FirstSolutionValue;
            }
        }

        private static double SphericalBesselJ_Series (int n, double x) {
            double xx = -0.5 * x * x;
            double df = PowerOverDoubleFactorial(x, n);
            //double df = MoreMath.Pow(x, n) / AdvancedIntegerMath.DoubleFactorial(2 * n + 1);
            //double df = Math.Exp(n * Math.Log(x) - AdvancedIntegerMath.LogDoubleFactorial(2 * n + 1));
            double f = df;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                df *= xx / (k * (2 * (n + k) + 1));
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselY_Series (int n, double x) {
            double xx = -0.5 * x * x;
            double df = -1.0 / PowerOverDoubleFactorial(x, n - 1) / (x * x);
            double f = df;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                df *= xx / (k * (2 * (k - n) - 1));
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        private static double SphericalBesselY_Zero (double x) {
            return -MoreMath.Cos(x) / x;
        }

        private static void SphericalBesselY_ZeroAndOne (double x, out double y0, out double y1) {
            y0 = -MoreMath.Cos(x) / x;
            y1 = (y0 - MoreMath.Sin(x)) / x;
        }

        // Miller's method assumes a value at some high N and recurrs downward
        // the result is then normalized using a sum relation or a known value

        private static double SphericalBesselJ_Miller (int n, double x) {

            int m = Bessel_Miller_Limit(n, x);
            Debug.Assert(m > n);

            double jp1 = 0.0;
            double j = 1.0E-150;

            // recur downward to order zero
            // the recurrence j_{k-1} = (2k+1)/x * j_k - j_{k+1} is stable in this direction
            for (int k = m; k > n; k--) {
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

}
