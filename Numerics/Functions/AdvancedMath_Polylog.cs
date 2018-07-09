using System;
using System.Diagnostics;
using Meta.Numerics;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the dilogarathm function, also called Spence's function.
        /// </summary>
        /// <param name="x">The argument, which must be less than or equal to unity.</param>
        /// <returns>The value Li<sub>2</sub>(x).</returns>
        /// <remarks>
        /// <para>The dilogarithm can be defined by an infinite sum.</para>
        /// <img src="../images/DilogSum.png" />
        /// <para>The function gets is name from the similarity of this series to the expansion of ln(1-x), the
        /// difference being that the integer in the denominator is raised to the second power.</para>
        /// <para>Li<sub>2</sub>(x) is real for -&#x221E; &lt; x &#x2264; 1; for values outside this range,
        /// use the complex version <see cref="AdvancedComplexMath.DiLog"/>.</para>
        /// </remarks>
        /// <seealso cref="AdvancedComplexMath.DiLog" />
        /// <seealso href="http://en.wikipedia.org/wiki/Dilogarithm" />
        public static double DiLog (double x) {
            if (x > 1.0) {
                throw new ArgumentOutOfRangeException(nameof(x));
            } else if (x > 0.625) {
                // use Li(x) + Li(1-x) = \frac{\pi^2}{6} - \log x \log (1-x)
                // to map x near 1 to x near 0
                if (x == 1.0) {
                    // Special case x = 1 exactly to avoid computing log(0).
                    return (Math.PI * Math.PI / 6.0);
                } else {
                    return (Math.PI * Math.PI / 6.0 - DiLog_Series(1.0 - x) - Math.Log(x) * MoreMath.LogOnePlus(-x));
                }
                // use series near 1
                //return (DiLog_Series_1(1.0 - x));
            } else if (x > -0.625) {
                // series near 0 (defining power series)
                return (DiLog_Series(x));
            } else if (x >= -1.0) {
                // Use Li(x) + Li(-x) = \frac{1}{2} Li(x^2)
                // to map negative x to positive x
                return (0.5 * DiLog(x * x) - DiLog(-x));
            } else if (x >= Double.NegativeInfinity) {
                // use formula for Li(1/x) to map to [-1,0]
                return (-Math.PI * Math.PI / 6.0 - 0.5 * MoreMath.Sqr(Math.Log(-x)) - DiLog(1.0 / x));
            } else {
                return (Double.NaN);
            }
        }

        // Series definition of DiLog
        //   Li_2(x) = \sum_{k=1}^{\infty} \frac{x^k}{k^2}
        // This converges relably to full accuracy within a few tens of iterations below x~1/2; by x~1 it no longer converges.

        private static double DiLog_Series (double x) {
            double xk = x; // tracks x^k
            double f = xk;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double f_old = f;
                xk *= x;
                f += xk / (k * k);
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        /*
        private static double DiLog_Series_1 (double e) {

            double f = Math.PI * Math.PI / 6.0;
            if (e == 0.0) return (f);

            double L = Math.Log(e);
            double ek = 1.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                ek *= e;
                double df = ek * (L - 1.0 / k) / k;
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }
        */

        /// <summary>
        /// Computes the polylogarithm function.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="x">The argument, which must be less than or equal to one.</param>
        /// <returns>The value of Li<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The nth polylog of x is defined via the series:</para>
        /// <img src="..\Images\PolyLogSeries.png" />
        /// <para>Its name comes from the fact that this is a generalization of the logarithm
        /// series. For n = 1 it reduces to -log(1-x). For n = 2 it reduces to the <see cref="AdvancedMath.DiLog"/> function.</para>
        /// <para>The polylogarithm function becomes complex for arguments larger than one.</para>
        /// </remarks>
        public static double PolyLog (int n, double x) {

            if (x > 1.0) throw new ArgumentOutOfRangeException(nameof(x));

            if (n == 0) {
                return (x / (1.0 - x));
            } else if (n == 1) {
                return (-MoreMath.LogOnePlus(-x));
            } else {

                if (x < -1.0) {
                    // For x < -1, reflect x -> 1/x
                    double w = Math.Log(-x);
                    double s = PolyLog_BernoulliSum(n, w);
                    double t = PolyLog(n, 1.0 / x);
                    if (n % 2 == 0) {
                        return (s - t);
                    } else {
                        return (s + t);
                    }
                } else if (x < -0.25) {
                    // For -1 <= x <= 0, reflect x -> -x, unless we are close enough to just directly use the series.
                    // The terms have opposite signs, so we should worry about cancelation. But since Li(x) increases
                    // monotonically, Li_n(x^2) will be smaller than Li(x) on 0 < x < 1, and that term is additionally
                    // suppressed by 2^{1-n} for all n > 1.
                    return (-PolyLog(n, -x) + MoreMath.Pow(2.0, 1 - n) * PolyLog(n, x * x));
                } else if (x < 0.25) {
                    // For |x| < 0.25, use the defining series.
                    return(PolyLog_Series(n, x));
                } else {
                    // For 0.25 < x <= 1.0, use log series
                    return(PolyLog_LogSeries(n, x));
                }

            }

        }

        private static double PolyLog_BernoulliSum (int n, double w) {

            double s = 0.0;
            for (int k = n; k > 1; k--) {

                double ds = (1.0 - MoreMath.Pow(2.0, 1 - k)) * Math.Abs(AdvancedIntegerMath.BernoulliNumber(k)) *
                    MoreMath.Pow(Global.TwoPI, k) / AdvancedIntegerMath.Factorial(k) *
                    MoreMath.Pow(w, n - k) / AdvancedIntegerMath.Factorial(n - k);
                s -= ds;
            }
            s -= MoreMath.Pow(w, n) / AdvancedIntegerMath.Factorial(n);
            return (s);

        }

        // The Taylor series for Li_s(z) around the origin, often taken as the definition is
        //   Li_s(z) = \sum_{k=1}^{\infty} \frac{z^k}{k^s} = z + \frac{z^2}{2^s} + \frac{z^3}{3^s} + \cdots
        // While technically convergent for all |z| < 1, convergence becomes slower further from the origin.

        private static double PolyLog_Series (int n, double x) {
            Debug.Assert(Math.Abs(x) < 1.0);
            double xk = x;
            double f = xk;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double f_old = f;
                xk *= x;
                f += xk / MoreMath.Pow(k, n);
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        // Li_n(x) = \sum_{k=0}^{\infty} \zeta(n-k) \frac{\log^k x}{k!} + \frac{\log^{n-1} x}{(n-1)!} \left( H_{n-1} -\log(-\log x) \right)

        private static double PolyLog_LogSeries (int n, double x) {
            
            double lnx = Math.Log(x);

            double f = AdvancedMath.RiemannZeta(n);

            if (lnx == 0.0) return (f);

            // c stores [log(x)]^k / k!
            double c = 1.0;

            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                c *= lnx / k;
                // argument of zeta
                int m = n - k;
                if (m < 0) {
                    // For negative arguments, use \zeta(-m) = \frac{B_{m+1}}{m+1} and that odd Bernoulli numbers vanish.
                    if (-m % 2 == 0) continue;
                    // This could theoretically overrun our stored Bernoulli values, but if we haven't converged after 32 negative terms, we are in trouble.
                    //f += c * AdvancedIntegerMath.Bernoulli[(-m + 1) / 2] / (-m + 1);
                    f += c * AdvancedMath.RiemannZeta(m);
                } else if (m == 1) {
                    // Special term in place of \zeta(1).
                    f += c * (AdvancedIntegerMath.HarmonicNumber(n - 1) - Math.Log(-lnx));
                } else {
                    // Otherwise just compute \zeta(m).
                    // We could reduce even m to Bernoulli references but then we would be in trouble for n > 32.
                    f += c * AdvancedMath.RiemannZeta(m);
                }
                if (f == f_old) return(f);
            }
            throw new NonconvergenceException();

        }

        /// <summary>
        /// Computes the Clausen integral.
        /// </summary>
        /// <param name="t">The argument.</param>
        /// <returns>The value of Cl<sub>2</sub>(t).</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Clausen%27s_function"/>
        public static double Clausen (double t) {

            // reduce t to [0,2\pi]
            if ((t < 0.0) || (t > Global.TwoPI)) {
                double z = t / Global.TwoPI;
                z = z - Math.Floor(z);
                t = z * Global.TwoPI;
            }

            // Pick either the expansion around 0 or the expansion around \pi.
            // The expansion around 0 converges faster, so we use it over a larger area.
            if (t < 2.0 * Math.PI / 3.0) {
                return (ClausenNearZero(t));
            } else if (t < 4.0 * Math.PI / 3.0) {
                return(ClausenNearPi(Math.PI - t));
            } else {
                return(ClausenNearZero(t - Global.TwoPI));
            }



        }

        // Abromowitz and Stegun 27.8.2
        //  Cl_2(\theta) = \theta - \theta \log |\theta| + \sum_{k=1}^{\infty} \frac{|B_{2k}| \theta^{2k+1}}{(2k)! (2k) (2k+1)}

        private static double ClausenNearZero (double t) {
            // avoid computation of log(0)
            if (t == 0.0) return (0.0);
            double f = 1.0 - Math.Log(Math.Abs(t));

            double t2 = t * t;
            double sk = 1.0; // tracks t^2k / k!
            for (int k = 1; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                double f_old = f;
                int tk = 2 * k;
                sk *= t2 / tk / (tk - 1);
                f += Math.Abs(AdvancedIntegerMath.Bernoulli[k]) * sk / tk / (tk + 1);
                if (f == f_old) return (t * f);
            }
            throw new NonconvergenceException();
        }

        // Abromowitz & Stegun 7.8.3
        // Cl_2(\pi - \theta) = \theta \log 2 - \sum_{k=1}^{\infty} \frac{|B_{2k}| (2^2k - 1) \theta^{2k+1}}{(2k!) (2k) (2k+1)}

        // Because of the additional factor 2^{2k} in numerator, this series converges much slower than the series around 0.
        // If we transition between the two of them at \pi/2, this series takes about 30 terms while the other takes only 10.
        // If we transition between the two of them at 2\pi/3, each takes about 15 terms.

        private static double ClausenNearPi (double t) {

            double f = Global.LogTwo;

            double t2 = t * t;
            double sk = 1.0; // tracks t^2k / k!
            for (int k = 1; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                double f_old = f;
                int tk = 2 * k;
                sk *= t2 / tk / (tk - 1);
                f -= Math.Abs(AdvancedIntegerMath.Bernoulli[k]) * (MoreMath.Pow(4.0, k) - 1.0) * sk / tk / (tk + 1);
                if (f == f_old) return (t * f);
            }
            throw new NonconvergenceException();
        }

    }


    public static partial class AdvancedComplexMath {

        /// <summary>
        /// Computes the complex dilogarithm function, also called Spence's function.
        /// </summary>
        /// <param name="z">The complex argument.</param>
        /// <returns>The value Li<sub>2</sub>(z).</returns>
        /// <remarks>
        /// <para>This function is the analytic continuation of the dilogarithm function (<see cref="AdvancedMath.DiLog"/>) into the complex plane.</para>
        /// <para>The image below shows the complex dilogarithm function near the origin, using domain coloring.</para>
        /// <img src="../images/ComplexDiLogPlot.png" />
        /// </remarks>
        /// <seealso cref="AdvancedMath.DiLog"/>
        /// <seealso href="http://mathworld.wolfram.com/Dilogarithm.html" />
        public static Complex DiLog (Complex z) {

            Complex f;
            double a0 = ComplexMath.Abs(z);
            if (a0 > 1.0) {
                // outside the unit disk, reflect into the unit disk
                Complex ln = ComplexMath.Log(-z);
                f = -Math.PI * Math.PI / 6.0 - ln * ln / 2.0 - DiLog(1.0 / z);
            } else {
                // inside the unit disk...
                if (a0 < 0.75) {
                    // close to 0, use the expansion about zero
                    f = DiLog_Series_0(z);
                } else {
                    // we are in the annulus near the edge of the unit disk
                    if (z.Re < 0.0) {
                        // reflect negative into positive half-disk
                        // this avoids problems with the log expansion near -1
                        f = DiLog(z * z) / 2.0 - DiLog(-z);
                    } else {
                        // figure out whether we are close to 1
                        Complex e = 1.0 - z;
                        if (ComplexMath.Abs(e) < 0.5) {
                            // close to 1, use the expansion about 1
                            f = DiLog_Series_1(e);
                        } else {
                            // otherwise, use the log expansion, which is good
                            // near the unit circle but not too close to 1 or -1
                            f = DiLog_Log_Series(z);
                        }
                    }
                }
            }

            if ((z.Re > 1.0) && (Math.Sign(f.Im) != Math.Sign(z.Im))) f = f.Conjugate;

            return (f);

        }

        private static Complex DiLog_Series_0 (Complex z) {
            Complex zz = z;
            Complex f = zz;
            for (int k = 2; k < Global.SeriesMax; k++) {
                Complex f_old = f;
                zz *= z;
                f += zz / (k * k);
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        private static Complex DiLog_Series_1 (Complex e) {

            Complex f = Math.PI * Math.PI / 6.0;
            if (e == 0.0) return (f);

            Complex L = ComplexMath.Log(e);
            Complex ek = 1.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                Complex f_old = f;
                ek *= e;
                Complex df = ek * (L - 1.0 / k) / k;
                f += df;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        private static Complex DiLog_Log_Series (Complex z) {

            Complex ln = ComplexMath.Log(z);
            Complex ln2 = ln * ln;

            Complex f = Math.PI * Math.PI / 6.0 + ln * (1.0 - ComplexMath.Log(-ln)) - ln2 / 4.0;

            Complex p = ln;
            for (int k = 1; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                Complex f_old = f;
                p *= ln2 / (2 * k + 1) / (2 * k);
                f += (-AdvancedIntegerMath.Bernoulli[k] / (2 * k)) * p;
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();

        }

    }

}