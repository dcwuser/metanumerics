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
        /// <para>Li<sub>2</sub>(x) is real for -&#x221E; &lt; x &#x2264; 1; for values outside this range,
        /// use the complex verion <see cref="AdvancedComplexMath.DiLog"/>.</para>
        /// </remarks>
        /// <seealso cref="AdvancedComplexMath.DiLog"/>
        public static double DiLog (double x) {
            if (x > 1.0) {
                throw new ArgumentOutOfRangeException("x");
            } else if (x > 0.7) {
                // use series near 1
                return (DiLog_Series_1(1.0 - x));
            } else if (x > -0.7) {
                // series near 0 (defining power series)
                return (DiLog_Series_0(x));
            } else if (x >= -1.0) {
                // use Li(-x) = 1/2 Li(x^2) - Li(-x) to map to [0,1]
                return (DiLog(x * x) / 2.0 - DiLog(-x));
            } else {
                // use formula for Li(1/x) to map to [-1,0]
                double ln = Math.Log(-x);
                return (-Math.PI * Math.PI / 6.0 - ln * ln / 2.0 - DiLog(1.0/x));
            }
        }

        // series definition of DiLog; this converges relably to full accuracy within a few tens of iterations below x~1/2; by x~1 it no longer converges

        private static double DiLog_Series_0 (double x) {
            double xx = x;
            double f = xx;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double f_old = f;
                xx *= x;
                f += xx / (k * k);
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

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

    }


    public static partial class AdvancedComplexMath {

        /// <summary>
        /// Computes the dilogarathm function, also called Spence's function.
        /// </summary>
        /// <param name="z">The complex argument.</param>
        /// <returns>The value Li<sub>2</sub>(z).</returns>
        /// <seealso cref="AdvancedMath.DiLog"/>
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
            for (int k = 1; k < DC.Length; k++) {
                Complex f_old = f;
                p *= ln2 / (2 * k + 1) / (2 * k);
                f += (-DC[k] / (2 * k)) * p;
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();

        }

        // Bernoulli numbers
        // these are coefficients in the log expansion for the dilog log series

        private static readonly double[] DC = new double[] {
            1.0, 1.0 / 6.0, -1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0, 5.0 / 66.0, - 691.0 / 2730.0, 7.0 / 6.0,
            -3617.0 / 510.0, 43867.0 / 798.0, -174611.0 / 330.0, 854513.0 / 138.0, -236364091.0 / 2730.0,
            -236364091.0 / 2730.0, 8553103.0 / 6.0, -23749461029.0 / 870.0, 8615841276005.0 / 14322.0
        };

    }

}