using System;


namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        /*
        private static double Maximum (params double[] values) {
            double max = 0.0;
            foreach (double value in values) {
                double abs = Math.Abs(value);
                if (abs > max) max = abs;
            }
            return (max);
        }
        */

        /// <summary>
        /// Computes the Carslon elliptic integral R<sub>F</sub>.
        /// </summary>
        /// <param name="x">The first parameter, which must be non-negative.</param>
        /// <param name="y">The second parameter, which must be non-negative.</param>
        /// <param name="z">The third parameter, which must be non-negative.</param>
        /// <returns>The value of R<sub>F</sub>(x,y,z).</returns>
        /// <remarks>
        /// <para>The Carlson F integral is:</para>
        /// <img src="../images/CarlsonFIntegral.png" />
        /// <para>As can be seen from the integral, all three parameters are equivilent, so R<sub>F</sub> is symmetric with respect to any
        /// permutation of the parameters.</para>
        /// <para>The Carlson integrals can be used to express integrals of rational functions. In that sense, they are replacements for
        /// the Legendre elliptic functions.</para>
        /// </remarks>
        /// <seealso cref="EllipticK"/>
        /// <seealso cref="EllipticF"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Carlson_symmetric_form"/>
        public static double CarlsonF (double x, double y, double z) {

            if (x < 0.0) throw new ArgumentOutOfRangeException("x");
            if (y < 0.0) throw new ArgumentOutOfRangeException("y");
            if (z < 0.0) throw new ArgumentOutOfRangeException("z");

            // if more than one is zero, the result diverges
            if (((x == 0.0) && ((y == 0.0) || (z == 0.0))) || ((y == 0.0) && (z == 0.0))) return (Double.NaN);

            for (int n = 0; n < Global.SeriesMax; n++) {

                // find out how close we are to the expansion point
                double m = (x + y + z) / 3.0;
                double dx = (x - m) / m;
                double dy = (y - m) / m;
                double dz = (z - m) / m;
                double e = Math.Max(Math.Abs(dx), Math.Max(Math.Abs(dy), Math.Abs(dz)));
                //double e = Maximum(dx, dy, dz);

                // if we are close enough, use the seventh order expansion

                if (e < 0.02) {
 
                    double E2 = dx * dy + dx * dz + dy * dz;
                    double E3 = dx * dy * dz;

                    //double F = 1.0 - E2 / 10.0 - E3 / 14.0 + E2 * E2 / 24.0 + 3.0 * E2 * E3 / 44.0;
                    double F = 1.0 - E2 / 10.0 - E3 / 14.0 + E2 * E2 / 24.0 + 3.0 * E2 * E3 / 44.0 - 5.0 / 208.0 * E2 * E2 * E2 + 3.0 / 104.0 * E3 * E3 - E2 * E2 * E3 / 16.0;

                    return (F / Math.Sqrt(m));

                }

                // if we are not close enough, use the duplication theory to move us closer

                double lambda = Math.Sqrt(x * y) + Math.Sqrt(x * z) + Math.Sqrt(y * z);

                x = (x + lambda) / 4.0;
                y = (y + lambda) / 4.0;
                z = (z + lambda) / 4.0;

            }

            throw new NonconvergenceException();

        }


        public static double CarlsonD (double x, double y, double z) {

            double s = 1.0;
            double t = 0.0;

            for (int n = 0; n < Global.SeriesMax; n++) {

                // find out how close we are to the expansion point
                double m = (x + y + 3.0 * z) / 5.0;
                double dx = (x - m) / m;
                double dy = (y - m) / m;
                double dz = (z - m) / m;
                double e = Math.Max(Math.Abs(dx), Math.Max(Math.Abs(dy), Math.Abs(dz)));

                if (e < 0.01) {
                }

                // we are not close enough; use the duplication theory to move us closer

                double lambda = Math.Sqrt(x * y) + Math.Sqrt(x * z) + Math.Sqrt(y * z);

                s *= 4.0;
                t += 3.0 / Math.Sqrt(z) / (z + lambda) / s;

                x = (x + lambda) / 4.0;
                y = (y + lambda) / 4.0;
                z = (z + lambda) / 4.0;

            }

            throw new NonconvergenceException();

        }

        private static double EllipticK_Series (double k) {

            // (Pi/2) ( 1 + (1/2 k)^2 + (1/2 k 3/4 k)^2 + (1/2 k 3/4 k 5/6 k)^2 + ... ) 

            double z = 1.0;
            double f = 1.0;
            for (int n = 1; n < Global.SeriesMax; n++) {
                double f_old = f;
                z = z * (2 * n - 1) / (2 * n) * k;
                f += z * z;
                if (f == f_old) return (Global.HalfPI * f);
            }

            throw new NonconvergenceException();

        }

        // Asymptotic expansion of Elliptic integral of first kind
        // K = \sum_n [ (1/2)_n k' / n! ]^2 [ ln(1/k') + Psi(k+1) - Psi(k+1/2) ]
        // This is good close to k=1 / k'=0 since it is a power series in k'
        // 8 terms at k' = 0.1, 12 terms at k' = 0.25, 23 terms at k' = 0.5, 53 terms at k' = 0.75
        // accurate even at k' = 0.75; because all terms are the same sign, no cancelations occur even when terms to not decrease rapidly

        private static double EllipticK_Asymptotic (double k1) {
            double p = 1.0;
            double q = Math.Log(1.0 / k1) + 2.0 * Global.LogTwo;
            double f = q;
            for (int m = 1; m < Global.SeriesMax; m++) {
                double f_old = f;
                p *= k1 / m * (m - 0.5);
                q -= 1.0 / m / (2 * m - 1);
                double df = p * p * q;
                f += df;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        // compute K via the arithmetic-geometric-mean algorithm

        private static double Elliptic_AGM (double k) {

            double a = 1.0 - k;
            double b = 1.0 + k;

            for (int n=0; n < Global.SeriesMax; n++) {

                double am = (a + b) / 2.0;

                // if a and b are seperated by a small e, the update only changes a and b by e^2
                // therefore we can stop as soon as a and b are within the square root of machine precision
                // in fact, we must not use (a == b) as a termination criterion, because we can get into a loop
                // where a and b dance around and never become precisely equal
                if (Math.Abs(a-b) < SqrtAccuracy) {
                    return (Global.HalfPI / am);
                }

                double gm = Math.Sqrt(a * b);

                a = am;
                b = gm;

            }

            throw new NonconvergenceException();

        }

        private static readonly double SqrtAccuracy = Math.Sqrt(Global.Accuracy);

        /// <summary>
        /// Computes the complete elliptic integral of the first kind.
        /// </summary>
        /// <param name="k">The elliptic modulus, which must lie between zero and one.</param>
        /// <returns>The value of the Legendre integral K(k).</returns>
        /// <remarks>
        /// <para>K(k) is defined as the complete elliptic integral:</para>
        /// <img src="../images/EllipticKIntegral.png" />
        /// <para>It appears in the Legendre reduction of integrals of rational funtions.</para>
        /// <para>Be aware that some authors use the the parameter m = k<sup>2</sup> instead of the modulus k.</para>
        /// </remarks>
        /// <seealso cref="EllipticF"/>
        public static double EllipticK (double k) {
            if ((k < 0) || (k > 1.0)) throw new ArgumentOutOfRangeException("k");
            if (k < 0.25) {
                return (EllipticK_Series(k));
            } if (k > 0.875) {
                double k1 = Math.Sqrt(1.0 - k * k);
                // k'=0.484 at k=0.875
                if (k1 == 0.0) return (Double.PositiveInfinity);
                return (EllipticK_Asymptotic(k1));
            } else {
                return (Elliptic_AGM(k));
            }
        }

        /// <summary>
        /// Computes the incomplete elliptic integral of the first kind.
        /// </summary>
        /// <param name="phi">The integration angle (in radians).</param>
        /// <param name="k">The elliptic modulus, which must lie between zero and one.</param>
        /// <returns>The value of the Legendre integral F(k).</returns>
        /// <remarks>
        /// <para>Legendre's first incomplete elliptic integral is:</para>
        /// <img src="../images/EllipticFIntegral.png" />
        /// <para>When the integral angle is &#x3C0;, this function reduces to the complete elliptic integral of the
        /// first kind (<see cref="EllipticK"/>.</para>
        /// <para>Be aware that some authors use the the parameter m = k<sup>2</sup> instead of the modulus k.</para>
        /// </remarks>
        /// <seealso cref="EllipticK"/>
        public static double EllipticF (double phi, double k) {
            if (Math.Abs(phi) > Global.HalfPI) throw new ArgumentOutOfRangeException("phi");
            if ((k < 0) || (k > 1.0)) throw new ArgumentOutOfRangeException("k");
            double s = AdvancedMath.Sin(phi, 0.0);
            double c = AdvancedMath.Cos(phi, 0.0);
            double z = s * k;
            return (s * CarlsonF(c * c, 1.0 - z * z, 1.0));
        }

        // Series for complete elliptic integral of the second kind
        // 7 terms for k=0.1, 12 for k=0.25, 22 for k=0.5, 49 for k=0.75; accurate for all these cases

        private static double EllipticE_Series (double k) {

            // (Pi/2) [ 1 - (1/2 k)^2 / 1 - (1/2 k 3/4 k)^2 / 3 - (1/2 k 3/4 k 5/6 k)^2 / 5 + ... ]

            double z = 1.0;
            double f = 1.0;
            for (int n = 1; n < Global.SeriesMax; n++) {
                double f_old = f;
                z = z * (n - 0.5) / n * k;
                f -= z * z / ( 2 * n - 1);
                if (f == f_old) return (Global.HalfPI * f);
            }

            throw new NonconvergenceException();

        }

        // 7 terms at k'=0.1, 12 at k'=0.25, 23 at k'=0.5, 52 at k'=0.75

        private static double EllipticE_Asymptotic (double k1) {

            double k12 = k1 * k1;
            double p = k12 / 2.0;
            double q = Math.Log(1.0 / k1) + 2.0 * Global.LogTwo - 0.5;
            double f = 1.0 + p * q;
            for (int m = 1; m < Global.SeriesMax; m++) {
                double f_old = f;
                p *= (m - 0.5) * (m + 0.5) / m / (m + 1) * k12;
                q -= 1.0 / (2 * m - 1) / (2 * m) + 1.0 / (2 * m + 1) / (2 * m + 2);
                f += p * q;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();

        }

        /// <summary>
        /// Computes the complete elliptic integral of the second kind.
        /// </summary>
        /// <param name="k">The elliptic modulus, which must lie between zero and one.</param>
        /// <returns>The value of the Legendre integral E(k).</returns>
        /// <remarks>
        /// <para>E(k) is defined as the complete elliptic integral:</para>
        /// <img src="../images/EllipticEIntegral.png" />
        /// <para>It appears in the Legendre reduction of integrals of rational funtions.</para>
        /// <para>The perimeter of an ellipse with major axis a and eccentricity e is 4 a E(e).</para>
        /// <para>Be aware that some authors use the the parameter m = k<sup>2</sup> instead of the modulus k.</para>
        /// </remarks>
        public static double EllipticE (double k) {
            if ((k < 0.0) || (k > 1.0)) throw new ArgumentOutOfRangeException("k");
            // these expansions are accurate in the intermediate region, but require many terms
            // it would be good to use a faster approach there, like we do for K
            if (k < 0.71) {
                return (EllipticE_Series(k));
            } else {
                double k1 = Math.Sqrt(1.0 - k * k);
                if (k1 == 0.0) return (1.0);
                return (EllipticE_Asymptotic(k1));
            }
        }

    }

}