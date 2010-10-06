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

            // check that at most one is zero

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
            } else {
                if (k == 1.0) return (Double.PositiveInfinity);
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

    }

}