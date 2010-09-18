using System;


namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        private static double Maximum (params double[] values) {
            double max = 0.0;
            foreach (double value in values) {
                double abs = Math.Abs(value);
                if (abs > max) max = abs;
            }
            return (max);
        }

        /// <summary>
        /// Computes the Carslon elliptic integral R<sub>F</sub>.
        /// </summary>
        /// <param name="x">The first parameter, which must be non-negative.</param>
        /// <param name="y">The second parameter, which must be non-negative.</param>
        /// <param name="z">The third parameter, which must be non-negative.</param>
        /// <returns>The value of R<sub>F</sub>(x,y,z).</returns>
        /// <remarks>
        /// <para>As can be seen from the integral, all three parameters are equivilent, so R<sub>F</sub> is symmetric with respect to any
        /// permutation of the parameters.</para>
        /// </remarks>
        /// <seealso cref="EllipticK"/>
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
                double e = Maximum(dx, dy, dz);

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
                if (f == f_old) return (Math.PI / 2.0 * f);
            }

            throw new NonconvergenceException();

        }

        private static readonly double SqrtAccuracy = Math.Sqrt(Global.Accuracy);

        private static double Elliptic_AGM (double k) {

            double a = 1.0 - k;
            double b = 1.0 + k;

            for (int n=0; n < Global.SeriesMax; n++) {

                double am = (a + b) / 2.0;

                if (Math.Abs(a-b) < SqrtAccuracy) {
                    return (Math.PI / 2.0 / am);
                }

                double gm = Math.Sqrt(a * b);

                a = am;
                b = gm;

            }

            throw new NonconvergenceException();

        }

        /// <summary>
        /// Computes the complete elliptic integral of the first kind.
        /// </summary>
        /// <param name="k">The elliptic modulus, which must lie between zero and one.</param>
        /// <returns>The value of the Legendre integral K(k).</returns>
        /// <remarks>
        /// <para>K(k) is defined via an elliptic integral.</para>
        /// <img src="../images/EllipticKIntegral.png" />
        /// <para>It appears in the Legendre reduction of integrals of rational funtions.</para>
        /// <para>Be aware that some authors use the the parameter m = k<sup>2</sup> instead of k.</para>
        /// </remarks>
        public static double EllipticK (double k) {
            if ((k < 0) || (k > 1.0)) throw new ArgumentOutOfRangeException("k");
            if (k < 0.25) {
                return (EllipticK_Series(k));
                //return (CarlsonF(0.0, 1.0 - k * k, 1.0));
            } else {
                return (Elliptic_AGM(k));
                //return (CarlsonF(0.0, 1.0 - k * k, 1.0));
            }
        }

        /// <summary>
        /// Computes the incomplete elliptic integral of the first kind.
        /// </summary>
        /// <param name="phi">The integration angle.</param>
        /// <param name="k">The elliptic modulus, which must lie between zero and one.</param>
        /// <returns>The value of the Legendre integral K(k).</returns>
        public static double EllipticF (double phi, double k) {
            if ((k < 0) || (k > 1.0)) throw new ArgumentOutOfRangeException("k");
            double s = AdvancedMath.Sin(phi, 0.0);
            double c = AdvancedMath.Cos(phi, 0.0);
            double z = s * k;
            return (s * CarlsonF(c * c, 1.0 - z * z, 1.0));
        }

    }

}