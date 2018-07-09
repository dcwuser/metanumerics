using System;
using System.Diagnostics;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the Carslon elliptic integral R<sub>F</sub>.
        /// </summary>
        /// <param name="x">The first parameter, which must be non-negative.</param>
        /// <param name="y">The second parameter, which must be non-negative.</param>
        /// <param name="z">The third parameter, which must be non-negative.</param>
        /// <returns>The value of R<sub>F</sub>(x,y,z).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/>, <paramref name="y"/>, or
        /// <paramref name="z"/> is negative.</exception>
        /// <remarks>
        /// <para>The Carlson F integral is:</para>
        /// <img src="../images/CarlsonFIntegral.png" />
        /// <para>As can be seen from the definition, R<sub>F</sub> is symmetric with respect to permutations among its parameters.</para>
        /// <para>The Carlson integrals can be used to express integrals of rational functions. In that sense, they are replacements for
        /// the Legendre elliptic functions.</para>
        /// </remarks>
        /// <seealso cref="EllipticK"/>
        /// <seealso cref="EllipticF"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Carlson_symmetric_form"/>
        public static double CarlsonF (double x, double y, double z) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if (y < 0.0) throw new ArgumentOutOfRangeException(nameof(y));
            if (z < 0.0) throw new ArgumentOutOfRangeException(nameof(z));

            // if more than one is zero, the result diverges
            if (((x == 0.0) && ((y == 0.0) || (z == 0.0))) || ((y == 0.0) && (z == 0.0))) return (Double.PositiveInfinity);

            // To compute CarlsonF, combine the shift identity
            //   R_F(x, y, z) = 2 R_F(x+\lambda, y+\lambda, z+\lambda)
            // where \lambda = \sqrt{xy} + \sqrt{yz} + \sqrt{xz}, and the homogeneity identity
            //   R_F(a x, a y, a z) = R_F(x, y, z) / \sqrt{a}
            // to obtain
            //   R_F(x, y, z) = R_F(\frac{x+\lambda}{4}, \frac{y+\lambda}{4}, \frac{z+\lambda}{4})
            // Applying this repeatedly moves the three arguments closer together. When they are all exactly equal (DLMF 19.20.1)
            //   R_F(x, x, x) = 1/\sqrt{x}
            // and when they are nearly equal, we can use an expansion around this point

            for (int n = 0; n < Global.SeriesMax; n++) {

                // find out how close we are to the expansion point
                double m = (x + y + z) / 3.0;
                double dx = (x - m) / m;
                double dy = (y - m) / m;
                double dz = (z - m) / m;
                double e = Math.Max(Math.Abs(dx), Math.Max(Math.Abs(dy), Math.Abs(dz)));

                // If we are close enough, use the seventh order expansion (DLMF 19.36.1)
                // Since the error is order e^8, this is good enough for e ~ (10^-16)^(1/8) = 0.01

                if (e < 0.01) {
 
                    double E2 = dx * dy + dx * dz + dy * dz;
                    double E3 = dx * dy * dz;

                    double F = 1.0 - E2 / 10.0 - E3 / 14.0 + E2 * E2 / 24.0 + 3.0 * E2 * E3 / 44.0 - 5.0 / 208.0 * E2 * E2 * E2 + 3.0 / 104.0 * E3 * E3 - E2 * E2 * E3 / 16.0;

                    // Rather unsymmetrically, Carlson defines dz = dx + dy which gives it the opposite sign, which gives E3 the opposite sign, which gives terms with odd powers
                    // of E3 the opposite sign. So don't get confused when you look at Carlson's DLMF article and 1994 paper.

                    return (F / Math.Sqrt(m));

                }

                // if we are not close enough, use the duplication theorem (DLMF 19.26.18) to move us closer

                double lambda = Math.Sqrt(x * y) + Math.Sqrt(x * z) + Math.Sqrt(y * z);

                x = (x + lambda) / 4.0;
                y = (y + lambda) / 4.0;
                z = (z + lambda) / 4.0;

            }

            throw new NonconvergenceException();

        }


        /// <summary>
        /// Computes the Carlson integral R<sub>D</sub>.
        /// </summary>
        /// <param name="x">The first parameter, which must be non-negative.</param>
        /// <param name="y">The second parameter, which must be non-negative.</param>
        /// <param name="z">The third parameter, which must be non-negative.</param>
        /// <returns>The value of R<sub>D</sub>(x, y, z)</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/>, <paramref name="y"/>, or
        /// <paramref name="z"/> is negative.</exception>
        /// <remarks>
        /// <para>The Carlson D integral is:</para>
        /// <img src="../images/CarlsonDIntegral.png" />
        /// <para>It is symmetric with respect to the interchange of the first two parameters, but not the third parameter.</para>
        /// <para>The Carlson integrals can be used to express integrals of rational functions. In that sense, they are replacements for
        /// the Legendre elliptic functions.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Carlson_symmetric_form"/>
        public static double CarlsonD (double x, double y, double z) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if (y < 0.0) throw new ArgumentOutOfRangeException(nameof(y));
            if (z < 0.0) throw new ArgumentOutOfRangeException(nameof(z));

            if ((x == 0.0) && (y == 0.0)) return (Double.PositiveInfinity);

            // To compute CarlsonD, we use a similar technique as we did for CarlsonF: use a shift identity
            // to move the arguments closer together, then expand about the equal-arguments point

            // The shift and homogeneity relationships are
            //   R_D(x, y, z) = 2 R_D(x+\lambda, y+\lambda, z+\lambda) + \frac{3}{\sqrt{z}(z+\lambda)}
            //   R_D(a x, a y, a z) = a^{-3/2} R_D(x, y, z)
            // which combine to give
            //   R_D(x, y, z) = R_D(\frac{x+\lambda}{c}, \frac{y+\lambda}{c}, \frac{z+\lambda}{c}) + \frac{3}{\sqrt{z}(z+\lambda)}
            // where c is the cube root of 4.

            // This is not quite as nice as the CarlsonF case, first because of the second term that we must track and sum,
            // and second because c~1.4 is not as large a factor as 4, so the arguments tend to increase. (The second problem
            // could be contained at the expense of keeping track of a factor, but don't do that here.)

            // variable to hold the sum of the second terms
            double t = 0.0;

            for (int n = 0; n < Global.SeriesMax; n++) {

                // find out how close we are to the expansion point
                double m = (x + y + 3.0 * z) / 5.0;
                double dx = (x - m) / m;
                double dy = (y - m) / m;
                double dz = (z - m) / m;
                double e = Math.Max(Math.Abs(dx), Math.Max(Math.Abs(dy), Math.Abs(dz)));

                // Our series development (DLMF 19.36.2) goes up to O(e^6). In order that the neglected term e^7 <~ 1.0E-16, we need e <~ 0.005.

                if (e < 0.005) {

                    double xy = dx * dy; double zz = dz * dz;
                    double E2 = xy - 6.0 * zz;
                    double E3 = (3.0 * xy - 8.0 * zz) * dz;
                    double E4 = 3.0 * (xy - zz) * zz;
                    double E5 = xy * zz * dz;

                    double F = 1.0 - 3.0 / 14.0 * E2 - E3 / 6.0 + 9.0 / 88.0 * E2 * E2 - 3.0 / 22.0 * E4 + 9.0 / 52.0 * E2 * E3 - 3.0 / 26.0 * E5
                        - E2 * E2 * E2 / 16.0 + 3.0 / 40.0 * E3 * E3 + 3.0 / 20.0 * E2 * E4;

                    // As in CarlsonF above, Carlson gives dz the opposite sign, which gives E3 and E5 the opposite sign,
                    // which changes the sign of some terms in the series relative to his paper.

                    return (F / Math.Pow(m, 3.0 / 2.0) + t);

                }

                // we are not close enough; use the duplication theory to move us closer

                double lambda = Math.Sqrt(x * y) + Math.Sqrt(x * z) + Math.Sqrt(y * z);

                t += 3.0 / Math.Sqrt(z) / (z + lambda);

                x = (x + lambda) / c4;
                y = (y + lambda) / c4;
                z = (z + lambda) / c4;

            }

            throw new NonconvergenceException();

        }

        // precompute the cube root of four, used by CarlsonD algorithm
        private static readonly double c4 = Math.Pow(2.0, 2.0 / 3.0);

        /// <summary>
        /// Computes the Carlson integral R<sub>G</sub>.
        /// </summary>
        /// <param name="x">The first argument, which must be non-negative.</param>
        /// <param name="y">The second argument, which must be non-negative.</param>
        /// <param name="z">The third argument, which must be non-negative.</param>
        /// <returns>The value of R<sub>G</sub>(x, y, z).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/>, <paramref name="y"/>, or
        /// <paramref name="z"/> is negative.</exception>
        /// <remarks>
        /// <para>The Carlson G integral is:</para>
        /// <img src="..\images\CarlsonGIntegral.png" />
        /// <para>As can be seen from the definition, it is symmetric with respect to interchanges of any of its arguments.</para>
        /// <para>The Carlson integrals can be used to express integrals of rational functions. In that sense, they are replacements for
        /// the Legendre elliptic functions.</para>
        /// </remarks>
        public static double CarlsonG (double x, double y, double z) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if (y < 0.0) throw new ArgumentOutOfRangeException(nameof(y));
            if (z < 0.0) throw new ArgumentOutOfRangeException(nameof(z));

            // We will use
            // 2 R_G = z R_F - 1/3 (x-z)(y-z) R_D + \sqrt{xy/z}

            // To avoid cancelation, we need (x-z)(y-z) < 0. Since R_G is symmetric,
            // we can permute arguments to ensure this. We need z to be the middle-sized
            // argument, e.g. x < z < y.
            if (x > y) Global.Swap(ref x, ref y);
            if (x > z) Global.Swap(ref x, ref z);
            if (z > y) Global.Swap(ref z, ref y);

            double s = (x - z) * (y - z);
            Debug.Assert(s <= 0.0);

            return ((z * CarlsonF(x, y, z) - s / 3.0 * CarlsonD(x, y, z) + Math.Sqrt(x * y / z)) / 2.0);

        }

        // Series for complete elliptic integral of the first kind (A&S 17.3.11):
        //   K(k) = (pi/2) ( 1 + (1/2 k)^2 + (1/2 k 3/4 k)^2 + (1/2 k 3/4 k 5/6 k)^2 + ... )

        private static double EllipticK_Series (double k) {
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

        // Asymptotic expansion of Elliptic integral of first kind (DLMF 19.21.1)
        //   K(k) = \sum_n [ (1/2)_n k' / n! ]^2 [ ln(1/k') + Psi(n+1) - Psi(n+1/2) ]
        // Since this is a power series in k', it is good close to k = 1
        // Converges in 8 terms at k' = 0.1, 12 terms at k' = 0.25, 23 terms at k' = 0.5, 53 terms at k' = 0.75
        // Accurate even at k' = 0.75; because all terms are the same sign, no cancelations occur even when terms do not decrease rapidly

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

        // K(k) can be expressed as an arithmetic-geometric mean (A&S 17.6)
        //   K(k) = = \frac{\pi}{2 AGM(1-k, 1+k)} = \frac{\pi}{2 AGM(1,k')}
        // AGM(a,b) is determined by taking the arithmetic mean (a+b)/2 and the geometric mean sqrt(ab) and re-inserting them as
        // arguments until convergence is achieved.

        private static double EllipticK_AGM (double k) {

            double a = 1.0;
            double b = Math.Sqrt(1.0 - k * k);

            // Starting from 1-k, 1+k, the first iteration will always take us to 1, k',
            // so although it looks prettier (more symmetric) to start with 1-k, 1-k,
            // it's computationally faster to start with 1, k'.

            for (int n=0; n < Global.SeriesMax; n++) {

                double am = (a + b) / 2.0;

                // if a and b are seperated by a small distance e, the update only changes a and b by e^2
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
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> lies outside [0, 1].</exception>
        /// <remarks>
        /// <para>K(k) is defined as the complete elliptic integral:</para>
        /// <img src="../images/EllipticKIntegral.png" />
        /// <para>It appears in the Legendre reduction of integrals of rational funtions.</para>
        /// <para>Be aware that some authors use the the parameter m = k<sup>2</sup> instead of the modulus k.</para>
        /// </remarks>
        /// <seealso cref="EllipticF"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Elliptic_integral"/>
        public static double EllipticK (double k) {
            if ((k < 0) || (k > 1.0)) throw new ArgumentOutOfRangeException(nameof(k));
            if (k < 0.25) {
                // For small k, use the series near k~0.
                return (EllipticK_Series(k));
            } else if (k < 0.875) {
                // For intermediate k, use the AGM method.
                return (EllipticK_AGM(k));
            } else {
                // For large k, use the asymptotic expansion. (k' = 0.484 at k=0.875)
                double k1 = Math.Sqrt(1.0 - k * k);
                if (k1 == 0.0) {
                    return (Double.PositiveInfinity);
                } else {
                    return (EllipticK_Asymptotic(k1));
                }

            }
        }

        /// <summary>
        /// Computes the incomplete elliptic integral of the first kind.
        /// </summary>
        /// <param name="phi">The amplitude (in radians).</param>
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
        /// <seealso href="http://mathworld.wolfram.com/EllipticIntegraloftheFirstKind.html"/>
        public static double EllipticF (double phi, double k) {
            if ((k < 0) || (k > 1.0)) throw new ArgumentOutOfRangeException(nameof(k));
            if (Math.Abs(phi) > Math.PI / 2.0) {
                RangeReduction.ReduceByPiHalves(0.5 * phi, out long m, out double m1);
                return (2 * m * EllipticK(k) + EllipticF(Math.PI * m1, k));
            }
            double s = MoreMath.Sin(phi);
            double c = MoreMath.Cos(phi);
            double z = s * k;
            return (s * CarlsonF(c * c, 1.0 - z * z, 1.0));
        }

        // Series for complete elliptic integral of the second kind
        // Converges in 7 terms for k = 0.1, 12 for k = 0.25, 22 for k = 0.5, 49 for k = 0.75; accurate for all these cases

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
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> lies outside [0, 1].</exception>
        /// <remarks>
        /// <para>E(k) is defined as the complete elliptic integral:</para>
        /// <img src="../images/EllipticEIntegral.png" />
        /// <para>It appears in the Legendre reduction of integrals of rational funtions.</para>
        /// <para>The perimeter of an ellipse with major axis a and eccentricity e is 4 a E(e).</para>
        /// <para>Be aware that some authors use the the parameter m = k<sup>2</sup> instead of the modulus k.</para>
        /// </remarks>
        /// <seealso cref="EllipticE(double,double)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Elliptic_integral"/>
        /// <seealso href="http://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html"/>
        public static double EllipticE (double k) {
            if ((k < 0.0) || (k > 1.0)) throw new ArgumentOutOfRangeException(nameof(k));
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


        /// <summary>
        /// Computes the incomplete elliptic integral of the second kind.
        /// </summary>
        /// <param name="phi">The amplitude (in radians).</param>
        /// <param name="k">The elliptic modulus, which must lie between zero and one.</param>
        /// <returns>The value of E(phi,k).</returns>
        /// <remarks>
        /// <para>The incomplete elliptic integral of the second kind is:</para>
        /// <img src="../images/EllipticEIncompleteIntegral.png" />
        /// <para>It appears in the Legendre reduction of integrals of rational funtions.</para>
        /// <para>Be aware that some authors use the the parameter m = k<sup>2</sup> instead of the modulus k.</para>
        /// </remarks>
        /// <seealso cref="EllipticE(double)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Elliptic_integral"/>
        /// <seealso href="http://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html"/>
        public static double EllipticE (double phi, double k) {

            if ((k < 0.0) || (k > 1.0)) throw new ArgumentOutOfRangeException(nameof(k));

            if (Math.Abs(phi) > Math.PI / 2.0) {
                RangeReduction.ReduceByPiHalves(0.5 * phi, out long m, out double m1);
                return (2 * m * EllipticE(k) + EllipticE(Math.PI * m1, k));
            }

            //  Arguments in Carlson F and D functions are x = \cos^2 \phi, y = 1 - k^2 \sin^2 \phi, z = 1, z = 1
            double s = MoreMath.Sin(phi);
            double x = 1.0 - s * s;
            double sk = s * k;
            double sk2 = sk * sk;
            double y = 1.0 - sk2;

            // I am a little worried that there could be cases where the cancelation between these two terms is significant
            return (s * (CarlsonF(x, y, 1.0) - sk2 / 3.0 * CarlsonD(x, y, 1.0)));

        }

    }

}