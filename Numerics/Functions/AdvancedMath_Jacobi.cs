using System;
using System.Diagnostics;

namespace Meta.Numerics.Functions {
    public static partial class AdvancedMath {

        // This implements Fukushima's method of computing the Jacobian elliptic functions
        // sn, cn, and dn. The method is essentially to use range reduction to move u into
        // the range -K/2 < u < K/2, then repeatedly halve that value of u until it gets
        // small enough to apply the Maclaurin series for small m to compute sn, then
        // using the doubling formula to get back to the value of u required.

        /// <summary>
        /// Computes the Jacobian elliptic function sn.
        /// </summary>
        /// <param name="u">The argument.</param>
        /// <param name="k">The modulus, which must be between 0 and 1.</param>
        /// <returns>The value of sn(u,k).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> lies outside [0, 1].</exception>
        /// <remarks>
        /// <para>Be aware that some authors use the parameter m = k^2 instead of the modulus k.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Jacobi_elliptic_functions"/>
        /// <seealso href="http://mathworld.wolfram.com/JacobiEllipticFunctions.html"/>
        /// <seealso href="http://dlmf.nist.gov/22"/>
        public static double JacobiSn (double u, double k) {

            if (u < 0.0) return (-JacobiSn(-u, k));
            if ((k < 0.0) || (k > 1.0)) throw new ArgumentOutOfRangeException(nameof(k));

            // For k=1, K is infinite and our reduction algorithm breaks down.
            // But even for k = 1-10^(-16), K~19.4 and the reduction works fine,
            // so we only need to do this at exactly k=1.
            if (k == 1.0) return (Math.Tanh(u));

            // Decompose u = u_0 K + u_1, where -K/2 < u_1 < K/2, and compute sn(u_1)
            long u0; double u1;
            double sn = JacobiSn_ViaRangeReduction(u, k, out u0, out u1);

            // Transform to the appropriate quadrant
            switch (MoreMath.Mod(u0, 4)) {
                case 0:
                    return (sn);
                case 1:
                    return (Math.Sqrt((1.0 - sn * sn) / (1.0 - MoreMath.Sqr(k * sn))));
                case 2:
                    return (-sn);
                case 3:
                    return (-Math.Sqrt((1.0 - sn * sn) / (1.0 - MoreMath.Sqr(k * sn))));
                default:
                    throw new InvalidOperationException();
            }

        }

        /// <summary>
        /// Compute the Jacobian elliptic function cn.
        /// </summary>
        /// <param name="u">The argument.</param>
        /// <param name="k">The modulus, which must be between 0 and 1.</param>
        /// <returns>The value of cn(u,k).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> lies outside [0, 1].</exception>
        /// <remarks>
        /// <para>Be aware that some authors use the parameter m = k^2 instead of the modulus k.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Jacobi_elliptic_functions"/>
        /// <seealso href="http://mathworld.wolfram.com/JacobiEllipticFunctions.html"/>
        /// <seealso href="http://dlmf.nist.gov/22"/>
        public static double JacobiCn (double u, double k) {

            if (u < 0.0) u = -u;
            if ((k < 0.0) || (k > 1.0)) throw new ArgumentOutOfRangeException(nameof(k));

            if (k == 1.0) return (1.0 / Math.Cosh(u));

            // Decompose u = u_0 K + u_1, where -K/2 < u_1 < K/2, and compute sn(u_1)
            long u0;
            double u1;
            double sn = JacobiSn_ViaRangeReduction(u, k, out u0, out u1);

            switch (MoreMath.Mod(u0, 4)) {
                case 0:
                    return (Math.Sqrt(1.0 - sn * sn));
                case 1:
                    return (-Math.Sqrt((1.0 - k * k) / (1.0 - MoreMath.Sqr(k * sn))) * sn);
                case 2:
                    return (-Math.Sqrt(1.0 - sn * sn));
                case 3:
                    return (Math.Sqrt((1.0 - k * k) / (1.0 - MoreMath.Sqr(k * sn))) * sn);
                default:
                    throw new InvalidOperationException();
            }

        }

        /// <summary>
        /// Compute the Jacobian elliptic function dn.
        /// </summary>
        /// <param name="u">The argument.</param>
        /// <param name="k">The modulus, which must be between 0 and 1.</param>
        /// <returns>The value of dn(u,k).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> lies outside [0, 1].</exception>
        /// <remarks>
        /// <para>Be aware that some authors use the parameter m = k^2 instead of the modulus k.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Jacobi_elliptic_functions"/>
        /// <seealso href="http://mathworld.wolfram.com/JacobiEllipticFunctions.html"/>
        /// <seealso href="http://dlmf.nist.gov/22"/>
        public static double JacobiDn (double u, double k) {

            if (u < 0.0) u = -u;
            if ((k < 0.0) || (k > 1.0)) throw new ArgumentOutOfRangeException(nameof(k));

            if (k == 1.0) return (1.0 / Math.Cosh(u));

            // Decompose u = u_0 (K / 2) + u_1, where -K/4 < u_1 < K/4, and compute sn(u_1)
            long u0; double u1;
            //double sn = JacobiSn_ViaMoreRangeReduction(u, k, out u0, out u1);
            double sn = JacobiSn_ViaRangeReduction(u, k, out u0, out u1);

            double t = k * sn;
            /*
            switch (MoreMath.Mod(u0, 4L)) {
                case 0:
                    return (Math.Sqrt(1.0 - t * t));
                case 1:
                    double k1 = Math.Sqrt(1.0 - k * k);
                    return (Math.Sqrt(k1) * ((1.0 + k1) * Math.Sqrt(1.0 - t * t) - k * t * Math.Sqrt(1.0 - sn * sn)) / (1.0 + k1 - t * t));
                case 2:
                    return (Math.Sqrt((1.0 - k * k) / (1.0 - t * t)));
                case 3:
                    k1 = Math.Sqrt(1.0 - k * k);
                    return (Math.Sqrt(k1) * ((1.0 + k1) * Math.Sqrt(1.0 - t * t) + k * t * Math.Sqrt(1.0 - sn * sn)) / (1.0 + k1 - t * t));
                default:
                    throw new InvalidOperationException();
            }
            */
            
            if (u0 % 2 == 0) {
                return (Math.Sqrt(1.0 - t * t));
                //return (Math.Sqrt((1.0 + t) * (1.0 - t)));
            } else {
                return (Math.Sqrt((1.0 - k * k) / (1.0 - t * t)));
                //return (Math.Sqrt((1.0 + k) * (1.0 - k)) / ((1.0 + t) * (1.0 - t)));
            }
            

        }

        private static double JacobiSn_Series (double u, double m) {

            double u2 = u * u;
            double u2k = 1.0;
            double sn = 1.0;

            for (int k = 1; k < 8; k++) {
                double sn_old = sn;
                u2k *= u2;
                double c;
                switch (k) {
                    case 0:
                        c = 1.0;
                        break;
                    case 1:
                        c = -(m + 1) / 6;
                        break;
                    case 2:
                        c = (m * (m + 14.0) + 1.0) / 120.0;
                        break;
                    case 3:
                        c = -(m * (m * (m + 135.0) + 135.0) + 1.0) / 5040.0;
                        break;
                    case 4:
                        c = (m * (m * (m * (m + 1228.0) + 5478.0) + 1228.0) + 1.0) / 362880.0;
                        break;
                    case 5:
                        c = -(m * (m * (m * (m * (m + 11068.0) + 165826.0) + 165826.0) + 11069.0) + 1.0) / 39916800.0;
                        break;
                    case 6:
                        c = (m * (m * (m * (m * (m * (m + 99642.0) + 4494351.0) + 13180268.0) + 4494351.0) + 99642.0) + 1.0) / 6227020800.0;
                        break;
                    case 7:
                        c = -(m * (m * (m * (m * (m * (m * (m + 896803.0) + 116294673.0) + 834687179.0) + 834687179.0) + 116294673.0) + 896803.0) + 1.0) / 1307674368000.0;
                        break;
                    default:
                        throw new NonconvergenceException();
                }
                sn += c * u2k;

                if (sn == sn_old) {
                    return (u * sn);
                }
            }

            return (u * sn);
        }

        private static double JacobiSn_SeriesLimit (double m) {
            // This is an approximate expression for the value of u
            // that makes the final term in our series 1 ulp.
            return (Math.Pow(1.0E-16 * 6227020800.0 / (3104954.0 * m + 1.0), 1.0 / 12.0));
            // The upshot is that the limit is about ~0.07 for m = 0 and ~0.30 for m = 1
        }

        private static double JacobiSn_ReduceToSeries (double u, double m) {

            // Halve u until it gets small enough to be computable via our series.
            // Then use the doubling formula to get back the original value.
            double um = JacobiSn_SeriesLimit(m);

            int n = 0;
            while (Math.Abs(u) > um) {
                n++;
                u = u / 2.0;
            }

            Debug.Assert(n < 8);

            double sn = JacobiSn_Series(u, m);

            for (int i = 0; i < n; i++) {
                double sn2 = sn * sn;
                sn = 2.0 * sn * Math.Sqrt((1.0 - sn2) * (1.0 - m * sn2)) / (1.0 - m * sn2 * sn2);
            }

            // For n >> 1, we would loose a lot of accuracy this way. But as long
            // as n is small, we will only loose about the last ~n bits.

            return (sn);

        }

        private static double JacobiSn_ViaRangeReduction (double u, double k, out long u0, out double u1) {

            double m = k * k;

            // First we compare u to a quick-to-compute lower bound for K / 2.
            // If it's below the bound, we can move directly to series for sn without having
            // to compute K or perform range reduction.
            double K = LowerBoundK(m);
            if (u <= K / 2.0) {
                u0 = 0;
                u1 = u;
            } else {
                // Too bad, we need to actually compute K and range reduce.
                K = AdvancedMath.EllipticK(k);
                double v = u / K;
                double v0 = Math.Round(v);
                if (v < Int64.MaxValue) {
                    u0 = (long) v0;
                    u1 = (v - v0) * K;
                } else {
                    u0 = Int64.MaxValue;
                    u1 = 0.0;
                }
            }
            Debug.Assert(Math.Abs(u1) <= K / 2.0);
            // Note that for u >> K, this has the same problem as naive trig function calculations
            // for x >> 2 \pi: we loose a lot of accuracy for u1 because it is computed via subtraction.
            // There is not much we can do about this, though, short of moving to arbitrary-precision
            // arithmetic, because K is different for each value of m.

            // Compute sn of the reduced -K/2 < u1 < K/2
            return (JacobiSn_ReduceToSeries(u1, m));

        }

        private static double JacobiSn_ViaMoreRangeReduction (double u, double k, out long u0, out double u1) {
            double m = k * k;

            // First we compare u to a quick-to-compute lower bound for K / 2.
            // If it's below the bound, we can move directly to series for sn without having
            // to compute K or perform range reduction.
            double K = LowerBoundK(m);
            if (Math.Abs(u) <= K / 4.0) {
                u0 = 0;
                u1 = u;
            } else {
                // Too bad, we need to actually compute K and range reduce.
                K = AdvancedMath.EllipticK(k);
                double v = u / (K / 2.0);
                double v0 = Math.Round(v);
                if (v < Int64.MaxValue) {
                    u0 = (long) v0;
                    u1 = (v - v0) * (K / 2.0);
                } else {
                    u0 = Int64.MaxValue;
                    u1 = 0.0;
                }
            }
            Debug.Assert(Math.Abs(u1) <= K / 4.0);
            // Note that for u >> K, this has the same problem as naive trig function calculations
            // for x >> 2 \pi: we loose a lot of accuracy for u1 because it is computed via subtraction.
            // There is not much we can do about this, though, short of moving to arbitrary-precision
            // arithmetic, because K is different for each value of m.

            // Compute sn of the reduced -K/2 < u1 < K/2
            return (JacobiSn_ReduceToSeries(u1, m));
        }

        private static double LowerBoundK (double m) {

            // This lower bound on K is computed by the first few terms of the series expansion
            // for small m (which is a lower bound because all terms are positive), and the
            // inequality DLMF 19.9.2 for large m.
            // It's always within about ~0.01 of the true K.

            if (m < 0.5) {
                return (Math.PI / 2.0 * (1.0 + m * (0.25 + m * (9.0 / 64.0 + m) * 25.0 / 256.0)));
            } else {
                double m1 = 1.0 - m;
                return ((1.0 + m1 / 8.0) * Math.Log(4.0 / Math.Sqrt(m1)));
            }

        }

    }
}
