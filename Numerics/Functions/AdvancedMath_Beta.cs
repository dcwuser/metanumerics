using System;
using System.Diagnostics;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the Beta function.
        /// </summary>
        /// <param name="x">The first parameter, which must be positive.</param>
        /// <param name="y">The second parameter, which must be positive.</param>
        /// <returns>The beta function B(x,y).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> or <paramref name="y"/> is non-positive.</exception>
        /// <remarks>
        /// <para>The Beta function can be defined by the integral:</para>
        /// <img src="..\images\BetaIntegral.png" />
        /// <para>Equivalently, it can be defined as a commonly occurring ratio of Gamma functions:</para>
        /// <img src="..\images\BetaGammaRelation.png" />
        /// <para>When evaluating such a ratio of Gamma functions, it is better to use this method than to call
        /// <see cref="Gamma(double)"/> three times and form the ratio explicitly. One reason is that calling this method
        /// will be faster. Another reason is that, for many values, the individual Gamma functions will overflow
        /// even though the Beta function does not; this method will not overflow in such cases. There are still
        /// other cases in which the value of Beta does overflow or underflow a <see cref="double"/>; in such cases,
        /// the method <see cref="LogBeta(double, double)"/> will still return an accurate value of its logarithm.</para>
        /// <para>One place this ratio occurs is in the expression for a binomial coefficient in terms of factorials, so the Beta
        /// function can used to generalize binomial coefficients (<see cref="AdvancedIntegerMath.BinomialCoefficient(int, int)"/>)
        /// to non-integer values.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Beta_function"/>
        /// <seealso href="http://mathworld.wolfram.com/BetaFunction.html"/>
        /// <seealso href="http://dlmf.nist.gov/5.12">DLMF on the Beta Function</seealso>
        public static double Beta(double x, double y) {

            if (Double.IsNaN(x) || Double.IsNaN(y)) return Double.NaN;

            if (x > y) Global.Swap(ref x, ref y);
            Debug.Assert(x <= y);

            if (x < 0.0) {
                return Math.PI / MoreMath.SinPi(x) / y / Beta(1.0 - x, x + y);
            }

            if (y < 0.0) {
                return Math.PI / MoreMath.SinPi(y) / x / Beta(x + y, 1.0 - y);
            }

            Debug.Assert(x >= 0.0);
            Debug.Assert(y >= 0.0);
            Debug.Assert(x <= y);

            double z = x + y;

            if (x < GammaSeries.Limit) {
                if (z < GammaSeries.Limit) {
                    // Close enough to the origin, transform both x and y to again get a beta near (1,1)
                    Debug.Assert(x <= GammaSeries.Limit && y <= GammaSeries.Limit && z < GammaSeries.Limit);
                    return Math.PI * MoreMath.SinPi(z) / (MoreMath.SinPi(x) * MoreMath.SinPi(y) * GammaSeries.BetaOnePlus(-x, -y) * (1.0 - z));
                } else {
                    // CLose enough to the axis (but not to the origin) 
                    // Note \pi, \sin(\pix), and \Gamma(1-x) are all O(1), so no need to worry about cancelling overflows
                    return Math.PI / (MoreMath.SinPi(x) * GammaSeries.GammaOnePlus(-x) * Pochhammer(y, x));
                }
            } else if (x >= 16.0) {
                // For both arguments large, use Stirling
                Debug.Assert(x >= 16.0 && y >= 16.0);
                return Stirling.Beta(x, y);
            } else if (Math.Abs(x - 1.0) < GammaSeries.Limit && Math.Abs(y - 1.0) < GammaSeries.Limit && Math.Abs(z - 2.0) < GammaSeries.Limit) {
                // Near enough to (1, 1), use the series in both x and y
                // We test for this last because the test itself is more involved
                return GammaSeries.BetaOnePlus(x - 1.0, y - 1.0);
            } else {
                // Everywhere else, fall back to Lanczos
                // This does leave some spots where one of the component Gamma functions has an argument that is near 1, 2, or in Stirling region,
                // so this logic does compute that component less accurately than it could. But other compomponents would still be Lanczos, so
                // there is not much to gain by trying to compute that one component extra-carefully.
                return Lanczos.Beta(x, y);
                // Here are some point that land here that I wish didn't: B(1,2), B(1,3) B(0.76,0.76)
            }

        }



        /// <summary>
        /// Computes the logarithm of the Beta function.
        /// </summary>
        /// <param name="x">The first parameter, which must be positive.</param>
        /// <param name="y">The second parameter, which must be positive.</param>
        /// <returns>The value of ln(B(x,y)).</returns>
        /// <remarks>
        /// <para>This function accurately computes ln(B(x,y)) even for values of x and y for which B(x,y) is
        /// too small or large to be represented by a double.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> or <paramref name="y"/> is negative or zero.</exception>
        /// <seealso cref="Beta(System.Double,System.Double)"/>
        public static double LogBeta(double x, double y) {

            if (x < 0) throw new ArgumentOutOfRangeException(nameof(x));
            if (y < 0) throw new ArgumentOutOfRangeException(nameof(y));

            if (Double.IsNaN(x) || Double.IsNaN(y)) return Double.NaN;

            if (y < x) Global.Swap<double>(ref x, ref y);
            Debug.Assert(x <= y);

            if (x <= 0.5) {
                // Close to the axis, the log term dominates as the result blows up.
                if (x == 0.0) return Double.PositiveInfinity;
                return LogBeta(1.0 + x, y) - Math.Log(x / (x + y));
            } else if (Math.Abs(x - 1.0) < GammaSeries.Limit && Math.Abs(y - 1.0) < GammaSeries.Limit && Math.Abs(x + y - 2.0) < GammaSeries.Limit) {
                // Use the series expansion in the box near (1, 1)
                return (GammaSeries.LogBetaOnePlus(x - 1.0, y - 1.0));
            } else if (x < 16.0) {
                return Lanczos.LogBeta(x, y);
            } else {
                // When both arguments are large, use the Stirling series.
                Debug.Assert((x >= 16.0) && (y >= 16.0));
                return Stirling.LogBeta(x, y);
            }

            /*
            if ((a > 16.0) && (b > 16.0)) {
                return (Stirling.LogBeta(a, b));
            } else if (a < 0.0) {
                throw new ArgumentOutOfRangeException(nameof(a));
            } else if (b < 0.0) {
                throw new ArgumentOutOfRangeException(nameof(b));
            } else if (Double.IsNaN(a) || Double.IsNaN(b)) {
                return (Double.NaN);
            } else if (a == 0.0 || b == 0.0) {
                return (Double.PositiveInfinity);
            } else {
                return (Lanczos.LogBeta(a, b));
            }
            */
        }

    }

}