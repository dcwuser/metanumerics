using System;

namespace Meta.Numerics.Extended {

    /// <summary>
    /// Contains method for computing advanced functions to double double accuracy.
    /// </summary>
    /// <remarks>
    /// <para>In order to take advantage of the increased accurary of these functions over the corresponding
    /// <see cref="Double"/>-based functions, you must ensure that the arguments you provide are also more
    /// accurate than <see cref="Double"/>. Because decimal numbers expressed in code are automatically
    /// intrepreted by the compiler as <see cref="Double"/>, it's easier than you might think do this wrong.
    /// For example, invoking <c>Gamma(0.2)</c> will <i>not</i> compute &#x393;(1/5) to <see cref="DoubleDouble"/>
    /// precision. The reason is that 0.2 is intrepreted by the compiler as a <see cref="Double"/>, and there
    /// is no <see cref="Double"/> value that is precisely 1/5. Instead, 0.2 parsed as a <see cref="Double"/>,
    /// is stored as 3602879701896397 X 2<sup>-54</sup> = 0.20000000000000001110223024625157... This equals
    /// 0.2 to within the 16 decimal-place accuracy of <see cref="Double"/>, but clearly not to within the 32
    /// decimal-place accuracy of <see cref="DoubleDouble"/>.</para>
    /// <para>There are a number of ways to ensure that you are providing an argument to <see cref="DoubleDouble"/> accuracy.
    /// One possibility is to use the <see cref="DoubleDouble"/> text parser, for example by invoking <c>new DoubleDouble("0.2")</c>.
    /// Another is to produce tha argument as the result of calculation from exact integers, <c>(DoubleDouble) 1 / 5</c>, which works
    /// because 1 and 5 (like all integers within range) are represented exactly and the because of the cast the division operation
    /// is <see cref="DoubleDouble.op_Division"/>.
    /// (The stored value in these cases is again not precisely 1/5, but is equal within the accuracy of <see cref="DoubleDouble"/>.)
    /// Finally, if you know that the argument you want is precisely represetable as a <see cref="Double"/>, you can safely
    /// use the compiler's parser. For example, invoking Gamma(0.25) does compute &#x393;(1/4) to <see cref="DoubleDouble"/>
    /// precision, because 1/4 is exactly representable by <see cref="Double"/>.</para>
    /// </remarks>
    public static partial class AdvancedDoubleDoubleMath {

        /// <summary>
        /// The double double value of Euler's gamma.
        /// </summary>
        /// <seealso href="https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant"/>
        public static readonly DoubleDouble EulerGamma = DoubleDouble.Parse("0.57721566490153286060651209008240243");

        /// <summary>
        /// Computes the error function with double double precision.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of erf(x).</returns>
        /// <seealso cref="Meta.Numerics.Functions.AdvancedMath.Erf(double)"/>
        public static DoubleDouble Erf (DoubleDouble x) {
            if (x < 0.0) {
                return (-Erf(-x));
            } else if (x < 4.0) {
                return (Erf_Series(x));
            } else {
                return (1.0 - Erfc_ContinuedFraction(x));
            }
        }

        /// <summary>
        /// Computes the complementary error function with double double precision.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of erfc(x).</returns>
        /// <seealso cref="Meta.Numerics.Functions.AdvancedMath.Erfc(double)"/>
        public static DoubleDouble Erfc (DoubleDouble x) {
            if (x < 4.0) {
                return (1.0 - Erf_Series(x));
            } else {
                return (Erfc_ContinuedFraction(x));
            }
        }

        private static DoubleDouble Erf_Series (DoubleDouble x) {
            DoubleDouble mx2 = -(x * x);
            DoubleDouble t = DoubleDouble.One;
            DoubleDouble f = t;
            for (int k = 1; k < 100; k++) {
                DoubleDouble f_old = f;
                t *= mx2 / k;
                f += t / (2 * k + 1);
                if (f == f_old) {
                    return (2.0 / DoubleDouble.Sqrt(DoubleDouble.Pi)) * x * f;
                }
            }
            throw new InvalidOperationException();
        }

        private static DoubleDouble Erfc_ContinuedFraction (DoubleDouble x) {
            DoubleDouble x2 = x * x;
            double aa = 1.0;
            DoubleDouble bb = x2 + 0.5;
            DoubleDouble D = 1.0 / bb;
            DoubleDouble Df = aa / bb;
            DoubleDouble f = Df;
            for (int k = 1; k < 100; k++) {
                DoubleDouble f_old = f;
                aa = -k * (k - 0.5);
                bb += 2.0;
                D = 1.0 / (bb + aa * D);
                Df = (bb * D - 1.0) * Df;
                f += Df;
                if (f == f_old) {
                    DoubleDouble e = DoubleDouble.Exp(-x2);
                    DoubleDouble g = x / DoubleDouble.Sqrt(DoubleDouble.Pi);
                    return (e * f * g);
                }
            }
            throw new InvalidOperationException();
        }
    }

}
