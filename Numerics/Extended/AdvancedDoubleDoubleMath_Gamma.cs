using System;
using System.Diagnostics;

namespace Meta.Numerics.Extended {


    // A key to our algorithm to compute Gamma for double arguments is the Lanczos algorithm.
    // This doesn't map easily to the DoubleDouble case, because we don't have the benefit of
    // a good set of coefficients computed by others. I also looked at Spouge's algorithm, but
    // it requires many more digits of precision for intermediates than for the result. So
    // I've fallen back to zeta series around 1 and 2, and the asymptotic approximation. For
    // results in other regiemes, I just use z \Gamma(z) = \Gamma(z+1) to step either inward or
    // outward. Surprisingly, even though this introduces cancellation from subtraction, we
    // don't seem to lose more than one decimal digit in the last place, when comparing to
    // Mathematica results.

    public static partial class AdvancedDoubleDoubleMath {

        /// <summary>
        /// Computes the logarithm of the Gamma function double double precision.
        /// </summary>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of Gamma(x).</returns>
        public static DoubleDouble LogGamma (DoubleDouble x) {

            if (x < DoubleDouble.Zero) {
                throw new ArgumentOutOfRangeException(nameof(x));
            } else if (x < 0.5) {
                // Use x G(x) = G(1+x) to get G(x) from G(1+x)
                return LogGammaOnePlus(x) - DoubleDouble.Log(x);
            } else if (x < 1.5) {
                DoubleDouble y = x - 1.0;
                Debug.Assert(DoubleDouble.Abs(y) <= 0.5);
                return LogGammaOnePlus(y);
            } else if (x < 16.0) {
                return LogGamma_ShiftToZetaSeries(x);
            } else {
                return LogGamma_ShiftToAsymptotic(x);
            }
        }

        // Series Regime

        // We need quite a few [zeta(k) - 1] values, up to k~60 for |z|=1/2. Rather
        // than store text and parse to get all of them, I do that for the first few,
        // then use the zeta series itself to compute [zeta(k) - 1] for higher k
        // on demand as required.

        private static readonly DoubleDouble[] zetaMinusOne = InitializeZetaMinusOne();

        private static DoubleDouble[] InitializeZetaMinusOne () {
            DoubleDouble[] zetaMinusOne = new DoubleDouble[64];
            zetaMinusOne[0] = -1.5;
            zetaMinusOne[1] = Double.PositiveInfinity;
            zetaMinusOne[2] = new DoubleDouble("0.64493406684822643647241516664602519");
            zetaMinusOne[3] = new DoubleDouble("0.20205690315959428539973816151144999");
            zetaMinusOne[4] = new DoubleDouble("0.082323233711138191516003696541167903");
            zetaMinusOne[5] = new DoubleDouble("0.036927755143369926331365486457034168");
            zetaMinusOne[6] = new DoubleDouble("0.017343061984449139714517929790920528");
            zetaMinusOne[7] = new DoubleDouble("8.3492773819228268397975498497967596E-3");
            zetaMinusOne[8] = new DoubleDouble("4.0773561979443393786852385086524653E-3");
            zetaMinusOne[9] = new DoubleDouble("2.0083928260822144178527692324120605E-3");
            zetaMinusOne[10] = new DoubleDouble("9.9457512781808533714595890031901701E-4");
            zetaMinusOne[11] = new DoubleDouble("4.9418860411946455870228252646993647E-4");
            zetaMinusOne[12] = new DoubleDouble("2.4608655330804829863799804773967096E-4");
            zetaMinusOne[13] = new DoubleDouble("1.2271334757848914675183652635739571E-4");
            zetaMinusOne[14] = new DoubleDouble("6.1248135058704829258545105135333747E-5");
            zetaMinusOne[15] = new DoubleDouble("3.0588236307020493551728510645062588E-5");
            return zetaMinusOne;
        }

        private static DoubleDouble ZetaMinusOne (int n) {
            // For n < ~16, we would needs more than 255 terms.
            // Look into using Euler-Maclauren to accelerate.
            DoubleDouble s = DoubleDouble.Zero;
            for (int k = 2; k < Global.SeriesMax; k++) {
                DoubleDouble s_old = s;
                s += DoubleDouble.Pow(k, -n);
                if (s == s_old) return s;
            }
            throw new NonconvergenceException();
        }

        private static DoubleDouble ZetaSeries (DoubleDouble x) {
            DoubleDouble s = 0.0;
            DoubleDouble xMinus = -x;
            DoubleDouble xPower = xMinus;
            for (int k = 2; k < zetaMinusOne.Length; k++) {
                DoubleDouble s_old = s;
                xPower *= xMinus;
                // If a [\zeta(k) - 1] value is not yet computed, compute it as needed.
                if (zetaMinusOne[k] == DoubleDouble.Zero) {
                    // Technically this is not thread-safe, because assignment is not atomic for non-native structs.
                    // But at worst we are filling in the same value from two different threads, so this would only
                    // be a problem if intermediate values are not directly replaced during assignment.
                    zetaMinusOne[k] = ZetaMinusOne(k);
                }
                s += zetaMinusOne[k] * xPower / k;
                if (s == s_old) return s;
            }
            throw new NonconvergenceException();
        }

        private static DoubleDouble LogGammaOnePlus (DoubleDouble y) {
            return LogGammaTwoPlus(y) - DoubleDouble.Log1P(y);
        }

        private static DoubleDouble LogGammaTwoPlus (DoubleDouble y) {
            return (DoubleDouble.One - AdvancedDoubleDoubleMath.EulerGamma) * y + ZetaSeries(y);
        }

        private static DoubleDouble LogGamma_ShiftToZetaSeries (DoubleDouble x) {

            Debug.Assert(x >= 1.5);

            DoubleDouble s = DoubleDouble.Zero;
            while (x > 2.5) {
                x -= DoubleDouble.One;
                s += DoubleDouble.Log(x);
            }

            DoubleDouble y = x - 2.0;
            Debug.Assert(DoubleDouble.Abs(y) <= 0.5);
            return (LogGammaTwoPlus(y) + s);
        }

        // Asymptotic Regime

        private static readonly DoubleDouble[] Bernoulli = new DoubleDouble[] {
            DoubleDouble.One, DoubleDouble.One / 6, -DoubleDouble.One / 30, DoubleDouble.One / 42, -DoubleDouble.One / 30,
            ((DoubleDouble) 5) / 66, -((DoubleDouble) 691)  / 2730, ((DoubleDouble) 7) / 6, -((DoubleDouble) 3617) / 510, ((DoubleDouble) 43867) / 798,
            -((DoubleDouble) 74611) / 330, ((DoubleDouble) 854513) / 138, -((DoubleDouble) 236364091) / 2730, ((DoubleDouble) 8553103) / 6, -((DoubleDouble) 23749461029) / 870
        };

        public static DoubleDouble BernoulliSum (DoubleDouble x) {
            DoubleDouble rxPower = 1.0 / x;
            DoubleDouble rxSquared = rxPower * rxPower;
            DoubleDouble f = 0.5 * Bernoulli[1] * rxPower;
            for (int k = 2; k < Bernoulli.Length; k++) {
                DoubleDouble f_old = f;
                rxPower *= rxSquared;
                f += Bernoulli[k] / ((2 * k) * (2 * k - 1)) * rxPower;
                if (f == f_old) {
                    return (f);
                }
            }
            throw new NonconvergenceException();
        }

        private static readonly DoubleDouble halfLogTwoPi = 0.5 * DoubleDouble.Log(2.0 * DoubleDouble.Pi);

        private static DoubleDouble LogGamma_Asymptotic (DoubleDouble x) {
            // Sum from smallest to largest terms to minimize error.
            return (BernoulliSum(x) + halfLogTwoPi - x + (x - 0.5) * DoubleDouble.Log(x));
        }

        private static DoubleDouble LogGamma_ShiftToAsymptotic (DoubleDouble x) {

            Debug.Assert(x > 0.0);

            DoubleDouble s = DoubleDouble.Zero;
            while (x < 34.0) {
                s += DoubleDouble.Log(x);
                x += DoubleDouble.One;
            }

            return (LogGamma_Asymptotic(x) - s);

        }
    }
}
