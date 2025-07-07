using System;
using System.Diagnostics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {


    /// <summary>
    /// Represents the asymptotic distribution of Kuiper's V statistic.
    /// </summary>
    public sealed class KuiperDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new Kuiper distribution.
        /// </summary>
        public KuiperDistribution () {
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return 0.0;
            } else if (x < 1.0) {
                return AsymptoticP(x);
            } else {
                return 1.0 - AsymptoticQ(x);
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return 1.0;
            } else if (x < 1.0) {
                return 1.0 - AsymptoticP(x);
            } else {
                return AsymptoticQ(x);
            }
        }

        // series \sqrt{2\pi}{x^3} \sum_{k=1}^{\infty} k^2 \pi^2 e^{-k^2 \pi^2 / 2 x^2}
        // useful for small x

        // For small x,
        //    P = \frac{\sqrt{2\pi}}{x} \sum_{k=1}^{\infty} \frac{k^2 \pi^2}{x^2} e^{-k^2 \pi^2 / 2 x^2}
        // W

        private static double AsymptoticP (double x) {

            Debug.Assert(x > 0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;
                double z = k * Math.PI / x;
                double z2 = z * z;
                double ds = z2 * Math.Exp(-0.5 * z2);
                s += ds;

                if (s == s_old) return Global.SqrtTwoPI / x * s;

            }

            throw new NonconvergenceException();

        }

        // series \sum_{k=1}^{\infty} (4 k^2 x^2 - 1) e^{-2 k^2 x^2}
        // useful for large x

        private static double AsymptoticQ (double x) {

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (4.0 * z2 - 1.0) * Math.Exp(-2.0 * z2);
                s += ds;

                if (s == s_old) return 2.0 * s;
            }

            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return 0.0;
            } else if (x < 1.0) {
                return AsymptoticPPrime(x);
            } else {
                return AsymptoticQPrime(x);
            }
        }

        private static double AsymptoticQPrime (double x) {

            double s = 0.0;
            for (int k = 1; k < 20; k++) {

                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (k * k) * (4.0 * z2 - 3.0) * Math.Exp(-2.0 * z2);
                s += ds;

                if (s == s_old) return 8.0 * x * s;
            }

            throw new NonconvergenceException();

        }

        private static double AsymptoticPPrime (double x) {

            Debug.Assert(x > 0.0);

            double s = 0.0;
            for (int k = 1; k < 20; k++) {

                double s_old = s;
                double z = Math.PI * k / x;
                double z2 = z * z;
                double ds = z2 * (z2 - 3.0) * Math.Exp(-z2 / 2.0);
                s += ds;

                if (s == s_old) return Global.SqrtTwoPI / (x * x) * s;

            }

            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return Interval.Semiinfinite;
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                // originally we got this via numerical integration, later recognized it as \sqrt{\frac{\pi}{2}}
                // this can be understood from moment formula by noting that \lim_{x \rightarrow 1} \zeta(x) (x-1) = 1
                return Global.SqrtHalfPI;
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return Math.PI / 2.0 * (Math.PI / 3.0 - 1.0);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (3.0 * AdvancedMath.Apery - Math.PI * (Math.PI - 2.0)) / (Math.PI * Math.Pow(Math.PI / 3.0 - 1.0, 3.0 / 2.0));
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else if (r == 1) {
                return Mean;
            } else if (r == 2) {
                return Math.PI * Math.PI / 6.0;
            } else {
                return AdvancedMath.RiemannZeta(r) * AdvancedMath.Gamma(1 + r / 2.0) * (r - 1) / Math.Pow(2.0, r / 2.0 - 1);
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else if (r == 1) {
                return 0.0;
            } else if (r == 2) {
                return Variance;
            } else {
                return CentralMomentFromRawMoment(r);
            }
        }

    }

    // Kuiper limiting distribution given by
    //   P_0 = \sum_{k=-\infty}^{+\infty} \left( 1 - 4 k^2 x^2 \right) e^{-2 k^2 x^2} = 1 - 2 \sum_{k=1}^{\infty} \left( 4 k^2 x^2 - 1 \right) e^{-2 k^2 x^2}
    // This form is useful for large x. To transform it into a form useful for small x, we apply the same technique as in Kolmogorov case. Note that
    //   \sum_{k=-\infty}^{+\infty} e^{-2 k^2 x^2} = \vartheta(0, 2 i x^2 / \pi)
    // where \vartheta is a theta function (http://en.wikipedia.org/wiki/Theta_function). Apply Jacobi transformation to get
    //   \sum_{k=-\infty}^{+\infty} e^{-2 k^2 x^2} = \frac{1}{x} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} e^{-\pi^2 k^2 / 2 x^2}
    // and take derivatives to get
    //   \sum_{k=-\infty}^{+\infty} k^2 e^{-2 k^2 x^2}
    //     = \frac{1}{4 x^3} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \left( 1 - \frac{\pi^2 k^2}{x^2} \right) e^{-\pi^2 k^2 / 2 x^2}
    //   \sum_{k=-\infty}^{+\infty} k^4 e^{-2 k^2 x^2}
    //     = \frac{1}{16 x^5} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \left( 3 - 6 \frac{\pi^2 k^2}{x^2} + \frac{\pi^4 k^4}{x^4} \right) e^{-\pi^2 k^2 / 2 x^2}
    // Plugging in the transformed expression for each required power of k gives
    //   P_0 = \frac{1}{x} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \frac{\pi^2 k^2}{x^2} e^{-\pi^2 k^2 / 2 x^2}
    // a form useful for small x.

    // The first correction is given by \frac{4}{3\sqrt{n}} P_1 where
    //   P_1 = \sum_{k=-\infty}^{+\infty} k^2 \left( 4 k^2 x^2 - 3 \right) e^{-2 k^2 x^2}
    // Applying the same transformation rules for each required power of k gives
    //   P_1 = \frac{1}{4 x^3} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \frac{\pi^2}{k^2}{x^2} \left( \frac{\pi^2}{k^2}{x^2} - 3 \right) e^{-\pi^2 k^2 / 2 x^2}

    // To compute the PDF instead of the CDF, just take derivatices. For the leading term
    //   p_0 = P_0' = 4 x \sum_{k=-\infty}^{+\infty} k^2 ( 4 k^2 x^2 - 3) e^{-2 k^2 x^2}
    //     = \frac{1}{x^2} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \left( \frac{\pi^2}{k^2}{x^2} - 3 \right) e^{-\pi^2 k^2 / 2 x^2}
    // (Note that P_1 = 4 x p_0.) And for the sub-leading term
    //   p_1 = P_1' = 4 x \sum_{k=-\infty}^{+\infty} k^4 ( 5 - 4 k^2 x^2 ) e^{-2 k^2 x^2}
    //     = \frac{1}{4 x^4} \sqrt{\frac{\pi}{2}}  \sum_{k=-\infty}^{+\infty} \frac{\pi^2}{k^2}{x^2} \left( \frac{\pi^4 k^4}{x^4} - 10 \frac{\pi^2 k^2}{x^2} + 15 \right) e^{-\pi^2 k^2 / 2 x^2}

    // To compute moments we need two facts
    //    \int_{0}^{\infty} x^m e^{-2 k^2 x^2}
    //    \sum_{k=1}^{\infty} \frac{1}{k^m} = \zeta(m)
    // This implies that the limiting moments are
    //    <x^m>_0 = \Gamma(m/2 + 1) (m - 1) / 2^{m/2 - 1} \zeta(m)
    // and the first correction is \frac{4}{3\sqrt{n}} times
    //    <x^m>_1 = - \Gamma(m/2 + 1) (m - 3) / 2^{m/2 - 1} \zeta(m - 2)


    internal class KuiperAsymptoticDistribution : ContinuousDistribution {

        public KuiperAsymptoticDistribution (int n) {
            if (n < 2) throw new ArgumentOutOfRangeException(nameof(n));
            this.n = n;
            this.sqrt_n = Math.Sqrt(n);
        }

        private readonly int n;
        private readonly double sqrt_n;

        public override Interval Support {
            get {
                return Interval.Semiinfinite;
            }
        }

        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x < 1.2) {
                return (P0(x) + P1(x) / sqrt_n);
            } else {
                return (1.0 - Q0(x) - Q1(x) / sqrt_n);
            }
        }

        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return(1.0);
            } else if (x < 1.2) {
                return(1.0 - P0(x) - P1(x) / sqrt_n);
            } else {
                return(Q0(x) + Q1(x) / sqrt_n);
            }
        }

        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x < 1.2) {
                return (P0Prime(x) + P1Prime(x) / sqrt_n);
            } else {
                return (-Q0Prime(x) - Q1Prime(x) / sqrt_n);
            }
        }

        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (Global.SqrtHalfPI * ( 1.0 - 1.0 / 6.0 / sqrt_n));
            } else if (r == 3) {
                return (Global.SqrtHalfPI * (3.0 / 2.0 * AdvancedMath.Apery - 3.0 / 4.0 / sqrt_n));
            }  else {
                // we needed to handle 1st and 3rd moments specially because \zeta divergences but multiplication by zero gives finite result
                return (AdvancedMath.Gamma(r / 2.0 + 1.0) / Math.Pow(2.0, r / 2.0 - 1.0) *
                    ((r - 1) * AdvancedMath.RiemannZeta(r) - (r - 3) * AdvancedMath.RiemannZeta(r - 2) / sqrt_n)
                );
            }
        }

        /*
        public override double Variance {
            get {
                return ((Math.PI / 3.0 - 1.0) * (Math.PI / 2.0 + 1.0 / 2.0 / sqrt_n));
                // the (1/36)(1/n) term is beyond our order of accuracy but we include it so that results are self-consistent
            }
        }
        */

        public static double Q0 (double x) {

            if (x >= Global.SqrtMax) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (4.0 * z2 - 1.0) * Math.Exp(-2.0 * z2);
                s += ds;
                if (s == s_old) { return (2.0 * s); }
            }
            throw new NonconvergenceException();

        }

        public static double P0 (double x) {

            if (x <= 1.0 / Global.SqrtMax) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                double z = Math.PI * k / x;
                double z2 = z * z;
                double ds = z2 * Math.Exp(-z2 / 2.0);
                s += ds;
                if (s == s_old) { return (Global.SqrtTwoPI / x * s); }
            }
            throw new NonconvergenceException();

        }

        public static double Q0Prime (double x) {

            return (4.0 * x * Q1(x));

        }

        public static double P0Prime (double x) {

            if (x <= 1.0 / Global.SqrtMax) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                double z = Math.PI * k / x;
                double z2 = z * z;
                double ds = z2 * (z2 - 3.0) * Math.Exp(-z2 / 2.0);
                s += ds;
                if (s == s_old) { return (Global.SqrtTwoPI / (x * x) * s); }
            }
            throw new NonconvergenceException();

        }

        public static double P1 (double x) {

            return (P0Prime(x) / (4.0 * x));

        }

        public static double P1Prime (double x) {

            if (x <= 1.0 / Global.SqrtMax) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                double z = Math.PI * k / x;
                double z2 = z * z;
                double ds = z2 * (z2 * z2 - 10.0 * z2 + 15.0) * Math.Exp(-z2 / 2.0);
                s += ds;
                if (s == s_old) { return (Global.SqrtTwoPI / MoreMath.Pow(x, 4) / 4.0 * s); }
            }
            throw new NonconvergenceException();


        }

        public static double Q1 (double x) {

            if (x >= Global.SqrtMax) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (k * k) * (4.0 * z2 - 3.0) * Math.Exp(-2.0 * z2);
                s += ds;
                if (s == s_old) return (-2.0 * s);
            }
            throw new NonconvergenceException();

        }

        public static double Q1Prime (double x) {

            if (x >= Global.SqrtMax) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (k * k * k * k) * (5.0 - 4.0 * z2) * Math.Exp(-2.0 * z2);
                s += ds;
                if (s == s_old) return (-8.0 * x * s);
            }
            throw new NonconvergenceException();

        }

    }

}
