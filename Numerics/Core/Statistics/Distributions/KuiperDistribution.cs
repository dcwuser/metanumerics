using System;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {


    /// <summary>
    /// Represents the asymptotic distribution of Kuiper's V statistic.
    /// </summary>
    public class KuiperDistribution : Distribution {

        /// <summary>
        /// Initializes a new Kuiper distribution.
        /// </summary>
        public KuiperDistribution () {
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x < 1.0) {
                return (AsymptoticP(x));
            } else {
                return (1.0 - AsymptoticQ(x));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else if (x < 1.0) {
                return (1.0 - AsymptoticP(x));
            } else {
                return (AsymptoticQ(x));
            }
        }

        // series \sqrt{2\pi}{x^3} \sum_{k=1}^{\infty} k^2 \pi^2 e^{-k^2 \pi^2 / 2 x^2}
        // useful for small x

        // For small x,
        //    P = \frac{\sqrt{2\pi}}{x} \sum_{k=1}^{\infty} \frac{k^2 \pi^2}{x^2} e^{-k^2 \pi^2 / 2 x^2}
        // W

        private static double AsymptoticP (double x) {

            if (x <= 0.0) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;
                double z = k * Math.PI / x;
                double z2 = z * z;
                double ds = z2 * Math.Exp(-z2 / 2.0);
                s += ds;

                if (s == s_old) return (Global.SqrtTwoPI / x * s);

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

                if (s == s_old) return (2.0 * s);
            }

            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x < 1.0) {
                return (AsymptoticPPrime(x));
            } else {
                return (AsymptoticQPrime(x));
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

                if (s == s_old) return (8.0 * x * s);
            }

            throw new NonconvergenceException();

        }

        private static double AsymptoticPPrime (double x) {

            if (x <= 0.0) return (0.0);

            double s = 0.0;
            for (int k = 1; k < 20; k++) {

                double s_old = s;
                double z = Math.PI * k / x;
                double z2 = z * z;
                double ds = z2 * (z2 - 3.0) * Math.Exp(-z2 / 2.0);
                s += ds;

                if (s == s_old) return (Global.SqrtTwoPI / (x * x) * s);

            }

            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                // originally we got this via numerical integration, later recognized it as \sqrt{\frac{\pi}{2}}
                // this can be understood from moment formula by noting that \lim_{x \rightarrow 1} \zeta(x) (x-1) = 1
                return (Global.SqrtHalfPI);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (Global.HalfPI * (Math.PI / 3.0 - 1.0));
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else if (n == 2) {
                return (Math.PI * Math.PI / 6.0);
            } else {
                return (AdvancedMath.RiemannZeta(n) * AdvancedMath.Gamma(1 + n / 2.0) * (n - 1) / Math.Pow(2.0, n / 2.0 - 1));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else if (n == 2) {
                return (Variance);
            } else {
                return (CentralMomentFromRawMoment(n));
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


    internal class KuiperAsymptoticDistribution : Distribution {

        public KuiperAsymptoticDistribution (int n) {
            if (n < 2) throw new ArgumentOutOfRangeException("n");
            this.n = n;
            this.sqrt_n = Math.Sqrt(n);
        }

        private readonly int n;
        private readonly double sqrt_n;

        public override Interval Support {
            get {
                return (Interval.FromEndpointAndWidth(0.0, Double.PositiveInfinity));
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

        public override double Moment (int m) {
            if (m < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (m == 0) {
                return (1.0);
            } else if (m == 1) {
                return (Global.SqrtHalfPI * ( 1.0 - 1.0 / 6.0 / sqrt_n));
            } else if (m == 3) {
                return (Global.SqrtHalfPI * (3.0 / 2.0 * AdvancedMath.RiemannZeta(3.0) - 3.0 / 4.0 / sqrt_n));
            }  else {
                // we needed to handle 1st and 3rd moments specially because \zeta divergences but multiplication by zero gives finite result
                return (AdvancedMath.Gamma(m / 2.0 + 1.0) / Math.Pow(2.0, m / 2.0 - 1.0) *
                    ((m - 1) * AdvancedMath.RiemannZeta(m) - (m - 3) * AdvancedMath.RiemannZeta(m - 2) / sqrt_n)
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

    // w = n V, support 1 < w < n


    // n = 2
    // 1 < w < 2    P = 1 (w - 1) = 1 - 1 (2 - w)

    // n = 3
    // 1 < w < 2    P = 2/3 (w - 1)^2
    // 2 < w < 3    P = 2/3 (-3 + 3 w - 1/2 w^2) = 1 - 1/3 (3 - w)^2

    // n = 4
    // 1 < w < 2    P = 3/8 (w - 1)^3
    // 2 < w < 3    P = 3/8 (1 - 4 w + 3 w^2 - 1/2 w^3)
    // 3 < w < 4    P = 3/8 (-8 + 8 w - 2 w^2 + 1/6 w^3) = 1 - 1/16 (4 - w)^3

    // n = 5
    // 1 < w < 2    P = 24/125 (w - 1)^4
    // 2 < 2 < 3    P = 24/125 (5 - 7 w + 3/2 w^2 + w^3 - 1/4 w^4)
    // 3 < w < 4    P = 24/125 (1/2 - 17/2 w + 31/4 w^2 - 2 w^3 + 1/6 w^4)
    // 4 < w < 5    P = 24/125 (-125/6 + 125/6 w - 25/4 w^2 + 5/6 w^3 - 1/24 w^4) = 1 - 1/125 (5 - w)^4

    // n = 6
    // 1 < w < 2    P = 5/54 (w - 1)^5
    // 2 < w < 3    P = 5/54 (-7 + 22 w - 23 w^2 + 19/2 w^3 - 5/4 w^4 + 0 w^5)
    // 3 < w < 4    P = 5/54 (35/4 - 65/4 w + 11/2 w^2 + 5/3 w^3 - 5/6 w^4 + 1/12 w^5)
    // 4 < w < 5    P = 5/54 (-23/12 - 227/12 w + 39/2 w^2 - 37/6 w^3 + 5/6 w^4 - 1/24 w^5)
    // 5 < w < 6    P = 5/54 (-54 + 54 w - 18 w^2 + 3 w^3 - 1/4 w^4 + 1/120 w^5) = 1 - 1/1296 (6 - w)^5

    // n = 7
    // 1 < w < 2    P = 720/16807 (w - 1)^6
    // 2 < w < 3    P = 720/16807 (-3 - 5 w + 51/2 w^2 - 28 w^3 + 25/2 w^4 - 9/4 w^5 + 1/8 w^6)
    // 3 < w < 4    P = 720/16807 (105/4 - 29 w + 5/2 w^2 + 9/2 w^3 - 5/8 w^4 - 1/6 w^5 + 1/36 w^6)
    // 4 < w < 5    P = 720/16807 (187/12 - 115/3 w + 97/6 w^2 + 91/36 w^3 - 35/16 w^4 + 3/8 w^5 - 1/48 w^6)
    // 5 < w < 6    P = 720/16807 (-251/24 - 1045/24 w + 2351/48 w^2 - 629/36 w^3 + 143/48 w^4 - 1/4 w^5 + 1/120 w^6)
    // 6 < w < 7    P = 1 - 1/16807 (7 - w)^6

    // At the ends P_n(1 < w < 2) = n! / n^(n-1) (w-1)^(n-1) and Q_n(n-1 < w < n) = (n-w)^(n-1) / n^(n-2)

    // Taking a derivative wrt w gives p(w) and integrating p(w) * w^m gives the mth moment.

    // n = 2: <w> = 3/2, <w^2> = 7/3, <w^3> = 15/4 => C2 = 1/12, C3 = 0
    // n = 3: <w> = 17/9 , <w^2> = 67/18, <w^3> = 229/30 => C2 = 25/162, C3 = 71/3645
    // n = 4: <w> = 71/32, <w^2> = 103/20, <w^3> = 1997/160 => C2 = 1163/5120, C3 = 3827/81920
    // n = 5: <w> = 1569/625, <w^2> = 2476/375, <w^3> = 11352/625 => C2 = 352217/1171875
    // n = 6: <w> = 899/324, <w^2> = 7847/972, <w^3> = 223099/9072 => C2 = 39275/104976
    // n = 7: <w> = 355081/117649, <w^2> = 1124366/117649, <w^3> = 11189840/352947 => C2 = 6198018973/13841287201

    // Moments provide an additional check. We should have <1> = 1 and <w> = \frac{n!}{n^n} \sum_{k=0}^{n-1} \frac{n^k}{k!}.
    // So far so good.

    internal class KuiperExactDistribution : Distribution {

        public KuiperExactDistribution (int n) {
            if (n < 2) throw new ArgumentOutOfRangeException("n");
            this.n = n;
        }

        private readonly int n;

        public override Interval Support {
            get {
                return (Interval.FromEndpoints(1.0, n));
            }
        }

        public override double Mean {
            get {
                return (ComputeMean());
            }
        }

        // As per Stephens, Biopmetrika (1965) 52, p. 309-, the mean
        //   <w> = \frac{n!}{n^n} \sum_{k=0}^{n-1} \frac{n^k}{k!}
        // This is derived by taking Birnbaum's expression of the mean of D_{\pm} and noting that <V> = <(D_+ + D_-)> = <D_+> + < D_->.
        // It's pretty cool that we can get this in simple analytic form. I'd love to do the same for <D> or <V^2>, but D=\max(D_+,D_-) does
        // not seperate into terms each involving only D_+ or D_- and V^2 involves the cross-term D_+ D_-, so both depend on the joint
        // distribution p(D_+,D_-).

        private double ComputeMean () {
            double t = 1.0;
            double s = t;
            for (int k = 1; k < n; k++) {
                t = t * n / k;
                s += t;
            }
            // Canceling factors of n, \frac{n!}{n^n} = \frac{(n-1)!}{n^{n-1}}, so the prefactor is actually just the inverse of the last term.
            return (s / t);
        }

        public override double Variance {
            get {
                switch (n) {
                    case 2:
                        // 1/12 = 0.041667 * 2
                        return (1.0 / 12.0);
                    case 3:
                        // 25/162 = 0.051440 * 3
                        return (25.0 / 162.0);
                    case 4:
                        // 1163/5120 = 0.056787 * 4
                        return (1163.0 / 5120.0);
                    case 5:
                        // 352217/1171875 = 0.060111 * 5
                        return (352217.0 / 1171875.0);
                    case 6:
                        // 39275/104976 = 0.062355 * 6
                        return (39275.0 / 104976.0);
                    case 7:
                        // 6198018973/13841287201 = 0.063970 * 7
                        return (6198018973.0 / 13841287201.0);
                    default:
                        // \frac{\pi}{2} \left( \frac{\pi}{3} - 1 \right) = 0.074138
                        // this will perform numerical integration, which will be very slow
                        return (base.MomentAboutMean(2));
                }
            }
        }

        public override double LeftProbability (double w) {
            if (w <= 1.0) {
                return (0.0);
            } else if (w < n) {
                return (DurbinMatrixP(w));
            } else {
                return (1.0);
            }
        }

        public override double RightProbability (double w) {
            if (w <= 1.0) {
                return (1.0);
            } else if (w < n) {
                return (1.0 - DurbinMatrixP(w));
            } else {
                return (1.0);
            }
        }

        public override double ProbabilityDensity (double w) {
            if (w <= 1.0) {
                return (0.0);
            } else if (w < n) {
                return (DurbinMatrixPPrime(w));
            } else {
                return (0.0);
            }
        }

        ///<inheritdoc />
        public override double Moment (int m) {
            if (m < 0) {
                throw new ArgumentOutOfRangeException("m");
            } else if (m == 0) {
                return (1.0);
            } else if (m == 1) {
                return (Mean);
            } else {
                return (base.Moment(m));
            }
        }

        ///<inheritdoc />
        public override double MomentAboutMean (int m) {
            if (m < 0) {
                throw new ArgumentOutOfRangeException("m");
            } else if (m == 0) {
                return (1.0);
            } else if (m == 1) {
                return (0.0);
            } else if (m == 2) {
                return (this.Variance);
            } else {
                return (base.MomentAboutMean(m));
            }
        }

        // Durbin derives a matrix result for the Kuiper statistic analogous to his result for the Kolmogorov-Smirnov that was used by Marsaglia.
        // In this case, if w = n V = k - h, H is k X k with entries of the form
        //       { 0  1           0           0           0        }
        //       { 0  1/1!        1           0           0        }
        //   H = { 0  1/2!        1/1!        1           0        }
        //       { 0  1/3!        1/2!        1/1!        1        }
        //       { 0  (1-h^4)/4!  (1-h^3)/3!  (1-h^2)/2!  (1-h)/1! }
        // and
        //   P = \frac{n!}{n^{n-1}} \left( H^n \right)_{1,2}
        // Note there is no discontinuity at half-integer values of w in the Kuiper case.
        // See Durbin, "Distribution Theory for Tests Based on the Sample Distribution Function", 1973, p. 12-13, 35.

        private double DurbinMatrixP (double w) {

            int k = (int) Math.Ceiling(w);
            double h = k - w;

            SquareMatrix H = GetDurbinMatrix(k, h);

            SquareMatrix Hn = H.Power(n);
            //SquareMatrix Hn = MatrixPower(H, n);

            double f = AdvancedIntegerMath.Factorial(n - 1) / MoreMath.Pow(n, n - 2);
            return (f * Hn[0, 1]);

        }

        private double DurbinMatrixPPrime (double w) {

            int k = (int) Math.Ceiling(w);
            double h = k - w;

            // compute derivative of H
            SquareMatrix DH = GetDurbinMatrixPrime(k, h);

            //PrintMatrix(DH);

            // compute powers of H
            SquareMatrix[] PowerH = new SquareMatrix[n];
            PowerH[1] = GetDurbinMatrix(k, h);
            for (int i = 2; i < n; i++) {
                PowerH[i] = PowerH[1] * PowerH[i - 1];
            }

            // use D(H^n) = (DH) H^(n-1) + H (DH) H^(n-2) + H^2 (DH) H^(n-3) + \cdots + H^(n-2) (DH) H + H^(n-1) (DH)
            SquareMatrix HnP = DH * PowerH[n - 1];
            for (int i = 1; i < (n - 1); i++) {
                HnP += PowerH[i] * DH * PowerH[n - 1 - i];
            }
            HnP += PowerH[n - 1] * DH;

            double f = AdvancedIntegerMath.Factorial(n - 1) / MoreMath.Pow(n, n - 2);
            return (f * HnP[0, 1]);


        }

        private static SquareMatrix GetDurbinMatrix (int k, double h) {

            // dimension of matrix
            int m = k;
            SquareMatrix H = new SquareMatrix(m);

            // populate the matrix along diagonals, since they share factorial factors

            // first superdiagonal is all 1s
            for (int j = 1; j < m; j++) {
                H[j - 1, j] = 1.0;
            }

            // bottom row and diagonals
            double Fi = 1.0; // variable for 1/i!
            double hi = h; // variable for h^i
            for (int i = 1; i < m; i++) {
                // bottom row
                H[m - 1, m - i] = Fi * (1.0 - hi);
                // diagonal
                for (int j = i + 1; j < m; j++) {
                    H[j - 1, j - i] = Fi;
                }
                // prepare for the next recurrsion
                hi = hi * h;
                Fi = Fi / (i + 1);
            }

            return (H);

        }

        private static SquareMatrix GetDurbinMatrixPrime (int k, double h) {

            // dimension of matrix
            int m = k;
            SquareMatrix DH = new SquareMatrix(m);

            // only lower row of H is non-constant, hence only lower row of H' is non-zero
            double Fi = 1.0;
            double hi = 1.0;
            for (int i = 1; i < m; i++) {
                DH[m - 1, m - i] = Fi * hi;
                Fi /= i;
                hi *= h;
            }

            return (DH);

        }

    }

}
