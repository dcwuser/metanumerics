using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {


    /// <summary>
    /// Represents the asymptotic distribution of Kuiper's V statistic.
    /// </summary>
    public class KuiperDistribution : Distribution {

        /// <summary>
        /// Instantiates a new Kuiper distribution.
        /// </summary>
        public KuiperDistribution () {
        }

        internal KuiperDistribution (double scale) {
            this.scale = scale;
        }

        private double scale = 1.0;

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x < scale) {
                return (AsymptoticP(x / scale));
            } else {
                return (1.0 - AsymptoticQ(x / scale));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < scale) {
                return (1.0 - AsymptoticP(x / scale));
            } else {
                return (AsymptoticQ(x / scale));
            }
        }

        // series \sqrt{2\pi}{x^3} \sum_{k=1}^{\infty} k^2 \pi^2 e^{-k^2 \pi^2 / 2 x^2}
        // useful for small x

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
            if (x < scale) {
                return (AsymptoticPPrime(x / scale) / scale);
            } else {
                return (AsymptoticQPrime(x / scale) / scale);
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
                // from numerical integral; would be nice to get an analytic result
                return (1.2533141373155 * scale);
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
                return (Math.PI * Math.PI / 6.0 * scale * scale);
            } else {
                return (AdvancedMath.RiemannZeta(n) * AdvancedMath.Gamma(1 + n / 2.0) * (n - 1) / Math.Pow(2.0, n / 2.0 - 1) * Math.Pow(scale, n));
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
            } else {
                return (CentralMomentFromRawMoment(n));
            }
        }

    }


}
