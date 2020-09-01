using System;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of the Kolmogorov-Smirnov D statistic.
    /// </summary>
    /// <remarks>
    /// <para>In the limit of large sample size, the D statistic of the Kolmogorov-Smirnov test (<see cref="Univariate.KolmogorovSmirnovTest(System.Collections.Generic.IReadOnlyList{double}, ContinuousDistribution)"/>)
    /// follows this distribution.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test"/>
    /// <seealse cref="Univariate.KolmogorovSmirnovTest(System.Collections.Generic.IReadOnlyList{double}, ContinuousDistribution)" />
    public sealed class KolmogorovDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new asymptotic Kolmogorov distribution.
        /// </summary>
        public KolmogorovDistribution () { }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return 0.0;
            } else if (x < 1.2) {
                return AsymptoticPPrime(x);
            } else {
                return AsymptoticQPrime(x);
            }

        }

        // the asymptotic PDF for x <~ 1

        private static double AsymptoticPPrime (double x) {

            if (x <= 0.0) return (0.0);

            double p = 0.0;
            for (int k = 1; k < Global.SeriesMax; k += 2) {
                double p_old = p;
                double z = k * Math.PI / x / 2.0;
                double dp = Math.Exp(-z * z / 2.0) * (z * z - 1.0);
                p += dp;
                if (p == p_old) return (Global.SqrtTwoPI / (x * x) * p);
            }

            throw new NonconvergenceException();
        }

        // the asymptotic PDF for x >~ 1

        private static double AsymptoticQPrime (double x) {

            double p = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double p_old = p;
                double kx = k * x;
                double dp = k * k * Math.Exp(-2.0 * kx * kx);
                if (k % 2 == 0) dp = -dp;
                p += dp;
                if (p == p_old) return (8.0 * p * x);
            }

            throw new NonconvergenceException();
        }


        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return 0.0;
            } else if (x < 1.2) {
                return AsymptoticP(x);
            } else {
                return 1.0 - AsymptoticQ(x);
            }

        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return 1.0;
            } else if (x < 1.2) {
                return 1.0 - AsymptoticP(x);
            } else {
                return AsymptoticQ(x);
            }

        }

        // implements \frac{\sqrt{2\pi}}{x1} \sum{k=0}{\infty} e^{ \frac{(2k+1)^2 \pi^2}{8 x1^2} }
        // convergence is rapid; 4 terms at x~1 and still just 10 terms at x~3

        private static double AsymptoticP (double x) {

            if (x <= 0.0) {
                return (0.0);
            } else {

                double p = 0.0;
                for (int k = 1; k < Global.SeriesMax; k += 2) {
                    double p_old = p;
                    double z = k * Math.PI / x / 2.0;
                    double dp = Math.Exp(-z * z / 2.0);
                    p += dp;
                    if (p == p_old) return (Global.SqrtTwoPI / x * p);
                }

                throw new NonconvergenceException();
            }
        }

        // implements \sum_{k=-\infty}^{\infty} (-1)^k e^{-2 k^2 x1^2}
        // convergence is very rapid; 5 terms at x~1 and just 2 terms at x~3

        private static double AsymptoticQ (double x) {
            double xx = x * x;
            double f = 0.0;
            int sign = -1;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                sign = -sign;
                double df = sign * Math.Exp(-(2 * k * k) * xx);
                f = f_old + df;
                if (f == f_old) {
                    return (2.0 * f);
                }
            }
            throw new NonconvergenceException();
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (Global.SqrtHalfPI * Global.LogTwo);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return Math.PI / 2.0 * (Math.PI / 6.0 - Global.LogTwo * Global.LogTwo);

            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                // this constant was determined empiricaly
                return (0.82757355518991);
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
            } else {
                return AdvancedMath.Gamma(r / 2.0 + 1.0) * AdvancedMath.DirichletEta(r) / Math.Pow(2.0, r / 2.0 - 1.0);
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
            } else {
                // Use integration; computation from raw moments suffers from cancelation.
                return base.CentralMoment(r);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return Interval.Semiinfinite;
            }
        }

    }

}
