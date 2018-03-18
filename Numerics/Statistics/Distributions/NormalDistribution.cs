using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a normal (Gaussian) distribution.
    /// </summary>
    /// <remarks>
    /// <para>A normal distribution is a bell-shaped curve centered at its mean and falling off symmetrically on each side. It is
    /// a two-parameter distribution determined by giving its mean and standard deviation, i.e. its center and width. Its range is the
    /// entire real number line, but the tails, i.e. points more than a few standard deviations from the means, fall off extremely rapidly.</para>
    /// <img src="../images/NormalPlot.png" />
    /// <para>A normal distribution with mean zero and standard deviation one is called a standard normal distribution. Any normal distribution
    /// can be converted to a standard normal distribution by re-parameterizing the data in terms of "standard deviations from the mean",
    /// i.e. z = (x - &#x3BC;) / &#x3C3;.</para>
    /// <para>Normal distribution appear in many contexts. In practical work, the normal distribution is often used as a crude
    /// model for the distribution of any continuous parameter that tends to cluster near its average, for example human height
    /// and weight. In more refined theoretical work, the normal distribution often emerges as a limiting distribution. For example,
    /// it can be shown that, if a large number of errors affect a measurement, then for nearly any underlying distribution
    /// of error terms, the distribution of total error tends to a normal distribution.</para>
    /// <para>The normal distribution is sometimes called a Gaussian, after the mathematician Friedrich Gauss.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Normal_distribution"/>
    public sealed class NormalDistribution : ContinuousDistribution {

        private readonly double mu, sigma;

        // a standard normal generator to be used by GetRandomValue()
        private readonly IDeviateGenerator normalRng;

        /// <summary>
        /// Initializes a new normal distribution with the given mean and standard deviation.
        /// </summary>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The standard deviation, which must be positive.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="sigma"/> is less than or equal to zero.</exception>
        public NormalDistribution (double mu, double sigma) {
            if (sigma <= 0.0) throw new ArgumentOutOfRangeException(nameof(sigma));
            this.mu = mu;
            this.sigma = sigma;
            this.normalRng = DeviateGeneratorFactory.GetNormalGenerator();
        }

        /// <summary>
        /// Initializes a new standard normal distribution.
        /// </summary>
        /// <remarks>The standard normal distribution has mean zero and standard deviation one.</remarks>
        public NormalDistribution () : this(0.0, 1.0) { }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            double z = (x - mu) / sigma;
            return ((1.0 / Global.SqrtTwoPI / sigma) * Math.Exp(-z * z / 2.0));
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double z = (x - mu) / sigma;
            return (Phi(z));
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double z = (x - mu) / sigma;
            return (Phi(-z));
        }

        // standard normal left CDF, usually denoted by capital Phi
        // this is used by several other distributions, so we expose it here

        internal static double Phi (double z) {
            if (z < 0.0) {
                return (0.5 * AdvancedMath.Erfc(-z / Global.SqrtTwo));
            } else {
                return (0.5 * (1.0 + AdvancedMath.Erf(z / Global.SqrtTwo)));
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (Mean);
            } else {
                double[] C = CentralMoments(r);
                return (MomentMath.CentralToRaw(Mean, C, r));
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if ((r % 2) == 0) {
                // (r-1)!! \sigma^r
                return (AdvancedIntegerMath.DoubleFactorial(r - 1) * MoreMath.Pow(sigma, r));
            } else {
                return (0.0);
            }
        }

        // when computing multiple central moments, it's more efficient to use recursion that to compute each separately
        // this is used when we compute raw moments

        internal override double[] CentralMoments (int rMax) {
            double[] C = new double[rMax + 1];
            C[0] = 1.0;
            double sigma2 = sigma * sigma;
            for (int r = 1; r < rMax; r += 2) {
                C[r + 1] = C[r - 1] * r * sigma2;
            }
            return (C);
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (0.0);
            } else if (r == 1) {
                return (mu);
            } else if (r == 2) {
                return (sigma * sigma);
            } else {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (sigma);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            double z = AdvancedMath.Probit(P, 1.0 - P);
            return (mu + sigma * z);
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            double z = AdvancedMath.Probit(1.0 - Q, Q);
            return (mu + sigma * z);
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException(nameof(rng));
            double z = normalRng.GetNext(rng);
            return (mu + sigma * z);
        }

        /// <summary>
        /// Computes the normal distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The result of the fit.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static NormalFitResult FitToSample (IReadOnlyList<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            double m, dm, s, ds;
            FitToSampleInternal(sample, out m, out dm, out s, out ds);

            NormalDistribution distribution = new NormalDistribution(m, s);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return (new NormalFitResult(new UncertainValue(m, dm), new UncertainValue(s, ds), distribution, test));

        }

        internal static void FitToSampleInternal (IEnumerable<double> sample,
            out double m, out double dm, out double s, out double ds) {

            // We factor out this method and make it internal because it is also
            // used by the Lognormal fit method.

            // Maximum likelihood estimate is straightforward.
            //   p_i = \frac{1}{\sqrt{2\pi}\sigma} \exp \left[ -\frac{1}{2} \left( \frac{x_i - \mu}{\sigma} \right)^2 \right]
            //   \ln p_i = -\ln (\sqrt{2\pi} \sigma) - \frac{1}{2} \left( \frac{x_i - \mu}{\sigma} \right)^2
            //   \ln L = \sum_i
            //  so
            //    \frac{\partial \ln L}{\partial \mu} = \sum_i \frac{x_i - \mu}{\sigma^2}
            //    \frac{\partial \ln L}{\partial \sigma} = -\frac{n}{\sigma} - \frac{1}{\sigma^2} \sum_i (x_i - \mu)^2
            //  Setting equal to zero and solving gives the unsurprising result
            //     \mu = n^{-1} \sum_i x_i
            //     \sigma^2 = n^{-1} \sum_i (x_i - \mu)^2
            //  that MLE says to estimate the model mean and variance by the sample mean and variance.

            // MLE estimators are guaranteed to be asymptotically unbiased, but they can be biased for finite n.
            // You can see that must be the case for \sigma because the denominator has n instead of n-1.

            // To un-bias our estimators, we will derive exact distributions for these quantities.

            // First the mean estimator. Start from x_i \sim N(\mu, \sigma). By the addition of normal deviates,
            //   \sum_i x_i \sim N(n \mu, \sqrt{n} \sigma). So
            //   m  = \frac{1}{n} \sum_i x_i \sim N(\mu, \sigma / \sqrt{n}).
            // which means the estimator m is normally distributed with mean \mu and standard deviation
            // \sigma / \sqrt{n}. Now we know that m is unbiased and we know its variance.

            int n;
            double ss;
            Univariate.ComputeMomentsUpToSecond(sample, out n, out m, out ss);
            dm = Math.Sqrt(ss) / n;

            // Next the variance estimator. By the definition of the chi squared distribution and a bit of algebra that
            // reduces the degrees of freedom by one, u^2 = \sum_i ( \frac{x_i - m}{\sigma} )^2 \sim \chi^2(n - 1), which has
            // mean n - 1 and variance 2(n-1). Therefore the estimator
            //   v = \sigma^2 u^2 / (n-1) = \frac{1}{n-1} \sum_i ( x_i - m )^2
            // has mean \sigma^2 and variance 2 \sigma^4 / (n-1). 

            // If we consider \sigma^2 the parameter, we are done -- we have derived an estimator that is unbiased and
            // know its variance. But we don't consider the parameter \sigma^2, we consider it \sigma.
            // The mean of the square root is not the square root of the mean, so the square root of an unbiased
            // estimator of \sigma^2 will not be an unbiased estimator of \sigma. If we want an unbiased estimator of \sigma
            // itself, we need to go a bit further. Since u^2 ~ \chi^2(n-1), u ~ \chi(n-1). It's mean is a complicated ratio
            // of Gamma functions and it's variance is an even more complicated difference whose evaluation can be delicate,
            // but our machinery in the ChiDistribution class handles that. To get an unbiased estimator of \sigma, we just
            // need to apply the same principal of dividing by the mean of this distribution.
            //   s = \sigma u / <u> = \sqrt{\sum_i (x_i - m)^2} / <u>
            // to get an estimator with mean \sigma and known variance.

            ChiDistribution d = new ChiDistribution(n - 1);
            s = Math.Sqrt(ss) / d.Mean;
            ds = d.StandardDeviation / d.Mean * s;

        }

    }

}