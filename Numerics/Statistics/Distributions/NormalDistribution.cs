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
    /// can be converted to a standard normal distribution by reparameterzing the data in terms of "standard deviations from the mean",
    /// i.e. z = (x - &#x3BC;) / &#x3C3;.</para>
    /// <para>Normal distribution appear in many contexts. In practical work, the normal distribution is often used as a crude
    /// model for the distribution of any continuous parameter that tends to cluster near its average, for example human height
    /// and weight. In more refined theoretical work, the normal distribution often emerges as a limiting distribution. For example,
    /// it can be shown that, if a large number of errors affect a measurement, then for nearly any underlying distribution
    /// of error terms, the distribution of total error tends to a normal distribution.</para>
    /// <para>The normal distribution is sometimes called a Gaussian distribtuion, after the mathematician Friedrich Gauss.</para>
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

        // when computing multiple central moments, it's more efficient to use recursion that to compute each seperately
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

            // Maximum likelyhood estimate is straightforward.
            //   p_i = \frac{1}{\sqrt{2\pi}\sigma} \exp \left[ -\frac{1}{2} \left( \frac{x_i - \mu}{\sigma} \right)^2 \right]
            //   \ln p_i = - \frac{1}{2} \ln (2 \pi) - \ln \sigma - \frac{1}{2} \left( \frac{x_i - \mu}{\sigma} \right)^2
            //   \ln L = \sum_i \ln p_i
            //  so
            //    \frac{\partial \ln L}{\partial \mu} =

            // These are guaranteed to be asymptotically unbiased, but they are biased for finite n. You can see that
            // must be the case because the denominator in the formula for \sigma is n instead of n-1.

            // To get unbiased formulas, we will derive exact distributions for the estimators.

            // First the mean estimator. Start from x_i \sim N(\mu, \sigma). By the addition of normal deviates,
            // \sum_i x_i \sim N(n \mu, \sqrt{n} \sigma). So
            //   m  = \frac{1}{n} \sum_i x_i \sim N(\mu, \sigma / \sqrt{n}).
            // which means the estimator m is normally distributed with mean \mu and standard deviation \sigma / \sqrt{n}.

            int n;
            double m, v;
            Univariate.ComputeMomentsUpToSecond(sample, out n, out m, out v);
            v = v / n;
            double dm = Math.Sqrt(v / n);
            /*
            int n = sample.Count;
            double m = sample.Mean;
            double v = sample.Variance;
            double dm = Math.Sqrt(v / n);
            */
            // Next the variance estimator. By the definition of the chi squared distribution and a bit of algebra that
            // reduces the degrees of freedom by one, \sum_i ( \frac{x_i - m}{\sigma} )^2 \sim \chi^2(n - 1), which has
            // mean n - 1 and variance 2(n-1). Therefore the estimator
            //   v = \frac{1}{n-1} \sum_i ( x_i - m )^2
            // has mean \sigma^2 and variance 

            // If we consider \sigma^2 the parameter, we are done. But we don't consider the parameter \sigma^2, we consider it
            // \sigma. The mean of the square root is not the square root of the mean, so the square root of an unibiased
            // estimator of \sigma^2 will not be an unbiased estimator of \sigma. If we want an unbiased estimator of \sigma
            // itself, we need to go a bit further.

            ChiDistribution d = new ChiDistribution(n - 1);
            double s = Math.Sqrt(n * v) / d.Mean;
            double ds = Math.Sqrt(d.Variance) / d.Mean * s;

            // maximum likelyhood estimates are guaranteed to be asymptotically unbiased, but not necessarily unbiased
            // this hits home for the maximum likelyhood estimate of the variance of a normal distribution, which fails
            // to include the N/(N-1) correction factor. since we know the bias, there is no reason for us not to correct
            // it, and we do so here

            //UncertainValue mu = sample.PopulationMean;
            //UncertainValue sigma = sample.PopulationStandardDeviation;
            //ContinuousDistribution distribution = new NormalDistribution(mu.Value, sigma.Value);
            //TestResult test = sample.KolmogorovSmirnovTest(distribution);
            TestResult test = sample.KolmogorovSmirnovTest(new NormalDistribution(m, s));

            // the best-fit sigma and mu are known to be uncorrelated
            // you can prove this by writing down the log likelyhood function and
            // computing its mixed second derivative, which you will see vanishes
            // at the minimum

            //return (new FitResult(mu.Value, mu.Uncertainty, sigma.Value, sigma.Uncertainty, 0.0, test));
            //return (new FitResult(m, dm, s, ds, 0.0, test));
            return (new NormalFitResult(new UncertainValue(m, dm), new UncertainValue(s, ds), test));

        }

    }

    /// <summary>
    /// Represents the result of a sample to a normal distribution.
    /// </summary>
    public sealed class NormalFitResult : DistributionFitResult<NormalDistribution> {

        internal NormalFitResult(UncertainValue mu, UncertainValue sigma, TestResult goodnessOfFit) : base(new NormalDistribution(mu.Value, sigma.Value), goodnessOfFit)
        {
            this.mu = mu;
            this.sigma = sigma;
        }

        private readonly UncertainValue mu;
        private readonly UncertainValue sigma;

        internal override ParameterCollection CreateParameters ()
        {
            ColumnVector p = new ColumnVector(mu.Value, sigma.Value);
            p.IsReadOnly = true;
            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = MoreMath.Sqr(mu.Uncertainty);
            C[1, 1] = MoreMath.Sqr(sigma.Uncertainty);
            C.IsReadOnly = true;
            ParameterCollection parameters = new ParameterCollection(
                new string[] { nameof(Mean), nameof(StandardDeviation) }, p, C
            );
            return (parameters);
        }

        /// <summary>
        /// Gets an estimate, with uncertainty, of the mean.
        /// </summary>
        public UncertainValue Mean
        {
            get
            {
                return (mu);
            }
        }

        /// <summary>
        /// Gets and estimate, with uncertainty, of the standard deviation.
        /// </summary>
        public UncertainValue StandardDeviation
        {
            get
            {
                return (sigma);
            }
        }

    }

}