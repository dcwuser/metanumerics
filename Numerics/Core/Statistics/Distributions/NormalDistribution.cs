using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics {


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
    public class NormalDistribution : Distribution, IParameterizedDistribution {

        private double mu, sigma;

        /// <summary>
        /// Initializes a new normal distribution with the given mean and standard deviation.
        /// </summary>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The standard deviation, which must be positive.</param>
        public NormalDistribution (double mu, double sigma) {
            if (sigma <= 0.0) throw new ArgumentOutOfRangeException("sigma");
            this.mu = mu;
            this.sigma = sigma;
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
            /*
            if (z < 0.0) {
                return (0.5 * (1.0 + AdvancedMath.Erf(-z / Global.SqrtTwo)));
            } else {
                return (0.5 * AdvancedMath.Erfc(z / Global.SqrtTwo));
            }
            */
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
        public override double Moment (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else {
                // use an expansion in powers of the mean and moments about the mean
                // i think this still has problems
                //double mu = Mean;
                double M = 0.0;
                for (int k = 0; k <= n; k = k + 2) {
                    M += AdvancedIntegerMath.BinomialCoefficient(n, k) * Math.Pow(mu, n - k) * MomentAboutMean(k);
                }
                return (M);
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if ((n % 2) == 0) {
                // compute (n-1)!! sigma^n
                double sigma2 = sigma * sigma;
                double m = 1.0;
                for (int j = 1; j < n; j = j + 2) {
                    m *= sigma2 * j;
                }
                return (m);
            } else {
                return (0.0);
            }
        }

        internal override double Cumulant (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (0.0);
            } else if (n == 1) {
                return (mu);
            } else if (n == 2) {
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
        public override double Median {
            get {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");

            double z = Global.SqrtTwo * AdvancedMath.InverseErf(2.0 * P - 1.0);

            return (mu + sigma * z);
        }

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { mu, sigma });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 2) throw new DimensionMismatchException();
            if (parameters[1] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            mu = parameters[0];
            sigma = parameters[1];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }

    }


}