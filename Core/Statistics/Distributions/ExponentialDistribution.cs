using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {


    /// <summary>
    /// Represents an exponential distribution.
    /// </summary>
    /// <remarks>
    /// <para>An exponential distribution falls off exponentially in the range from zero to infinity. It is a one-parameter
    /// distribution, determined entirely by its rate of fall-off.</para>
    /// <img src="../images/ExponentialPlot.png" />
    /// <para>The exponential distribution describes the distribution of decay times of radioactive particles.</para>
    /// <para>An exponential distribution with mean one is called a standard exponential distribution. Any exponential distribution
    /// can be converted to a standard exponential by reparameterizing the data into "fractions of the mean,"
    /// i.e. z = x / &#x3BC;.</para>
    /// <para>Processes resulting in events that are exponentially distributed in time are said to be "ageless" because the hazard function
    /// of the exponential distribution is constant. The Weibull distribution is a generalization of the exponential distribution which the
    /// hazard function changes (typically by increasing) with time.</para>
    /// </remarks>
    /// <seealso href="WeibullDistribution"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Exponential_distribution"/>
    public sealed class ExponentialDistribution : Distribution, IParameterizedDistribution {

        private double mu;

        /// <summary>
        /// Initializes a new exponential distribution with the given mean.
        /// </summary>
        /// <param name="mu">The mean, which must be positive.</param>
        public ExponentialDistribution (double mu) {
            if (mu <= 0.0) throw new ArgumentOutOfRangeException("mu");
            this.mu = mu;
        }

        /// <summary>
        /// Initializes a new standard exponential distribution.
        /// </summary>
        public ExponentialDistribution () {
            this.mu = 1.0;
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x < 0.0) {
                return (0.0);
            } else {
                return (Math.Exp(-x / mu) / mu);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x < 0.0) {
                return (0.0);
            } else {
                // 1 - e^{-x/mu}
                return (-MoreMath.ExpMinusOne(-x / mu));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < 0.0) {
                return (1.0);
            } else {
                return (Math.Exp(-x / mu));
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            // if P is very small, use the inverted series
            return (-mu * Math.Log(1.0 - P));
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                return (MoreMath.Pow(mu, n) * AdvancedIntegerMath.Factorial(n));
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
                double f = 0.0;
                double df = -1.0;
                for (int k = 2; k <= n; k++) {
                    df = -df / k;
                    f += df;
                }
                return (MoreMath.Pow(mu, n) * AdvancedIntegerMath.Factorial(n) * f);
            }
        }

        internal override double Cumulant (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (0.0);
            } else {
                return (MoreMath.Pow(mu, n) * AdvancedIntegerMath.Factorial(n - 1));
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
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (mu * Global.LogTwo);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (2.0);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { mu });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 1) throw new DimensionMismatchException();
            if (parameters[0] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            mu = parameters[0];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }

        /// <summary>
        /// Computes the exponential distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The fit result.</returns>
        /// <remarks>
        /// <para>The fit result contains one parameter, which is the mean of the best-fit exponential distribution.</para>
        /// <para>The goodness-of-fit test is a Kolmogorov-Smirnov test (<see cref="Sample.KolmogorovSmirnovTest(Distribution)"/>)
        /// which is not adjusted to reflect the fit. It is therefore unreliable for very small sample sizes.</para>
        /// </remarks>
        public static FitResult FitToSample (Sample sample) {

            if (sample == null) throw new ArgumentNullException("sample");
            if (sample.Count < 2) throw new InsufficientDataException();

            // none of the data is allowed to be negative
            foreach (double value in sample) {
                if (value < 0.0) throw new InvalidOperationException();
            }

            // the best-fit exponential's mean sample mean, with corresponding uncertainly

            double lambda = sample.Mean;
            double dLambda = lambda / Math.Sqrt(sample.Count);

            Distribution distribution = new ExponentialDistribution(lambda);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return (new FitResult(lambda, dLambda, test));

        }

    }

}