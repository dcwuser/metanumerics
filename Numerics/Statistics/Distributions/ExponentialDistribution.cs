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
    /// of the exponential distribution is constant. The Weibull distribution (<see cref="WeibullDistribution"/>) is a generalization
    /// of the exponential distribution which the hazard function changes (typically by increasing) with time.</para>
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
            return (-mu * MoreMath.LogOnePlus(-P));
        }

        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException("Q");
            return (-mu * Math.Log(Q));
        }

        /// <inheritdoc />
        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else {
                return (AdvancedIntegerMath.Factorial(r) * MoreMath.Pow(mu, r));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                // Subfactorial !r, see http://mathworld.wolfram.com/Subfactorial.html for properties of subfactorial.
                // Most relevent for fast computation is (!r) = round[ r! / e ].
                return (Math.Round(AdvancedIntegerMath.Factorial(r) / Math.E) * MoreMath.Pow(mu, r));
            }
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (0.0);
            } else {
                return (AdvancedIntegerMath.Factorial(r - 1) * MoreMath.Pow(mu, r));
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
        public override double ExcessKurtosis {
            get {
                return (6.0);
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
        /// <returns>The best fit parameter.</returns>
        /// <remarks>
        /// <para>The returned fit parameter is &#x3BC; (the <see cref="Mean"/>).
        /// This is the same parameter that is required by the <see cref="ExponentialDistribution(double)"/> constructor to
        /// specify a new exponential distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than two values.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        public static FitResult FitToSample (Sample sample) {

            if (sample == null) throw new ArgumentNullException("sample");
            if (sample.Count < 2) throw new InsufficientDataException();

            // none of the data is allowed to be negative
            foreach (double value in sample) {
                if (value < 0.0) throw new InvalidOperationException();
            }

            // the best-fit exponential's mean is the sample mean, with corresponding uncertainly

            double lambda = sample.Mean;
            double dLambda = lambda / Math.Sqrt(sample.Count);

            Distribution distribution = new ExponentialDistribution(lambda);
            TestResult test = sample.KolmogorovSmirnovTest(distribution);

            return (new FitResult(lambda, dLambda, test));

        }

    }

}