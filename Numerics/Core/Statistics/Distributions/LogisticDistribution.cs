using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a logistic distribution.
    /// </summary>
    /// <remarks>
    /// <para>Like the normal distribution, the logistic distribution is a symmetric, unimodal distribution
    /// distribution with exponentially supressed tails.</para>
    /// <para>A logistic distribution with mean zero and standard deviation one is called a standard logistic distribution. Any logistic distribution
    /// can be converted to a standard logistic distribution by reparameterzing into z = (x-m)/s.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Logistic_distribution" />
    public sealed class LogisticDistribution : Distribution, IParameterizedDistribution {

        /// <summary>
        /// Initializes a new standard logistic distribution.
        /// </summary>
        /// <remarks>
        /// <para>The standard logistic distribution has zero mean and a unit width parameter.</para>
        /// </remarks>
        public LogisticDistribution () {
            m = 0.0;
            s = 1.0;
        }

        /// <summary>
        /// Initializes a new logistic distribution with the given mean and width parameters.
        /// </summary>
        /// <param name="m">The mean.</param>
        /// <param name="s">The width parameter.</param>
        public LogisticDistribution (double m, double s) {
            if (s <= 0.0) throw new ArgumentOutOfRangeException("s");
            this.m = m;
            this.s = s;
        }

        // mean and width parameter

        private double m, s;

        /// <inheritdoc />
        public override double Mean {
            get {
                return (m);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (m);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.PI * s / Global.SqrtThree);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (Math.PI * Math.PI * s * s / 3.0);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            double z = (x - m) / s;
            double cosh = Math.Cosh(z / 2.0);
            return (1.0 / (4.0 * s * cosh * cosh));
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double z = (x - m) / s;
            return (1.0 / (1.0 + Math.Exp(-z)));
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double z = (x - m) / s;
            return (1.0 / (1.0 + Math.Exp(z)));
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (m + s * Math.Log(P / (1.0 - P)));
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else {
                return (RawMomentFromCentralMoments(n));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                if (n % 2 == 0) {
                    return (2.0 * AdvancedIntegerMath.Factorial(n) * AdvancedMath.DirichletEta(n) * Math.Pow(s, n));
                } else {
                    return (0.0);
                }
            }
        }

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { m, s });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 2) throw new DimensionMismatchException();
            if (parameters[1] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            m = parameters[0];
            s = parameters[1];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }

    }

}