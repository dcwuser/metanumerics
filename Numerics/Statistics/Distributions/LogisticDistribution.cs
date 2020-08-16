using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a logistic distribution.
    /// </summary>
    /// <remarks>
    /// <para>Like the normal distribution, the logistic distribution is a bell-shaped (symmetric, unimodal) distribution
    /// distribution with exponentially supressed tails.</para>
    /// <para>A logistic distribution with mean zero and standard deviation one is called a standard logistic distribution. Any logistic distribution
    /// can be converted to a standard logistic distribution by reparameterzing into z = (x-m)/s.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Logistic_distribution" />
    /// <seealso href="https://mathworld.wolfram.com/LogisticDistribution.html"/>
    public sealed class LogisticDistribution : ContinuousDistribution {

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
        /// <param name="s">The scale parameter.</param>
        public LogisticDistribution (double m, double s) {
            if (s <= 0.0) throw new ArgumentOutOfRangeException(nameof(s));
            this.m = m;
            this.s = s;
        }

        // mean and width parameter

        private readonly double m, s;

        /// <inheritdoc />
        public override double Mean {
            get {
                return (m);
            }
        }

        /// <summary>
        /// Gets the scale parameter for the distribution.
        /// </summary>
        public double Scale {
            get {
                return (s);
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
        public override double ExcessKurtosis {
            get {
                return (6.0 / 5.0);
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
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (InverseProbability(P, 1.0 - P));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (InverseProbability(1.0 - Q, Q));
        }

        private double InverseProbability (double P, double Q) {
            if (P == 0.0) {
                return (Double.NegativeInfinity);
            } else if (Q == 0.0) {
                return (Double.PositiveInfinity);
            } else {
                return (m + s * Math.Log(P / Q));
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
            } else if (r % 2 != 0) {
                return (0.0);
            } else  {
                return (2.0 * AdvancedIntegerMath.Factorial(r) * AdvancedMath.DirichletEta(r) * MoreMath.Pow(s, r));
            }
        }

        internal override double[] CentralMoments (int rMax) {
            // This is just a recursive formulation of the direct formula.
            double[] C = new double[rMax + 1];
            C[0] = 1.0;
            double t = 1.0; // Keep track of r! s^r
            for (int r = 2; r <= rMax; r += 2) {
                t *= r * (r - 1) * (s * s);
                C[r] = 2.0 * t * AdvancedMath.DirichletEta(r);
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
                return (m);
            } else if (r % 2 != 0) {
                return (0.0);
            } else {
                return (2.0 * AdvancedIntegerMath.Factorial(r - 1) * AdvancedMath.RiemannZeta(r) * MoreMath.Pow(s, r));
            }
        }

    }

}