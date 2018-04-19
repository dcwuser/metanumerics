using System;
using System.Diagnostics;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes a test statistic with a continuous distribution.
    /// </summary>
    public sealed class ContinuousTestStatistic {

        internal ContinuousTestStatistic (string name, double value, ContinuousDistribution distribution) {
            Debug.Assert(name != null);
            Debug.Assert(distribution != null);
            this.name = name;
            this.value = value;
            this.distribution = distribution;
        }

        private readonly string name;
        private readonly double value;
        private readonly ContinuousDistribution distribution;

        /// <summary>
        /// Gets the name of the test statistic.
        /// </summary>
        public string Name { get { return (name); } }

        /// <summary>
        /// Gets the value of the test statistic.
        /// </summary>
        public double Value { get { return (value); } }

        /// <summary>
        /// Gets the distribution of the test statistic under the null hypothesis.
        /// </summary>
        public ContinuousDistribution Distribution { get { return (distribution); } }

        internal double Probability (TestType type) {
            switch (type) {
                case TestType.LeftTailed:
                    return (distribution.LeftProbability(value));
                case TestType.RightTailed:
                    return (distribution.RightProbability(value));
                case TestType.TwoTailed:
                    // This certainly works for symmetric distributions.
                    // It also works for F-test, which is the one asymmetric case we care about
                    // for the distributions we cover. Not clear that it is true in general, though.
                    if (value < distribution.Median) {
                        return (2.0 * distribution.LeftProbability(value));
                    } else {
                        return (2.0 * distribution.RightProbability(value));
                    }
                default:
                    throw new InvalidOperationException();
            }
        }

        /// <summary>
        /// Converts a continuous test statistic to its value.
        /// </summary>
        /// <param name="statistic">The continuous test statistic.</param>
        /// <returns>The <see cref="Value"/> of the test statistic.</returns>
        public static implicit operator double (ContinuousTestStatistic statistic) {
            return (statistic.Value);
        }

        /// <inheritdoc/>
        public override string ToString () {
            return ($"{this.Name} = {this.Value}");
        }

    }

    /// <summary>
    /// Describes a test statistic with a discrete distribution.
    /// </summary>
    public sealed class DiscreteTestStatistic {

        internal DiscreteTestStatistic (string name, int value, DiscreteDistribution distribution) {
            Debug.Assert(name != null);
            Debug.Assert(distribution != null);
            this.name = name;
            this.value = value;
            this.distribution = distribution;
        }

        private readonly string name;
        private readonly int value;
        private readonly DiscreteDistribution distribution;

        /// <summary>
        /// Gets the name of the test statistic.
        /// </summary>
        public string Name { get { return (name); } }

        /// <summary>
        /// Gets the value of the test statistic.
        /// </summary>
        public int Value { get { return (value); } }

        /// <summary>
        /// Gets the distribution of the test statistic under the null hypothesis.
        /// </summary>
        public DiscreteDistribution Distribution { get { return (distribution); } }

        internal double Probability (TestType type) {
            switch (type) {
                case TestType.LeftTailed:
                    return (distribution.LeftExclusiveProbability(value) + 0.5 * distribution.ProbabilityMass(value));
                case TestType.RightTailed:
                    return (distribution.RightExclusiveProbability(value) + 0.5 * distribution.ProbabilityMass(value));
                case TestType.TwoTailed:
                    if (distribution.Skewness == 0.0) {
                        if (value < distribution.Mean) {
                            return (2.0 * (distribution.LeftExclusiveProbability(value) + 0.5 * distribution.ProbabilityMass(value)));
                        } else {
                            return (2.0 * (distribution.RightExclusiveProbability(value) + 0.5 * distribution.ProbabilityMass(value)));
                        }
                    } else {
                        throw new NotImplementedException();
                    }
                default:
                    throw new NotSupportedException();
            }
        }

        /// <inheritdoc/>
        public override string ToString () {
            return ($"{this.Name} = {this.Value}");
        }

    }
}
