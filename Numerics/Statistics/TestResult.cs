using System;
using System.Collections.Generic;
using System.Text;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents the result of a statistical test.
    /// </summary>
    /// <remarks>
    /// <para>A statistical test compares a data set to a model (or to another data set) and computes a single, real
    /// number, called the test statistic, which measures how much the data set differs from model (or the other data set).
    /// The key to a useful statistical test is that the distribution of the test statistic, under the assumption that
    /// the model actually explains the data (or that the other data set is drawn from the same distribution) is known.
    /// This assumption is called the null hypothesis.</para>
    /// </remarks>
    public class TestResult {

        private readonly double statistic;
        private readonly ContinuousDistribution distribution;

        internal TestResult (double statistic, ContinuousDistribution distribution) {
            this.statistic = statistic;
            this.distribution = distribution;
        }

        internal TestResult (string name, double statistic, TestType type, ContinuousDistribution distribution) : this(statistic, distribution) {
            this.name = name;
            this.type = type;
        }

        /// <summary>
        /// Gets the value of the test statistic.
        /// </summary>
        public double Statistic {
            get {
                return (statistic);
            }
        }

        /// <summary>
        /// Gets the distribution of the test statistic under the null hypothesis.
        /// </summary>
        public ContinuousDistribution Distribution {
            get {
                return (distribution);
            }
        }

        /// <summary>
        /// Get the probability, under the null hypothesis, of obtaining a test statistic value as small or smaller than the one actually obtained. 
        /// </summary>
        public double LeftProbability {
            get {
                return (distribution.LeftProbability(statistic));
            }
        }

        /// <summary>
        /// Get the probability, under the null hypothesis, of obtaining a test statistic value as large as or larger than the one actually obtained. 
        /// </summary>
        public double RightProbability {
            get {
                return (distribution.RightProbability(statistic));
            }
        }

        private readonly string name;

        private TestType type;

        private string Name {
            get {
                return (name);
            }
        }

        /// <summary>
        /// Gets a value indicating the type of statistical test.
        /// </summary>
        public TestType Type {
            get {
                return (type);
            }
        }

        /// <summary>
        /// Gets the probability of such an extreme value of the satistic.
        /// </summary>
        /// <value>The P-value of the test-statistic.</value>
        public double Probability {
            get {
                switch (type) {
                    case TestType.LeftTailed:
                        return (distribution.LeftProbability(statistic));
                    case TestType.RightTailed:
                        return (distribution.RightProbability(statistic));
                    case TestType.TwoTailed:
                        // This implementation is only good for distributions
                        // symmetric about the mean, so if we ever implement
                        // a two-tailed test that does not fit this requirement,
                        // we will need to revisit it.
                        if (statistic < distribution.Mean) {
                            return (2.0 * distribution.LeftProbability(statistic));
                        } else {
                            return (2.0 * distribution.RightProbability(statistic));
                        }
                    default:
                        throw new NotSupportedException();
                }
            }
        }

    }

    /// <summary>
    /// Describes the sidedness of a statistical test.
    /// </summary>
    public enum TestType {
        
        /// <summary>
        /// The P-value gives the probability of obtaining a value as low as observed.
        /// </summary>
        LeftTailed,

        /// <summary>
        /// The P-value gives the probability of obtaining a value as high as observed.
        /// </summary>
        RightTailed,

        /// <summary>
        /// The P-values gives the probabily of obtain a value as large in absolute value as observed.
        /// </summary>
        TwoTailed
    
    }

}
