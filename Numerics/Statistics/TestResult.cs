using System;
using System.Diagnostics;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a statistical test.
    /// </summary>
    /// <remarks>
    /// <para>A statistical test compares data to a prediction, called the null hypothesis, and computes a number, 
    /// called the test statistic, which measures how much the data differs from the prediction.
    /// A key aspect of a useful statistical test is that the distribution of the test statistic under the null
    /// hypothesis can be computed. Given this distribution, one can say how likely it is for the test statistic
    /// value to take on as extreme a value was observed. If this probability is low, the null hypothesis is
    /// likely false. If this probability is not too low, the data is compatible with the null hypothesis.</para>
    /// </remarks>
    public sealed class TestResult {

        internal TestResult (string name, double value, ContinuousDistribution distribution, TestType type) {
            this.continuousStatistic = new ContinuousTestStatistic(name, value, distribution);
            this.discreteStatistic = null;
            this.type = type;
        }

        internal TestResult (string name, int value, DiscreteDistribution distribution, TestType type) {
            this.continuousStatistic = new ContinuousTestStatistic(name, value, new DiscreteAsContinuousDistribution(distribution));
            this.discreteStatistic = new DiscreteTestStatistic(name, value, distribution);
            this.type = type;
        }

        internal TestResult (
            string discreteName, int discreteValue, DiscreteDistribution discreteDistribution,
            string continuousName, double continuousValue, ContinuousDistribution continuousDistribution,
            TestType type) {
            this.discreteStatistic = new DiscreteTestStatistic(discreteName, discreteValue, discreteDistribution);
            this.continuousStatistic = new ContinuousTestStatistic(continuousName, continuousValue, continuousDistribution);
            this.type = type;
        }

        private readonly ContinuousTestStatistic continuousStatistic;
        private readonly DiscreteTestStatistic discreteStatistic;
        private TestType type;

        /// <summary>
        /// Gets or sets the sided-ness of the statistical test.
        /// </summary>
        public TestType Type { get { return (type); }  set { type = value; } }

        /// <summary>
        /// Gets the P-value of the statistical test.
        /// </summary>
        /// <value>The probability, under the test's null hypothesis,
        /// of obtaining such an extreme value of the statistic.</value>
        public double Probability {
            get {
                if (discreteStatistic != null) {
                    return (discreteStatistic.Probability(type));
                } else {
                    return (continuousStatistic.Probability(type));
                }
            }
        }

        /// <summary>
        /// Gets the test statistic.
        /// </summary>
        public ContinuousTestStatistic Statistic { get { return (continuousStatistic); } }

        /// <summary>
        /// Gets the underlying discrete test statistic.
        /// </summary>
        public DiscreteTestStatistic UnderlyingStatistic { get { return (discreteStatistic); } }

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
        /// The P-values gives the probability of obtaining a value as large in absolute value as observed.
        /// </summary>
        TwoTailed
    
    }

}
