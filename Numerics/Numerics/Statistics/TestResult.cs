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

        private double statistic;
        private Distribution distribution;

        internal TestResult (double statistic, Distribution distribution) {
            this.statistic = statistic;
            this.distribution = distribution;
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
        public Distribution Distribution {
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

    }


}
