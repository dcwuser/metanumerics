using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Meta.Numerics.Statistics {

    
    /// <summary>
    /// Tracks summary statistics for a stream of data points.
    /// </summary>
    /// <remarks>
    /// <para>Unlike the methods of the <see cref="Univariate"/> class, which by and large expect
    /// to work with sample data that is stored entirely within memory, this class is designed
    /// to compute the most important summary statistics on streamed data, i.e. individual
    /// data points which need only be presented once and may afterwards be discarded.</para>
    /// <para>The <see cref="Combine(SummaryStatistics, SummaryStatistics)"/> method allows you
    /// to combine two SummaryStatistics objects formed by processing disjoint data streams to
    /// obtain a combined SummaryStatistics object which has accurate summary statistics
    /// for the combined data set.</para>
    /// </remarks>
    public sealed class SummaryStatistics {

        /// <summary>
        /// Initializes a new summary for an empty sample.
        /// </summary>
        public SummaryStatistics () { }


        /// <summary>
        /// Initializes a new summary of the given sample values.
        /// </summary>
        /// <param name="values">An enumerator which can iterate over the values to summarize.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is <see langword="null"/>.</exception>
        public SummaryStatistics (IEnumerable<double> values) : this() {
            Add(values);
        }

        private int n = 0;
        private double m = 0.0;
        private double s2 = 0.0;
        private double s3 = 0.0;
        private double s4 = 0.0;
        private double min = Double.MaxValue;
        private double max = Double.MinValue;

        /// <summary>
        /// Adds a new value to the sample.
        /// </summary>
        /// <param name="value">The value to add.</param>
        public void Add (double value) {

            double s2_old = s2;
            double s3_old = s3;

            n++;
            double d = (value - m) / n;
            m += d;
            double e = d * d * n * (n - 1);
            s2 += e;
            s3 += d * (e * (n - 2) - 3.0 * s2_old);
            s4 += d * (d * (e * (n * n - 3 * n + 3) + 6.0 * s2_old) - 4.0 * s3_old);

            if (value < min) min = value;
            if (value > max) max = value;
        }

        /// <summary>
        /// Adds new values to the sample.
        /// </summary>
        /// <param name="values">An enumerator which can iterate over the values to add.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is <see langword="null"/>.</exception>
        public void Add (IEnumerable<double> values) {
            if (values == null) throw new ArgumentNullException(nameof(values));
            foreach (double value in values) Add(value);
        }

        /// <summary>
        /// Gets the number of values observed.
        /// </summary>
        public int Count {
            get {
                return (n);
            }
        }

        /// <summary>
        /// Gets the mean of the observed data.
        /// </summary>
        public double Mean {
            get {
                return (m);
            }
        }

        /// <summary>
        /// Get the variance of the observed data.
        /// </summary>
        public double Variance {
            get {
                return (s2 / n);
            }
        }

        /// <summary>
        /// Gets the standard deviation of the observed data.
        /// </summary>
        public double StandardDeviation {
            get {
                return (Math.Sqrt(Variance));
            }
        }

        /// <summary>
        /// Gets the skewness of the observed data.
        /// </summary>
        public double Skewness {
            get {
                return (s3 / n / Math.Pow(s2 / n, 3.0 / 2.0));
            }
        }

        /// <summary>
        /// Gets the smallest value observed.
        /// </summary>
        public double Minimum {
            get {
                return (min);
            }
        }

        /// <summary>
        /// Gets the largest value observed.
        /// </summary>
        public double Maximum {
            get {
                return (max);
            }
        }

        /// <summary>
        /// Estimates the mean of the underlying population.
        /// </summary>
        /// <value>As estimate, with uncertainty, of the mean of the distribution from which the sample was drawn.</value>
        public UncertainValue PopulationMean {
            get {
                return (new UncertainValue(m, Math.Sqrt(s2 / n / (n - 1))));
            }
        }

        /// <summary>
        /// Estimates the variance of the underlying population.
        /// </summary>
        /// <value>As estimate, with uncertainty, of the variance of the distribution from which the sample was drawn.</value>

        public UncertainValue PopulationVariance {
            get {
                return (new UncertainValue(s2 / (n - 1), Math.Sqrt((s4 / n - MoreMath.Sqr(s2 / n)) / (n - 1))));
            }
        }

        /// <summary>
        /// Estimates the standard deviation of the underlying population.
        /// </summary>
        /// <value>An estimate, with uncertainty, of the standard deviation of the
        /// distribution from which the sample was drawn.</value>
        public UncertainValue PopulationStandardDeviation {
            get {
                double c2 = s2 / n;
                double c4 = s4 / n;
                double v = (c4 - c2 * c2) / (n - 1);
                double s = Math.Sqrt(s2 / (n - 1)) * (1.0 + v / (8.0 * c2 * c2));
                double ds = 0.5 * Math.Sqrt(v) / s;

                return (new UncertainValue(s, ds));
            }
        }


        /// <summary>
        /// Combines two sample summaries.
        /// </summary>
        /// <param name="a">One sample summary.</param>
        /// <param name="b">Another sample summary.</param>
        /// <returns>Summary statistics of the combined sample.</returns>
        public static SummaryStatistics Combine (SummaryStatistics a, SummaryStatistics b) {

            double d = b.m - a.m;
            double d2 = d * d;
            double d3 = d2 * d;
            double d4 = d2 * d2;

            SummaryStatistics ab = new SummaryStatistics();
            ab.n = a.n + b.n;
            ab.m = (a.n * a.m + b.n * b.m) / ab.n;
            ab.s2 = a.s2 + b.s2 + d2 * a.n * b.n / ab.n;
            ab.s3 = a.s3 + b.s3 + (d / ab.n) * (d2 * a.n * b.n / ab.n * (a.n - b.n) + 3.0 * (a.n * b.s2 - b.n * a.s2));
            ab.s4 = a.s4 + b.s4 + (d / ab.n) * (d2 * a.n * b.n / ab.n * (d / ab.n) * (a.n * a.n - a.n * b.n + b.n * b.n) + 6.0 * (d / ab.n) * (a.n * a.n * b.s2 + b.n * b.n * a.s2) + 4.0 * (a.n * b.s3 - b.n * a.s3));

            ab.min = Math.Min(a.min, b.min);
            ab.max = Math.Max(a.max, b.max);

            return (ab);
        }

    }

}
