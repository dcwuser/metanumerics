using System;
using System.Collections;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a multivariate sample.
    /// </summary>
    /// <remarks>
    /// <para>A multivariate sample is simply a sample in which more than one number is associated with each
    /// data point. A study which records only the height of each subject could use the <see cref="Sample"/>
    /// class to store its data, but a study which records the income and height of each subject should use
    /// MutlivariateSample class. In addition to descriptive statistics, this class offers tests for studying
    /// the associations between the recorded variables, and routines for fitting the sample to a model.
    /// </para>
    /// </remarks>
    public class MultivariateSample : ICollection<double[]>, IEnumerable<double[]>, IEnumerable {

        /// <summary>
        /// Initializes a multivariate sample.
        /// </summary>
        /// <param name="dimension">The dimension of the sample space, that is, the number of variables
        /// recorded for each sample entry.</param>
        public MultivariateSample (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException("dimension");
            n = dimension;
            means = new double[n];
        }

        private int n;
        private List<double[]> data = new List<double[]>();
        private double[] means;

        /// <summary>
        /// Gets the dimension of the sample.
        /// </summary>
        public int Dimension {
            get {
                return (n);
            }
        }

        /// <summary>
        /// Adds an entry to the sample.
        /// </summary>
        /// <param name="value">The values associated with the entry.</param>
        public void Add (params double[] value) {
            Add((IList<double>) value);
        }

        /// <summary>
        /// Adds an entry to the sample.
        /// </summary>
        /// <param name="value">The values associated with the entry.</param>
        public void Add (IList<double> value) {

            if (value.Count != n) throw new DimensionMismatchException();

            double[] datum = new double[n];
            for (int i = 0; i < n; i++) {
                datum[i] = value[i];
                means[i] += (value[i] - means[i]) / (data.Count + 1);
            }
            data.Add(datum);

        }

        /// <summary>
        /// Removes an entry from the sample.
        /// </summary>
        /// <param name="value">The values associated with the entry to remove.</param>
        /// <returns>Whether the entry was found and removed.</returns>
        public bool Remove (params double[] value) {
            return (Remove((IList<double>) value));
        }

        /// <summary>
        /// Removes an entry from the sample.
        /// </summary>
        /// <param name="value">The values associated with the entry to remove.</param>
        /// <returns>Whether the entry was found and removed.</returns>
        public bool Remove (IList<double> value) {

            if (value.Count != n) throw new DimensionMismatchException();

            for (int i = 0; i < Count; i++) {
                if (Matches(data[i], value)) {
                    // fix mean
                    data.RemoveAt(i);
                    return (true);
                }
            }

            return (false);

        }

        /// <summary>
        /// Determines whether the sample contains a given entry.
        /// </summary>
        /// <param name="value">The values associated with the entry to search for.</param>
        /// <returns>Whether the sample contains the given entry.</returns>
        public bool Contains (params double[] value) {
            return (Contains((IList<double>) value));
        }

        /// <summary>
        /// Determines whether the sample contains a given entry.
        /// </summary>
        /// <param name="value">The values associated with the entry to search for.</param>
        /// <returns>Whether the sample contains the given entry.</returns>
        public bool Contains (IList<double> value) {
            for (int i = 0; i < Count; i++) {
                if (Matches(data[i], value)) return (true);
            }
            return (false);
        }

        private bool Matches (IList<double> a, IList<double> b) {
            for (int i = 0; i < Dimension; i++) {
                if (a[i] != b[i]) return (false);
            }
            return (true);
        }

        /// <summary>
        /// Removes all entries from the sample.
        /// </summary>
        public void Clear () {
            data.Clear();
            means = new double[n];
        }

        /// <summary>
        /// Gets the number of enties in the sample.
        /// </summary>
        public int Count {
            get {
                return (data.Count);
            }
        }

        /// <summary>
        /// Computes the mean of one of the sample variables.
        /// </summary>
        /// <param name="d">The (zero-based) variable index.</param>
        /// <returns>The mean of the variable.</returns>
        public double Mean (int d) {
            if ((d < 0) || (d >= n)) throw new ArgumentOutOfRangeException("d");
            return (means[d]);
        }

        /// <summary>
        /// Computes the standard deviation of a sample variable.
        /// </summary>
        /// <param name="d">The (zero-based) index of the variable.</param>
        /// <returns>The standard deviation of the variable.</returns>
        public double StandardDeviation (int d) {
            return (Math.Sqrt(Variance(d)));
        }

        /// <summary>
        /// Computes the variance of a sample variable.
        /// </summary>
        /// <param name="d">The (zero-based) index of the variable.</param>
        /// <returns>The variance of the variable.</returns>
        public double Variance (int d) {
            if ((d < 0) || (d >= n)) throw new ArgumentOutOfRangeException("d");

            double M = Mean(d);
            double V = 0.0;
            for (int i = 0; i < Count; i++) {
                double z = data[i][d] - M;
                V += z * z;
            }
            V = V / Count;

            return (V);
        }

        /// <summary>
        /// Compute the covariance between two sample variables.
        /// </summary>
        /// <param name="d1">The (zero-based) index of the first variable.</param>
        /// <param name="d2">The (zero-based) index of the second variable.</param>
        /// <returns>The covariance of the two variables.</returns>
        public double Covariance (int d1, int d2) {
            if ((d1 < 0) || (d1 >= n)) throw new ArgumentOutOfRangeException("d1");
            if ((d2 < 0) || (d2 >= n)) throw new ArgumentOutOfRangeException("d2");

            double m1 = Mean(d1);
            double m2 = Mean(d2);

            double C = 0.0;
            foreach (double[] datum in data) {
                C += (datum[d1] - m1) * (datum[d2] - m2);
            }
            C = C / data.Count;

            return (C);
        }

        public double Moment (params int[] powers) {
            return (Moment((IList<int>) powers));

        }

        public double Moment (IList<int> powers) {
            if (powers == null) throw new ArgumentNullException("powers");
            if (powers.Count != Dimension) throw new DimensionMismatchException();

            double M = 0.0;
            for (int i = 0; i < data.Count; i++) {
                double t = 1.0;
                for (int j = 0; j < powers.Count; j++) {
                    t *= Math.Pow(data[i][j], powers[j]);
                }
                M += t;
            }
            M = M / Count;

            return (M);

        }

        public double MomentAboutMean (params int[] powers) {
            return (MomentAboutMean((IList<int>) powers));

        }

        public double MomentAboutMean (IList<int> powers) {
            if (powers == null) throw new ArgumentNullException("powers");
            if (powers.Count != Dimension) throw new DimensionMismatchException();

            // if any power is one, the moment vanishes by symmetry
            for (int i = 0; i < powers.Count; i++) {
                if (powers[i] == 1) return (0.0);
            }

            double C = 0.0;
            for (int i = 0; i < data.Count; i++) {
                double t = 1.0;
                for (int j = 0; j < powers.Count; j++) {
                    t *= Math.Pow(data[i][j] - means[j], powers[j]);
                }
                C += t;
            }
            C = C / Count;

            return (C);

        }

        /// <summary>
        /// Estimates of the population mean of one sample variable.
        /// </summary>
        /// <param name="d">The (zero-based) index of the variable.</param>
        /// <returns>An estimate, with associated uncertainty, of the population mean of the variable.</returns>
        public UncertainValue PopulationMean (int d) {
            if ((d < 0) || (d >= n)) throw new ArgumentOutOfRangeException("d");
            return (new UncertainValue(Mean(d), Math.Sqrt(Variance(d) / (Count-1))));
        }

        /// <summary>
        /// Estimates of the population covariance of two sample variable.
        /// </summary>
        /// <param name="d1">The (zero-based) index of a variable.</param>
        /// <param name="d2">The (zero-based) index of a variable.</param>
        /// <returns>An estimate, with associated uncertainty, of the population covariance of the variable.</returns>
        public UncertainValue PopulationCovariance (int d1, int d2) {

            if ((d1 < 0) || (d1 >= n)) throw new ArgumentOutOfRangeException("d1");
            if ((d2 < 0) || (d2 >= n)) throw new ArgumentOutOfRangeException("d2");

            // get the sample covariance
            double C = Covariance(d1, d2);

            // estimate the population covariance
            double PC = C / (1.0 - 1.0 / Count);

            // estimate the error in the population covariance
            // this is a large-N approximation
            int[] p4 = new int[Dimension];
            p4[d1] += 2;
            p4[d2] += 2;
            double dPC = Math.Sqrt((MomentAboutMean(p4) - C*C) / Count);

            return (new UncertainValue(PC, dPC));

        }

        /// <summary>
        /// Performs a Pearson correlation test of association between two variables.
        /// </summary>
        /// <param name="d1">The (zero-based) index of the first variable.</param>
        /// <param name="d2">The (zero-based) index of the second variable.</param>
        /// <returns></returns>
        /// <remarks>
        /// <para>This test measures the strength of the linear correlation between two variables. The
        /// test statistic r is simply the covariance of the two variables, scaled by their respective
        /// standard deviations so as to obtain a number between -1 (perfect linear anti-correlation)
        /// and +1 (perfect linear correlation).</para>
        /// <para>The Pearson test cannot reliably detect or rule out non-linear correlations.</para>
        /// <para>The Pearson correlation test requires O(N) operations.</para>
        /// </remarks>
        /// <seealso cref="SpearmanRhoTest"/>
        /// <seealso cref="KendallTauTest"/>
        public TestResult PearsonRTest (int d1, int d2) {
            double r = Covariance(d1, d2) / StandardDeviation(d1) / StandardDeviation(d2);
            Distribution p = new NormalDistribution(0.0, 1.0 / Math.Sqrt(Count));
            return (new TestResult(r, p));
        }

        /// <summary>
        /// Performs a Spearman rank-order test of association between two variables.
        /// </summary>
        /// <param name="d1">The (zero-based) index of the first variable.</param>
        /// <param name="d2">The (zero-based) index of the second variable.</param>
        /// <returns></returns>
        /// <remarks>
        /// <para>The Spearman rank-order test of association is a non-parametric test for association between
        /// two variables. The test statistic rho is the correlation coefficient of the <em>rank</em> of
        /// each entry in the sample. It is thus invariant over monotonic reparameterizations of the data,
        /// and will, for example, detect a quadratic or exponential association just as well as a linear
        /// association.</para>
        /// <para>The Spearman rank-order test requires O(N log N) operations.</para>
        /// </remarks>
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="KendallTauTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient"/>
        public TestResult SpearmanRhoTest (int d1, int d2) {

            if ((d1 < 0) || (d1 >= n)) throw new ArgumentOutOfRangeException("d1");
            if ((d2 < 0) || (d2 >= n)) throw new ArgumentOutOfRangeException("d2");

            // compute the correlated ranks
            double[] x1 = new double[Count];
            double[] x2 = new double[Count];
            for (int i = 0; i < Count; i++) {
                x1[i] = data[i][d1];
                x2[i] = data[i][d2];
            }
            Array.Sort(x1, x2);
            for (int i = 0; i < Count; i++) x1[i] = i;
            Array.Sort(x2, x1);
            for (int i = 0; i < Count; i++) x2[i] = i;

            // analytic expressions for the mean and variance of ranks
            double M = (Count - 1) / 2.0;
            double V = (Count + 1) * (Count - 1) / 12.0;

            // compute the covariance of ranks
            double C = 0.0;
            for (int i = 0; i < Count; i++) {
                C += (x1[i] - M) * (x2[i] - M);
            }
            C = C / Count;

            // compute rho
            double rho = C / V;

            return (new TestResult(rho, new NormalDistribution(0.0, 1.0 / Math.Sqrt(Count))));

        }

        /// <summary>
        /// Performs a Kendall concordance test for association between two variables.
        /// </summary>
        /// <param name="d1">The (zero-based) index of the first variable.</param>
        /// <param name="d2">The (zero-based) index of the first variable.</param>
        /// <returns></returns>
        /// <remarks>
        /// <para>Kendall's &#x3C4; is a non-parameteric and robust test of association
        /// between two variables. It simply measures the number of cases where an increase
        /// in one variable is associated with an increase in the other (corcordant pairs),
        /// compared with the number of cases where an increase in one variable is associated
        /// with a decrease in the other (discordant pairs).</para>
        /// <para>Because &#x3C4; depends only on the sign
        /// of a change and not its magnitude, it is not skewed by outliers exhibiting very large
        /// changes, nor by cases where the degree of change in one variable associated with
        /// a given change in the other changes over the range of the varibles. Of course, it may
        /// still miss an association whoose sign changes over the range of the variables. For example,
        /// if data points lie along a semi-circle in the plane, an increase in the first variable
        /// is associated with an increase in the second variable along the rising arc and and decrease in
        /// the second variable along the falling arc. No test that looks for single-signed correlation
        /// will catch this association.
        /// </para>
        /// <para>Because it examine all pairs of data points, the Kendall test requirest
        /// O(N<sup>2</sup>) operations. It is thus impractical for very large data sets. While
        /// not quite as robust as the Kendall test, the Spearman test is a good fall-back in such cases.</para>
        /// </remarks>
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="SpearmanRhoTest"/>
        public TestResult KendallTauTest (int d1, int d2) {

            if ((d1 < 0) || (d1 >= n)) throw new ArgumentOutOfRangeException("d1");
            if ((d2 < 0) || (d2 >= n)) throw new ArgumentOutOfRangeException("d2");

            // loop over all pairs, counting concordant and discordant
            int C = 0;
            int D = 0;
            for (int i = 0; i < Count; i++) {
                for (int j = 0; j < i; j++) {

                    // note the way each variable varies in the pair
                    int s1 = Math.Sign( data[i][d1] - data[j][d1] );
                    int s2 = Math.Sign( data[i][d2] - data[j][d2] );

                    // if they vary in the same way, they are concordant, otherwise they are discordant
                    if (s1 == s2) {
                        C++;
                    } else {
                        D++;
                    }
                    // note this does not count ties specially, as is sometimes done

                }
            }

            // compute tau
            double t = 1.0 * (C - D) / (C + D);

            // analytic expression for variance of tau
            double dt = Math.Sqrt( (4 * Count + 10) / 9.0 / Count / (Count - 1));

            return (new TestResult(t, new NormalDistribution(0.0, dt)));

        }

        // interface implementations

        IEnumerator<double[]> IEnumerable<double[]>.GetEnumerator () {
            return (data.GetEnumerator());
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (data.GetEnumerator());
        }

        void ICollection<double[]>.CopyTo (double[][] array, int start) {
            data.CopyTo(array, start);
        }

        bool ICollection<double[]>.IsReadOnly {
            get {
                return (false);
            }
        }

    }

}
