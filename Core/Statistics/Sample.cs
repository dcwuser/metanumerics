using System;
using System.Collections;
using System.Collections.Generic;
#if !SILVERLIGHT
using System.Data;
#endif
using System.Diagnostics;
using System.Globalization;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    // the SampleStorage class is the internal representation of a sample
    // it is used internally by the Sample, BivariateSample, and MultivariateSample classes

    internal class SampleStorage  {

        public SampleStorage () {
            data = new List<double>();
            M = 0.0; SS = 0.0;
            order = null;
        }

        private string name;

        // the actual values
        private List<double> data;

        // the mean and standard deviation
        private double M, SS;

        // the order array
        // order[0] is index of lowest value, so data[order[0]] is lowest value
        // order[1] is index of next lowest value, etc.
        // order is null is the order has not been computed
        private int[] order;

        public string Name {
            get {
                return (name);
            }
            set {
                name = value;
            }
        }

        public int Count {
            get {
                return (data.Count);
            }
        }

        public double Mean {
            get {
                return (M);
            }
        }

        public double SumOfSquaredDeviations {
            get {
                return (SS);
            }
        }

        public double Variance {
            get {
                return (SS / data.Count);
            }
        }

        public void Add (double value) {
            data.Add(value);
            int n = data.Count;
            double dM = (value - M);
            M += dM / n;
            SS += (n - 1.0) / n * dM * dM;
            order = null;
        }

        public bool Remove (double value) {
            int i = data.IndexOf(value);
            if (i < 0) {
                return (false);
            } else {
                RemoveAt(i);
                return (true);
            }
        }

        public void RemoveAt (int i) {
            double value = data[i];
            data.RemoveAt(i);
            int n = data.Count;
            if (n > 0) {
                double dM = (value - M);
                M -= dM / n;
                SS -= (n + 1.0) / n * dM * dM;
            } else {
                M = 0.0;
                SS = 0.0;
            }
            order = null;
        }

        public void Clear () {
            data.Clear();
            M = 0.0;
            SS = 0.0;
            order = null;
        }

        public bool Contains (double value) {
            return (data.Contains(value));
        }

        public double this[int index] {
            get {
                return (data[index]);
            }
            set {
                data[index] = value;
            }
        }

        public void Transform (Func<double, double> transformFunction) {
            M = 0.0; SS = 0.0;
            for (int i = 0; i < data.Count; i++) {
                double value = transformFunction(data[i]);
                double dM = (value - M);
                M += dM / (i + 1);
                SS += dM * dM * i / (i + 1);
                data[i] = value;
            }
            order = null;
        }

        public int[] GetSortOrder () {
            if (order == null) {
                order = new int[data.Count];
                for (int i = 0; i < order.Length; i++) {
                    order[i] = i;
                }
                Array.Sort<int> (order, delegate (int xi, int yi) {
                    double x = data[xi];
                    double y = data[yi];
                    if (x < y) {
                        return(-1);
                    } else if (x > y) {
                        return(+1);
                    } else {
                        return(0);
                    }
                } );
                return (order);
            } else {
                return (order);
            }
        }

        public int[] GetRanks () {
            int[] order = GetSortOrder();
            int[] ranks = new int[order.Length];
            for (int i = 0; i < order.Length; i++) {
                ranks[order[i]] = i;
            }
            return (ranks);
        }

        public SampleStorage Copy () {
            SampleStorage copy = new SampleStorage();
            copy.data = new List<double>(this.data);
            copy.M = this.M;
            copy.SS = this.SS;
            copy.order = this.order;
            return (copy);
        }

        public IEnumerator<double> GetEnumerator () {
            return (data.GetEnumerator());
        }

    }

    // Summaries the count, mean, and variance of a sample.
    // These values are updatable, so that as we add values
    // we can know how the summary statistics change without
    // referencing past values. This enables "one pass" computation
    // of the mean and standard deviation.
    // See http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance and Knuth

    internal class SampleSummary {

        public SampleSummary () { }

        public SampleSummary (IEnumerable<double> values) : this() {
            Add(values);
        }

        int N;
        double M1;
        double M2;

        public int Count {
            get {
                return (N);
            }
        }

        public double Mean {
            get {
                return (M1);
            }
        }

        public double SumOfSquaredDeviations {
            get {
                return (M2);
            }
        }

        public double Variance {
            get {
                // Note this is the sample variance, not the population variance
                return (M2 / N);
            }
        }

        public void Add (double value) {
            N++;
            double dM = value - M1;
            M1 += dM / N;
            M2 += (value - M1) * dM;
            // note in M2 update, one factor is delta from old mean, the other the delta from new mean
            // M2 can also be updated as dM * dM * (n-1) / n, but this requires an additional flop
        }

        public void Add (IEnumerable<double> values) {
            if (values == null) throw new ArgumentNullException("values");
            foreach (double value in values) {
                Add(value);
            }
        }

        public void Remove (double value) {
            if (N == 0) {
                throw new InvalidOperationException();
            } else if (N == 1) {
                N = 0;
                M1 = 0.0;
                M2 = 0.0;
            } else {
                N--;
                double dM = value - M1;
                M1 -= dM / N;
                M2 -= (value - M1) * dM;
            }
        }

        public void Clear () {
            N = 0;
            M1 = 0.0;
            M2 = 0.0;
        }

    }


    /// <summary>
    /// Represents a set of data points, where each data point consists of a single real number.
    /// </summary>
    /// <remarks>
    /// <para>A univariate sample is a data set which records one number for each independent
    /// observation. For example, data from a study which measured the weight of each subject could be
    /// stored in the Sample class. The class offers descriptive statistics for the sample, estimates
    /// of descriptive statistics of the underlying population distribution, and statistical
    /// tests to compare the sample distribution to other sample distributions or theoretical models.</para>
    /// </remarks>
    public sealed class Sample : ICollection<double>, IEnumerable<double>, IEnumerable {

        private SampleStorage data;

        // isReadOnly is true for samples returned as columns of larger data sets
        // values cannot be added to or deleted from such samples because that would
        // distrub the correspondence with other columns
        private bool isReadOnly;

        /// <summary>
        /// Initializes a new, empty sample.
        /// </summary>
        public Sample () {
            this.data = new SampleStorage();
            this.isReadOnly = false;
        }

        /// <summary>
        /// Initializes a new, empty sample with the given name.
        /// </summary>
        /// <param name="name">The name of the sample.</param>
        public Sample (string name) : this () {
            this.data.Name = name;
        }

        /// <summary>
        /// Initializes a new sample from a list of values.
        /// </summary>
        /// <param name="values">Values to add to the sample.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        public Sample (IEnumerable<double> values) : this() {
            Add(values);
        }

        /// <summary>
        /// Initializes a new sample from a list of values.
        /// </summary>
        /// <param name="values">Values to add to the sample.</param>
        public Sample (params double[] values) : this() {
            Add((IEnumerable<double>)values);
        }

        internal Sample (SampleStorage data, bool isReadOnly) {
            this.data = data;
            this.isReadOnly = isReadOnly;
        }

        /// <summary>
        /// Gets or sets the name of the sample.
        /// </summary>
        public string Name {
            get {
                return (data.Name);
            }
            set {
                data.Name = value;
            }
        }

        /// <summary>
        /// Adds a value to the sample.
        /// </summary>
        /// <param name="value">The value to add.</param>
        public void Add (double value) {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                data.Add(value);
            }
        }

        /// <summary>
        /// Adds multiple values to the sample.
        /// </summary>
        /// <param name="values">An enumerable set of the values to add.</param>
        public void Add (IEnumerable<double> values) {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                if (values == null) throw new ArgumentNullException("values");
                foreach (double value in values) {
                    data.Add(value);
                }
            }
        }

        /// <summary>
        /// Adds multiple values to the sample.
        /// </summary>
        /// <param name="values">An arbitrary number of values.</param>
        public void Add (params double[] values) {
            Add((IEnumerable<double>) values);
        }

        /// <summary>
        /// Removes a given value from the sample.
        /// </summary>
        /// <param name="value">The value to remove.</param>
        /// <returns>True if the value was found and removed, otherwise false.</returns>
        public bool Remove (double value) {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                return (data.Remove(value));
            }
        }

        /// <summary>
        /// Remove all values from the sample.
        /// </summary>
        public void Clear () {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                data.Clear();
            }
        }

        /// <summary>
        /// Transforms all values using a user-supplied function.
        /// </summary>
        /// <param name="transformFunction">The function used to transform the values, which must not be null.</param>
        /// <remarks>
        /// <para>For example, to replace all values with their logarithms, apply a transform using <see cref="Math.Log"/>.</para>
        /// <para>If the supplied transform function throws an excaption, or returns infinite or NaN values, the transformation
        /// may be incomplete or the data corrupted.</para>
        /// </remarks>
        public void Transform (Func<double, double> transformFunction) {
            if (transformFunction == null) throw new ArgumentNullException("transformFunction");
            data.Transform(transformFunction);
            // it's okay to apply a transform even to a "read-only" sample because read-only is there to
            // keep us from adding or removing entries, not from updating them
        }

        /// <summary>
        /// Gets a value indicating whether the sample is read-only.
        /// </summary>
        public bool IsReadOnly {
            get {
                return (isReadOnly);
            }
        }

        /// <summary>
        /// Determines whether the sample contains the given value.
        /// </summary>
        /// <param name="value">The value to check for.</param>
        /// <returns>True if the sample contains <paramref name="value"/>, otherwise false.</returns>
        public bool Contains (double value) {
            return (data.Contains(value));
        }

        /// <summary>
        /// Copies the sample.
        /// </summary>
        /// <returns>An independent copy of the sample.</returns>
        public Sample Copy () {
            return (new Sample(data.Copy(), false));
        }

        /// <summary>
        /// Gets the number of values in the sample.
        /// </summary>
        public int Count {
            get {
                return (data.Count);
            }
        }

        /// <summary>
        /// Gets the sample mean.
        /// </summary>
        /// <remarks>
        /// <para>The mean is the average of all values in the sample.</para>
        /// </remarks>
        /// <seealso cref="PopulationMean"/>
        public double Mean {
            get {
                return (data.Mean);
            }
        }

        // in some cases (e.g. t-tests) its helpful to access the sum of squared deviations directly
        // without this property, we would need to multiply by N or N-1 just to get back to the 
        // number we divided by N or N-1 in order to return the sample or population variance

        internal double SumOfSquareDeviations {
            get {
                return (data.SumOfSquaredDeviations);
            }
        }


        /// <summary>
        /// Gets the sample variance.
        /// </summary>
        /// <remarks>
        /// <para>Note this is the actual variance of the sample values, not the infered variance
        /// of the underlying population; to obtain the latter use <see cref="PopulationVariance" />.</para>
        /// </remarks>
        public double Variance {
            get {
                return (data.Variance);
            }
        }

        /// <summary>
        /// Gets the sample standard deviation.
        /// </summary>
        /// <remarks>
        /// <para>Note this is the actual standard deviation of the sample values, not the infered standard
        /// deviation of the underlying population; to obtain the latter use
        /// <see cref="PopulationStandardDeviation" />.</para>
        /// </remarks>        
        public double StandardDeviation {
            get {
                return (Math.Sqrt(Variance));
            }
        }

        /// <summary>
        /// Computes the given sample moment.
        /// </summary>
        /// <param name="n">The order of the moment to compute.</param>
        /// <returns>The <paramref name="n"/>th moment of the sample.</returns>
        public double Moment (int n) {
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else if (n == 2) {
                return (Variance + MoreMath.Pow2(Mean));
            } else {
                double M = 0.0;
                for (int i = 0; i < data.Count; i++) {
                    M += MoreMath.Pow(data[i], n);
                }
                return (M / data.Count);
            }
        }

        /// <summary>
        /// Computes the given sample moment about its mean.
        /// </summary>
        /// <param name="n">The order of the moment to compute.</param>
        /// <returns>The <paramref name="n"/>th moment about its mean of the sample.</returns>
        public double MomentAboutMean (int n) {
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else if (n == 2) {
                return (Variance);
            } else {
                double mu = Mean;
                double C = 0.0;
                for (int i = 0; i < data.Count; i++) {
                    C += MoreMath.Pow(data[i] - mu, n);
                }
                return (C / data.Count);
            }
        }

        // median and other measures based on quantile functions

        /// <summary>
        /// Gets the sample median.
        /// </summary>
        public double Median {
            get {
                return (InverseLeftProbability(0.5));
            }
        }

        /// <summary>
        /// Gets the interquartile range of sample measurmements.
        /// </summary>
        /// <remarks>The interquartile range is the interval between the 25th and the 75th percentile.</remarks>
        /// <seealso cref="InverseLeftProbability"/>
        public Interval InterquartileRange {
            get {
                return (Interval.FromEndpoints(InverseLeftProbability(0.25), InverseLeftProbability(0.75)));
            }
        }

        /// <summary>
        /// Gets the smallest value in the sample.
        /// </summary>
        public double Minimum {
            get {
                if (data.Count < 1) throw new InsufficientDataException();
                int[] order = data.GetSortOrder();
                return (data[order[0]]);
            }
        }

        /// <summary>
        /// Gets the largest value in the sample.
        /// </summary>
        public double Maximum {
            get {
                if (data.Count < 1) throw new InsufficientDataException();
                int[] order = data.GetSortOrder();
                return (data[order[data.Count - 1]]);
            }
        }

        /// <summary>
        /// Gets the fraction of values equal to or less than the given value.
        /// </summary>
        /// <param name="value">The reference value.</param>
        /// <returns>The fraction of values in the sample that are less than or equal to the given reference value.</returns>
        public double LeftProbability (double value) {
            int[] order = data.GetSortOrder();
            for (int i = 0; i < data.Count; i++) {
                if (data[order[i]] > value) return (((double) i) / data.Count);
            }
            return (1.0);
        }

        /// <summary>
        /// Gets the sample value corresponding to a given percentile score.
        /// </summary>
        /// <param name="P">The percentile, which must lie between zero and one.</param>
        /// <returns>The corresponding value.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="P"/> lies outside [0,1].</exception>
        /// <exception cref="InsufficientDataException"><see cref="Sample.Count"/> is less than two.</exception>
        public double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            if (data.Count < 2) throw new InsufficientDataException();
            int[] order = data.GetSortOrder();
            double n = P * (data.Count - 1);
            int n1 = (int) Math.Floor(n);
            int n2 = (int) Math.Ceiling(n);
            double p1 = n2 - n;
            double p2 = 1.0 - p1;
            return (p1 * data[order[n1]] + p2 * data[order[n2]]);
        }

        // infered population-level statistics

        /// <summary>
        /// Gets an estimate of the population mean from the sample.
        /// </summary>
        public UncertainValue PopulationMean {
            get {
                // population mean is estimated by sample mean
                // population standard deviation is estimated by Sqrt(N/(N-1)) * sample standard deviation
                // standard deviation of sample mean is population standard deviation / Sqrt(N)
                if (Count < 2) throw new InsufficientDataException();
                return (new UncertainValue(Mean, Math.Sqrt(Variance / (Count - 1))));
            }
        }

        /// <summary>
        /// Gets an estimate of the population variance from the sample.
        /// </summary>
        public UncertainValue PopulationVariance {
            get {
                // the best estimate of the population variance is the sample variance C2 increased by the standard (N-1)/N correction factor
                // the uncertainty in the estimate can be computed exactly in terms of the C2 and C4 sample moments 
                int n = data.Count;
                double C2 = MomentAboutMean(2);
                double C4 = MomentAboutMean(4);
                return (new UncertainValue(C2 * n / (n - 1.0), Math.Sqrt((C4 - C2 * C2) / n)));
            }
        }

        /// <summary>
        /// Gets an estimate of the population standard deviation from the sample.
        /// </summary>
        public UncertainValue PopulationStandardDeviation {
            get {
                return (UncertainMath.Sqrt(PopulationVariance));
            }
        }

        /// <summary>
        /// Estimates the given population moment using the sample.
        /// </summary>
        /// <param name="n">The order of the moment.</param>
        /// <returns>An estimate of the <paramref name="n"/>th moment of the population.</returns>
        public UncertainValue PopulationMoment (int n) {
            if (n < 0) {
                // we don't do negative moments
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                // the zeroth moment is exactly one for any distribution
                return (new UncertainValue(1.0, 0.0));
            } else if (n == 1) {
                // the first momemt is just the mean
                return (PopulationMean);
            } else {
                // moments of order two and higher
                double M_n = Moment(n);
                double M_2n = Moment(2 * n);
                return (new UncertainValue(M_n, Math.Sqrt((M_2n - M_n * M_n) / data.Count)));
            }
        }

        /// <summary>
        /// Estimates the given population moment about the mean using the sample.
        /// </summary>
        /// <param name="n">The order of the moment.</param>
        /// <returns>An estimate, with uncertainty, of the <paramref name="n"/>th moment about the mean
        /// of the underlying population.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        public UncertainValue PopulationMomentAboutMean (int n) {
            if (n < 0) {
                // we don't do negative moments
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                // the zeroth moment is exactly one for any distribution
                return (new UncertainValue(1.0, 0.0));
            } else if (n == 1) {
                // the first moment about the mean is exactly zero by the defintion of the mean
                return (new UncertainValue(0.0, 0.0));
            } else if (n == 2) {
                return (PopulationVariance);
            } else {
                // moments of order three and higher
                // in these cases we use formulas that include only the leading and first subleading term in a 1/N expansion

                // the formula for the best estimate of the n'th population moment involves C_n, C_{n-2}, and C_2
                // the formula for the uncertainty in this estimate involves C_{2n}, C_{n+1}, C_{n-1}, and C_2
                // we need all these sample moments
                int N = data.Count;
                double C_2 = MomentAboutMean(2);
                double C_nm2 = MomentAboutMean(n - 2);
                double C_nm1 = MomentAboutMean(n - 1);
                double C_n = MomentAboutMean(n);
                double C_np1 = MomentAboutMean(n + 1);
                double C_2n = MomentAboutMean(2 * n);

                // now that we have the moments, compute the estimate and its uncertainty.
                double C = C_n - (1.0 * n / N) * (0.5 * (n - 1) * C_2 * C_nm2 - C_n);
                double V = (C_2n - 2 * n * C_nm1 * C_np1 - C_n * C_n + C_2 * C_nm1 * C_nm1 * n * n) / N;
                return (new UncertainValue(C, Math.Sqrt(V)));
            }

        }

        // statistical tests

        /// <summary>
        /// Performs a z-test.
        /// </summary>
        /// <param name="referenceMean">The mean of the comparison population.</param>
        /// <param name="referenceStandardDeviation">The standard deviation of the comparison population.</param>
        /// <returns>A test result indicating whether the sample mean is significantly different from that of the comparison population.</returns>
        /// <remarks>
        /// <para>A z-test determines whether the sample is compatible with a normal population with known mean and standard deviation.
        /// In most cases, Student's t-test (<see cref="StudentTTest(double)"/>), which does not assume a known population standard deviation,
        /// is more appropriate.</para>
        /// </remarks>
        /// <example>
        /// <para>Suppose a standardized test exists, for which it is known that the mean score is 100 and the standard deviation is 15
        /// across the entire population. The test is administered to a small sample of a subpopulation, who obtain a mean sample score of 95.
        /// You can use the z-test to determine how likely it is that the subpopulation mean really is lower than the population mean,
        /// that is that their slightly lower mean score in your sample is not merely a fluke.</para>
        /// </example>
        /// <exception cref="InsufficientDataException"><see cref="Sample.Count"/> is zero.</exception>
        /// <seealso cref="StudentTTest(double)"/>
        public TestResult ZTest (double referenceMean, double referenceStandardDeviation) {
            if (this.Count < 1) throw new InsufficientDataException();
            double z = (this.Mean - referenceMean) / (referenceStandardDeviation / Math.Sqrt(this.Count));
            return (new TestResult(z, new NormalDistribution()));

        }

        /// <summary>
        /// Tests whether the sample mean is compatible with the reference mean.
        /// </summary>
        /// <param name="referenceMean">The reference mean.</param>
        /// <returns>The result of the test. The test statistic is a t-value. If t &gt; 0, the one-sided likelyhood
        /// to obtain a greater value under the null hypothesis is the (right) propability of that value. If t &lt; 0, the
        /// corresponding one-sided likelyhood is the (left) probability of that value. The two-sided likelyhood to obtain
        /// a t-value as far or farther from zero as the value obtained is just twice the one-sided likelyhood.</returns>
        /// <remarks>
        /// <para>The test statistic of Student's t-test is the difference between the sample mean and the reference mean,
        /// measured in units of the sample mean uncertainty. For samples drawn from a normally
        /// distributed population with the given reference mean, this statistic can be shown to follow a Student distribution
        /// (<see cref="StudentDistribution"/>). If t is
        /// far from zero, then the sample is unlikely to have been drawn from a population with the given reference
        /// mean.</para>
        /// <para>Because the distribution of a t-statistic assumes a normally distributed population, this
        /// test should only be used only on sample data compatible with a normal distribution. The sign test (<see cref="SignTest"/>)
        /// is a non-parametric alternative that can be used to test the compatibility of the sample median with an assumed population median.</para>
        /// </remarks>
        /// <example>
        /// <para>In some country, the legal limit blood alcohol limit for drivers is 80 on some scale. Because they
        /// have noticed that the results given by their measuring device fluctuate, the police perform three
        /// seperate measurements on a suspected drunk driver.
        /// They obtain the results 81, 84, and 93. They argue that, because all three results exceed the
        /// limit, the court should be very confident that the driver's blood alcohol level did, in fact, exceed
        /// the legal limit. You are the driver's lawyer. Can you make an argument to that the court shouldn't be so
        /// sure?</para>
        /// <para>Here is some code that computes the probability of obtaining such high measured values,
        /// assuming that the true level is exactly 80.</para>
        /// <code lang="c#">
        /// Sample values = new Sample();
        /// values.Add(81, 84, 93);
        /// TestResult result = values.StudentTTest(80);
        /// return(result.RightProbability);
        /// </code>
        /// <para>What level of statistical confidence do you think should a court require in order to pronounce a defendant guilty?</para>
        /// </example>
        /// <exception cref="InsufficientDataException">There are fewer than two data points in the sample.</exception>
        /// <seealso cref="StudentDistribution" />
        public TestResult StudentTTest (double referenceMean) {

            // we need to be able to compute a mean and standard deviation in order to do this test; the standard deviation requires at least 2 data points
            if (this.Count < 2) throw new InsufficientDataException();
            double sigma = Math.Sqrt(this.SumOfSquareDeviations / (this.Count - 1.0));
            //double sigma = Math.Sqrt(this.Count / (this.Count - 1.0)) * this.StandardDeviation;
            double se = sigma / Math.Sqrt(this.Count);

            double t = (this.Mean - referenceMean) / se;
            int dof = this.Count - 1;

            return (new TestResult(t, new StudentDistribution(dof)));
        }

        /// <summary>
        /// Tests whether the sample median is compatible with the given reference value.
        /// </summary>
        /// <param name="referenceMedian">The reference median.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>The sign test is a non-parametric alternative to the Student t-test (<see cref="StudentTTest(double)"/>).
        /// It tests whether the sample is consistent with the given refernce median.</para>
        /// <para>The null hypothesis for the test is that the median of the underlying population from which the sample is
        /// drawn is the reference median. The test statistic is simply number of sample values that lie above the median. Since
        /// each sample value is equally likely to be below or above the population median, each draw is an independent Bernoulli
        /// trial, and the total number of values above the population median is distributed accordng to a binomial distribution
        /// (<see cref="BinomialDistribution"/>).</para>
        /// <para>The left probability of the test result is the chance of the sample median being so low, assuming the sample to have been
        /// drawn from a population with the reference median. The right probability of the test result is the chance of the sample median
        /// being so high, assuming the sample to have been drawn from a population with the reference median.</para>
        /// </remarks>
        /// <seealso cref="StudentTTest(double)"/>
        public TestResult SignTest (double referenceMedian) {

            // we need some data to do a sign test
            if (this.Count < 1) throw new InsufficientDataException();

            // count the number of entries that exceed the reference median
            int W = 0;
            foreach (double value in data) {
                if (value > referenceMedian) W++;
            }

            // W should be distributed binomially
            return (new TestResult(W, new DiscreteAsContinuousDistribution(new BinomialDistribution(0.5, data.Count))));

        }

        /// <summary>
        /// Tests whether one sample mean is compatible with another sample mean.
        /// </summary>
        /// <param name="a">The first sample, which must contain at least two entries.</param>
        /// <param name="b">The second sample, which must contain at least two entries.</param>
        /// <returns>The result of an (equal-variance) Student t-test.</returns>
        public static TestResult StudentTTest (Sample a, Sample b) {

            if (a == null) throw new ArgumentNullException("a");
            if (b == null) throw new ArgumentNullException("b");

            // get counts
            int na = a.Count;
            int nb = b.Count;

            if (na < 2) throw new InsufficientDataException();
            if (nb < 2) throw new InsufficientDataException();

            // get means
            double ma = a.Mean;
            double mb = b.Mean;

            // get pooled variance and count
            double v = (a.SumOfSquareDeviations + b.SumOfSquareDeviations) / (na + nb - 2);
            double n = 1.0 / (1.0 / na + 1.0 / nb);

            // evaluate t
            double t = (ma - mb) / Math.Sqrt(v / n);

            return (new TestResult(t, new StudentDistribution(na + nb - 2)));

        }

        /// <summary>
        /// Tests whether the sample median is compatible with the mean of another sample.
        /// </summary>
        /// <param name="a">The fisrt sample.</param>
        /// <param name="b">The second sample.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>The Mann-Whitney test is a non-parametric alternative to Student's t-test.
        /// Essentially, it supposes that the medians of the two samples are equal and tests
        /// the likelihood of this null hypothesis. Unlike the t-test, it does not assume that the sample distributions are normal.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Mann-Whitney_U_test"/>
        public static TestResult MannWhitneyTest (Sample a, Sample b) {

            if (a == null) throw new ArgumentNullException("a");
            if (b == null) throw new ArgumentNullException("b");

            // essentially, we want to order the entries from the two samples and find out how
            // many times "a's beat b's". In the ordering ababb, a beats b 5 times.

            // implementing this naively would be O(N^2), so instead we use a formula that
            // relates this quantity to the sum of ranks of a's and b's; to get those ranks
            // we do seperate sorts O(N ln N) and a merge sort O(N).

            // sort the two samples
            int[] aOrder = a.data.GetSortOrder();
            int[] bOrder = b.data.GetSortOrder();

            // now we essentially do a merge sort, but instead of actually forming the merged list,
            // we just keep track of the ranks that elements from each sample would have in the merged list

            // variables to track the sum of ranks for each sample
            int r1 = 0;
            int r2 = 0;

            // pointers to the current position in each list
            int p1 = 0;
            int p2 = 0;

            // the current rank
            int r = 1;

            while (true) {

                if (a.data[aOrder[p1]] < b.data[bOrder[p2]]) {
                    r1 += r;
                    p1++;
                    r++;
                    if (p1 >= a.Count) {
                        while (p2 < b.Count) {
                            r2 += r;
                            p2++;
                            r++;
                        }
                        break;
                    }

                } else {
                    r2 += r;
                    p2++;
                    r++;
                    if (p2 >= b.Count) {
                        while (p1 < a.Count) {
                            r1 += r;
                            p1++;
                            r++;
                        }
                        break;
                    }

                }
            }

            // relate u's to r's
            int u1 = r1 - a.Count * (a.Count + 1) / 2;
            int u2 = r2 - b.Count * (b.Count + 1) / 2;
            Debug.Assert(u1 + u2 == a.Count * b.Count);

            // return the result

            // if possible, use the exact distribution
            // to compute it, we need to do exact integer arithmetic on numbers of order the total number of possible orderings
            // since decimal is the built-in type that can hold the largest exact integers, we use it for the computation
            // therefore, to generate the exact distribution, the total number of possible orderings must be less than the capacity of a decimal 
            Distribution uDistribution;
            double lnTotal = AdvancedIntegerMath.LogFactorial(a.Count + b.Count) - AdvancedIntegerMath.LogFactorial(a.Count) - AdvancedIntegerMath.LogFactorial(b.Count);
            if (lnTotal > Math.Log((double)Decimal.MaxValue)) {
                double mu = a.Count * b.Count / 2.0;
                double sigma = Math.Sqrt(mu * (a.Count + b.Count + 1) / 6.0);
                uDistribution = new NormalDistribution(mu, sigma);
            } else {
                uDistribution = new DiscreteAsContinuousDistribution(new MannWhitneyExactDistribution(a.Count, b.Count));
            }

            return (new TestResult(u1, uDistribution));

        }

        /// <summary>
        /// Performs a one-way ANOVA.
        /// </summary>
        /// <param name="samples">The samples to compare.</param>
        /// <returns>The result of an F-test comparing the between-group variance to
        /// the within-group variance.</returns>
        /// <remarks>
        /// <para>The one-way ANOVA is an extension of the Student t-test (<see cref="StudentTTest(Sample,Sample)"/>)
        /// to more than two groups. The test's null hypothesis is that all the groups' data are drawn
        /// from the same distribution. If the null hypothesis is rejected, it indicates
        /// that at least one of the groups differs significantly from the others.</para>
        /// <para>Given more than two groups, you should use an ANOVA to test for differences
        /// in the means of the groups rather than perform multiple t-tests. The reason
        /// is that each t-test incurs a small risk of a false positive, so multiple t-tests increase
        /// the total risk of a false positive. For example, given a 95% confidence requirement,
        /// there is only a 5% chance that a t-test will incorrectly diagnose a significant
        /// difference. But given 5 samples, there are 5 * 4 /2 = 10 t-tests to be
        /// performed, giving about a 40% chance that at least one of them will incorrectly
        /// diagnose a significant difference! The ANOVA avoids the accumulation of risk
        /// by performing a single test at the required confidence level to test for
        /// any significant differences between the groups.</para>
        /// <para>A one-way ANOVA performed on just two samples is equivilent to a
        /// t-test (<see cref="Sample.StudentTTest(Sample,Sample)" />).</para>
        /// <para>ANOVA is an acronym for "Analysis of Variance". Do not be confused
        /// by the name and by the use of a ratio-of-variances test statistic: an
        /// ANOVA is primarily (although not exclusively) sensitive to changes in the
        /// <i>mean</i> between samples. The variances being compared by the test are not the
        /// variances of the individual samples; instead the test is comparing the variance of
        /// all samples considered together as one single sample to the variances of the samples
        /// considered individually. If the means of some groups differ significantly,
        /// then the variance of the unified sample will be much larger than the vairiances of the
        /// individual samples, and the test will signal a significant difference. Thus the
        /// test uses variance as a tool to detect shifts in mean, not because it
        /// is interesed in the individual sample variances per se.</para>
        /// <para>ANOVA is most appropriate when the sample data are approximately normal
        /// and the samples are distinguished by a nominal variable. For example, given
        /// a random sampling of the ages of members of five different political parties,
        /// a one-way ANOVA would be an appropriate test of the whether the different
        /// parties tend to attract different-aged memberships.</para>
        /// <para>On the other hand, given
        /// data on the incomes and vacation lengths of a large number of people,
        /// dividing the people into five income quintiles and performing a one-way ANOVA
        /// to compare the vacation day distribution of each quintile would <i>not</i> be an
        /// appropriate way to test the hypothesis that richer people take longer
        /// vacations. Since income is a cardinal variable, it would be better to
        /// in this case of put the data into a <see cref="BivariateSample"/> and
        /// perform a test of association, such as a <see cref="BivariateSample.PearsonRTest" />,
        /// <see cref="BivariateSample.SpearmanRhoTest" />, or <see cref="BivariateSample.KendallTauTest" />
        /// between the two variables. If you have measurements
        /// of additional variables for each indiviual, a <see cref="MultivariateSample.LinearRegression(int)" />
        /// analysis would allow you to adjust for confounding effects of the other variables.</para>
        /// </remarks>
        public static OneWayAnovaResult OneWayAnovaTest (params Sample[] samples) {
            return (OneWayAnovaTest((IList<Sample>)samples));
        }

        /// <summary>
        /// Performs a one-way ANOVA.
        /// </summary>
        /// <param name="samples">The samples to compare.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>For detailed information, see the variable argument overload.</para>
        /// </remarks>
        public static OneWayAnovaResult OneWayAnovaTest (IList<Sample> samples) {

            if (samples == null) throw new ArgumentNullException("samples");
            if (samples.Count < 2) throw new InvalidOperationException();

            // determine total count, mean, and within-group sum-of-squares
            int n = 0;
            double mean = 0.0;
            double SSW = 0.0;
            for (int i = 0; i < samples.Count; i++) {
                n += samples[i].Count;
                mean += samples[i].Count * samples[i].Mean;
                SSW += samples[i].Count * samples[i].Variance;
            }
            mean = mean / n;

            // determine between-group sum-of-squares
            double SSB = 0.0;
            for (int i = 0; i < samples.Count; i++) {
                SSB += samples[i].Count * MoreMath.Pow2(samples[i].Mean - mean);
            }

            // determine degrees of freedom associated with each sum-of-squares
            int dB = samples.Count - 1;
            int dW = n - 1 - dB;

            // determine F statistic
            double F = (SSB / dB) / (SSW / dW);

            AnovaRow factor = new AnovaRow(SSB, dB);
            AnovaRow residual = new AnovaRow(SSW, dW);
            AnovaRow total = new AnovaRow(SSB + SSW, n - 1);
            return (new OneWayAnovaResult(factor, residual, total));

        }

        /// <summary>
        /// Performs a Kruskal-Wallis test on the given samples.
        /// </summary>
        /// <param name="samples">The set of samples to compare.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>Kruskal-Wallis tests for differences between the samples. It is a non-parametric alternative to the
        /// one-way ANOVA (<see cref="OneWayAnovaTest(Sample[])"/>).</para>
        /// <para>The test is essentially a one-way ANOVA performed on the <i>ranks</i> of sample values instead of the sample
        /// values themselves.</para>
        /// <para>A Kruskal-Wallis test on two samples is equivilent to a Mann-Whitney test (see <see cref="MannWhitneyTest"/>).</para>
        /// </remarks>
        public static TestResult KruskalWallisTest (IList<Sample> samples) {
            if (samples == null) throw new ArgumentNullException("samples");
            if (samples.Count < 2) throw new InvalidOperationException();

            // sort each sample individually and compute count total from all samples
            int N = 0;
            int[][] orders = new int[samples.Count][];
            for (int i = 0; i < samples.Count; i++) {
                N += samples[i].data.Count;
                orders[i] = samples[i].data.GetSortOrder();
            }

            // do a multi-merge sort to determine ranks sums

            // initialize a pointer to the current active index in each ordered list
            int[] p = new int[samples.Count];

            // keep track the current rank to be assigned
            int r = 0;

            // keep track of rank sums
            // this is all that we need for KW
            // the ranks of individual entries can be made to disappear from the final formula using sum identities
            int[] rs = new int[samples.Count];

            while (true) {

                // increment the rank
                // (programmers may think ranks start from 0, but in the definition of the KW test they start from 1)
                r++;

                // determine the smallest of current entries
                int j = -1;
                double f = Double.PositiveInfinity;
                for (int k = 0; k < samples.Count; k++) {
                    if ((p[k] < orders[k].Length) && (samples[k].data[orders[k][p[k]]] < f)) {
                        j = k;
                        f = samples[k].data[orders[k][p[k]]];
                    }
                }

                // test for all lists complete
                if (j < 0) break;

                // increment the pointer and the rank sum for that column
                p[j]++;
                rs[j] += r;

            }

            // compute the KW statistic
            double H = 0.0;
            for (int i = 0; i < samples.Count; i++) {
                double z = ((double)rs[i]) / orders[i].Length - (N + 1) / 2.0;
                H += orders[i].Length * (z * z);
            }
            H = 12.0 / N / (N + 1) * H;

            // use the chi-squared approximation to the null distribution
            return (new TestResult(H, new ChiSquaredDistribution(samples.Count - 1)));

        }

        /// <summary>
        /// Performs a Kruskal-Wallis test on the given samples.
        /// </summary>
        /// <param name="samples">The set of samples to compare.</param>
        /// <returns>The result of the test.</returns>
        public static TestResult KruskalWallisTest (params Sample[] samples) {
            return (KruskalWallisTest((IList<Sample>)samples));
        }


        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="distribution">The test distribution.</param>
        /// <returns>The test result. The test statistic is the D statistic and the likelyhood is the right probability
        /// to obtain a value of D as large or larger than the one obtained.</returns>
        /// <remarks><para>The Kolmogorov-Smirnov test measures the departure of a sample from a hypothesized population
        /// distribution by comparing the cumulative probability function of the data to the cumulative probability function
        /// of the hypothesized population distribution. The test statistic D is the maximum seperation between the two curves.</para>
        /// <para>Under the null hypothesis N<sup>1/2</sup>D is known to be distributed according to the Kolomogorov distribution
        /// in the large-N limit. Because of the large-N assumption, this test should not be used with small (less than ~50)
        /// data sets.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> is null.</exception>
        /// <seealso cref="KolmogorovDistribution"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public TestResult KolmogorovSmirnovTest (Distribution distribution) {
            if (distribution == null) throw new ArgumentNullException("distribution");
            if (this.Count < 1) throw new InsufficientDataException();

            return (KolmogorovSmirnovTest(distribution, 0));
        }

        private TestResult KolmogorovSmirnovTest (Distribution distribution, int count) {

            // compute the D statistic, which is the maximum of the D+ and D- statistics
            double D, DP, DM;
            ComputeDStatistics(distribution, out DP, out DM);
            if (DP > DM) {
                D = DP;
            } else {
                D = DM;
            }

            // reduce n by the number of fit parameters, if any
            // this is heuristic and needs to be changed to deal with specific fits
            int n = data.Count - count;

            //return (new TestResult(D, new FiniteKolmogorovDistribution(n)));
            return (new TestResult(D, new KolmogorovDistribution(1.0 / Math.Sqrt(n))));
        }

        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="distribution">The test distribution.</param>
        /// <returns>The test result. The test statistic is the V statistic and the likelyhood is the right probability
        /// to obtain a value of V as large or larger than the one obtained.</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Kuiper%27s_test"/>
        public TestResult KuiperTest (Distribution distribution) {

            if (distribution == null) throw new ArgumentNullException("distribution");
            if (data.Count < 1) throw new InsufficientDataException();

            // compute the V statistic, which is the sum of the D+ and D- statistics
            double DP, DM;
            ComputeDStatistics(distribution, out DP, out DM);
            double V = DP + DM;

            return (new TestResult(V, new KuiperDistribution(1.0 / Math.Sqrt(data.Count))));

        }


        private void ComputeDStatistics (Distribution distribution, out double D1, out double D2) {

            int[] order = data.GetSortOrder();

            D1 = 0.0;
            D2 = 0.0;

            double P1 = 0.0;
            for (int i = 0; i < data.Count; i++) {

                // compute the theoretical CDF value at the data-point
                double x = data[order[i]];
                double P = distribution.LeftProbability(x);

                // compute the experimental CDF value at x+dx
                double P2 = (i + 1.0) / data.Count;

                // compare the the theoretical and experimental CDF values at x-dx and x+dx
                // update D if this is the largest seperation yet observed

                double DU = P2 - P;
                if (DU > D1) D1 = DU;
                double DL = P - P1;
                if (DL > D2) D2 = DL;

                //double D1 = Math.Abs(P - P1);
                //if (D1 > D) D = D1;
                //double D2 = Math.Abs(P - P2);
                //if (D2 > D) D = D2;

                // remember the experimental CDF value at x+dx;
                // this will be the experimental CDF value at x-dx for the next step
                P1 = P2;

            }

        }

        /// <summary>
        /// Tests whether the sample is compatible with another sample.
        /// </summary>
        /// <param name="sample">The other sample.</param>
        /// <returns>The test result. The test statistic is the D statistic and the likelyhood is the right probability
        /// to obtain a value of D as large or larger than the one obtained.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public TestResult KolmogorovSmirnovTest (Sample sample) {
            if (sample == null) throw new ArgumentNullException("sample");

            // we must have data to do the test
            if (this.Count < 1) throw new InsufficientDataException();
            if (sample.Count < 1) throw new InsufficientDataException();

            // the samples must be sorted
            int[] aOrder = this.data.GetSortOrder();
            int[] bOrder = sample.data.GetSortOrder();

            // keep track of where we are in each sample
            int c0 = 0;
            int c1 = 0;

            // keep track of the cdf of each sample
            double d0 = 0.0;
            double d1 = 0.0;

            // find the maximum cdf seperation
            double d = 0.0;
            while ((c0 < this.Count) && (c1 < sample.Count)) {
                if (this.data[aOrder[c0]] < sample.data[bOrder[c1]]) {
                    // the next data point is in sample 0
                    d0 = 1.0 * (c0 + 1) / this.Count;
                    c0++;
                } else {
                    // the next data point is in sample 1
                    d1 = 1.0 * (c1 + 1) / sample.Count;
                    c1++;
                }
                double dd = Math.Abs(d1 - d0);
                if (dd > d) d = dd;
            }

            // compute the effective degrees of freedom;
            double ne = 1.0 / (1.0 / this.Count + 1.0 / sample.Count);

            // return the result
            return (new TestResult(d, new KolmogorovDistribution(1.0 / Math.Sqrt(ne))));

        }

        /// <summary>
        /// Tests whether the variance of two samples is compatible.
        /// </summary>
        /// <param name="a">The first sample.</param>
        /// <param name="b">The second sample.</param>
        /// <returns>The result of the test.</returns>
        public static TestResult FisherFTest (Sample a, Sample b) {
            if (a == null) throw new ArgumentNullException("a");
            if (b == null) throw new ArgumentNullException("b");

            // compute population variances
            double v1 = a.Count / (a.Count - 1.0) * a.Variance;
            double v2 = b.Count / (b.Count - 1.0) * b.Variance;

            // compute the ratio
            double F = v1 / v2;

            return (new TestResult(F, new FisherDistribution(a.Count - 1, b.Count - 1)));

        }

        // enumeration

        /// <summary>
        /// Gets an enumerator of sample values.
        /// </summary>
        /// <returns>An enumerator of sample values.</returns>
        public IEnumerator<double> GetEnumerator () {
            return (data.GetEnumerator());
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (data.GetEnumerator());
        }

        void ICollection<double>.CopyTo (double[] array, int start) {
            if (array == null) throw new ArgumentNullException("array");
            // this is not fast; move implementation to SampleStorage and use CopyTo for underlying List
            for (int i = 0; i < data.Count; i++) {
                array[start + i] = data[i];
            }
        }

        /// <summary>
        /// Performs a maximum likelihood fit.
        /// </summary>
        /// <param name="distribution">The distribution to fit to the data.</param>
        /// <returns>The result of the fit, containg the parameters that result in the best fit,
        /// covariance matrix among those parameters, and a test of the goodness of fit.</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Maximum_likelihood"/>
        public FitResult MaximumLikelihoodFit (IParameterizedDistribution distribution) {

            // define a log likeyhood function
            Func<double[], double> L = delegate(double[] parameters) {
                distribution.SetParameters(parameters);
                double lnP = 0.0;
                foreach (double value in data) {
                    double P = distribution.Likelihood(value);
                    if (P == 0.0) throw new InvalidOperationException();
                    lnP += Math.Log(P);
                }
                return (-lnP);
            };

            // maximize it
            double[] v0 = distribution.GetParameters();
            SpaceExtremum min = FunctionMath.FindMinimum(L, v0);

            // turn this into a fit result
            FitResult result = new FitResult(min.Location(), min.Curvature().CholeskyDecomposition().Inverse(), null);

            return (result);

        }

#if !SILVERLIGHT
        /// <summary>
        /// Loads values from a data reader.
        /// </summary>
        /// <param name="reader">The data reader.</param>
        /// <param name="dbIndex">The column number.</param>
        public void Load (IDataReader reader, int dbIndex) {
            if (reader == null) throw new ArgumentNullException("reader");
            if (isReadOnly) throw new InvalidOperationException();
            while (reader.Read()) {
                if (reader.IsDBNull(dbIndex)) continue;
                object value = reader.GetValue(dbIndex);
                Add(Convert.ToDouble(value, CultureInfo.InvariantCulture));
            }
        }
#endif

    }

}
