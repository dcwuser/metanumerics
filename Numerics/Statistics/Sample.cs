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
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    // the SampleStorage class is the internal representation of a sample
    // it is used internally by the Sample, BivariateSample, and MultivariateSample classes

    internal class SampleStorage  {

        public SampleStorage () {
            data = new List<double>();
            //M = 0.0; SS = 0.0;
            order = null;
        }

        private string name;

        // the actual values
        private List<double> data;

        // the mean and standard deviation
        //private double M, SS;
        private SampleSummary summary = new SampleSummary();

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
                //return (M);
                return (summary.Mean);
            }
        }

        public double SumOfSquaredDeviations {
            get {
                //return (SS);
                return (summary.SumOfSquaredDeviations);
            }
        }

        public double Variance {
            get {
                //return (SS / data.Count);
                return (summary.Variance);
            }
        }

        public void Add (double value) {
            data.Add(value);
            int n = data.Count;
            //double dM = (value - M);
            //M += dM / n;
            //SS += (n - 1.0) / n * dM * dM;
            summary.Add(value);
            order = null;
        }

        public int IndexOf (double value) {
            return (data.IndexOf(value));
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
            /*
            int n = data.Count;
            if (n > 0) {
                double dM = (value - M);
                M -= dM / n;
                SS -= (n + 1.0) / n * dM * dM;
            } else {
                M = 0.0;
                SS = 0.0;
            }
            */
            summary.Remove(value);
            order = null;
        }

        public void Clear () {
            data.Clear();
            //M = 0.0;
            //SS = 0.0;
            summary.Clear();
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
                summary.Change(data[index], value);
                data[index] = value;
            }
        }

        public void CopyTo (double[] array, int startIndex) {
            data.CopyTo(array, startIndex);
        }

        public void Transform (Func<double, double> transformFunction) {
            //M = 0.0; SS = 0.0;
            summary = new SampleSummary();
            for (int i = 0; i < data.Count; i++) {
                double value = transformFunction(data[i]);
                //double dM = (value - M);
                //M += dM / (i + 1);
                //SS += dM * dM * i / (i + 1);
                summary.Add(value);
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
            //copy.M = this.M;
            //copy.SS = this.SS;
            copy.summary = this.summary;
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

    internal struct SampleSummary {

        /*
        public SampleSummary () {

        }
        */

        public SampleSummary (IEnumerable<double> values) : this() {
            Add(values);
        }

        private int N;
        private double M1;
        private double M2;

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

        public void Change (double oldValue, double newValue) {
            Remove(oldValue);
            Add(newValue);
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
        public Sample (string name) : this() {
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
            Add((IEnumerable<double>) values);
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
        /// <exception cref="InvalidOperationException">The sample is read-only.</exception>
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
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is <see langword="null"/>.</exception>
        /// <exception cref="InvalidOperationException">The sample is read-only.</exception>
        public void Add (IEnumerable<double> values) {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                if (values == null) throw new ArgumentNullException(nameof(values));
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
        /// <exception cref="InvalidOperationException">The sample is read-only.</exception>
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
        /// <exception cref="InvalidOperationException">The sample is read-only.</exception>
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
        /// <param name="transformFunction">The function used to transform the values.</param>
        /// <remarks>
        /// <para>For example, to replace all values with their logarithms, apply a transform using <see cref="Math.Log(double)"/>.</para>
        /// <para>If the supplied transform function throws an excaption, or returns infinite or NaN values, the transformation
        /// may be incomplete or the data corrupted.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="transformFunction"/> is <see langword="null"/>.</exception>
        public void Transform (Func<double, double> transformFunction) {
            if (transformFunction == null) throw new ArgumentNullException("transformFunction");
            data.Transform(transformFunction);
            // it's okay to apply a transform even to a "read-only" sample because read-only is there to
            // keep us from adding or removing entries, not from updating them
        }

        /// <summary>
        /// Gets a value indicating whether the sample is read-only.
        /// </summary>
        /// <value><see langword="true"/> if the sample is read-only, otherwise <see langword="false"/>.</value>
        /// <remarks>
        /// <para>If a sample is read-only and you need to make changes to it, you can use <see cref="Copy"/> to
        /// obtain a modifyable copy.</para>
        /// </remarks>
        public bool IsReadOnly {
            get {
                return (isReadOnly);
            }
        }

        /// <summary>
        /// Determines whether the sample contains the given value.
        /// </summary>
        /// <param name="value">The value to check for.</param>
        /// <returns>True if the sample contains <paramref name="value"/>, otherwise <see langword="false"/>.</returns>
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
        /// <para>This is the actual variance of the sample values, not the infered variance
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
        /// <para>This is the actual standard deviation of the sample values, not the infered standard
        /// deviation of the underlying population; to obtain the latter use
        /// <see cref="PopulationStandardDeviation" />.</para>
        /// </remarks>        
        public double StandardDeviation {
            get {
                return (Math.Sqrt(Variance));
            }
        }

        /// <summary>
        /// Gets the sample skewness.
        /// </summary>
        /// <remarks>
        /// <para>Skewness is the third central moment, measured in units of the appropriate power of the standard deviation.</para>
        /// </remarks>
        public double Skewness {
            get {
                return (CentralMoment(3) / Math.Pow(Variance, 3.0 / 2.0));
            }
        }

        /// <summary>
        /// Computes the given sample raw moment.
        /// </summary>
        /// <param name="r">The order of the moment to compute.</param>
        /// <returns>The <paramref name="r"/>th raw moment of the sample.</returns>
        public double RawMoment (int r) {
            if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (Mean);
            } else if (r == 2) {
                return (Variance + MoreMath.Sqr(Mean));
            } else {
                // Negative moments are okay.
                double M = 0.0;
                for (int i = 0; i < data.Count; i++) {
                    M += MoreMath.Pow(data[i], r);
                }
                return (M / data.Count);
            }
        }

        /// <summary>
        /// Computes the given sample central moment.
        /// </summary>
        /// <param name="r">The order of the moment to compute.</param>
        /// <returns>The <paramref name="r"/>th central moment.</returns>
        /// <remarks>
        /// <para>This method computes the central momements of the sample data, not the estiamted
        /// central moments of the underlying population; to obtain the latter, use <see cref="PopulationCentralMoment(int)"/>.</para>
        /// </remarks>
        public double CentralMoment (int r) {
            if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else if (r == 2) {
                return (Variance);
            } else {
                // Negative moments are okay
                double mu = Mean;
                double C = 0.0;
                for (int i = 0; i < data.Count; i++) {
                    C += MoreMath.Pow(data[i] - mu, r);
                }
                return (C / data.Count);
            }
        }

        // median and other measures based on quantile functions

        /// <summary>
        /// Gets the sample median.
        /// </summary>
        /// <seealso href="https://en.wikipedia.org/wiki/Median"/>
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
        /// <seealso href="https://en.wikipedia.org/wiki/Interquartile_range"/>
        public Interval InterquartileRange {
            get {
                return (Interval.FromEndpoints(InverseLeftProbability(0.25), InverseLeftProbability(0.75)));
            }
        }

        /// <summary>
        /// Gets the smallest value in the sample.
        /// </summary>
        /// <exception cref="InsufficientDataException">The sample contains no data.</exception>
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
        /// <exception cref="InsufficientDataException">The sample contains no data.</exception>
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
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
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
                return (EstimateFirstCumulant());
            }
        }

        /// <summary>
        /// Gets an estimate of the population variance from the sample.
        /// </summary>
        public UncertainValue PopulationVariance {
            get {
                return (EstimateSecondCumulant());
            }
        }

        /// <summary>
        /// Gets an estimate of the population standard deviation from the sample.
        /// </summary>
        /// <remarks>
        /// <para>Note that the value returned by this method is not exactly the square root
        /// of the value returned by <see cref="PopulationVariance"/>. This is not an error.</para>
        /// </remarks>
        public UncertainValue PopulationStandardDeviation {
            get {
                return (EstimateStandardDeviation());
            }
        }

        /// <summary>
        /// Gets an estimate of the population skewness from the sample.
        /// </summary>
        public UncertainValue PopulationSkewness {
            get {
                return (EstimateSkewness());
            }
        }

        /// <summary>
        /// Estimates the given population raw moment from the sample.
        /// </summary>
        /// <param name="r">The order of the moment.</param>
        /// <returns>An estimate of the <paramref name="r"/>th raw moment of the population.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> is negative.</exception>
        public UncertainValue PopulationRawMoment (int r) {
            if (r < 0) {
                // we don't do negative moments
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                // the zeroth moment is exactly one for any distribution
                return (new UncertainValue(1.0, 0.0));
            } else if (r == 1) {
                // the first momemt is just the mean
                return (EstimateFirstCumulant());
            } else {
                // moments of order two and higher
                int n = data.Count;
                double M_r = RawMoment(r);
                double M_2r = RawMoment(2 * r);
                return (new UncertainValue(M_r, Math.Sqrt((M_2r - M_r * M_r) / n)));
            }
        }

        /// <summary>
        /// Estimates the given population central moment from the sample.
        /// </summary>
        /// <param name="r">The order of the moment.</param>
        /// <returns>An estimate, with uncertainty, of the <paramref name="r"/>th moment about the mean
        /// of the underlying population.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> is negative.</exception>
        public UncertainValue PopulationCentralMoment (int r) {
            if (r < 0) {
                // we don't do negative moments
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                // the zeroth moment is exactly one for any distribution
                return (new UncertainValue(1.0, 0.0));
            } else if (r == 1) {
                // the first moment about the mean is exactly zero by the defintion of the mean
                return (new UncertainValue(0.0, 0.0));
            } else if (r == 2) {
                return (EstimateSecondCumulant());
            } else if (r == 3) {
                return (EstimateThirdCumulant());
            } else {
                return (EstimateCentralMoment(r));
            }

        }

        // First three cumulants are M1, C2, and C3. Unbiased estimators of these are called k-statistics.
        // See http://mathworld.wolfram.com/k-Statistic.html for k-statistic formulas in terms of sample central
        // moments and expressions for their variances in terms of population cumulants.

        // These formulas are also found in Dodge & Rousson, "The Complications of the Fourth Central Moment",
        // The American Statistician 53 (1999) 267. Also Stuart & Ord, "Kendall's Advanced Theory of Statistics",
        // Volume I, Chapter 12

        // The first cumulant, the mean:
        //   k_1 = m_1
        //   V(k_1) = \frac{K_2}{n}
        // To get V(k_1), we need K_2. We have an unbiased estimate k_2 = \frac{n}{n-1} c_2, and if we plug that in, we get
        //   V(k_1) = \frac{c_2}{n-1}
        // None of this is surprising. We already knew that the m_1 is an unbiased estimator of M_1, with variance C_2 / n.

        private UncertainValue EstimateFirstCumulant () {
            int n = this.Count;
            if (n < 2) throw new InsufficientDataException();
            double m1 = this.Mean;
            double c2 = this.Variance;
            return (new UncertainValue(m1, Math.Sqrt(c2 / (n - 1))));
        }

        // The second cumulant, the variance:
        //   k_2 = \frac{n}{n-1} c_2
        //   V(k_2) = \frac{K_4}{n} + \frac{2 K_2^2}{n-1}
        // Evaluating V(k_2) is not trivial. Just as we used k_2 for K_2 above, we could just use k_2 and k_4 for K_2 and K_4.
        // But that isn't exactly right, because an even though k_2 is an unbiased estimator of K_2, k_2^2 isn't an unbiased
        // estimator of K_2^2. Using polykays, Stuart & Ord derive an unbiased estimator (formula 12.131), which is also
        // quoted by Mathworld:
        //   V(k_2) = \frac{1}{n+1} \left[ \frac{(n-1)}{n} k_4 + 2 k_2^2 \right]
        // Plugging in the formulas for k-statistics in terms of central moments and doing a lot of algebra gets us to
        //   V(k_2) = \frac{n}{(n-2)(n-3)} \left[ c_4 - \frac{(n^2 - 3)}{(n - 1)^2} c_2^2 \right]
        // Numerical simulations show this to indeed be unbiased, but it has a problem. The coefficient of c_2^2 is
        // greater than unity, so the quantity in brackets can be negative. Therefore we fall back on the asymptotic
        // limit
        //   V(k_2) = \frac{c_4 - c_2^2}{n}
        // which cannot be negative because of the inequality (verified in numerical simulations) c_4 > c_2^2. 

        // See if we can do better, either by adjusting denominator and/or making exact for normal distributions.
        // Being a little off here isn't terrible, since this is just an error estimate.

        // For normal, V(k_2) = \frac{2 \sigma^4}{n-1}, so change denominator to n-1?

        private UncertainValue EstimateSecondCumulant () {
            int n = this.Count;
            if (n < 3) throw new InsufficientDataException();
            double c2 = this.CentralMoment(2);
            double c4 = this.CentralMoment(4);
            Debug.Assert(c2 >= 0.0);
            Debug.Assert(c4 >= 0.0);
            double k2 = c2 * n / (n - 1);
            double v = c4 - c2 * c2;
            Debug.Assert(v >= 0.0);
            return (new UncertainValue(k2, Math.Sqrt(v / (n - 1))));
        }

        // The third cumulant, which characterizes skew:
        //   k_3 = \frac{n^2}{(n-1)(n-2)} c_3
        //   V(k_3) = \frac{K_6}{n} + \frac{9 K_4 K_2}{n-1} + \frac{9 K_3^2}{n-1} + \frac{6 n K_2^3}{(n-1)(n-2)}    
        // Evaluating V(k_3) has the same issues as V(k_2) above. In this case, we won't even attempt an unbiased
        // estimator, but instead just go directly to the asymptotic form. In terms of central moments, this is:
        //   V(k_3) = \frac{c_6 - 6 c_4 c_2 - c_3^2 + 9 c_2^3}{n}
        // Numerical simulations suggest this is guaranteed positive, but it would be nice to find a proof.

        private UncertainValue EstimateThirdCumulant () {
            int n = this.Count;
            if (n < 4) throw new InsufficientDataException();
            double c2 = this.Variance;
            double c3 = this.CentralMoment(3);
            double c4 = this.CentralMoment(4);
            double c6 = this.CentralMoment(6);
            double k3 = c3 * n * n / (n - 1) / (n - 2);
            double v = c6 - c3 * c3 - 3.0 * (2.0 * c4 - 3.0 * c2 * c2) * c2;
            Debug.Assert(v >= 0.0);
            return (new UncertainValue(k3, Math.Sqrt(v / n)));
        }

        // Stuart & Ord, "Kendall's Advanced Theory of Statistics", Volume I, Chapter 10 derives
        // formulas for the the expectation and variance of arbitrary c_r in terms of C_r as a series in 1/n.
        //   E(c_r) = C_r - \frac{r C_r - r(r-1) / 2 C_{r-2} C_2}{n} + \cdots
        //   V(c_r) = \frac{C_{2r} - C_{r}^2 + r^2 C_{r-1}^2 C_2 - 2 r C_{r+1} C_{r-1}}{n} + \cdots
        // We can turn these around to get an estimator for C_r and its variance in terms of c_r.
        // We have verified that these agree, to these orders, with the exact expressions for M1, C2, and C3 above.

        private UncertainValue EstimateCentralMoment (int r) {

            Debug.Assert(r >= 2);

            // the formula for the best estimate of the r'th population moment involves C_r, C_{r-2}, and C_2
            // the formula for the uncertainty in this estimate involves C_{2r}, C_{r+1}, C_{r-1}, and C_2
            // we need all these sample moments
            int n = data.Count;
            if (n < 3) throw new InsufficientDataException();
            double c_2 = CentralMoment(2);
            double c_rm2 = CentralMoment(r - 2);
            double c_rm1 = CentralMoment(r - 1);
            double c_r = CentralMoment(r);
            double c_rp1 = CentralMoment(r + 1);
            double c_2r = CentralMoment(2 * r);

            double c = c_r + r * (c_r - 0.5 * (r - 1) * c_rm2 * c_2) / n;
            double v = c_2r - c_r * c_r + r * r * c_rm1 * c_rm1 * c_2 - 2 * r * c_rp1 * c_rm1;
            Debug.Assert(v >= 0.0);
            return (new UncertainValue(c, Math.Sqrt(v / n)));
        }

        // An obvious estimate of the standard deviation is just the square root of the variance, computed
        // with our usual error propagation methods. This is right as a first approximation, and is what
        // we did in previous versions, but we can do better.

        // The reason it's not exactly right is that the square root of the mean of a distributed quantity
        // is not the mean of the square root of the quantity. This means that, even though we have found
        // in k_2 an unbiased estimator of K_2 = C_2, it's square root is not an unbiased estimator of
        // \sqrt{C_2}. To discover it's bias, use Taylor expansion
        //   E(f(x)) = \int \! dp \, p(x) \, f(x)
        //           = \int \! dp \, p(x) \left[ f(x_0) + f'(x_0) (x - x_0) + \half f''(x_0) (x - x_0)^2 + \cdots \right]
        //           = f(x_0) \int \! dp \, p(x) + f'(x_0) \int \! dp \, p(x) (x - x_0) + \half f''(x_0) \int \! dp \, p(x) (x - x_0)^2 + \cdots
        //           = f(x_0) + 0 + \half f''(x_0) V(x) + \cdots
        //  In our case f(x) = \sqrt{k_2} so
        //    E(\sqrt{k_2}) = C_2^{1/2} - \frac{1}{8} C_2^{-3/2} V(k_2) + \cdots
        //  So less unbiased estimator of \sqrt{C_2} is
        //    \hat{\sigma} = k_2^{1/2} + \frac{1}{8} k_2^{-3/2} V(k_2)
        //                 = \sqrt{k_2} \left[ 1 + \frac{c_4 - c_2^2}{8 n c_2^2} 
        // This is covered in Stuart & Ord in exercise 10.20.

        private UncertainValue EstimateStandardDeviation () {
            int n = this.Count;
            if (n < 3) throw new InsufficientDataException();
            double c2 = this.CentralMoment(2);
            double c4 = this.CentralMoment(4);
            double v = (c4 - c2 * c2) / (n - 1);
            double s = Math.Sqrt(c2 * n / (n - 1)) * (1.0 + v / (8.0 * c2 * c2));
            double ds = 0.5 * Math.Sqrt(v) / s;
            return (new UncertainValue(s, ds));
        }

        private UncertainValue EstimateSkewness () {
            int n = this.Count;
            if (n < 4) throw new InsufficientDataException();

            // Compute the moments we need
            double c2 = this.CentralMoment(2);
            double c3 = this.CentralMoment(3);
            double c4 = this.CentralMoment(4);
            double c5 = this.CentralMoment(5);
            double c6 = this.CentralMoment(6);

            // Compute the vairance of c3 and c2 and c2, c3 covariance to leading order
            double v2 = (c4 - c2 * c2) / (n - 1);
            double v3 = (c6 - c3 * c3 + 9.0 * c2 * c2 * c2 - 6.0 * c4 * c2) / (n - 1);
            double v32 = (c5 - 4.0 * c3 * c2) / (n - 1);

            // Use these to compute the leading-order bias correction for skewness, and estimated variance
            double g0 = c3 * n * n / (n - 1) / (n - 2) / Math.Pow(c2 * n / (n - 1), 1.5);
            double g1 = g0 * (1.0 - 15.0 / 8.0 * v2 / (c2 * c2) + 3.0 / 2.0 * v32 / (c3 * c2));
            double vg = g0 * g0 * (v3 / (c3 * c3) - 3.0 * v32 / (c3 * c2) + 9.0 / 4.0 * v2 / (c2 * c2));
            Debug.Assert(vg >= 0.0);

            return (new UncertainValue(g1, Math.Sqrt(vg)));

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
        /// In most cases, Student's t-test (<see cref="StudentTTest(System.Double)"/>), which does not assume a known population standard deviation,
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
            return (ZTest(referenceMean, referenceStandardDeviation, TestType.TwoTailed));
        }

        /// <summary>
        /// Performs a z-test with the given sidedness.
        /// </summary>
        /// <param name="referenceMean">The mean of the comparison population.</param>
        /// <param name="referenceStandardDeviation">he standard deviation of the comparison population.</param>
        /// <param name="type">The sidedness of the test to perform.</param>
        /// <returns>A test result indicating whether the sample mean is significantly different from that of the comparison population
        /// in the direction indicated by <paramref name="type"/>.</returns>
        /// <seealso cref="ZTest(double, double)"/>
        public TestResult ZTest (double referenceMean, double referenceStandardDeviation, TestType type) {
            if (this.Count < 1) throw new InsufficientDataException();
            double z = (this.Mean - referenceMean) / (referenceStandardDeviation / Math.Sqrt(this.Count));
            return (new TestResult("z", z, type, new NormalDistribution()));
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
        /// measured in units of the sample mean uncertainty. Under the null hypothesis that the sample was drawn from a normally
        /// distributed population with the given reference mean, this statistic can be shown to follow a Student distribution
        /// (<see cref="StudentDistribution"/>). If t is far from zero, with correspondingly small left or right tail probability,
        /// then the sample is unlikely to have been drawn from a population with the given reference mean.</para>
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
        /// <seealso href="https://en.wikipedia.org/wiki/Student%27s_t-test"/>
        public TestResult StudentTTest (double referenceMean) {
            return (StudentTTest(referenceMean, TestType.TwoTailed));
        }

        /// <summary>
        /// Tests whether the sample mean differs from the reference mean in the specified direction.
        /// </summary>
        /// <param name="referenceMean">The reference mean.</param>
        /// <param name="type">The sidedness of the test to perform.</param>
        /// <returns>A test result indicating whether the sample mean is significantly different from the reference mean
        /// in the direction indicated by <paramref name="type"/>.</returns>
        public TestResult StudentTTest (double referenceMean, TestType type) {

            // We need to be able to compute a mean and standard deviation in order to do this test; the standard deviation requires at least 2 data points.

            if (this.Count < 2) throw new InsufficientDataException();
            double sigma = Math.Sqrt(this.SumOfSquareDeviations / (this.Count - 1.0));
            double se = sigma / Math.Sqrt(this.Count);

            double t = (this.Mean - referenceMean) / se;
            int dof = this.Count - 1;

            return (new TestResult("t", t, type, new StudentDistribution(dof)));
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
            return (new TestResult("W", W, TestType.TwoTailed, new DiscreteAsContinuousDistribution(new BinomialDistribution(0.5, data.Count))));

        }

        /// <summary>
        /// Tests whether one sample mean is compatible with another sample mean.
        /// </summary>
        /// <param name="a">The first sample, which must contain at least two entries.</param>
        /// <param name="b">The second sample, which must contain at least two entries.</param>
        /// <returns>The result of the test. The statistic is the Student's t and the probability
        /// is the chance of obtaining such an extreme value of t if the two samples are drawn
        /// from the same distribution.</returns>
        /// <remarks>
        /// <para>Given two samples, a back-of-the-envelope way to determine whether their means differ in a statistically
        /// significant way is to compare their <see cref="PopulationMean"/> values. If their error bars overlap, they
        /// are probably statistically compatible; if they do not, the difference in means is probably statistically
        /// significant. Student's t-test is a way to refine this back-of-the-envelope procedure into a statistical test
        /// that can determine exactly how likely a given seperation of means is under the null hypothesis that the
        /// two samples are drawn from the same distribution.</para>
        /// <para>The t-statistic is proportional to the mean of <paramref name="a"/> minus the mean of <paramref name="b"/>,
        /// so t > 0 indicates that <paramref name="a"/> has a greater mean.</para>
        /// <para>Student's t-test was one of the first statistical tests. It was described by William Sealy Gosset,
        /// a chemist who worked for the Guiness brewing company. Since Guiness was concerned that other breweries might take
        /// advantage of a technique published by one of its chemists, Gosset published his work under the pseudonym Student.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="a"/> or <paramref name="b"/> contains less than two values.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Student's_t-test"/>
        public static TestResult StudentTTest (Sample a, Sample b) {

            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

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

            return (new TestResult("t", t, TestType.TwoTailed, new StudentDistribution(na + nb - 2)));

        }

        /// <summary>
        /// Tests whether one sample median is compatible with another sample median.
        /// </summary>
        /// <param name="a">The fisrt sample.</param>
        /// <param name="b">The second sample.</param>
        /// <returns>The result of the test. The statistic is the Mann-Whitney U value and the probability
        /// is the chance of obtaining such an extreme value of U if the two samples are drawn from the
        /// same distribution.</returns>
        /// <remarks>
        /// <para>The Mann-Whitney test is a non-parametric alternative to Student's t-test (<see cref="StudentTTest(Sample, Sample)"/>).
        /// Essentially, it supposes that the medians of the two samples are equal and tests
        /// the likelihood of this null hypothesis. Unlike the t-test, it does not assume that the sample distributions are normal.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Mann-Whitney_U_test"/>
        public static TestResult MannWhitneyTest (Sample a, Sample b) {

            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            // Essentially, we want to order the entries from both samples together and find out how
            // many times "a's beat b's". In the ordering ababb, a beats b 5 times.

            // Implementing this naively would be O(N^2), so instead we use a formula that
            // relates this quantity to the sum of ranks of a's and b's; to get those ranks
            // we do seperate sorts O(N ln N) and a merge sort O(N).

            // Sort the two samples.
            int[] aOrder = a.data.GetSortOrder();
            int[] bOrder = b.data.GetSortOrder();

            // Now we essentially do a merge sort, but instead of actually forming the merged list,
            // we just keep track of the ranks that elements from each sample would have in the merged list

            // Variables to track the sum of ranks for each sample
            int r1 = 0;
            int r2 = 0;

            // Pointers to the current position in each list
            int p1 = 0;
            int p2 = 0;

            // The current rank
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

            // Relate u's to r's,
            int u1 = r1 - a.Count * (a.Count + 1) / 2;
            int u2 = r2 - b.Count * (b.Count + 1) / 2;
            Debug.Assert(u1 + u2 == a.Count * b.Count);

            // Return the result

            // If possible, we want to use the exact distribution of U.
            // To compute it, we need to do exact integer arithmetic on numbers of order the total number of possible orderings,
            // which is (a.Count + b.Count!).
            // Since decimal is the built-in type that can hold the largest exact integers, we use it for the computation.
            // Therefore, to generate the exact distribution, the total number of possible orderings must be less than the capacity of a decimal. 
            ContinuousDistribution uDistribution;
            double lnTotal = AdvancedIntegerMath.LogFactorial(a.Count + b.Count) - AdvancedIntegerMath.LogFactorial(a.Count) - AdvancedIntegerMath.LogFactorial(b.Count);
            if (lnTotal > Math.Log((double)Decimal.MaxValue)) {
                double mu = a.Count * b.Count / 2.0;
                double sigma = Math.Sqrt(mu * (a.Count + b.Count + 1) / 6.0);
                uDistribution = new NormalDistribution(mu, sigma);
            } else {
                uDistribution = new DiscreteAsContinuousDistribution(new MannWhitneyExactDistribution(a.Count, b.Count));
            }

            return (new TestResult("U", u1, TestType.TwoTailed, uDistribution));

        }

        /// <summary>
        /// Performs a one-way analysis of variance (ANOVA).
        /// </summary>
        /// <param name="samples">The samples to compare.</param>
        /// <returns>ANOVA data, including an F-test comparing the between-group variance to
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
        /// there is only a 5% chance that an individual t-test will incorrectly diagnose a significant
        /// difference. But given 5 samples, there are 5 * 4 / 2 = 10 t-tests to be
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
        /// variances of the individual samples; instead the test compares the variance of
        /// all samples considered together as one single, large sample to the variances of the samples
        /// considered individually. If the means of some groups differ significantly,
        /// then the variance of the unified sample will be much larger than the vairiances of the
        /// individual samples, and the test will signal a significant difference. Thus the
        /// test uses variance as a tool to detect shifts in mean, not because it
        /// is interesed in the individual sample variances per se.</para>
        /// <para>ANOVA is most appropriate when the sample data are continuous and approximately normal,
        /// and the samples are distinguished by a nominal variable. For example, given
        /// a random sampling of the heights of members of five different political parties,
        /// a one-way ANOVA would be an appropriate test of the whether the different
        /// parties tend to attract people with different heights.</para>
        /// <para>Given a continuous independent variable, binning in order to define
        /// groups and perform an ANOVA is generally not appropriate.
        /// For exapmple, given the incomes and heights of a large number of people,
        /// dividing these people into low-height, medium-height, and high-height groups
        /// and performing an ANOVA of the income of people in each group is not a
        /// good way to test whether height influences income.
        /// In a case like this, it would be better to put the data into a <see cref="BivariateSample"/> and
        /// perform a test of association, such as a <see cref="BivariateSample.PearsonRTest" />,
        /// <see cref="BivariateSample.SpearmanRhoTest" />, or <see cref="BivariateSample.KendallTauTest" />
        /// between the two variables. If you have measurements
        /// of additional variables for each indiviual, a <see cref="MultivariateSample.LinearRegression(int)" />
        /// analysis would allow you to adjust for the confounding effects of the other variables. If you
        /// define arbitrary bins of continuously variable data in order to form groups, then your
        /// ANOVA results will depend on your choice of bins.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="samples"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="samples"/> contains fewer than two samples.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Analysis_of_variance"/>
        /// <seealso href="https://en.wikipedia.org/wiki/One-way_analysis_of_variance"/>
        public static OneWayAnovaResult OneWayAnovaTest (params Sample[] samples) {
            return (OneWayAnovaTest((IList<Sample>)samples));
        }

        /// <summary>
        /// Performs a one-way analysis of variance (ANOVA).
        /// </summary>
        /// <param name="samples">The samples to compare.</param>
        /// <returns>ANOVA data, including an F-test comparing the between-group variance to
        /// the within-group variance.</returns>
        /// <remarks>
        /// <para>For detailed information, see <see cref="OneWayAnovaTest(Sample[])"/>.</para>
        /// </remarks>
        public static OneWayAnovaResult OneWayAnovaTest (ICollection<Sample> samples) {

            if (samples == null) throw new ArgumentNullException(nameof(samples));
            if (samples.Count < 2) throw new ArgumentException("There must be at least two samples in the sample list.", "samples");

            // determine total count, mean, and within-group sum-of-squares
            int n = 0;
            double mean = 0.0;
            double SSW = 0.0;
            foreach (Sample sample in samples) {
                if (sample == null) throw new ArgumentNullException("sample");
                n += sample.Count;
                mean += sample.Count * sample.Mean;
                SSW += sample.Count * sample.Variance;
            }
            mean = mean / n;

            // determine between-group sum-of-squares
            double SSB = 0.0;
            foreach (Sample sample in samples) {
                SSB += sample.Count * MoreMath.Sqr(sample.Mean - mean);
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
        /// Performs a two-way analysis of variance.
        /// </summary>
        /// <param name="samples">A two-dimensional array of samples, all of which must have equal counts.</param>
        /// <returns>The result of the analysis.</returns>
        /// <remarks>
        /// <para>A two-way ANOVA analyzes the effects of two seperate input factors, each with
        /// two or more nominal values, on a continuous output variable.</para>
        /// <para>The only design supported is complete and balanced: samples must exist
        /// for all combinations of treatment factors, and each of those samples must
        /// contain the same number of data points. These samples are passed to the method
        /// as a two-dimensional array <paramref name="samples"/>, whose (i, j)th entry
        /// contains the sample with ith ith row factor value and jth column factor value.</para>
        /// <para>For more information on ANOVA tests and when to use them, see the remarks for
        /// <see cref="OneWayAnovaTest(Sample[])"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="samples"/> is null, or one of its entries is null.</exception>
        /// <exception cref="InvalidOperationException">The design is not complete or balanced.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Two-way_analysis_of_variance"/>
        public static TwoWayAnovaResult TwoWayAnovaTest (Sample[,] samples) {

            if (samples == null) throw new ArgumentNullException(nameof(samples));

            int rCount = samples.GetLength(0);
            int cCount = samples.GetLength(1);
            if (rCount < 2) throw new InvalidOperationException();
            if (cCount < 2) throw new InvalidOperationException();

            // Determine mean, within group variance
            int m = -1;
            int n = 0;
            double mean = 0.0;
            double SSE = 0.0;
            for (int r = 0; r < rCount; r++) {
                for (int c = 0; c < cCount; c++) {
                    Sample sample = samples[r, c];
                    if (sample == null) throw new ArgumentNullException(String.Format("{0}[{1},{2}]", nameof(samples), r, c));
                    if (sample.Count != m) {
                        if (m < 0) {
                            m = sample.Count;
                        } else {
                            throw new InvalidOperationException();
                        }
                    }
                    n += sample.Count;
                    mean += sample.Count * sample.Mean;
                    SSE += sample.Count * sample.Variance;
                }
            }
            mean = mean / n;

            // Determine between group variance
            double SSF = 0.0;
            for (int r = 0; r < rCount; r++) {
                for (int c = 0; c < cCount; c++) {
                    SSF += MoreMath.Sqr(samples[r, c].Mean - mean);
                }
            }
            SSF = SSF * samples[0, 0].Count;

            // Determine row-wise variance
            double SSA = 0.0;
            for (int r = 0; r < rCount; r++) {
                double rMean = 0.0;
                for (int c = 0; c < cCount; c++) {
                    rMean += samples[r, c].Mean;
                }
                rMean = rMean / samples.GetLength(1);

                SSA += MoreMath.Sqr(rMean - mean);
            }
            SSA = SSA * samples[0, 0].Count * cCount;

            // Determine column-wise variance
            double SSB = 0.0;
            for (int c = 0; c < cCount; c++) {
                double cMean = 0.0;
                for (int r = 0; r < rCount; r++) {
                    cMean += samples[r, c].Mean;
                }
                cMean = cMean / rCount;

                SSB += MoreMath.Sqr(cMean - mean);
            }
            SSB = SSB * samples[0, 0].Count * rCount;

            // Finding the interaction effect by subtraction allows us to determine it
            // without storing multiple row and column means. But it introduces a risk
            // of getting a tiny negative value due to roundoff error, so we should
            // probably just go ahead and store the row and column means.
            double SSI = SSF - SSA - SSB;
            Debug.Assert(SSI >= 0.0);

            AnovaRow row = new AnovaRow(SSA, rCount - 1);
            AnovaRow column = new AnovaRow(SSB, cCount - 1);
            AnovaRow interaction = new AnovaRow(SSI, (rCount - 1) * (cCount - 1));
            AnovaRow error = new AnovaRow(SSE, n - rCount * cCount);

            TwoWayAnovaResult result = new TwoWayAnovaResult(row, column, interaction, error);
            return (result);

        }

        /// <summary>
        /// Performs a Kruskal-Wallis test on the given samples.
        /// </summary>
        /// <param name="samples">The set of samples to compare.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>For detailed information, see <see cref="KruskalWallisTest(Sample[])"/>.</para>
        /// </remarks>
        public static TestResult KruskalWallisTest (IList<Sample> samples) {
            if (samples == null) throw new ArgumentNullException(nameof(samples));
            if (samples.Count < 2) throw new ArgumentException("There must be at least two samples in the sample list.", "samples");

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
        /// <remarks>
        /// <para>Kruskal-Wallis tests for differences between the samples. It is a non-parametric alternative to the
        /// one-way ANOVA (<see cref="OneWayAnovaTest(Sample[])"/>), which is more appropriate when the data far from
        /// normally distributed.</para>
        /// <para>The test is essentially a one-way ANOVA performed on the <i>ranks</i> of sample values instead of the sample
        /// values themselves.</para>
        /// <para>A Kruskal-Wallis test on two samples is equivilent to a Mann-Whitney test (see <see cref="MannWhitneyTest(Sample, Sample)"/>).</para>
        /// <para>As with a normal ANOVA, it is not appropriate to bin a continuous independent variable in order to form
        /// groups for a Kruskal-Wallis test. Kruskal-Wallis addresses the non-normality of the dependent variable, not
        /// the non-discreteness of the independent variable.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="samples"/> is <see langword="null"/>.</exception>
        /// <see cref="OneWayAnovaTest(Sample[])"/>
        /// <see href="https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance"/>
        public static TestResult KruskalWallisTest (params Sample[] samples) {
            return (KruskalWallisTest((IList<Sample>)samples));
        }


        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="distribution">The test distribution.</param>
        /// <returns>The test result. The test statistic is the D statistic and the probability is the chance of
        /// obtaining such a large value of D under the assumption that the sample is drawn from the given distribution.</returns>
        /// <remarks>
        /// <para>The null hypothesis of the Kolmogorov-Smirnov (KS) test is that the sample is drawn from the given continuous distribution.
        /// The test statsitic D is the maximum deviation of the sample's empirical distribution function (EDF) from
        /// the distribution's cumulative distribution function (CDF). A high value of the test statistic, corresponding
        /// to a low right tail probability, indicates that the sample distribution disagrees with the given distribution
        /// to a degree unlikely to arise from statistical fluctuations.</para>
        /// <para>For small sample sizes, we compute the null distribution of D exactly. For large sample sizes, we use an accurate
        /// asympototic approximation. Therefore it is safe to use this method for all sample sizes.</para>
        /// <para>A variant of this test, <see cref="KolmogorovSmirnovTest(Sample, Sample)"/>, allows you to non-parametrically
        /// test whether two samples are drawn from the same underlying distribution, without having to specify that distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">There is no data in the sample.</exception>
        /// <seealso cref="KolmogorovDistribution"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public TestResult KolmogorovSmirnovTest (ContinuousDistribution distribution) {
            if (distribution == null) throw new ArgumentNullException(nameof(distribution));
            if (this.Count < 1) throw new InsufficientDataException();

            return (KolmogorovSmirnovTest(distribution, 0));
        }

        private TestResult KolmogorovSmirnovTest (ContinuousDistribution distribution, int count) {

            // compute the D statistic, which is the maximum of the D+ and D- statistics
            double DP, DM;
            ComputeDStatistics(distribution, out DP, out DM);
            double D = Math.Max(DP, DM);

            // reduce n by the number of fit parameters, if any
            // this is heuristic and needs to be changed to deal with specific fits
            int n = data.Count - count;
            
            ContinuousDistribution DDistribution;
            if (n < 32) {
                DDistribution = new TransformedDistribution(new KolmogorovExactDistribution(n), 0.0, 1.0 / n);
            } else {
                DDistribution = new TransformedDistribution(new KolmogorovAsymptoticDistribution(n), 0.0, 1.0 / Math.Sqrt(n));
            }
            return (new TestResult("D", D, TestType.RightTailed, DDistribution));
            
        }

        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="distribution">The test distribution.</param>
        /// <returns>The test result. The test statistic is the V statistic and the chance to obtain such a large
        /// value of V under the assumption that the sample is drawn from the given distribution.</returns>
        /// <remarks>
        /// <para>Like the Kolmogorov-Smirnov test ((<see cref="KolmogorovSmirnovTest(ContinuousDistribution)"/>,
        /// Kuiper's test compares the EDF of the sample to the CDF of the given distribution.</para>
        /// <para>For small sample sizes, we compute the null distribution of V exactly. For large sample sizes, we use an accurate
        /// asympototic approximation. Therefore it is safe to use this method for all sample sizes.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">There is no data in the sample.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Kuiper%27s_test"/>
        public TestResult KuiperTest (ContinuousDistribution distribution) {

            if (distribution == null) throw new ArgumentNullException(nameof(distribution));
            if (data.Count < 1) throw new InsufficientDataException();

            // compute the V statistic, which is the sum of the D+ and D- statistics
            double DP, DM;
            ComputeDStatistics(distribution, out DP, out DM);
            double V = DP + DM;

            ContinuousDistribution VDistribution;
            if (this.Count < 32) {
                VDistribution = new TransformedDistribution(new KuiperExactDistribution(this.Count), 0.0, 1.0 / this.Count);
            } else {
                VDistribution = new TransformedDistribution(new KuiperAsymptoticDistribution(this.Count), 0.0, 1.0 / Math.Sqrt(this.Count));
            }
            return (new TestResult("V", V, TestType.RightTailed, VDistribution));

        }


        private void ComputeDStatistics (ContinuousDistribution distribution, out double D1, out double D2) {

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
        /// <param name="a">One sample.</param>
        /// <param name="b">The other sample.</param>
        /// <returns>The test result. The test statistic is the D statistic and the likelyhood is the right probability
        /// to obtain a value of D as large or larger than the one obtained.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">One or both of the samples is empty.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public static TestResult KolmogorovSmirnovTest (Sample a, Sample b) {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            // we must have data to do the test
            if (a.Count < 1) throw new InsufficientDataException("Messages.InsufficientData");
            if (b.Count < 1) throw new InsufficientDataException("Messages.InsufficientData");

            // the samples must be sorted
            int[] aOrder = a.data.GetSortOrder();
            int[] bOrder = b.data.GetSortOrder();

            // keep track of where we are in each sample
            int ia = 0;
            int ib = 0;

            // keep track of the cdf of each sample
            double d0 = 0.0;
            double d1 = 0.0;

            // find the maximum cdf seperation
            double d = 0.0;
            while ((ia < a.Count) && (ib < b.Count)) {
                if (a.data[aOrder[ia]] < b.data[bOrder[ib]]) {
                    // the next data point is in sample 0
                    d0 = 1.0 * (ia + 1) / a.Count;
                    ia++;
                } else {
                    // the next data point is in sample 1
                    d1 = 1.0 * (ib + 1) / b.Count;
                    ib++;
                }
                double dd = Math.Abs(d1 - d0);
                if (dd > d) d = dd;
            }

            // return the result
            // currently we use the fast that, asymptotically, \sqrt{\frac{n m}{n + m}} D goes to the limiting Kolmogorov distribution
            // but this limit is known to be approached slowly; we should look into the exact distribution
            // the exact distribution is actually discrete (since steps are guaranteed to be 1/n or 1/m, D will always be an integer times 1/LCM(n,m))
            // the exact distribution is known and fairly simple in the n=m case (and in that case D will always be an integer times 1/n)
            ContinuousDistribution nullDistribution;
            if (AdvancedIntegerMath.BinomialCoefficient(a.Count + b.Count, a.Count) < Int64.MaxValue) {
                nullDistribution = new DiscreteAsContinuousDistribution(new KolmogorovTwoSampleExactDistribution(a.Count, b.Count), Interval.FromEndpoints(0.0, 1.0));
            } else {
                nullDistribution = new TransformedDistribution(new KolmogorovDistribution(), 0.0, Math.Sqrt(1.0 * (a.Count + b.Count) / a.Count / b.Count));
            }

            return (new TestResult("D", d, TestType.RightTailed, nullDistribution));

        }

        /// <summary>
        /// Tests whether the variances of two samples are compatible.
        /// </summary>
        /// <param name="a">The first sample.</param>
        /// <param name="b">The second sample.</param>
        /// <returns>The result of the test.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is <see langword="null"/>.</exception>
        public static TestResult FisherFTest (Sample a, Sample b) {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if ((a.Count < 2) || (b.Count < 2)) throw new InsufficientDataException();

            // compute population variances
            double v1 = a.Count / (a.Count - 1.0) * a.Variance;
            double v2 = b.Count / (b.Count - 1.0) * b.Variance;

            // compute the ratio
            double F = v1 / v2;

            return (new TestResult("F", F, TestType.RightTailed, new FisherDistribution(a.Count - 1, b.Count - 1)));

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
            if (array == null) throw new ArgumentNullException(nameof(array));
            data.CopyTo(array, start);
        }

        /// <summary>
        /// Performs a maximum likelihood fit.
        /// </summary>
        /// <param name="factory">A function that returns a distribution when given its defining parameters.</param>
        /// <param name="start">An initial guess for the defining parameters.</param>
        /// <returns>The result of the fit, containg the parameters that result in the best fit,
        /// covariance matrix among those parameters, and a test of the goodness of fit.</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Maximum_likelihood"/>
        public FitResult MaximumLikelihoodFit (Func<IList<double>, ContinuousDistribution> factory, IList<double> start) {

            if (factory == null) throw new ArgumentNullException(nameof(factory));
            if (start == null) throw new ArgumentNullException(nameof(start));

            // Define a log likelyhood function
            Func<IList<double>, double> L = (IList<double> parameters) => {
                ContinuousDistribution distribution = factory(parameters);
                double lnP = 0.0;
                foreach (double value in data) {
                    double P = distribution.ProbabilityDensity(value);
                    if (P == 0.0) throw new InvalidOperationException();
                    lnP += Math.Log(P);
                }
                return (-lnP);
            };

            // Maximize it
            MultiExtremum min = MultiFunctionMath.FindLocalMinimum(L, start);

            FitResult result = new FitResult(min.Location, min.HessianMatrix.CholeskyDecomposition().Inverse(), null);
            return (result);

        }

        /// <summary>
        /// Performs a Shapiro-Francia test of normality on the sample.
        /// </summary>
        /// <returns>The result of the test.</returns>
        /// <exception cref="InsufficientDataException">There are fewer than 16 values in the sample.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Shapiro%E2%80%93Francia_test"/>
        public TestResult ShapiroFranciaTest () {

            int n = this.Count;
            if (n < 16) throw new InsufficientDataException();

            // Determine V = 1 - r^2
            double mx = this.Mean;
            double cxx = 0.0;
            double cmm = 0.0;
            double cmx = 0.0;
            int[] order = data.GetSortOrder();
            for (int i = 0; i < n; i++) {
                double zx = data[order[i]] - mx;
                // We could use symmetry to make half as many calls to NormalOrderStatistic,
                // but the code would be less straightforward, because we can't proceed linearly
                // through the data.
                double m = NormalOrderStatistic(i + 1, n);
                cxx += zx * zx;
                cmm += m * m;
                cmx += zx * m;
            }
            // W = r^2 = cmx / sqrt(cxx * cmm). But if we compute r^2 like this and then 1 - r^2, we will
            // get significant cancellation errors. So re-form V = 1 - r^2 analytically.
            double q = cxx * cmm;
            double p = Math.Sqrt(cxx * cmm);
            double lnV = Math.Log((p - cmx) * (p + cmx) / q);

            // The observed mean and variance of log(1 - W) can be fit to six digits over n=16-8192
            // by simple polynomials in log(n). To get these results, I ran simulations for n = 16, 24, 32, 48, 64, ..., 8192.
            // The number of simulated samples was 10^8 for n=10-100, 10^7 for n=100-1000, and 10^6 for n=1000-10000.
            // Folloying Royston, I transformed to log(1 - W) and did polynomial fits to K1 and log(K2) as
            // functions of log(n), increasing the order until I was able to reproduce both cumulants to
            // within errors for all n=16-8192. Note that K3 is detectably non-zero, so we really should
            // introduce a correction for the remaining skewness, perhaps via an Edgeworth-style expansion
            // or by fitting to a skew-normal distribution.
            double t = Math.Log(n);
            double mu = -2.41348986 + t * (0.34785948 + t * (-0.31462118 + t * (0.04298010 + t * (-0.00311513 + t * 0.00009230))));
            double sigma = Math.Exp((0.56767288 + t * (-1.18529022 + t * (0.31916472 + t * (-0.04912452 + t * (0.00383692 + t * -0.00011891))))) / 2.0);
            NormalDistribution nullDistribution = new NormalDistribution(mu, sigma);

            // Need to handle 3-15 seperately.

            // N = 3: W' = W, 3/4 <= W <= 1, P = (6 / Pi)(Asin(Sqrt(W)) - Asin(Sqrt(3/4)))

            return (new TestResult("ln(1-W')", lnV, TestType.RightTailed, nullDistribution));

        }

        private static double NormalOrderStatistic (int i, int n) {

            // This uses the David and Johnson asymptotic expansion, as presented
            // in Arnold, Balakrishnan, and Nagaraja, "A First Course in Order Statistics",
            // Section 5.5, specialized to the normal distribution.
            
            // It is more accurate than simpler approximations, but still looses
            // accuracy in the tails. I have verified that dropping the last term
            // doesn't change the moments in my simulations to within their accuracy,
            // so it's not necessary to add more terms. If we ever do want to add
            // another term, see Childs and Balakrisnan, "Series approximations
            // for moments of order statistics using MAPLE", Computations Statistics
            // and Data Analysis 38 (2002) 331, which gives the lengthy next term.

            double p = (double) i / (n + 1);
            double q = (double) (n - i + 1) / (n + 1);

            double f0 = AdvancedMath.Probit(p, q);
            double p0 = Math.Exp(-f0 * f0 / 2.0) / Global.SqrtTwoPI;

            double f2 = f0 / (p0 * p0);
            double f3 = (1.0 + 2.0 * f0 * f0) / (p0 * p0 * p0);
            double f4 = f0 * (7.0 + 6.0 * f0 * f0) / (p0 * p0 * p0 * p0);
            double f5 = (7.0 + 46.0 * f0 * f0 + 24.0 * f0 * f0 * f0 * f0) / (p0 * p0 * p0 * p0 * p0);
            double f6 = f0 * (127.0 + 326.0 * f0 * f0 + 120.0 * f0 * f0 * f0 * f0) / (p0 * p0 * p0 * p0 * p0 * p0);

            double qmp = q - p;

            double m = f0;

            m += p * q / (n + 2) * f2 / 2.0;

            m += p * q / (n + 2) / (n + 2) * (qmp / 3.0 * f3 + p * q / 8.0 * f4);

            m += p * q / (n + 2) / (n + 2) / (n + 2) * (-qmp / 3.0 * f3 + (qmp * qmp - p * q) / 4.0 * f4 + p * q * qmp / 6.0 * f5 + p * p * q * q / 48.0 * f6);

            return (m);

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
