using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
#if !SILVERLIGHT
using System.Data;
#endif
using System.Globalization;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

using System.Diagnostics;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a set of data points, where each data point is described by a pair of real numbers.
    /// </summary>
    /// <remarks>
    /// <para>A bivariate sample consists of pairs of real numbers, where each pair is an independent
    /// measurement. For example, if you measure the height and weight of a sample of people, the data
    /// could be stored as a bivariate sample. The class can compute various descriptive statistics
    /// for the sample, perform appropriate statistical tests on the sample data, and fit the sample
    /// data to various models.</para>
    /// </remarks>
    public sealed class BivariateSample : ICollection<XY> {

        /// <summary>
        /// Initializes a new bivariate sample.
        /// </summary>
        public BivariateSample () {
            xData = new SampleStorage();
            yData = new SampleStorage();
            isReadOnly = false;
        }

        /// <summary>
        /// Initializes a new bivariate sample with the given variable names.
        /// </summary>
        /// <param name="xName">The name of the x-variable.</param>
        /// <param name="yName">The name of the y-variable.</param>
        public BivariateSample (string xName, string yName) : this () {
            xData.Name = xName;
            yData.Name = yName;
        }

        internal BivariateSample (SampleStorage xData, SampleStorage yData, bool isReadOnly) {
            this.xData = xData;
            this.yData = yData;
            this.isReadOnly = isReadOnly;
        }

        private SampleStorage xData;
        private SampleStorage yData;
        private bool isReadOnly;

        // manipulation methods

        /// <summary>
        /// Adds a data point to the sample.
        /// </summary>
        /// <param name="x">The x-value of the data point.</param>
        /// <param name="y">The y-value of the data point.</param>
        public void Add (double x, double y) {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                xData.Add(x);
                yData.Add(y);
            }
        }

        /// <summary>
        /// Adds a data point to the sample.
        /// </summary>
        /// <param name="point">The data point.</param>
        public void Add (XY point) {
            Add(point.X, point.Y);
        }

        /// <summary>
        /// Adds multiple data points to the sample.
        /// </summary>
        /// <param name="points">The data points.</param>
        public void Add (IEnumerable<XY> points) {
            if (points == null) throw new ArgumentNullException(nameof(points));
            foreach (XY point in points) Add(point);
        }

        /// <summary>
        /// Adds points from two lists to the sample.
        /// </summary>
        /// <param name="x">The x values of the data points.</param>
        /// <param name="y">The y values of the data points.</param>
        /// <exception cref="DimensionMismatchException">The lengths of the two lists are not equal.</exception>
        public void Add (IReadOnlyList<double> x, IReadOnlyList<double> y) {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();
            for (int i = 0; i < x.Count; i++) {
                Add(x[i], y[i]);
            }
        }

        private int IndexOf (double x, double y) {
            for (int i = 0; i < Count; i++) {
                if ((xData[i] == x) && (yData[i] == y)) {
                    return (i);
                }
            }
            return (-1);
        }

        /// <summary>
        /// Removes a data point from the sample.
        /// </summary>
        /// <param name="x">The x-value of the data point.</param>
        /// <param name="y">The y-value of the data point.</param>
        /// <returns>True if the given data point was found and removed, otherwise false.</returns>
        public bool Remove (double x, double y) {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                int index = IndexOf(x, y);
                if (index < 0) {
                    return (false);
                } else {
                    xData.RemoveAt(index);
                    yData.RemoveAt(index);
                    return (true);
                }
            }
        }

        /// <summary>
        /// Removes a data point from the sample.
        /// </summary>
        /// <param name="point">The point to remove.</param>
        /// <returns>True if the given data point was found and removed, otherwise false.</returns>
        public bool Remove (XY point) {
            return(Remove(point.X, point.Y));
        }

        /// <summary>
        /// Removes all data points from the sample.
        /// </summary>
        public void Clear () {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                xData.Clear();
                yData.Clear();
            }
        }

        /// <summary>
        /// Determines whether the sample contains a given data point.
        /// </summary>
        /// <param name="x">The x-value of the data point.</param>
        /// <param name="y">The y-value of the data point.</param>
        /// <returns>True if the sample contains the given data point, otherwise false.</returns>
        public bool Contains (double x, double y) {
            int index = IndexOf(x, y);
            return ((index >= 0));
        }

        /// <summary>
        /// Determines whether the sample contains a given data point.
        /// </summary>
        /// <param name="xy">The data point.</param>
        /// <returns>True if the sample contains the given data point, otherwise false.</returns>
        public bool Contains (XY xy) {
            return (Contains(xy.X, xy.Y));
        }

        /// <summary>
        /// Gets a value indicating whether the bivariate sample is read-only.
        /// </summary>
        public bool IsReadOnly {
            get {
                return (isReadOnly);
            }
        }

        /// <summary>
        /// Copies the bivariate sample.
        /// </summary>
        /// <returns>An independent copy of the bivariate sample.</returns>
        public BivariateSample Copy () {
            return (new BivariateSample(xData.Copy(), yData.Copy(), false));
        }


        /// <summary>
        /// Gets an enumerator of sample values.
        /// </summary>
        /// <returns>An enumerator of sample values.</returns>
        public IEnumerator<XY> GetEnumerator () {
            int i = 0;
            while (i < this.Count) {
                yield return (new XY(xData[i], yData[i]));
                i++;
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (this.GetEnumerator());
        }

        void ICollection<XY>.CopyTo (XY[] array, int start) {
            if (array == null) throw new ArgumentNullException(nameof(array));
            // this is not fast; move implementation to SampleStorage and use CopyTo for underlying List
            for (int i = 0; i < this.Count; i++) {
                array[start + i] = new XY(xData[i], yData[i]);
            }
        }

        /// <summary>
        /// Swaps the X and Y variables in the bivariate sample.
        /// </summary>
        public void TransposeXY () {
            SampleStorage t = yData;
            yData = xData;
            xData = t;
        }

        /// <summary>
        /// Gets a read-only univariate sample consisting of the x-values of the data points.
        /// </summary>
        /// <remarks>
        /// <para>Use this method to obtain information specific to the x-vales, such as their <see cref="Sample.Median"/> or
        /// <see cref="Sample.Variance"/>.</para>
        /// <para>Note that this is a fast, O(1) operation, which does not create an independent copy of the data.
        /// The advantage of this is that you can access x-data as a <see cref="Sample"/> as often as you like without
        /// worrying about performance. The disadvantage of this is that the returned sample cannot be altered. If you
        /// need to alter x-data independent of the bivariate sample, use the <see cref="Sample.Copy"/>
        /// method to obtain an independent copy.</para>
        /// </remarks>
        public Sample X {
            get {
                return (new Sample(xData, false));
            }
        }

        /// <summary>
        /// Gets a read-only univariate sample consisting of the y-values of the data points.
        /// </summary>
        /// <remarks>
        /// <para>Use this method to obtain information specific to the y-vales, such as their <see cref="Sample.Median"/> or
        /// <see cref="Sample.Variance"/>.</para>
        /// <para>Note that this is a fast, O(1) operation, which does not create an independent copy of the data.
        /// The advantage of this is that you can access y-data as a <see cref="Sample"/> as often as you like without
        /// worrying about performance. The disadvantage of this is that the returned sample cannot be altered. If you
        /// need to alter y-data independent of the bivariate sample, use the <see cref="Sample.Copy"/>
        /// method to obtain an independent copy.</para>
        /// </remarks>
        public Sample Y {
            get {
                return (new Sample(yData, false));
            }
        }

        // basic properties

        /// <summary>
        /// Gets the number of data points.
        /// </summary>
        public int Count {
            get {
                return (xData.Count);
            }
        }

        /// <summary>
        /// Gets the covariance of the two variables.
        /// </summary>
        public double Covariance {
            get {
                return (Bivariate.Covariance(xData, yData));
            }
        }

        /// <summary>
        /// Gets the correlation coefficient between the two variables.
        /// </summary>
        public double CorrelationCoefficient {
            get {
                return (Bivariate.CorrelationCoefficient(xData, yData));
            }
        }

        /// <summary>
        /// Estimates of the population covariance of two variables.
        /// </summary>
        /// <returns>An estimate, with associated uncertainty, of the population covariance.</returns>
        public UncertainValue PopulationCovariance {
            get {
                return(Bivariate.PopulationCovariance(xData, yData));
            }
        }

        /// <summary>
        /// Performs a Pearson correlation test for association.
        /// </summary>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>This test measures the strength of the linear correlation between two variables. The
        /// test statistic r is simply the covariance of the two variables, scaled by their respective
        /// standard deviations so as to obtain a number between -1 (perfect linear anti-correlation)
        /// and +1 (perfect linear correlation).</para>
        /// <para>The Pearson test cannot reliably detect or rule out non-linear associations. For example,
        /// variables with a perfect quadratic association may have only a weak linear correlation. If
        /// you wish to test for associations that may not be linear, consider using the Spearman or
        /// Kendall tests instead.</para>
        /// <para>The Pearson correlation test requires O(N) operations.</para>
        /// <para>The Pearson test requires at least three bivariate values.</para>
        /// </remarks>
        /// <exception cref="InsufficientDataException"><see cref="Count"/> is less than three.</exception>
        /// <seealso cref="SpearmanRhoTest"/>
        /// <seealso cref="KendallTauTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Pearson_correlation_coefficient" />
        public TestResult PearsonRTest () {
            return (Bivariate.PearsonRTest(xData, yData));
        }

        /// <summary>
        /// Performs a Spearman rank-order test of association between the two variables.
        /// </summary>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>The Spearman rank-order test of association is a non-parametric test for association between
        /// two variables. The test statistic rho is the correlation coefficient of the <em>rank</em> of
        /// each entry in the sample. It is thus invariant over monotonic re-parameterizations of the data,
        /// and will, for example, detect a quadratic or exponential association just as well as a linear
        /// association.</para>
        /// <para>The Spearman rank-order test requires O(N log N) operations.</para>
        /// </remarks>
        /// <exception cref="InsufficientDataException">There are fewer than three data points.</exception>
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="KendallTauTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient"/>
        public TestResult SpearmanRhoTest () {
            return (Bivariate.SpearmanRhoTest(xData, yData));
        }

        /// <summary>
        /// Performs a Kendall concordance test for association.
        /// </summary>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>Kendall's &#x3C4; is a non-parametric and robust test of association
        /// between two variables. It simply measures the number of cases where an increase
        /// in one variable is associated with an increase in the other (concordant pairs),
        /// compared with the number of cases where an increase in one variable is associated
        /// with a decrease in the other (discordant pairs).</para>
        /// <para>Because &#x3C4; depends only on the sign
        /// of a change and not its magnitude, it is not skewed by outliers exhibiting very large
        /// changes, nor by cases where the degree of change in one variable associated with
        /// a given change in the other changes over the range of the variables. Of course, it may
        /// still miss an association whose sign changes over the range of the variables. For example,
        /// if data points lie along a semi-circle in the plane, an increase in the first variable
        /// is associated with an increase in the second variable along the rising arc and and decrease in
        /// the second variable along the falling arc. No test that looks for single-signed correlation
        /// will catch this association.
        /// </para>
        /// <para>Because it examine all pairs of data points, the Kendall test requires
        /// O(N<sup>2</sup>) operations. It is thus impractical for very large data sets. While
        /// not quite as robust as the Kendall test, the Spearman test is a good fall-back in such cases.</para>
        /// </remarks>
        /// <exception cref="InsufficientDataException"><see cref="Count"/> is less than two.</exception>
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="SpearmanRhoTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kendall_tau_test" />
        public TestResult KendallTauTest () {
            return (Bivariate.KendallTauTest(xData, yData));
        }

        /// <summary>
        /// Performs a paired Student t-test.
        /// </summary>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>Like a two-sample, unpaired t-test (<see cref="Sample.StudentTTest(Sample,Sample)" />),
        /// a paired t-test compares two samples to detect a difference in means.
        /// Unlike the unpaired version, the paired version assumes that each </para>
        /// </remarks>
        /// <exception cref="InsufficientDataException">There are fewer than two data points.</exception>
        public TestResult PairedStudentTTest () {
            return (Bivariate.PairedStudentTTest(xData, yData));
        }

        /// <summary>
        /// Performs a Wilcoxon signed rank test.
        /// </summary>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>The Wilcoxon signed rank test is a non-parametric alternative to the
        /// paired t-test (<see cref="PairedStudentTTest"/>). Given two measurements on
        /// the same subjects, this method tests for changes in the distribution between
        /// the two measurements. It is sensitive primarily to shifts in the median.
        /// Note that the distributions of the individual measurements
        /// may be far from normal, and may be different for each subject.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test"/>
        public TestResult WilcoxonSignedRankTest () {
            return (Bivariate.WilcoxonSignedRankTest(xData, yData));
        }

        /// <summary>
        /// Computes the best-fit linear regression from the data.
        /// </summary>
        /// <returns>The result of the fit.</returns>
        /// <remarks>
        /// <para>Linear regression assumes that the data have been generated by a function y = a + b x + e, where e is
        /// normally distributed noise, and determines the values of a and b that best fit the data. It also
        /// determines an error matrix on the parameters a and b, and does an F-test to</para>
        /// <para>The fit result is two-dimensional. The first parameter is the intercept a, the second is the slope b.
        /// The goodness-of-fit test is a F-test comparing the variance accounted for by the model to the remaining,
        /// unexplained variance.</para>
        /// </remarks>
        /// <exception cref="InsufficientDataException">There are fewer than three data points.</exception>
        public LinearRegressionResult LinearRegression () {
            return (Bivariate.LinearRegression(yData, xData));
        }

        /// <summary>
        /// Computes the polynomial of given degree which best fits the data.
        /// </summary>
        /// <param name="m">The degree, which must be non-negative.</param>
        /// <returns>The fit result.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> is negative.</exception>
        /// <exception cref="InsufficientDataException">There are fewer data points than coefficients to be fit.</exception>
        public PolynomialRegressionResult PolynomialRegression (int m) {
            return (Bivariate.PolynomialRegression(yData, xData, m));
        }


        /// <summary>
        /// Finds the parameterized function that best fits the data.
        /// </summary>
        /// <param name="f">The parameterized function.</param>
        /// <param name="start">An initial guess for the parameters.</param>
        /// <returns>The fit result.</returns>
        /// <remarks>
        /// <para>
        /// In the returned <see cref="FitResult"/>, the parameters appear in the same order as in
        /// the supplied fit function and initial guess vector. No goodness-of-fit test is returned.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> or <paramref name="start"/> is null.</exception>
        /// <exception cref="InsufficientDataException">There are not more data points than fit parameters.</exception>
        /// <exception cref="DivideByZeroException">The curvature matrix is singular, indicating that the data is independent of
        /// one or more parameters, or that two or more parameters are linearly dependent.</exception>
        public NonlinearRegressionResult NonlinearRegression (Func<IReadOnlyList<double>, double, double> f, IReadOnlyList<double> start) {
            if (start == null) throw new ArgumentNullException(nameof(start));
            return (Bivariate.NonlinearRegression(yData, xData, f, start));
        }

        /// <summary>
        /// Computes the best-fit linear logistic regression from the data.
        /// </summary>
        /// <returns>The fit result.</returns>
        /// <remarks>
        /// <para>Linear logistic regression is a way to fit binary outcome data to a linear model.</para>
        /// <para>The method assumes that binary outcomes are encoded as 0 and 1. If any y-values other than
        /// 0 and 1 are encountered, it throws an <see cref="InvalidOperationException"/>.</para>
        /// <para>The fit result is two-dimensional. The first parameter is a, the second b.</para>
        /// </remarks>
        /// <exception cref="InsufficientDataException">There are fewer than three data points.</exception>
        /// <exception cref="InvalidOperationException">There is a y-value other than 0 or 1.</exception>
        public LinearLogisticRegressionResult LinearLogisticRegression () {
            List<bool> y = yData.Select(v => { if (v == 0.0) { return (false); } else if (v == 1.0) { return (true); } else { throw new InvalidOperationException(); } }).ToList();
            return (Bivariate.LinearLogisticRegression(y, xData));
        }

#if !SILVERLIGHT
        /// <summary>
        /// Loads values from a data reader.
        /// </summary>
        /// <param name="reader">The data reader.</param>
        /// <param name="xIndex">The column number of the x-variable.</param>
        /// <param name="yIndex">The column number of the y-variable.</param>
        public void Load (IDataReader reader, int xIndex, int yIndex) {
            if (reader == null) throw new ArgumentNullException("reader");
            if (isReadOnly) throw new InvalidOperationException();
            while (reader.Read()) {
                if (reader.IsDBNull(xIndex) || reader.IsDBNull(yIndex)) continue;
                object xValue = reader.GetValue(xIndex);
                object yValue = reader.GetValue(yIndex);
                Add(Convert.ToDouble(xValue, CultureInfo.InvariantCulture), Convert.ToDouble(yValue, CultureInfo.InvariantCulture));
            }
        }
#endif

    }

}
