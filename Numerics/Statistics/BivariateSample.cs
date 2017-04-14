using System;
using System.Collections;
using System.Collections.Generic;
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
        public void Add (IList<double> x, IList<double> y) {
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
            internal set {
                isReadOnly = value;
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

        // marginals

        /// <summary>
        /// Gets a read-only univariate sample consisting of the x-values of the data points.
        /// </summary>
        /// <remarks>
        /// <para>Use this method to obtain sinformation specific to the x-vales, such as their <see cref="Sample.Median"/> or
        /// <see cref="Sample.Variance"/>.</para>
        /// <para>Note that this is a fast, O(1) operation, which does not create an independent copy of the data.
        /// The advantage of this is that you can access x-data as a <see cref="Sample"/> as often as you like without
        /// worying about performance. The disadvantage of this is that the returned sample cannot be altered. If you
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
        /// <para>Use this method to obtain sinformation specific to the y-vales, such as their <see cref="Sample.Median"/> or
        /// <see cref="Sample.Variance"/>.</para>
        /// <para>Note that this is a fast, O(1) operation, which does not create an independent copy of the data.
        /// The advantage of this is that you can access y-data as a <see cref="Sample"/> as often as you like without
        /// worying about performance. The disadvantage of this is that the returned sample cannot be altered. If you
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
                int n = Count;
                if (n < 2) throw new InsufficientDataException();
                double mx = xData.Mean;
                double my = yData.Mean;

                double C = 0.0;
                for (int i = 0; i < n; i++) {
                    C += (xData[i] - mx) * (yData[i] - my);
                }
                C = C / n;

                return (C);

                // this is an O(N) operation, which really shouldn't be done in a property
                // in the Sample class we maintain mean and variance as each value is
                // added or removed; ideally we should do that with covariance in the
                // bivariate sample, but Transforms, and maybe other operations, would
                // screw that up
            }
        }

        /// <summary>
        /// Gets the correlation coefficient between the two variables.
        /// </summary>
        public double CorrelationCoefficient {
            get {
                return (this.Covariance / Math.Sqrt(xData.Variance * yData.Variance));
            }
        }

        private double MomentAboutMean (int nx, int ny) {
            double mx = xData.Mean;
            double my = yData.Mean;

            double C = 0.0;
            for (int i = 0; i < Count; i++) {
                C += MoreMath.Pow(xData[i] - mx, nx) * MoreMath.Pow(yData[i] - my, ny);
            }
            return (C / Count);
        }

        /// <summary>
        /// Estimates of the population covariance of two variables.
        /// </summary>
        /// <returns>An estimate, with associated uncertainty, of the population covariance.</returns>
        public UncertainValue PopulationCovariance {
            get {
                // two moments are involved for best value and error
                int n = Count;
                double C_xy = MomentAboutMean(1, 1);
                double C_xxyy = MomentAboutMean(2, 2);

                return (new UncertainValue(C_xy * n / (n - 1), Math.Sqrt((C_xxyy - C_xy * C_xy) / n)));
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
            if (Count < 3) throw new InsufficientDataException();
            double r = this.Covariance / Math.Sqrt(xData.Variance * yData.Variance);
            ContinuousDistribution p = new PearsonRDistribution(Count);
            return (new TestResult(r, p));
        }

        /// <summary>
        /// Performs a Spearman rank-order test of association between the two variables.
        /// </summary>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>The Spearman rank-order test of association is a non-parametric test for association between
        /// two variables. The test statistic rho is the correlation coefficient of the <em>rank</em> of
        /// each entry in the sample. It is thus invariant over monotonic reparameterizations of the data,
        /// and will, for example, detect a quadratic or exponential association just as well as a linear
        /// association.</para>
        /// <para>The Spearman rank-order test requires O(N log N) operations.</para>
        /// </remarks>
        /// <exception cref="InsufficientDataException">There are fewer than three data points.</exception>
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="KendallTauTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient"/>
        public TestResult SpearmanRhoTest () {

            if (Count < 3) throw new InsufficientDataException();

            // analytic expressions for the mean and variance of ranks
            double M = (Count - 1) / 2.0;
            double V = (Count + 1) * (Count - 1) / 12.0;

            // compute the covariance of ranks
            int[] rx = xData.GetRanks();
            int[] ry = yData.GetRanks();
            double C = 0.0;
            for (int i = 0; i < Count; i++) {
                C += (rx[i] - M) * (ry[i] - M);
            }
            C = C / Count;

            // compute rho
            double rho = C / V;

            // for small enough sample, use the exact distribution
            ContinuousDistribution rhoDistribution;
            if (Count <= 10) {
                // for small enough sample, use the exact distribution
                // it would be nice to do this for at least slightly higher n, but computation time grows dramatically
                // would like to ensure return in less than 100ms; current timings n=10 35ms, n=11 72ms, n=12 190ms
                rhoDistribution = new DiscreteAsContinuousDistribution(new SpearmanExactDistribution(Count), Interval.FromEndpoints(-1.0, 1.0));
            } else {
                // for larger samples, use the normal approximation
                // would like to fit support and C_4 too; look into logit-normal
                // i was not happy with Edgeworth expansion, which can fit C_4 but screws up tails badly, even giving negative probabilities
                rhoDistribution = new NormalDistribution(0.0, 1.0 / Math.Sqrt(Count - 1));
            }

            return (new TestResult(rho, rhoDistribution));

        }

        /// <summary>
        /// Performs a Kendall concordance test for association.
        /// </summary>
        /// <returns>The result of the test.</returns>
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
        /// <para>Because it examine all pairs of data points, the Kendall test requires
        /// O(N<sup>2</sup>) operations. It is thus impractical for very large data sets. While
        /// not quite as robust as the Kendall test, the Spearman test is a good fall-back in such cases.</para>
        /// </remarks>
        /// <exception cref="InsufficientDataException"><see cref="Count"/> is less than two.</exception>
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="SpearmanRhoTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kendall_tau_test" />
        public TestResult KendallTauTest () {

            int n = xData.Count;
            if (n < 2) throw new InsufficientDataException();

            // loop over all pairs, counting concordant and discordant
            int C = 0;
            int D = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {

                    // note the way each variable varies in the pair
                    int sx = Math.Sign(xData[i] - xData[j]);
                    int sy = Math.Sign(yData[i] - yData[j]);

                    // if they vary in the same way, they are concordant, otherwise they are discordant
                    if (sx == sy) {
                        C++;
                    } else {
                        D++;
                    }
                    // note this does not count ties specially, as is sometimes done

                }
            }

            // compute tau
            double t = 1.0 * (C - D) / (C + D);

            // compute tau distribution
            ContinuousDistribution tauDistribution;
            if (n <= 20) {
                tauDistribution = new DiscreteAsContinuousDistribution(new KendallExactDistribution(n), Interval.FromEndpoints(-1.0, 1.0));
            } else {
                double dt = Math.Sqrt((4 * n + 10) / 9.0 / n / (n - 1));
                tauDistribution = new NormalDistribution(0.0, dt);
            }

            return (new TestResult(t, tauDistribution));

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

            if (Count < 2) throw new InsufficientDataException();

            // the paired t-test is just a normal t-test of one sample against a mean,
            // but where the sample consists of the differences between the paired measurements,
            // and the mean being tested against is (usually) zero

            // loop over pairs, computing mean and standard deviation of differences
            double m = 0.0;
            double v = 0.0;
            for (int i = 0; i < Count; i++) {
                double z = xData[i] - yData[i];
                v += MoreMath.Pow2(z - m) * i / (i + 1);
                m += (z - m) / (i + 1);
            }
            v = v / (Count - 1);

            // compute standard error
            double s = Math.Sqrt(v / Count);

            // t is the mean deviation as a fraction of standard error
            double t = m / s;

            return (new TestResult("t", t, TestType.TwoTailed, new StudentDistribution(Count - 1)));

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

            int n = this.Count;
            if (n < 3) throw new InsufficientDataException();

            // The means and covariances are the inputs to most of the regression formulas.
            double mx = xData.Mean;
            double my = yData.Mean;
            double cxx = xData.Variance;
            double cyy = yData.Variance;
            double cxy = this.Covariance;

            Debug.Assert(cxx >= 0.0);
            Debug.Assert(cyy >= 0.0);

            // Compute the best-fit parameters
            double b = cxy / cxx;
            double a = my - b * mx;
            // Since cov(x,y) = (n S_xy - S_x S_y)/n^2 and var(x) = (n S_xx - S_x^2) / n^2,
            // these formulas are equivilent to the 
            // to the usual formulas for a and b involving sums, but it is more stable against round-off
            ColumnVector v = new ColumnVector(a, b);
            v.IsReadOnly = true;

            // Compute Pearson r value
            double r = cxy / Math.Sqrt(cxx * cyy);
            TestResult rTest = new TestResult("r", r, TestType.TwoTailed, new Distributions.PearsonRDistribution(n));

            // Compute residuals and other sum-of-squares
            double SSR = 0.0;
            double SSF = 0.0;
            Sample residuals = new Sample();
            foreach (XY point in this) {
                double y = a + b * point.X;
                double z = point.Y - y;
                SSR += z * z;
                residuals.Add(z);
                SSF += MoreMath.Sqr(y - my);
            }
            double SST = cyy * n;
            // Note SST = SSF + SSR because \sum_{i} ( y_i - \bar{y})^2 = \sum_i (y_i - f_i)^2 + \sum_i (f_i - \bar{y})^2

            // Use sums-of-squares to do ANOVA
            AnovaRow fit = new AnovaRow(SSF, 1);
            AnovaRow residual = new AnovaRow(SSR, n - 2);
            AnovaRow total = new AnovaRow(SST, n - 1);
            OneWayAnovaResult anova = new OneWayAnovaResult(fit, residual, total);

            // Compute covariance of parameters matrix
            double s2 = SSR / (n - 2);
            double cbb = s2 / cxx / n;
            double cab = -mx * cbb;
            double caa = (cxx + mx * mx) * cbb;

            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = caa;
            C[1, 1] = cbb;
            C[0, 1] = cab;
            C.IsReadOnly = true;

            // Prepare the prediction function
            Func<double, UncertainValue> predict = (double x) => {
                double y = a + b * x;
                return (new UncertainValue(y, Math.Sqrt(s2 * (1.0 + (1.0 + MoreMath.Sqr(x - mx) / cxx) / n))));
            };

            return (new LinearRegressionResult(v, C, rTest, anova, residuals, predict));

        }

        /// <summary>
        /// Computes the polynomial of given degree which best fits the data.
        /// </summary>
        /// <param name="m">The degree, which must be non-negative.</param>
        /// <returns>The fit result.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> is negative.</exception>
        /// <exception cref="InsufficientDataException">There are fewer data points than coefficients to be fit.</exception>
        public PolynomialRegressionResult PolynomialRegression (int m) {

            if (m < 0) throw new ArgumentOutOfRangeException(nameof(m));

            int n = Count;
            if (n < (m + 1)) throw new InsufficientDataException();

            // Construct the n X m design matrix X_{ij} = x_{i}^{j}
            RectangularMatrix X = new RectangularMatrix(n, m + 1);
            ColumnVector y = new ColumnVector(n);
            for (int i = 0; i < n; i++) {
                double x = xData[i];
                X[i, 0] = 1.0;
                for (int j = 1; j <= m; j++) {
                    X[i, j] = X[i, j - 1] * x;
                }
                y[i] = yData[i];
            }

            // Use X = QR to solve X b = y and compute C
            ColumnVector b;
            SymmetricMatrix C;
            QRDecomposition.SolveLinearSystem(X, y, out b, out C);

            // Compute residuals
            double SSR = 0.0;
            double SSF = 0.0;
            ColumnVector yHat = X * b;
            Sample residuals = new Sample();
            for (int i = 0; i < n; i++) {
                double z = yData[i] - yHat[i];
                residuals.Add(z);
                SSR += z * z;
                SSF += MoreMath.Sqr(yHat[i] - yData.Mean);
            }
            double sigma2 = SSR / (n - (m + 1));

            // Scale up C by \sigma^2
            // (It sure would be great to be able to overload *=.)
            for (int i = 0; i <= m; i++) {
                for (int j = i; j <= m; j++) {
                    C[i, j] = C[i, j] * sigma2;
                }
            }

            // Compute remaing sums-of-squares
            double SST = yData.Variance * n;

            // Use sums-of-squares to do ANOVA
            AnovaRow fit = new AnovaRow(SSF, m);
            AnovaRow residual = new AnovaRow(SSR, n - (m + 1));
            AnovaRow total = new AnovaRow(SST, n - 1);
            OneWayAnovaResult anova = new OneWayAnovaResult(fit, residual, total);

            return (new PolynomialRegressionResult(b, C, anova, residuals));

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
        public FitResult NonlinearRegression (Func<IList<double>, double, double> f, IList<double> start) {

            if (f == null) throw new ArgumentNullException(nameof(f));
            if (start == null) throw new ArgumentNullException(nameof(start));

            int n = this.Count;
            int d = start.Count;
            if (n <= d) throw new InsufficientDataException();

            MultiExtremum min = MultiFunctionMath.FindLocalMinimum((IList<double> a) => {
                double ss = 0.0;
                for (int i = 0; i < n; i++) {
                    double r = yData[i] - f(a, xData[i]);
                    ss += r * r;
                }
                return (ss);
            }, start);

            CholeskyDecomposition cholesky = min.HessianMatrix.CholeskyDecomposition();
            if (cholesky == null) throw new DivideByZeroException();
            SymmetricMatrix curvature = cholesky.Inverse();
            curvature = (2.0 * min.Value / (n - d)) * curvature;

            FitResult result = new FitResult(min.Location, curvature, null);

            return (result);
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
        public FitResult LinearLogisticRegression () {

            // check size of data set
            if (Count < 3) throw new InsufficientDataException();

            // check limits on Y
            double p = Y.Mean; double q = 1.0 - p;
            if ((p <= 0.0) || (q <= 0.0)) throw new InvalidOperationException();

            // Define the log-likelyhood as a function of the parameters
            Func<IList<double>,double> f = delegate (IList<double> a) {
                double L = 0.0;
                for (int i = 0; i < Count; i++) {
                    double x = xData[i];
                    double z = a[0] + a[1] * x;
                    double ez = Math.Exp(z);
                    double y = yData[i];
                    if (y == 0.0) {
                        L += Math.Log(1.0 + ez);
                    } else if (y == 1.0) {
                        L += Math.Log(1.0 + 1.0 / ez);
                    } else {
                        throw new InvalidOperationException();
                    }
                }
                return(L);
            };

            // We need an initial guess at the parameters.

            // This particular initial guess is driven by the fact that, in a logistic model
            //    \frac{\partial p}{\partial x} = \beta p ( 1 - p)
            // Evaluating at means, and noting that p (1 - p) = var(y) and that, in a development around the means,
            //    cov(p, x) = \frac{\partial p}{\partial x} var(x)
            // we get
            //    \beta = \frac{cov(y, x)}{var(x) var(y)}
            // This approximation gets the sign right, but it looks like it usually gets the magnitude quite wrong.

            double b0 = Covariance / X.Variance / Y.Variance;
            double a0 = Math.Log(p / q) - b0 * X.Mean;

            MultiExtremum m = MultiFunctionMath.FindLocalMinimum(f, new double[2] { a0, b0 });
            return (new FitResult(m.Location, m.HessianMatrix.Inverse(), null));
             
            //SpaceExtremum m = FunctionMath.FindMinimum(f, new double[2] { a0, b0 });
            //return (new FitResult(m.Location(), m.Curvature().Inverse(), null));

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
