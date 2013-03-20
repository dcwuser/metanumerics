using System;
#if !SILVERLIGHT
using System.Data;
#endif
using System.Globalization;

using Meta.Numerics.Functions;
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
    public sealed class BivariateSample {

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
        /// Copies the bivariate sample.
        /// </summary>
        /// <returns>An independent copy of the bivariate sample.</returns>
        public BivariateSample Copy () {
            return (new BivariateSample(xData.Copy(), yData.Copy(), false));
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
            Distribution p = new PearsonRDistribution(Count);
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
            Distribution rhoDistribution;
            if (Count <= 10) {
                // for small enough sample, use the exact distribution
                // it would be nice to do this for at least slightly higher n, but computation time grows dramatically
                // would like to ensure return in less than 100ms; current timings n=10 35ms, n=11 72ms, n=12 190ms
                rhoDistribution = new SpearmanDistribution(Count);
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

            // analytic expression for variance of tau
            double dt = Math.Sqrt((4 * n + 10) / 9.0 / n / (n - 1));

            return (new TestResult(t, new NormalDistribution(0.0, dt)));

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

            return (new TestResult(t, new StudentDistribution(Count - 1)));

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
        public FitResult LinearRegression () {

            int n = Count;
            if (n < 3) throw new InsufficientDataException();

            double b = this.Covariance / xData.Variance;
            double a = yData.Mean - b * xData.Mean;

            // since cov(x,y) = (n S_xy - S_x S_y)/n^2 and var(x) = (n S_xx - S_x^2) / n^2 this is equivilent
            // to the usual formulas for a and b involving sums, but it is more stable against round-off

            double SSR = 0.0;
            double SXX = 0.0;
            for (int i = 0; i < n; i++) {
                double z = a + b * xData[i] - yData[i];
                SSR += z * z;
                SXX += MoreMath.Pow2(xData[i]);
            }

            double cbb = SSR / xData.Variance / n / (n - 2);
            double cab = -xData.Mean * cbb;
            double caa = SXX / n * cbb;

            // compute F-statistic
            // total variance dof = n - 1, explained variance dof = 1, unexplained variance dof = n - 2
            double totalVarianceSum = yData.Variance * n;
            double unexplainedVarianceSum = SSR;
            double explainedVarianceSum = totalVarianceSum - unexplainedVarianceSum;
            double F = (explainedVarianceSum / 1) / (unexplainedVarianceSum / (n - 2));
            TestResult test = new TestResult(F, new FisherDistribution(1, n - 2));

            return (new FitResult(a, Math.Sqrt(caa), b, Math.Sqrt(cbb), cab, test));
        }

        /// <summary>
        /// Computes the polynomial of given degree which best fits the data.
        /// </summary>
        /// <param name="m">The degree, which must be non-negative.</param>
        /// <returns>The fit result.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> is negative.</exception>
        /// <exception cref="InsufficientDataException">There are fewer data points than coefficients to be fit.</exception>
        public FitResult PolynomialRegression (int m) {

            if (m < 0) throw new ArgumentOutOfRangeException("m");

            int n = Count;
            if (n < m + 1) throw new InsufficientDataException();

            // Construct the n X m design matrix A_{ij} = x_{i}^{j}
            RectangularMatrix A = new RectangularMatrix(n, m + 1);
            ColumnVector y = new ColumnVector(n);
            for (int i = 0; i < n; i++) {
                double x = xData[i];
                A[i, 0] = 1.0;
                for (int j = 1; j <= m; j++) {
                    A[i, j] = A[i, j - 1] * x;
                }
                y[i] = yData[i];
            }

            // Our problem is to solve A a = y, where a is the vector of coefficients. QR gives the least squares
            // solution that minimizes | A a - y |.
            QRDecomposition QR = A.QRDecomposition();
            RectangularMatrix R = QR.RMatrix();
            ColumnVector a = QR.Solve(y);

            // Compute the residual vector r and s^2 = r^2 / dof
            ColumnVector r = A * a - y;
            double V = r.Transpose() * r;
            double ss2 = V / (n - (m + 1));

            // Ake Bjorck, "Numerical Methods for Least Squares Problems", pp. 118-120

            // The covariance matrix V = s^2 C, with C = (R^T R)^{-1}. This is a matrix multiplication plus an inversion.
            
            // A faster way is to use the direct solution
            //   for k = n ... 1
            //     C_{kk} = R_{kk}^{-1} \left[ R_{kk}^{-1} - \sum_{j=k+1}^{n} R_{kj} C_{kj} \right]
            //     for i = k-1 ... 1
            //       C_{ik} = -R_{ii}^{-1} \left[ \sum_{j=i+1}^{k} R_{ij} C_{jk} + \sum_{j=k+1}^{n} R_{ij} C_{kj} \right]
            //     end
            //   end
            // This is detailed in Ake Bjorck, "Numerical Methods for Least Squares Problems", pp. 118-120

            //Stopwatch s4 = Stopwatch.StartNew();
            SymmetricMatrix C2 = new SymmetricMatrix(m + 1);
            double c; // used for storage so we don't call bounds-checking accessors each time
            for (int k = m; k >= 0; k--) {
                c = 1.0 / R[k, k];;
                for (int j = k + 1; j <= m; j++) {
                    c -= R[k, j] * C2[k, j];
                }
                C2[k, k] = c / R[k, k];
                for (int i = k - 1; i >= 0; i--) {
                    c = 0.0;
                    for (int j = i + 1; j <= k; j++) {
                        c += R[i, j] * C2[j, k];
                    }
                    for (int j = k + 1; j <= m; j++) {
                        c += R[i, j] * C2[k, j];
                    }
                    C2[i, k] = -c / R[i, i];
                }
            }
            //s4.Stop();
            //Console.WriteLine("C new {0}", s4.ElapsedMilliseconds);
            for (int i = 0; i <= m; i++) {
                for (int j = i; j <= m; j++) {
                    C2[i, j] = C2[i, j] * ss2;
                }
            }


            // compute F-statistic
            // total variance dof = n - 1, explained variance dof = m, unexplained variance dof = n - (m + 1)
            double totalVarianceSum = yData.Variance * n;
            double unexplainedVarianceSum = V;
            double explainedVarianceSum = totalVarianceSum - unexplainedVarianceSum;
            double unexplainedVarianceDof = n - (m + 1);
            double explainedVarianceDof = m;
            double F = (explainedVarianceSum / explainedVarianceDof) / (unexplainedVarianceSum / unexplainedVarianceDof);
            TestResult test = new TestResult(F, new FisherDistribution(explainedVarianceDof, unexplainedVarianceDof));

            return (new FitResult(a, C2, test));
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

            // make an initial guess at the parameters
            double b0 = Covariance / X.Variance / Y.Variance;
            double a0 = Math.Log(p / q) - b0 * X.Mean;

            Func<double[],double> f = delegate (double[] a) {

                // define a logistic log-likelyhood function
                double L = 0.0;
                for (int i = 0; i < Count; i++) {
                    double x = xData[i];
                    double z = a[0] + a[1] * x;
                    double ez = Math.Exp(z);
                    double y = yData[i];
                    if (y == 0.0) {
                        L -= Math.Log(1.0 / (1.0 + ez));
                    } else if (y == 1.0) {
                        L -= Math.Log(ez / (1.0 + ez));
                    } else {
                        throw new InvalidOperationException();
                    }
                }
                return(L);
            };

            SpaceExtremum m = FunctionMath.FindMinimum(f, new double[2] { a0, b0 });

            return (new FitResult(m.Location(), m.Curvature().Inverse(), null));

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
