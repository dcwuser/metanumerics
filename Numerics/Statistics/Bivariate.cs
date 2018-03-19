using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Contains methods for performing statistics on bivariate samples.
    /// </summary>
    public static class Bivariate {

        internal static void ComputeBivariateMomentsUpToTwo (IEnumerable<double> x, IEnumerable<double> y, out int n, out double xMean, out double yMean, out double xxSum, out double yySum, out double xySum) {

            Debug.Assert(x != null);
            Debug.Assert(y != null);

            IEnumerator<double> xEnumerator = x.GetEnumerator();
            IEnumerator<double> yEnumerator = y.GetEnumerator();

            n = 0;
            xMean = 0.0;
            yMean = 0.0;
            xxSum = 0.0;
            yySum = 0.0;
            xySum = 0.0;

            while (true) {

                bool xFlag = xEnumerator.MoveNext();
                bool yFlag = yEnumerator.MoveNext();
                Debug.Assert(xFlag == yFlag);
                if (!xFlag) {
                    break;
                }

                n++;
                double e = 1.0 / n;

                double xValue = xEnumerator.Current;
                double xDelta = xValue - xMean;
                xMean += xDelta * e;

                double yValue = yEnumerator.Current;
                double yDelta = yValue - yMean;
                yMean += yDelta * e;

                e = 1.0 - e;

                xxSum += xDelta * xDelta * e;
                yySum += yDelta * yDelta * e;
                xySum += xDelta * yDelta * e;

            }

        }


        /// <summary>
        /// Computes the covariance of the two variables.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <returns>The covariance between the two variables.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than two entries in the sample.</exception>
        public static double Covariance (this IReadOnlyList<double> x, IReadOnlyList<double> y) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();
            if (x.Count < 2) throw new InsufficientDataException();

            int n;
            double xMean, yMean, xxSum, yySum, xySum;
            ComputeBivariateMomentsUpToTwo(x, y, out n, out xMean, out yMean, out xxSum, out yySum, out xySum);
            Debug.Assert(n == x.Count);

            return (xySum / n);
        }

        /// <summary>
        /// Computes the correlation coefficient between the two variables.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <returns>The correlation coefficient between the two variables.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">Are fewer than two entries in the sample.</exception>
        public static double CorrelationCoefficient (this IReadOnlyList<double> x, IReadOnlyList<double> y) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();
            if (x.Count < 2) throw new InsufficientDataException();

            int n;
            double xMean, yMean, xxSum, yySum, xySum;
            ComputeBivariateMomentsUpToTwo(x, y, out n, out xMean, out yMean, out xxSum, out yySum, out xySum);

            Debug.Assert(xxSum >= 0.0);
            Debug.Assert(yySum >= 0.0);

            return (xySum / Math.Sqrt(xxSum * yySum));

        }

        /// <summary>
        /// Estimates the covariance of the two variables in the population.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <returns>An estimate, with uncertainty, of the covariance of the two variables in the distribution from which the sample was drawn.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">Are fewer than three entries in the sample.</exception>
        public static UncertainValue PopulationCovariance (this IReadOnlyList<double> x, IReadOnlyList<double> y) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();
            if (x.Count < 3) throw new InsufficientDataException();

            double xm = x.Mean();
            double ym = y.Mean();

            int n = x.Count;
            double c_xy = 0.0;
            double c_xxyy = 0.0;
            for (int i = 0; i < n; i++) {
                double xz = x[i] - xm;
                double yz = y[i] - ym;
                double xy = xz * yz;
                c_xy += xy;
                c_xxyy += xy * xy;
            }
            c_xy /= n;
            c_xxyy /= n;

            return (new UncertainValue(c_xy * n / (n - 1), Math.Sqrt((c_xxyy - c_xy * c_xy) / n)));
        }

        /// <summary>
        /// Computes the best-fit linear regression from the data.
        /// </summary>
        /// <param name="x">The values of the independent variable (x-coordinates).</param>
        /// <param name="y">The values of the dependent variables (y-coordinates).</param>
        /// <returns>The result of the fit.</returns>
        /// <remarks>
        /// <para>Linear regression assumes that the data have been generated by a function y = a + b x + e, where e is
        /// normally distributed noise, and determines the values of a and b that best fit the data. It also
        /// determines a covariance matrix on the parameters a and b, and computes an ANOVA analysis of the fit.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than three data points.</exception>
        public static LinearRegressionResult LinearRegression(this IReadOnlyList<double> y, IReadOnlyList<double> x) {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();
            if (x.Count < 3) throw new InsufficientDataException();

            return (new LinearRegressionResult(x, y));
        }

        /// <summary>
        /// Computes the polynomial of given degree which best fits the data.
        /// </summary>
        /// <param name="x">The values of the independent variable (x-coordinates).</param>
        /// <param name="y">The values of the dependent variables (y-coordinates).</param>
        /// <param name="m">The degree, which must be non-negative.</param>
        /// <returns>The result of the fit.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> is negative.</exception>
        /// <exception cref="InsufficientDataException">There are fewer data points than coefficients to be fit.</exception>
        public static PolynomialRegressionResult PolynomialRegression (this IReadOnlyList<double> y, IReadOnlyList<double> x, int m) {
            return (new PolynomialRegressionResult(x, y, m));
        }


        public static NonlinearRegressionResult NonlinearRegression (
            this IReadOnlyList<double> y,
            IReadOnlyList<double> x,
            Func<IReadOnlyDictionary<string, double>, double, double> function,
            IReadOnlyDictionary<string, double> start
        ) {
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (start == null) throw new ArgumentNullException(nameof(start));
            if (x.Count != y.Count) throw new DimensionMismatchException();

            List<string> names = start.Keys.ToList();
            List<double> startValues = start.Values.ToList();

            Func<IReadOnlyList<double>, double, double> f2 = (a, x2) => {
                Debug.Assert(a.Count == names.Count);
                Dictionary<string, double> parameters = new Dictionary<string, double>(names.Count);
                for (int i = 0; i < a.Count; i++) {
                    string name = names[i];
                    parameters.Add(name, a[i]);
                };
                return (function(parameters, x2));
            };

            return new NonlinearRegressionResult(x, y, f2, startValues, names);
        }

        /// <summary>
        /// Finds the parameterized function that best fits the data.
        /// </summary>
        /// <param name="x">The ordinate values.</param>
        /// <param name="y">The abscissa values.</param>
        /// <param name="function">The parameterized function.</param>
        /// <param name="start">An initial guess for the parameters.</param>
        /// <returns>The fit result.</returns>
        /// <remarks>
        /// <para>
        /// In the returned result, the fit parameters appear in the same order as in
        /// the supplied fit function and initial guess vector.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> or <paramref name="function"/> or <paramref name="start"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The sizes of <paramref name="x"/> and <paramref name="y"/> do not match.</exception>
        /// <exception cref="InsufficientDataException">There are not more data points than fit parameters.</exception>
        /// <exception cref="DivideByZeroException">The curvature matrix is singular, indicating that the data is independent of
        /// one or more parameters, or that two or more parameters are linearly dependent.</exception>
        public static NonlinearRegressionResult NonlinearRegression (
            this IReadOnlyList<double> y,
            IReadOnlyList<double> x,
            Func<IReadOnlyList<double>, double, double> function,
            IReadOnlyList<double> start) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (start == null) throw new ArgumentNullException(nameof(start));
            if (x.Count != y.Count) throw new DimensionMismatchException();

            List<string> names = new List<string>(start.Count);
            for (int i = 0; i < start.Count; i++) names.Add(i.ToString());

            return (new NonlinearRegressionResult(x, y, function, start, names));
        }

        /// <summary>
        /// Computes the best-fit linear logistic regression from the data.
        /// </summary>
        /// <param name="y">The values of the dependent variable.</param>
        /// <param name="x">The values of the independent variable.</param>
        /// <returns>The fit result.</returns>
        /// <remarks>
        /// <para>Linear logistic regression is a way to fit binary outcome data to a linear model.</para>
        /// <para>The method assumes that binary outcomes are encoded as 0 and 1. If any y-values other than
        /// 0 and 1 are encountered, it throws an <see cref="InvalidOperationException"/>.</para>
        /// <para>The fit result is two-dimensional. The first parameter is a, the second b.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The sizes of <paramref name="x"/> and <paramref name="y"/> do not match.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than three data points.</exception>
        public static LinearLogisticRegressionResult LinearLogisticRegression (this IReadOnlyList<bool> y, IReadOnlyList<double> x) {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();
            return (new LinearLogisticRegressionResult(x, y));
        }

        /// <summary>
        /// Performs a Pearson correlation test for association.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
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
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than 3 entries in the sample.</exception>
        /// <seealso cref="SpearmanRhoTest"/>
        /// <seealso cref="KendallTauTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Pearson_correlation_coefficient" />
        public static TestResult PearsonRTest (IReadOnlyList<double> x, IReadOnlyList<double> y) {
            return (PearsonRTest(x, y, TestType.TwoTailed));
        }

        /// <summary>
        /// Performs a Pearson correlation test for association with the given sided-ness.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <param name="type">The sided-ness of the test.</param>
        /// <returns>The result of the test.</returns>
        public static TestResult PearsonRTest (IReadOnlyList<double> x, IReadOnlyList<double> y, TestType type) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();

            int n = x.Count;
            if (n < 3) throw new InsufficientDataException();

            double xMean, yMean, xxSum, yySum, xySum;
            ComputeBivariateMomentsUpToTwo(x, y, out n, out xMean, out yMean, out xxSum, out yySum, out xySum);
            double r = xySum / Math.Sqrt(xxSum * yySum);
            ContinuousDistribution p = new PearsonRDistribution(n);
            return (new TestResult("r", r, type, p));
        }

        /// <summary>
        /// Performs a Spearman rank-order test of association between the two variables.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>The Spearman rank-order test of association is a non-parametric test for association between
        /// two variables. The test statistic rho is the correlation coefficient of the <em>rank</em> of
        /// each entry in the sample. It is thus invariant over monotonic re-parameterizations of the data,
        /// and will, for example, detect a quadratic or exponential association just as well as a linear
        /// association.</para>
        /// <para>The Spearman rank-order test requires O(N log N) operations.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than three data points.</exception>
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="KendallTauTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient"/>
        public static TestResult SpearmanRhoTest (IReadOnlyList<double> x, IReadOnlyList<double> y) {
            return (SpearmanRhoTest(x, y, TestType.TwoTailed));
        }

        /// <summary>
        /// Performs a Spearman rank-order test of association between the two variables with the given sided-ness.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <param name="type">The sided-ness of the test.</param>
        /// <returns>The result of the test.</returns>
        public static TestResult SpearmanRhoTest (IReadOnlyList<double> x, IReadOnlyList<double> y, TestType type) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();

            int n = x.Count;
            if (n < 3) throw new InsufficientDataException();

            // analytic expressions for the mean and variance of ranks
            double M = (n - 1) / 2.0;
            double V = (n + 1) * (n - 1) / 12.0;

            // compute the covariance of ranks
            int[] rx = Univariate.GetRanks(x);
            int[] ry = Univariate.GetRanks(y);
            double C = 0.0;
            for (int i = 0; i < n; i++) {
                C += (rx[i] - M) * (ry[i] - M);
            }
            C = C / n;

            // compute rho
            double rho = C / V;

            // for small enough sample, use the exact distribution
            ContinuousDistribution rhoDistribution;
            if (n < 12) {
                // for small enough sample, use the exact distribution
                // it would be nice to do this for at least slightly higher n, but computation time grows dramatically
                // would like to ensure return in less than 100ms; current timings n=10 35ms, n=11 72ms, n=12 190ms
                rhoDistribution = new DiscreteAsContinuousDistribution(new SpearmanExactDistribution(n), Interval.FromEndpoints(-1.0, 1.0));
            } else {
                // for larger samples, use the normal approximation
                // would like to fit support and C_4 too; look into logit-normal
                // i was not happy with Edgeworth expansion, which can fit C_4 but screws up tails badly, even giving negative probabilities
                rhoDistribution = new NormalDistribution(0.0, 1.0 / Math.Sqrt(n - 1));
            }

            return (new TestResult("ρ", rho, type, rhoDistribution));

        }

        /// <summary>
        /// Performs a Kendall concordance test for association.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>Kendall's &#x3C4; is a non-parametric and robust test of association
        /// between two variables. It simply measures the number of cases where an increase
        /// in one variable is associated with an increase in the other (concordant pairs),
        /// compared with the number of cases where an increase in one variable is associated
        /// with a decrease in the other (discordant pairs).</para>
        /// <para>Because &#x3C4; depends only on the sign
        /// of the difference and not its magnitude, it is not skewed by outliers exhibiting very large
        /// changes, nor by cases where the degree of difference
        /// changes over the ranges of the variables. Of course, it may
        /// still miss an association whose sign changes over the range of the variables. For example,
        /// if data points lie along a semi-circle in the plane, an increase in the first variable
        /// is associated with an increase in the second variable along the rising arc and and decrease in
        /// the second variable along the falling arc.
        /// </para>
        /// <para>Because it examines all pairs of data points, the Kendall test requires
        /// O(N<sup>2</sup>) operations. It is thus impractical for very large data sets. While
        /// not quite as robust as the Kendall test, the Spearman test is a good fall-back in such cases.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than two entries in the sample.</exception>
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="SpearmanRhoTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kendall_tau_test" />
        public static TestResult KendallTauTest (IReadOnlyList<double> x, IReadOnlyList<double> y) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();

            int n = x.Count;
            if (n < 2) throw new InsufficientDataException();

            // loop over all pairs, counting concordant and discordant
            int C = 0;
            int D = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {

                    // note the way each variable varies in the pair
                    int sx = Math.Sign(x[i] - x[j]);
                    int sy = Math.Sign(y[i] - y[j]);

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

            return (new TestResult("τ", t, TestType.TwoTailed, tauDistribution));

        }

        /// <summary>
        /// Performs a paired Student t-test.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>Like a two-sample, unpaired t-test (<see cref="Sample.StudentTTest(Sample,Sample)" />),
        /// a paired t-test compares two samples to detect a difference in means.
        /// Unlike the unpaired version, the paired version assumes that each </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than two data points.</exception>
        public static TestResult PairedStudentTTest (IReadOnlyList<double> x, IReadOnlyList<double> y) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();

            int n = x.Count;
            if (n < 2) throw new InsufficientDataException();

            // the paired t-test is just a normal t-test of one sample against a mean,
            // but where the sample consists of the differences between the paired measurements,
            // and the mean being tested against is (usually) zero

            // loop over pairs, computing mean and standard deviation of differences
            double m = 0.0;
            double v = 0.0;
            for (int i = 0; i < n; i++) {
                double z = x[i] - y[i];
                v += MoreMath.Sqr(z - m) * i / (i + 1);
                m += (z - m) / (i + 1);
            }
            v = v / (n - 1);

            // compute standard error
            double s = Math.Sqrt(v / n);

            // t is the mean deviation as a fraction of standard error
            double t = m / s;

            return (new TestResult("t", t, TestType.TwoTailed, new StudentDistribution(n - 1)));

        }


        /// <summary>
        /// Performs a Wilcoxon signed rank test.
        /// </summary>
        /// <param name="x">The values of the first variable.</param>
        /// <param name="y">The values of the second variable.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>The Wilcoxon signed rank test is a non-parametric alternative to the
        /// paired t-test (<see cref="PairedStudentTTest"/>). Given two measurements on
        /// the same subjects, this method tests for changes in the distribution between
        /// the two measurements. It is sensitive primarily to shifts in the median.
        /// Note that the distributions of the individual measurements
        /// may be far from normal, and may be different for each subject.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException"><paramref name="x"/> and <paramref name="y"/> do not contain the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than two data points.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test"/>
        public static TestResult WilcoxonSignedRankTest (IReadOnlyList<double> x, IReadOnlyList<double> y) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new DimensionMismatchException();

            int n = x.Count;
            if (n < 2) throw new InsufficientDataException();

            double[] z = new double[n];
            for (int i = 0; i < z.Length; i++) {
                z[i] = x[i] - y[i];
            }

            Array.Sort(z, (xValue, yValue) => Math.Abs(xValue).CompareTo(Math.Abs(yValue)));

            int W = 0;
            for (int i = 0; i < z.Length; i++) {
                if (z[i] > 0.0) W += (i + 1);
            }

            ContinuousDistribution nullDistribution;
            if (n < 32) {
                DiscreteDistribution wilcoxon = new WilcoxonDistribution(n);
                nullDistribution = new DiscreteAsContinuousDistribution(wilcoxon);
            } else {
                double mu = n * (n + 1.0) / 4.0;
                double sigma = Math.Sqrt(mu * (2.0 * n + 1.0) / 6.0);
                nullDistribution = new NormalDistribution(mu, sigma);
            }

            return (new TestResult("W", W, TestType.TwoTailed, nullDistribution));

        }

        /// <summary>
        /// Produces a cross-tabulation.
        /// </summary>
        /// <typeparam name="R">The type of row data.</typeparam>
        /// <typeparam name="C">The type of column data</typeparam>
        /// <param name="rowValues">The data values for rows.</param>
        /// <param name="columnValues">The data values for columns.</param>
        /// <returns>A cross-tabular summary of the number of joint occurrences of each row and column value.</returns>
        public static ContingencyTable<R, C> Crosstabs<R, C> (IReadOnlyList<R> rowValues, IReadOnlyList<C> columnValues) {

            if (rowValues == null) throw new ArgumentNullException(nameof(rowValues));
            if (columnValues == null) throw new ArgumentNullException(nameof(columnValues));
            if (rowValues.Count != columnValues.Count) throw new DimensionMismatchException();

            // This is coded as a two passes over the data. The first pass determines
            // the distinct row and column values. The second pass determines the cell counts.
            // It's possible to do this in one pass, but we need auxiliary memory and/or
            // dynamic storage and the details become messier.

            // Determine the distinct row and column values
            NullableDictionary<R, int> rowMap = new NullableDictionary<R, int>();
            NullableDictionary<C, int> columnMap = new NullableDictionary<C, int>();
            for (int i = 0; i < rowValues.Count; i++) {
                R rowValue = rowValues[i];
                C columnValue = columnValues[i];

                if (!rowMap.ContainsKey(rowValue)) rowMap.Add(rowValue, rowMap.Count);
                if (!columnMap.ContainsKey(columnValue)) columnMap.Add(columnValue, columnMap.Count);
            }

            // Fill out the cells and marginal totals
            int[,] counts = new int[rowMap.Count, columnMap.Count];
            int[] rowCounts = new int[rowMap.Count];
            int[] columnCounts = new int[columnMap.Count];
            int totalCount = 0;
            for (int i = 0; i < rowValues.Count; i++) {
                R rowValue = rowValues[i];
                C columnValue = columnValues[i];

                int rowIndex = rowMap[rowValue];
                int columnIndex = columnMap[columnValue];

                counts[rowIndex, columnIndex]++;
                rowCounts[rowIndex]++;
                columnCounts[columnIndex]++;
                totalCount++;
            }

            return new ContingencyTable<R, C>(rowMap, columnMap, counts, rowCounts, columnCounts, totalCount);
        }

    }
}
