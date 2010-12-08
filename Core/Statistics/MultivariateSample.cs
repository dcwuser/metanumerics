using System;
using System.Collections;
using System.Collections.Generic;

using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics.Distributions;

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

            if (value == null) throw new ArgumentNullException("value");
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
        /// <param name="values">The values associated with the entry to remove.</param>
        /// <returns>Whether the entry was found and removed.</returns>
        public bool Remove (IList<double> values) {

            if (values == null) throw new ArgumentNullException("values");
            if (values.Count != n) throw new DimensionMismatchException();

            for (int i = 0; i < Count; i++) {
                if (Matches(data[i], values)) {
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

        /// <summary>
        /// Computes the given sample raw moment.
        /// </summary>
        /// <param name="powers">The power to which each component should be raised.</param>
        /// <returns>The specified moment.</returns>
        public double Moment (params int[] powers) {
            return (Moment((IList<int>) powers));

        }

        /// <summary>
        /// Computes the given sample raw moment.
        /// </summary>
        /// <param name="powers">The power to which each component should be raised.</param>
        /// <returns>The specified moment.</returns>
        public double Moment (IList<int> powers) {
            if (powers == null) throw new ArgumentNullException("powers");
            if (powers.Count != Dimension) throw new DimensionMismatchException();

            double M = 0.0;
            for (int i = 0; i < data.Count; i++) {
                double t = 1.0;
                for (int j = 0; j < powers.Count; j++) {
                    t *= MoreMath.Pow(data[i][j], powers[j]);
                }
                M += t;
            }
            M = M / Count;

            return (M);

        }

        /// <summary>
        /// Computes the given sample central moment.
        /// </summary>
        /// <param name="powers">The power to which each component should be raised.</param>
        /// <returns>The specified moment.</returns>
        public double MomentAboutMean (params int[] powers) {
            return (MomentAboutMean((IList<int>) powers));

        }

        /// <summary>
        /// Computes the given sample central moment.
        /// </summary>
        /// <param name="powers">The power to which each component should be raised.</param>
        /// <returns>The specified moment.</returns>
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
                    t *= MoreMath.Pow(data[i][j] - means[j], powers[j]);
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
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>This test measures the strength of the linear correlation between two variables. The
        /// test statistic r is simply the covariance of the two variables, scaled by their respective
        /// standard deviations so as to obtain a number between -1 (perfect linear anti-correlation)
        /// and +1 (perfect linear correlation).</para>
        /// <para>The Pearson test cannot reliably detect or rule out non-linear correlations.</para>
        /// <para>The Pearson correlation test requires O(N) operations.</para>
        /// <para>The Pearson test requires at least three bivariate values.</para>
        /// </remarks>
        /// <seealso cref="SpearmanRhoTest"/>
        /// <seealso cref="KendallTauTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Pearson_correlation_coefficient" />
        public TestResult PearsonRTest (int d1, int d2) {
            if (Count < 3) throw new InvalidOperationException();
            double r = Covariance(d1, d2) / StandardDeviation(d1) / StandardDeviation(d2);
            Distribution p = new PearsonRDistribution(Count);
            //Distribution p = new NormalDistribution(0.0, 1.0 / Math.Sqrt(Count));
            return (new TestResult(r, p));
        }

        /// <summary>
        /// Performs a Spearman rank-order test of association between two variables.
        /// </summary>
        /// <param name="d1">The (zero-based) index of the first variable.</param>
        /// <param name="d2">The (zero-based) index of the second variable.</param>
        /// <returns>The result of the test.</returns>
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
            if (Count < 2) throw new InvalidOperationException();

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
        /// <seealso cref="PearsonRTest"/>
        /// <seealso cref="SpearmanRhoTest"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kendall_tau_test" />
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

        /// <summary>
        /// Performs a paired Student t-test.
        /// </summary>
        /// <param name="d1">The column representing the first group.</param>
        /// <param name="d2">The column representing the second group.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>Like a two-sample, unpaired t-test (<see cref="Sample.StudentTTest(Sample,Sample)" />),
        /// a paired t-test compares two samples to detect a difference in means.
        /// Unlike the unpaired version, the paired version assumes that each </para>
        /// </remarks>
        public TestResult PairedStudentTTest (int d1, int d2) {

            if ((d1 < 0) || (d1 >= n)) throw new ArgumentOutOfRangeException("d1");
            if ((d2 < 0) || (d2 >= n)) throw new ArgumentOutOfRangeException("d2");
            if (Count < 2) throw new InvalidOperationException();

            // the paired t-test is just a normal t-test of one sample against a mean,
            // but where the sample consists of the differences between the paired measurements,
            // and the mean being tested against is (usually) zero

            // loop over pairs, computing mean and standard deviation of differences
            double m = 0.0;
            double v = 0.0;
            for (int i = 0; i < Count; i++) {
                double z = data[i][d1] - data[i][d2];
                v += MoreMath.Pow2(z - m) * i / (i + 1);
                m += (z - m) / (i + 1);
            }
            v = v / (Count-1);
            //Console.WriteLine("m={0} v={1}", m, v);

            // compute standard error
            double s = Math.Sqrt(v / Count);

            // t is the mean deviation as a fraction of standard error
            double t = m / s;

            return (new TestResult(t, new StudentDistribution(Count-1)));

        }

        /// <summary>
        /// Performs a linear regression analysis using the input variables to predict the output variable.
        /// </summary>
        /// <param name="inputIndexes">The indices of the input variables to use.</param>
        /// <param name="outputIndex">The index of the variable to predict.</param>
        /// <returns>The result of the regression.</returns>
        /// <remarks>
        /// <para>Linear regression finds the linear combination of the given input variables
        /// that best predicts the given output variable.</para>
        /// <img src="../images/LinearRegressionEquation.png" />
        /// <para>The noise term epsilon is assumed to be drawn from the same normal distribution
        /// for each data point. Note that the model makes no assumptions about the distribution
        /// of the x's; it merely asserts a particular underlying relationship between the x's and
        /// the y.</para>
        /// <h4>Inputs and Outputs</h4>
        /// <para>The <paramref name="inputIndexes"/> array specifies the input variables
        /// as particular columns of the multivariate sample. The <paramref name="outputIndex" />
        /// specifies column containing the output variable to predict as a linear combination
        /// of input variables. The <paramref name="outputIndex" /> may not appear in the
        /// <paramref name="inputIndexes"/> array. By leaving other indices out of the
        /// <paramref name="inputIndexes"/> array, you can perform a regression fit using a
        /// reduced set of input variables. Of course, all indices must correspond to an
        /// existing column and indices may not be repeated in the <paramref name="inputIndexes"/>
        /// array.</para>
        /// <para>The parameters of the returned fit result are the liklihood-maximizing values of
        /// the coefficients, in the order specified by the <paramref name="inputIndexes"/> array,
        /// followed by the intercept parameter. For example, given
        /// <paramref name="inputIndexes"/> = (1, 4, 3), the fit result parameter list would
        /// be (&#x3B2;<sub>1</sub>, &#x3B2;<sub>4</sub>, &#x3B2;<sub>3</sub>, &#x3B1;).</para>
        /// <para>The correlation matrix among fit parameters is also returned with the fit
        /// result, as is an F-test for the goodness of the fit. If the result of the F-test
        /// is not significant, no conclusions should be drawn from the regression
        /// coefficients.</para>
        /// <h4>Regression vs. Correlation</h4>
        /// <para>If a given coefficient is significantly positive, then a change in the value
        /// of the corresponding input variable, <i>holding all other input variables constant</i>,
        /// will tend to increase the output variable. Note that italicized condition means
        /// that a linear regression coefficient measures something different than a linear
        /// correlation coefficient.</para>
        /// <para>Suppose, for example, we take a large number of measurements
        /// of water temperature, plankton concentration, and fish density in a large number of
        /// different locations. A correlation analysis might indicate that fish density is
        /// positively correlated with both water temperature and plankton concentration. But
        /// a regression analysis might reveal that increasing water temperature actually
        /// decreases the fish density. This seeming paradoxical situation might occur
        /// because fish do much better with more plankton, and plankton do much better at higher
        /// temperatures, and this positive knock-on effect of temperature on fish is larger than the
        /// negative direct effect.</para>
        /// <para>If we are in a situation where we can control the input
        /// variables independently -- for example we are running an aquarium -- we would ceratainly
        /// want to know the specific effect of one variable -- that our fishes would actually prefer us to turn
        /// down the temperature while maintaining a high plankton level -- rather than the observed
        /// effect as a variable changes along with all the others that tend to change with it.
        /// This does not mean that the correlation analysis is wrong -- higher temperatures
        /// are indeed associated with higher fish densitites in our hypothetical data set. It
        /// simply means that you need to be careful to ask the right question for your purpose.</para>
        /// <para>In most cases, it is indeed the specific effect of one variable when others are
        /// held constant that we seek. In a controlled experiment, the confounding effects of
        /// other variables are removed by the experimental design, either by random
        /// assignment or specific controls. In an observational experiment, though,
        /// confounding effects can be, and often are, large, and correlation analysis is
        /// not sufficient. It is worthwhile keeping this in find in politically charged
        /// debates in which easily observed correlations are likely to be bandied about as
        /// evidence, while a more difficult regression analysis that would actually be required
        /// to support an assertion is left undone.</para>
        /// <h4>Cavets</h4>
        /// <para>It can occur that two theoretically independent variables are so closely
        /// correlated in the observational data that a regression analsysis cannot reliably
        /// tease out the independent effect of each. In the case, a fit using only one
        /// of the variables will be as good or nearly as good as a fit using both, and
        /// the covariance between their corresponding linear regression coefficients will
        /// be large. In a situation like this, you should be wary of drawing any conclusions
        /// about their seperate effects.</para>
        /// <para>It can also occur that an input variable or a set of input variables is indeed good
        /// predictor of an output variable, but via a complex and non-linear relationship that
        /// a linear regression analysis will completely miss.</para>
        /// </remarks>
        /// <seealso cref="LinearRegression(int)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Linear_regression"/>
        public FitResult LinearRegression (IList<int> inputIndexes, int outputIndex) {
            // check inputs
            if (inputIndexes == null) throw new ArgumentNullException("inputIndexes");
            for (int i = 0; i < inputIndexes.Count; i++) {
                if ((inputIndexes[i] < 0) || (inputIndexes[i] >= n)) throw new ArgumentOutOfRangeException("inputIndexes");
                if (inputIndexes[i] == outputIndex) throw new InvalidOperationException();
                for (int j = 0; j < i; j++) {
                    if (inputIndexes[j] == inputIndexes[i]) throw new InvalidOperationException();
                }
            }
            if ((outputIndex < 0) || (outputIndex >= n)) throw new ArgumentOutOfRangeException("outputIndex");
            // do the fit
            return (LinearRegression_Internal(inputIndexes, outputIndex));
        }

        // the internal linear regression routine, which assumes inputs are entirely valid

        private FitResult LinearRegression_Internal (IList<int> inputIndices, int outputIndex) {

            // to do a fit, we need more data than parameters
            if (data.Count <= (inputIndices.Count + 1)) throw new InvalidOperationException();

            // construct the design matrix
            SymmetricMatrix D = new SymmetricMatrix(inputIndices.Count + 1);
            for (int i = 0; i < inputIndices.Count; i++) {
                for (int j = 0; j <= i; j++) {
                    double Dij = 0.0;
                    for (int k = 0; k < data.Count; k++) {
                        Dij += data[k][inputIndices[i]] * data[k][inputIndices[j]];
                    }
                    D[i,j] = Dij;
                }
                D[inputIndices.Count, i] = means[inputIndices[i]] * data.Count;
            }
            D[inputIndices.Count, inputIndices.Count] = data.Count;

            // construct the right hand side
            ColumnVector b = new ColumnVector(n);
            for (int i = 0; i < inputIndices.Count; i++) {
                double bi = 0.0;
                for (int k = 0; k < data.Count; k++) {
                    bi += data[k][outputIndex] * data[k][inputIndices[i]];
                }
                b[i] = bi;
            }
            b[inputIndices.Count] = means[outputIndex] * data.Count;

            // solve the system for the linear model parameters
            CholeskyDecomposition CD = D.CholeskyDecomposition();
            ColumnVector parameters = CD.Solve(b);

            // determine total variance (variance before fit)
            double m = means[outputIndex];
            double totalVarianceSum = 0.0;
            for (int k = 0; k < data.Count; k++) {
                double z = (data[k][outputIndex] - means[outputIndex]);
                totalVarianceSum += z * z;
            }
            // associated dof = # points - 1 (the one is for the variance-minimizing mean)

            // determine unexmplained remaining variance (variance after fit)
            double unexplainedVarianceSum = 0.0;
            for (int k = 0; k < data.Count; k++) {
                double y = parameters[inputIndices.Count];
                for (int i = 0; i < inputIndices.Count; i++) {
                    y += parameters[i] * data[k][inputIndices[i]];
                }
                double z = data[k][outputIndex] - y;
                unexplainedVarianceSum += z * z;
            }
            // associated dof = # points - # paramaters
            int unexplainedVarianceDof = data.Count - (inputIndices.Count + 1);
            // sigma-squared for the model error term is given by the unexplained variance
            double unexplainedVariance = unexplainedVarianceSum / unexplainedVarianceDof;


            // explained variance = total variance - unexplained variance
            double explainedVarianceSum = totalVarianceSum - unexplainedVarianceSum;
            // associated dof = (# points - 1) - (# points - # parameters) = # parameters - 1
            int explainedVarianceDof = inputIndices.Count;
            double explainedVariance = explainedVarianceSum / explainedVarianceDof;

            // F statistic = explainedVariance / unexplainedVariance
            double F = explainedVariance / unexplainedVariance;
            Distribution FDistribution = new FisherDistribution(explainedVarianceDof, unexplainedVarianceDof);
            TestResult test = new TestResult(F, FDistribution);

            // covariance matrix is proportional to inverse of design matrix
            SymmetricMatrix covariance = unexplainedVariance * CD.Inverse();

            return (new FitResult(parameters, covariance, test));

        }

        /// <summary>
        /// Performs a linear regression analysis.
        /// </summary>
        /// <param name="outputIndex">The index of the variable to be predicted.</param>
        /// <returns>The result of the analysis.</returns>
        /// <remarks>
        /// <para>This method returns the parameter of a linear model which uses the values
        /// of all other columns to predict the output column.</para>
        /// <para>In the returned fit result, the best-fit linear regression coefficients
        /// are returned in the order of the columns, followed by the intercept parameter.
        /// For example, if <paramref name="outputIndex"/>=2 in a 4-column multivariate
        /// sample, the returned parameters are (&#x3B2;<sub>0</sub>, &#x3B2;<sub>1</sub>,
        /// &#x3B2;<sub>3</sub>, &#x3B1;).</para>
        /// <para>To limit the set of columns used for input variables, and for moreextensive
        /// information on linear regression analysis, see <see cref="LinearRegression(IList&lt;int&gt;,int)" />.</para>
        /// </remarks>
        /// <seealso cref="LinearRegression(IList&lt;int&gt;,int)" />
        public FitResult LinearRegression (int outputIndex) {

            if ((outputIndex < 0) || (outputIndex >= n)) throw new ArgumentOutOfRangeException("outputIndex");

            List<int> inputIndices = new List<int>(n - 1);
            for (int i = 0; i < n; i++) {
                if (i != outputIndex) inputIndices.Add(i);
            }

            return (LinearRegression_Internal(inputIndices, outputIndex));

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
