using System;
using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.Globalization;

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
    public sealed class MultivariateSample : ICollection<double[]>, IEnumerable<double[]>, IEnumerable {

        /// <summary>
        /// Initializes a multivariate sample.
        /// </summary>
        /// <param name="dimension">The dimension of the sample space, that is, the number of variables
        /// recorded for each sample entry.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="dimension"/> is less than one.</exception>
        public MultivariateSample (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException("dimension");
            storage = new SampleStorage[dimension];
            for (int j = 0; j < storage.Length; j++) {
                storage[j] = new SampleStorage();
            }
            isReadOnly = false;
        }

        internal MultivariateSample (SampleStorage[] columns, bool isReadOnly) {
            this.storage = columns;
            this.isReadOnly = isReadOnly;
        }

        private List<double[]> data = new List<double[]>();

        private SampleStorage[] storage;
        private bool isReadOnly;

        private int IndexOf (IList<double> value) {
            for (int i = 0; i < Count; i++) {
                if (Matches(value, i)) return (i);
            }
            return (-1);
        }

        private bool Matches (IList<double> value, int i) {
            for (int j = 0; j < Dimension; j++) {
                if (storage[j][i] != value[j]) return (false);
            }
            return (true);
        }

        /// <summary>
        /// Gets the dimension of the sample.
        /// </summary>
        public int Dimension {
            get {
                return (storage.Length);
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
            if (value.Count != Dimension) throw new DimensionMismatchException();
            if (isReadOnly) throw new InvalidOperationException();

            for (int i = 0; i < Dimension; i++) {
                storage[i].Add(value[i]);
            }

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
            if (values.Count != Dimension) throw new DimensionMismatchException();
            if (isReadOnly) throw new InvalidOperationException();

            int i = IndexOf(values);
            if (i < 0) {
                return (false);
            } else {
                for (int j = 0; j < Dimension; j++) {
                    storage[j].RemoveAt(i);
                }
                return (true);
            }

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
            if (value == null) throw new ArgumentNullException("value");
            if (value.Count != Dimension) throw new DimensionMismatchException();
            return(IndexOf(value) >= 0); 
        }

        /// <summary>
        /// Removes all entries from the sample.
        /// </summary>
        public void Clear () {
            if (isReadOnly) {
                throw new InvalidOperationException();
            } else {
                for (int j = 0; j < Dimension; j++) {
                    storage[j].Clear();
                }
            }
        }

        /// <summary>
        /// Gets the number of enties in the sample.
        /// </summary>
        public int Count {
            get {
                return (storage[0].Count);
            }
        }

        /// <summary>
        /// Gets the indicated column as a univariate <see cref="Sample"/>.
        /// </summary>
        /// <param name="c">The (zero-based) column index.</param>
        /// <returns>A read-only <see cref="Sample"/> containing all values in the indicated column.</returns>
        /// <remarks>
        /// <para>Use this method to obtain column-specific information, such as the <see cref="Sample.Median"/> or
        /// <see cref="Sample.Variance"/> of the column.</para>
        /// <para>Note that this is a fast, O(1) operation, which does not create an independent copy of the column.
        /// The advantage of this is that you can access columns as independent samples as often as you like without
        /// worying about performance. The disadvantage of this is that the returned sample cannot be altered. If you
        /// need to alter values in a column independent of the multi-variate sample, use the <see cref="Sample.Copy"/>
        /// method to obtain an independent copy of the column.</para>
        /// </remarks>
        public Sample Column (int c) {
            if ((c < 0) || (c >= Dimension)) throw new ArgumentOutOfRangeException("c");
            return (new Sample(storage[c], true));
        }

        /// <summary>
        /// Gets the indicated columns as a <see cref="BivariateSample"/>.
        /// </summary>
        /// <param name="cx">The (zero-based) column index of the X variable.</param>
        /// <param name="cy">The (zero-based) column index of the Y variable.</param>
        /// <returns>A read-only <see cref="BivariateSample"/> consisting of the indicated columns..</returns>
        /// <remarks>
        /// <para>Use this method to obtain information specific to the two columns, such as the <see cref="BivariateSample.Covariance"/>,
        /// or to perform tests specific to the two columns, such as a <see cref="BivariateSample.PearsonRTest"/>.</para>
        /// <para>Note that this is a fast, O(1) operation, which does not create independent copies of the columns.
        /// The advantage of this is that you can access pairs of columns as bivariate samples as often as you like without
        /// worying about performance. The disadvantage of this is that the returned bivariate sample cannot be altered. If you
        /// need to alter values independent of the multi-variate sample, use the <see cref="BivariateSample.Copy"/>
        /// method to obtain an independent copy of the bivariate sample.</para>
        /// </remarks>
        public BivariateSample TwoColumns (int cx, int cy) {
            if ((cx < 0) || (cx >= Dimension)) throw new ArgumentOutOfRangeException("cx");
            if ((cy < 0) || (cy >= Dimension)) throw new ArgumentOutOfRangeException("cy");
            return (new BivariateSample(storage[cx], storage[cy], true));
        }

        /// <summary>
        /// Gets the indicated columns as a multivariate sample.
        /// </summary>
        /// <param name="columnIndexes">A list of column indexes.</param>
        /// <returns>A read-only <see cref="MultivariateSample"/> consisting of the indicated columns.</returns>
        /// <remarks>
        /// <para>Use this method to perform multivariate analyses, such as regression and principal component analyis, using
        /// only a subset of the variables in the original multivariate sample.</para>
        /// <para>Note that this is a fast operation, which does not create independent copies of the columns.</para>
        /// </remarks>
        public MultivariateSample Columns (IList<int> columnIndexes) {
            if (columnIndexes == null) throw new ArgumentNullException("columnIndexes");
            SampleStorage[] columns = new SampleStorage[columnIndexes.Count];
            for (int i = 0; i < columns.Length; i++) {
                int ci = columnIndexes[i];
                if ((ci < 0) || (ci >= Dimension)) throw new ArgumentOutOfRangeException();
                columns[i] = storage[ci];
            }
            return new MultivariateSample(columns, true);
        }

        /// <summary>
        /// Gets the indicated columns as a multivariate sample.
        /// </summary>
        /// <param name="columnIndexes">A list of columns indexes.</param>
        /// <returns>A read-only <see cref="MultivariateSample"/> consisting of the indicated columns.</returns>
        public MultivariateSample Columns (params int[] columnIndexes) {
            return (Columns((IList<int>)columnIndexes));
        }

        /// <summary>
        /// Copies the multivariate sample.
        /// </summary>
        /// <returns>An independent copy of the multivariate sample.</returns>
        public MultivariateSample Copy () {
            SampleStorage[] columnCopies = new SampleStorage[storage.Length];
            for (int i = 0; i < columnCopies.Length; i++) {
                columnCopies[i] = storage[i].Copy();
            }
            return (new MultivariateSample(columnCopies, false));
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
        /// <exception cref="ArgumentNullException"><paramref name="powers"/> is null.</exception>
        /// <exception cref="DimensionMismatchException">The length of <paramref name="powers"/> is not
        /// equal to the <see cref="Dimension"/> of the multivariate sample.</exception>
        public double Moment (IList<int> powers) {
            if (powers == null) throw new ArgumentNullException("powers");
            if (powers.Count != Dimension) throw new DimensionMismatchException();

            double M = 0.0;
            for (int i = 0; i < Count; i++) {
                double t = 1.0;
                for (int j = 0; j < Dimension; j++) {
                    t *= MoreMath.Pow(storage[j][i], powers[j]);
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
        /// <exception cref="ArgumentNullException"><paramref name="powers"/> is null.</exception>
        /// <exception cref="DimensionMismatchException">The length of <paramref name="powers"/> is not
        /// equal to the <see cref="Dimension"/> of the multivariate sample.</exception>
        public double MomentAboutMean (IList<int> powers) {
            if (powers == null) throw new ArgumentNullException("powers");
            if (powers.Count != Dimension) throw new DimensionMismatchException();

            double C = 0.0;
            for (int i = 0; i < Count; i++) {
                double t = 1.0;
                for (int j = 0; j < Dimension; j++) {
                    t *= MoreMath.Pow(storage[j][i] - storage[j].Mean, powers[j]);
                }
                C += t;
            }
            C = C / Count;

            return (C);

        }

        /// <summary>
        /// Computes the given sample central moment.
        /// </summary>
        /// <param name="powers">The power to which each component should be raised.</param>
        /// <returns>The specified moment.</returns>
        public double MomentAboutMean (params int[] powers) {
            return (MomentAboutMean((IList<int>) powers));
        }

#if PAST

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
#endif

        // the internal linear regression routine, which assumes inputs are entirely valid

        private FitResult LinearRegression_Internal (int outputIndex) {

            // to do a fit, we need more data than parameters
            if (Count < Dimension) throw new InsufficientDataException();

            // construct the design matrix
            SymmetricMatrix D = new SymmetricMatrix(Dimension);
            for (int i = 0; i < Dimension; i++) {
                for (int j = 0; j <= i; j++) {
                    if (i == outputIndex) {
                        if (j == outputIndex) {
                            D[i, j] = Count;
                        } else {
                            D[i, j] = storage[j].Mean * Count;
                        }
                    } else {
                        if (j == outputIndex) {
                            D[i, j] = storage[i].Mean * Count;
                        } else {
                            double Dij = 0.0;
                            for (int k = 0; k < Count; k++) {
                                Dij += storage[i][k] * storage[j][k];
                            }
                            D[i, j] = Dij;
                        }
                    }
                }
            }

            // construct the right hand side
            ColumnVector b = new ColumnVector(Dimension);
            for (int i = 0; i < Dimension; i++) {
                if (i == outputIndex) {
                    b[i] = storage[i].Mean * Count;
                } else {
                    double bi = 0.0;
                    for (int k = 0; k < Count; k++) {
                        bi += storage[outputIndex][k] * storage[i][k];
                    }
                    b[i] = bi;
                }
            }

            // solve the system for the linear model parameters
            CholeskyDecomposition CD = D.CholeskyDecomposition();
            ColumnVector parameters = CD.Solve(b);

            // find total sum of squares, with dof = # points - 1 (minus one for the variance-minimizing mean)
            double totalSumOfSquares = storage[outputIndex].Variance * Count;

            // find remaining unexplained sum of squares, with dof = # points - # parameters
            double unexplainedSumOfSquares = 0.0;
            for (int r = 0; r < Count; r++) {
                double y = 0.0;
                for (int c = 0; c < Dimension; c++) {
                    if (c == outputIndex) {
                        y += parameters[c];
                    } else {
                        y += parameters[c] * storage[c][r];
                    }
                }
                unexplainedSumOfSquares += MoreMath.Pow2(y - storage[outputIndex][r]);
            }
            int unexplainedDegreesOfFreedom = Count - Dimension;
            double unexplainedVariance = unexplainedSumOfSquares / unexplainedDegreesOfFreedom;

            // find explained sum of squares, with dof = # parameters - 1
            double explainedSumOfSquares = totalSumOfSquares - unexplainedSumOfSquares;
            int explainedDegreesOfFreedom = Dimension - 1;
            double explainedVariance = explainedSumOfSquares / explainedDegreesOfFreedom;

            // compute F statistic from sums of squares
            double F = explainedVariance / unexplainedVariance;
            Distribution fDistribution = new FisherDistribution(explainedDegreesOfFreedom, unexplainedDegreesOfFreedom);

            SymmetricMatrix covariance = unexplainedVariance * CD.Inverse();

            return (new FitResult(parameters, covariance, new TestResult(F, fDistribution)));

        }


        /// <summary>
        /// Performs a linear regression analysis.
        /// </summary>
        /// <param name="outputIndex">The index of the variable to be predicted.</param>
        /// <returns>The result of the regression.</returns>
        /// <remarks>
        /// <para>Linear regression finds the linear combination of the other variables
        /// that best predicts the output variable.</para>
        /// <img src="../images/LinearRegressionEquation.png" />
        /// <para>The noise term epsilon is assumed to be drawn from the same normal distribution
        /// for each data point. Note that the model makes no assumptions about the distribution
        /// of the x's; it merely asserts a particular underlying relationship between the x's and
        /// the y.</para>
        /// <h4>Inputs and Outputs</h4>
        /// <para>In the returned fit result, the indices of the parameters correspond to indices
        /// of the coefficients. The intercept parameter has the index of the output variable.
        /// Thus if a linear regression analaysis is done on a 4-dimensional multivariate sample
        /// to predict variable number 2, the coefficients of variables 0, 1, and 3 will be
        /// parameters 0, 1, and 3, of the returned fit result, and the intercept will be
        /// parameter 2.</para>
        /// <para>If you want to include fewer input variables in your regression, use the
        /// <see cref="Columns(IList{Int32})"/> method to create a multivariate sample that includes only
        /// the variables you want to use in your regression.</para>
        /// <para>The correlation matrix among fit parameters is also returned with the fit
        /// result, as is an F-test for the goodness of the fit. If the result of the F-test
        /// is not significant, no conclusions should be drawn from the regression
        /// coefficients.</para>
        /// <h4>Regression vs. Correlation</h4>
        /// <para>If a given coefficient is significantly positive, then a change in the value
        /// of the corresponding input variable, <i>holding all other input variables constant</i>,
        /// will tend to increase the output variable. Note that italicized condition means
        /// that, when there is more than one input variable, a linear regression coefficient measures
        /// something different than a linear correlation coefficient.</para>
        /// <para>Suppose, for example, we take a large number of measurements
        /// of water temperature, plankton concentration, and fish density in a large number of
        /// different locations. A simple correlation analysis might indicate that fish density is
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
        /// tease out the independent effect of each. In that case, a fit using only one
        /// of the variables will be as good or nearly as good as a fit using both, and
        /// the covariance between their corresponding linear regression coefficients will
        /// be large. In a situation like this, you should be wary of drawing any conclusions
        /// about their seperate effects.</para>
        /// <para>It can also occur that an input variable or a set of input variables is indeed good
        /// predictor of an output variable, but via a complex and non-linear relationship that
        /// a linear regression analysis will completely miss.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="outputIndex"/> is outside the range of allowed indexes.</exception>
        /// <exception cref="InsufficientDataException">There are fewer entries than the dimension of the multivariate sample.</exception>
        /// <seealso cref="LinearRegression(int)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Linear_regression"/>
        public FitResult LinearRegression (int outputIndex) {

            if ((outputIndex < 0) || (outputIndex >= Dimension)) throw new ArgumentOutOfRangeException("outputIndex");

            return (LinearRegression_Internal(outputIndex));

        }

        // interface implementations

        /// <summary>
        /// Gets an enumerator over the sample entries.
        /// </summary>
        /// <returns>An iterator over the sample entries.</returns>
        public IEnumerator<double[]> GetEnumerator () {
            for (int r = 0; r < Count; r++) {
                double[] entry = new double[Dimension];
                for (int c = 0; c < Dimension; c++) {
                    entry[c] = storage[c][r];
                }
                yield return(entry);
            }
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return(GetEnumerator());
        }

        void ICollection<double[]>.CopyTo (double[][] array, int start) {
            data.CopyTo(array, start);
        }

        /// <summary>
        /// Gets a value indicating whether the multivariate sample can be modified.
        /// </summary>
        public bool IsReadOnly {
            get {
                return (isReadOnly);
            }
        }

        /// <summary>
        /// Performs a principal component analysis of the data.
        /// </summary>
        /// <returns>The result of the principal component analysis.</returns>
        /// <exception cref="InsufficientDataException">The number of data entries (<see cref="Count"/>) is
        /// less than the number of variables (<see cref="Dimension"/>).</exception>
        /// <seealso cref="PrincipalComponentAnalysis"/>
        public PrincipalComponentAnalysis PrincipalComponentAnalysis () {

            if (Count < Dimension) throw new InsufficientDataException();

            // construct a (Count X Dimension) matrix of mean-centered data
            double[] store = MatrixAlgorithms.AllocateStorage(Count, Dimension);
            int i = 0;
            for (int c = 0; c < Dimension; c++) {
                double mu = storage[c].Mean;
                for (int r = 0; r < Count; r++) {
                    store[i] = storage[c][r] - mu;
                    i++;
                }
            }

            // bidiagonalize it
            double[] a, b;
            MatrixAlgorithms.Bidiagonalize(store, Count, Dimension, out a, out b);

            // form the U and V matrices
            double[] left = MatrixAlgorithms.AccumulateBidiagonalU(store, Count, Dimension);
            double[] right = MatrixAlgorithms.AccumulateBidiagonalV(store, Count, Dimension);

            // find the singular values of the bidiagonal matrix
            MatrixAlgorithms.ExtractSingularValues(a, b, left, right, Count, Dimension);

            // sort them
            MatrixAlgorithms.SortValues(a, left, right, Count, Dimension);

            PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis(left, a, right, Count, Dimension);

            return (pca);
            
        }

        /// <summary>
        /// Loads values from a data reader.
        /// </summary>
        /// <param name="reader">The data reader.</param>
        /// <param name="dbIndexes">The database column indexes of the sample columns.</param>
        public void Load (IDataReader reader, IList<int> dbIndexes) {
            if (reader == null) throw new ArgumentNullException("reader");
            if (dbIndexes == null) throw new ArgumentNullException("dbIndexes");
            if (dbIndexes.Count != Dimension) throw new DimensionMismatchException();
            if (isReadOnly) throw new InvalidOperationException();

                // create an array to store values, which we will re-use as we move through the data
                double[] entry = new double[Dimension];
                // move through the data
                while (reader.Read()) {
                    // check each entry and, if value, add it to the sample
                    if (ReadValues(reader, dbIndexes, entry)) Add(entry);
                }

        }

        private bool ReadValues (IDataReader reader, IList<int> dbIndexes, double[] entry) {
            for (int c = 0; c < Dimension; c++) {
                int i = dbIndexes[c];
                if (reader.IsDBNull(i)) {
                    return (false);
                } else {
                    entry[c] = Convert.ToDouble(reader.GetValue(i), CultureInfo.InvariantCulture);
                }
            }
            return (true);
        }

        /// <summary>
        /// Loads values from a data reader.
        /// </summary>
        /// <param name="reader">The data reader.</param>
        /// <param name="dbIndexes">The database column indexes of the sample columns.</param>
        public void Load (IDataReader reader, params int[] dbIndexes) {
            Load(reader, (IList<int>)dbIndexes);
        }

    }

}
