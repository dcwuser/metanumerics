﻿using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;

using Meta.Numerics.Analysis;
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
    /// the associations between the recorded variables, and routines for fitting the sample to a model.</para>
    /// <para>NOTE: This class will be retired in a future release. Its functionality is replaced by the
    /// <see cref="Multivariate"/> class, which can operatate on arbitrary data lists.</para>
    /// </remarks>
    public sealed class MultivariateSample : ICollection<double[]>, IEnumerable<double[]>, IEnumerable {

        /// <summary>
        /// Initializes a new multivariate sample.
        /// </summary>
        /// <param name="dimension">The dimension of the sample space, that is, the number of variables
        /// recorded for each sample entry.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="dimension"/> is less than one.</exception>
        public MultivariateSample (int dimension) {
            if (dimension < 1) throw new ArgumentOutOfRangeException(nameof(dimension));
            storage = new SampleStorage[dimension];
            for (int j = 0; j < storage.Length; j++) {
                storage[j] = new SampleStorage();
            }
            isReadOnly = false;
        }

        /// <summary>
        /// Initializes a new multivariate sample with the given variable names.
        /// </summary>
        /// <param name="names">The names of the variables.</param>
        public MultivariateSample (params string[] names) : this(names.Length) {
            for (int j = 0; j < storage.Length; j++) {
                storage[j].Name = names[j];
            }
        }

        internal MultivariateSample (SampleStorage[] columns, bool isReadOnly) {
            this.storage = columns;
            this.isReadOnly = isReadOnly;
        }

        private List<double[]> data = new List<double[]>();

        private SampleStorage[] storage;
        private bool isReadOnly;

        private int IndexOf (IReadOnlyList<double> value) {
            for (int i = 0; i < Count; i++) {
                if (Matches(value, i)) return (i);
            }
            return (-1);
        }

        private bool Matches (IReadOnlyList<double> value, int i) {
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
        /// <param name="values">The values associated with the entry.</param>
        public void Add (params double[] values) {
            Add((IReadOnlyList<double>) values);
        }

        /// <summary>
        /// Adds an entry to the sample.
        /// </summary>
        /// <param name="values">The values associated with the entry.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="values"/> does not match the
        /// dimension of the multivariate sample.</exception>
        /// <exception cref="InvalidOperationException">The multivariate sample is read-only.</exception>
        public void Add (IReadOnlyList<double> values) {

            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Count != Dimension) throw new DimensionMismatchException();
            if (isReadOnly) throw new InvalidOperationException();

            for (int i = 0; i < Dimension; i++) {
                storage[i].Add(values[i]);
            }

        }

        /// <summary>
        /// Removes an entry from the sample.
        /// </summary>
        /// <param name="value">The values associated with the entry to remove.</param>
        /// <returns>Whether the entry was found and removed.</returns>
        public bool Remove (params double[] value) {
            return (Remove((IReadOnlyList<double>) value));
        }

        /// <summary>
        /// Removes an entry from the sample.
        /// </summary>
        /// <param name="values">The values associated with the entry to remove.</param>
        /// <returns>Whether the entry was found and removed.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="values"/> does not match the
        /// dimension of the multivariate sample.</exception>
        /// <exception cref="InvalidOperationException">The multivariate sample is read-only.</exception>
        public bool Remove (IReadOnlyList<double> values) {

            if (values == null) throw new ArgumentNullException(nameof(values));
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
        /// <param name="values">The values associated with the entry to search for.</param>
        /// <returns>Whether the sample contains the given entry.</returns>
        public bool Contains (params double[] values) {
            return (Contains((IReadOnlyList<double>) values));
        }

        /// <summary>
        /// Determines whether the sample contains a given entry.
        /// </summary>
        /// <param name="values">The values associated with the entry to search for.</param>
        /// <returns>Whether the sample contains the given entry.</returns>
        public bool Contains (IReadOnlyList<double> values) {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Count != Dimension) throw new DimensionMismatchException();
            return(IndexOf(values) >= 0); 
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
        /// Gets the number of entries in the sample.
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
        /// worrying about performance. The disadvantage of this is that the returned sample cannot be altered. If you
        /// need to alter values in a column independent of the multi-variate sample, use the <see cref="Sample.Copy"/>
        /// method to obtain an independent copy of the column.</para>
        /// </remarks>
        public Sample Column (int c) {
            if ((c < 0) || (c >= Dimension)) throw new ArgumentOutOfRangeException(nameof(c));
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
        /// worrying about performance. The disadvantage of this is that the returned bivariate sample cannot be altered. If you
        /// need to alter values independent of the multi-variate sample, use the <see cref="BivariateSample.Copy"/>
        /// method to obtain an independent copy of the bivariate sample.</para>
        /// </remarks>
        public BivariateSample TwoColumns (int cx, int cy) {
            if ((cx < 0) || (cx >= Dimension)) throw new ArgumentOutOfRangeException(nameof(cx));
            if ((cy < 0) || (cy >= Dimension)) throw new ArgumentOutOfRangeException(nameof(cy));
            return (new BivariateSample(storage[cx], storage[cy], true));
        }

        /// <summary>
        /// Gets the indicated columns as a multivariate sample.
        /// </summary>
        /// <param name="columnIndexes">A list of column indexes.</param>
        /// <returns>A read-only <see cref="MultivariateSample"/> consisting of the indicated columns.</returns>
        /// <remarks>
        /// <para>Use this method to perform multivariate analyses, such as regression and principal component analysis, using
        /// only a subset of the variables in the original multivariate sample.</para>
        /// <para>Note that this is a fast operation, which does not create independent copies of the columns.</para>
        /// </remarks>
        public MultivariateSample Columns (IReadOnlyList<int> columnIndexes) {
            if (columnIndexes == null) throw new ArgumentNullException(nameof(columnIndexes));
            SampleStorage[] columns = new SampleStorage[columnIndexes.Count];
            for (int i = 0; i < columns.Length; i++) {
                int ci = columnIndexes[i];
                if ((ci < 0) || (ci >= Dimension)) throw new ArgumentOutOfRangeException(nameof(columnIndexes));
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
            return (Columns((IReadOnlyList<int>)columnIndexes));
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
        public double RawMoment (params int[] powers) {
            return (Multivariate.RawMoment(storage, powers));

        }

        /// <summary>
        /// Computes the given sample raw moment.
        /// </summary>
        /// <param name="powers">The power to which each component should be raised.</param>
        /// <returns>The specified moment.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="powers"/> is null.</exception>
        /// <exception cref="DimensionMismatchException">The length of <paramref name="powers"/> is not
        /// equal to the <see cref="Dimension"/> of the multivariate sample.</exception>
        public double RawMoment (IReadOnlyList<int> powers) {
            return (Multivariate.RawMoment(storage, powers));
        }

        /// <summary>
        /// Computes the given sample central moment.
        /// </summary>
        /// <param name="powers">The power to which each component should be raised.</param>
        /// <returns>The specified moment.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="powers"/> is null.</exception>
        /// <exception cref="DimensionMismatchException">The length of <paramref name="powers"/> is not
        /// equal to the <see cref="Dimension"/> of the multivariate sample.</exception>
        public double CentralMoment (IReadOnlyList<int> powers) {
            if (powers == null) throw new ArgumentNullException(nameof(powers));
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
        public double CentralMoment (params int[] powers) {
            return (CentralMoment((IReadOnlyList<int>) powers));
        }

        // the internal linear regression routine, which assumes inputs are entirely valid
        
        private MultiLinearRegressionResult LinearRegression_Internal (int outputIndex) {

            IReadOnlyList<double> yColumn = storage[outputIndex];
            IReadOnlyList<double>[] xColumns = new IReadOnlyList<double>[Dimension];
            string[] xNames = new string[Dimension];
            for (int c = 0; c < Dimension; c++) {
                if (c == outputIndex) {
                    xColumns[c] = null;
                    xNames[c] = "Intercept";
                } else {
                    SampleStorage xColumn = storage[c];
                    xColumns[c] = xColumn;
                    if (String.IsNullOrWhiteSpace(xColumn.Name)) {
                        xNames[c] = c.ToString();
                    } else {
                        xNames[c] = xColumn.Name;
                    }
                }
            }

            return (new MultiLinearRegressionResult(yColumn, xColumns, xNames));
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
        /// Thus if a linear regression analysis is done on a 4-dimensional multivariate sample
        /// to predict variable number 2, the coefficients of variables 0, 1, and 3 will be
        /// parameters 0, 1, and 3, of the returned fit result, and the intercept will be
        /// parameter 2.</para>
        /// <para>If you want to include fewer input variables in your regression, use the
        /// <see cref="Columns(IReadOnlyList{Int32})"/> method to create a multivariate sample that includes only
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
        /// variables independently -- for example we are running an aquarium -- we would certainly
        /// want to know the specific effect of one variable -- that our fishes would actually prefer us to turn
        /// down the temperature while maintaining a high plankton level -- rather than the observed
        /// effect as a variable changes along with all the others that tend to change with it.
        /// This does not mean that the correlation analysis is wrong -- higher temperatures
        /// are indeed associated with higher fish densities in our hypothetical data set. It
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
        /// correlated in the observational data that a regression analysis cannot reliably
        /// tease out the independent effect of each. In that case, a fit using only one
        /// of the variables will be as good or nearly as good as a fit using both, and
        /// the covariance between their corresponding linear regression coefficients will
        /// be large. In a situation like this, you should be wary of drawing any conclusions
        /// about their separate effects.</para>
        /// <para>It can also occur that an input variable or a set of input variables is indeed good
        /// predictor of an output variable, but via a complex and non-linear relationship that
        /// a linear regression analysis will completely miss.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="outputIndex"/> is outside the range of allowed indexes.</exception>
        /// <exception cref="InsufficientDataException">There are fewer entries than the dimension of the multivariate sample.</exception>
        /// <exception cref="DivideByZeroException">The curvature matrix is singular, indicating that the data is independent of
        /// one or more parameters, or that two or more parameters are linearly dependent.</exception>
        /// <seealso cref="LinearRegression(int)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Linear_regression"/>
        public MultiLinearRegressionResult LinearRegression (int outputIndex) {

            if ((outputIndex < 0) || (outputIndex >= Dimension)) throw new ArgumentOutOfRangeException(nameof(outputIndex));

            return (LinearRegression_Internal(outputIndex));

        }

        /// <summary>
        /// Performs a linear logistic regression analysis.
        /// </summary>
        /// <param name="outputIndex">The index of the column to predict.</param>
        /// <returns>A logistic multi-linear model fit. The kth parameter is the slope of the multi-linear model with respect to
        /// the kth column, except for k equal to the <paramref name="outputIndex"/>, for which it is the intercept.</returns>
        /// <remarks>Logistic linear regression is suited to situations where multiple input variables, either continuous or binary indicators, are used to predict
        /// the value of a binary output variable. Like a linear regression, a logistic linear regression tries to find a model that predicts the output variable using
        /// a linear combination of input variables. Unlike a simple linear regression, the model does not assume that this linear
        /// function predicts the output directly; instead it assumes that this function value is then fed into a logit link function, which
        /// maps the real numbers into the interval (0, 1), and interprets the value of this link function as the probability of obtaining success value
        /// for the output variable.</remarks>
        /// <exception cref="InvalidOperationException">The column to be predicted contains values other than 0 and 1.</exception>
        /// <exception cref="InsufficientDataException">There are not more rows in the sample than columns.</exception>
        /// <exception cref="DivideByZeroException">The curvature matrix is singular, indicating that the data is independent of
        /// one or more parameters, or that two or more parameters are linearly dependent.</exception>
        public MultiLinearLogisticRegressionResult LogisticLinearRegression (int outputIndex) {

            if ((outputIndex < 0) || (outputIndex >= this.Dimension)) throw new ArgumentOutOfRangeException(nameof(outputIndex));
            if (this.Count <= this.Dimension) throw new InsufficientDataException();

            List<bool> yColumn = new List<bool>(this.Count);
            SampleStorage yStorage = storage[outputIndex];
            foreach (double yValue in yStorage) {
                if (yValue == 0.0) {
                    yColumn.Add(false);
                } else if (yValue == 1.0) {
                    yColumn.Add(true);
                } else {
                    throw new InvalidOperationException();
                }
            }

            IReadOnlyList<double>[] xColumns = new IReadOnlyList<double>[Dimension];
            string[] xNames = new string[Dimension];
            for (int c = 0; c < Dimension; c++) {
                if (c == outputIndex) {
                    xColumns[c] = null;
                    xNames[c] = "Intercept";
                } else {
                    SampleStorage xColumn = storage[c];
                    xColumns[c] = xColumn;
                    if (String.IsNullOrWhiteSpace(xColumn.Name)) {
                        xNames[c] = c.ToString();
                    } else {
                        xNames[c] = xColumn.Name;
                    }
                }
            }

            return (new MultiLinearLogisticRegressionResult(yColumn, xColumns, xNames));
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
            return(Multivariate.PrincipalComponentAnalysis(storage));
        }

        /// <summary>
        /// Compute k-means clusters of the data.
        /// </summary>
        /// <param name="m">The number of clusters to compute.</param>
        /// <returns>A description of the identified clusters.</returns>
        public MeansClusteringResult MeansClustering (int m) {
            return (Multivariate.MeansClustering(storage, m));
        }
        
    }

}
