using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Meta.Numerics.Analysis;
using Meta.Numerics.Data;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Contains methods for analyzing multivariate samples.
    /// </summary>
    /// <remarks>
    /// <para>A multivariate sample is one in which each element of the sample consists of measurements of
    /// multiple quantities. For example, a survey in which we record the sex, income, education level,
    /// and race of each participant is a multivariate sample. Note that the types of each quantity
    /// measured need not be the same.</para>
    /// <para>One common form of statistical analysis of a multivariate sample is to try to predict
    /// one variable on the basis of some of the others. To do this you can use
    /// <see cref="MultiLinearRegression(IReadOnlyList{double}, IReadOnlyList{double}[])"/>
    /// or, if the dependent variable is Boolean,
    /// <see cref="MultiLinearLogisticRegression(IReadOnlyList{bool}, IReadOnlyList{double}[])"/>.</para>
    /// <para>Another common operation of multivariate data is to try to organize the observations
    /// into clusters, which you can do using
    /// <see cref="MeansClustering(IReadOnlyList{IReadOnlyList{double}}, int)"/>.</para>
    /// <para>Instead of trying to reduce all the observations to a few important representatives,
    /// you can try to reduce all the variables to a few important directions using
    /// <see cref="PrincipalComponentAnalysis(IReadOnlyList{IReadOnlyList{double}})"/>.</para>
    /// </remarks>
    public static class Multivariate {

        /// <summary>
        /// Performs a multivariate linear regression on the named columns.
        /// </summary>
        /// <param name="yColumn">A column of dependent variable observations to be predicted.</param>
        /// <param name="xColumnDictionary">A dictionary of columns of independent variable observations
        /// that serve as inputs to the regression function.</param>
        /// <returns>A multi-linear fit.</returns>
        /// <exception cref="NullReferenceException"><paramref name="yColumn"/> is <see langword="null"/>,
        /// or <paramref name="xColumnDictionary"/> is <see langword="null"/>,
        /// or one of the x-columns it contains is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all of the columns contain the same number of observations.</exception>
        /// <exception cref="InsufficientDataException">There are too few observations to perform the fit.</exception>
        public static MultiLinearRegressionResult MultiLinearRegression(this IReadOnlyList<double> yColumn, IReadOnlyDictionary<string, IReadOnlyList<double>> xColumnDictionary) {
            if (yColumn == null) throw new ArgumentNullException(nameof(yColumn));
            if (xColumnDictionary == null) throw new ArgumentNullException(nameof(xColumnDictionary));

            List<IReadOnlyList<double>> xColumns;
            List<string> xNames;
            PrepareColumnsFromDictionary(xColumnDictionary, yColumn.Count, out xColumns, out xNames);

            return (new MultiLinearRegressionResult(yColumn, xColumns, xNames));
        }

        /// <summary>
        /// Performs a multivariate linear regression.
        /// </summary>
        /// <param name="yColumn">A column of dependent variable observations to be predicted.</param>
        /// <param name="xColumns">Columns of independent variable observations
        /// that serve as inputs to the regression function.</param>
        /// <returns>A multi-linear fit.</returns>
        /// <exception cref="NullReferenceException"><paramref name="yColumn"/> is <see langword="null"/>,
        /// or <paramref name="xColumns"/> is <see langword="null"/>,
        /// or one of the x-columns is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all of the columns contain the same number of observations.</exception>
        /// <exception cref="InsufficientDataException">There are too few observations to perform the fit.</exception>
        public static MultiLinearRegressionResult MultiLinearRegression(this IReadOnlyList<double> yColumn, params IReadOnlyList<double>[] xColumns) {
            if (xColumns == null) throw new ArgumentNullException(nameof(xColumns));
            if (yColumn == null) throw new ArgumentNullException(nameof(yColumn));
            return (MultiLinearRegression(yColumn, (IReadOnlyList<IReadOnlyList<double>>) xColumns));
        }

        /// <summary>
        /// Performs a multivariate linear regression on the listed columns.
        /// </summary>
        /// <param name="yColumn">A column of dependent variable observations to be predicted.</param>
        /// <param name="xColumns">A list of columns of independent variable observations
        /// that serve as inputs to the regression function.</param>
        /// <returns>A multi-linear fit.</returns>
        /// <exception cref="NullReferenceException"><paramref name="yColumn"/> is <see langword="null"/>,
        /// or <paramref name="xColumns"/> is <see langword="null"/>,
        /// or one of the x-columns in the list is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all of the columns contain the same number of observations.</exception>
        /// <exception cref="InsufficientDataException">There are too few observations to perform the fit.</exception>
        /// <remarks>
        /// <para>Note that the <paramref name="xColumns"/> argument is column-oriented, i.e. it is a list of columns,
        /// not row-oriented, i.e. a list of rows. Either of these would match the method signature, but only
        /// column-oriented inputs will produce correct results.</para>
        /// </remarks>
        public static MultiLinearRegressionResult MultiLinearRegression (IReadOnlyList<double> yColumn, IReadOnlyList<IReadOnlyList<double>> xColumns) {
            if (yColumn == null) throw new ArgumentNullException(nameof(yColumn));
            if (xColumns == null) throw new ArgumentNullException(nameof(xColumns));

            List<string> xNames;
            List<IReadOnlyList<double>> xColumnsCopy;
            PrepareColumns(xColumns, yColumn.Count, out xColumnsCopy, out xNames);

            return (new MultiLinearRegressionResult(yColumn, xColumnsCopy, xNames));
        }

        // This method is used to prepare sets of x-columns in the format expected by the multivariate fitting routines.

        private static void PrepareColumns (
            IReadOnlyList<IReadOnlyList<double>> xColumns, int expectedLength,
            out List<IReadOnlyList<double>> xColumnsCopy, out List<string> xNames) {

            Debug.Assert(xColumns != null);

            xNames = new List<string>();
            xColumnsCopy = new List<IReadOnlyList<double>>();
            foreach (IReadOnlyList<double> xColumn in xColumns) {
                if (xColumn == null) throw new ArgumentNullException(nameof(xColumns));
                if (xColumn.Count != expectedLength) throw new DimensionMismatchException();
                INamed named = xColumn as INamed;
                if (named == null) {
                    xNames.Add(xNames.Count.ToString());
                } else {
                    xNames.Add(named.Name);
                }
                xColumnsCopy.Add(xColumn);
            }
            xNames.Add("Intercept");
            xColumnsCopy.Add(null);
        }

        private static void PrepareColumnsFromDictionary (
            IReadOnlyDictionary<string, IReadOnlyList<double>> xColumnDictionary, int expectedLength,
            out List<IReadOnlyList<double>> xColumns, out List<string> xNames) {

            Debug.Assert(xColumnDictionary != null);

            xNames = new List<string>();
            xColumns = new List<IReadOnlyList<double>>();
            foreach (KeyValuePair<string, IReadOnlyList<double>> xEntry in xColumnDictionary) {
                IReadOnlyList<double> xColumn = xEntry.Value;
                if (xColumn == null) throw new ArgumentNullException(nameof(xColumns));
                if (xColumn.Count != expectedLength) throw new DimensionMismatchException();
                xNames.Add(xEntry.Key);
                xColumns.Add(xEntry.Value);
            }
            xNames.Add("Intercept");
            xColumns.Add(null);
        }

        /// <summary>
        /// Performs a multivariate linear logistic regression.
        /// </summary>
        /// <param name="yColumn">A column of dependent Boolean variable observations to be predicted.</param>
        /// <param name="xColumns">Columns of independent variable observations
        /// that serve as inputs to the regression function.</param>
        /// <returns>A multi-linear fit.</returns>
        /// <exception cref="NullReferenceException"><paramref name="yColumn"/> is <see langword="null"/>,
        /// or <paramref name="xColumns"/> is <see langword="null"/>,
        /// or one of the x-columns is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all of the columns contain the same number of observations.</exception>
        /// <exception cref="InsufficientDataException">There are too few observations to perform the fit.</exception>
        /// <exception cref="DivideByZeroException">The curvature matrix is singular, indicating that the data is independent of
        /// one or more of the independent variables, or that two or more variables are linearly dependent.</exception>
        public static MultiLinearLogisticRegressionResult MultiLinearLogisticRegression (this IReadOnlyList<bool> yColumn, params IReadOnlyList<double>[] xColumns) {
            return(MultiLinearLogisticRegression(yColumn, (IReadOnlyList<IReadOnlyList<double>>) xColumns));
        }

        /// <summary>
        /// Performs a multivariate linear logistic regression on the listed columns.
        /// </summary>
        /// <param name="yColumn">A column of dependent Boolean variable observations to be predicted.</param>
        /// <param name="xColumns">Columns of independent variable observations
        /// that serve as inputs to the regression function.</param>
        /// <returns>A multi-linear fit.</returns>
        /// <exception cref="NullReferenceException"><paramref name="yColumn"/> is <see langword="null"/>,
        /// or <paramref name="xColumns"/> is <see langword="null"/>,
        /// or one of the x-columns is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all of the columns contain the same number of observations.</exception>
        /// <exception cref="InsufficientDataException">There are too few observations to perform the fit.</exception>
        /// <exception cref="DivideByZeroException">The curvature matrix is singular, indicating that the data is independent of
        /// one or more of the independent variables, or that two or more variables are linearly dependent.</exception>
        /// <remarks>
        /// <para>Note that the <paramref name="xColumns"/> argument is column-oriented, i.e. it is a list of columns,
        /// not row-oriented, i.e. a list of rows. Either of these would match the method signature, but only
        /// column-oriented inputs will produce correct results.</para>
        /// </remarks>
        public static MultiLinearLogisticRegressionResult MultiLinearLogisticRegression (this IReadOnlyList<bool> yColumn, IReadOnlyList<IReadOnlyList<double>> xColumns) {
            if (yColumn == null) throw new ArgumentNullException(nameof(yColumn));
            if (xColumns == null) throw new ArgumentNullException(nameof(xColumns));

            List<string> xNames;
            List<IReadOnlyList<double>> xColumnsCopy;
            PrepareColumns(xColumns, yColumn.Count, out xColumnsCopy, out xNames);

            return (new MultiLinearLogisticRegressionResult(yColumn, xColumnsCopy, xNames));
        }

        /// <summary>
        /// Performs a multivariate linear logistic regression on the named columns.
        /// </summary>
        /// <param name="yColumn">A column of dependent Boolean variable observations to be predicted.</param>
        /// <param name="xColumnDictionary">A dictionary of columns of independent variable observations
        /// that serve as inputs to the regression function.</param>
        /// <returns>A multi-linear fit.</returns>
        /// <exception cref="NullReferenceException"><paramref name="yColumn"/> is <see langword="null"/>,
        /// or <paramref name="xColumnDictionary"/> is <see langword="null"/>,
        /// or one of the x-columns is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all of the columns contain the same number of observations.</exception>
        /// <exception cref="InsufficientDataException">There are too few observations to perform the fit.</exception>
        /// <exception cref="DivideByZeroException">The curvature matrix is singular, indicating that the data is independent of
        /// one or more of the independent variables, or that two or more variables are linearly dependent.</exception>
        public static MultiLinearLogisticRegressionResult MultiLinearLogisticRegression (this IReadOnlyList<bool> yColumn, IReadOnlyDictionary<string, IReadOnlyList<double>> xColumnDictionary) {
            if (yColumn == null) throw new ArgumentNullException(nameof(yColumn));
            if (xColumnDictionary == null) throw new ArgumentNullException(nameof(xColumnDictionary));

            List<IReadOnlyList<double>> xColumns;
            List<string> xNames;
            PrepareColumnsFromDictionary(xColumnDictionary, yColumn.Count, out xColumns, out xNames);

            return (new MultiLinearLogisticRegressionResult(yColumn, xColumns, xNames));
        }

        /// <summary>
        /// Computes the specified multivariate raw moment.
        /// </summary>
        /// <param name="columns">The multivariate sample in column-oriented form.</param>
        /// <param name="powers">The power to which each column variable should be raised.</param>
        /// <returns>The requested multivariate moment.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="columns"/> or <paramref name="powers"/> is
        /// <see langword="null"/>, or one of the columns in <paramref name="columns"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The length of <paramref name="powers"/> is not
        /// equal to the length of <paramref name="columns"/>, or lengths of all the component
        /// columns of <paramref name="columns"/> are not equal.</exception>
        public static double RawMoment (this IReadOnlyList<IReadOnlyList<double>> columns, params int[] powers) {
            return (RawMoment(columns, (IReadOnlyList<int>) powers));
        }

        /// <summary>
        /// Computes the multivariate raw moment specified by the listed powers.
        /// </summary>
        /// <param name="columns">The multivariate sample in column-oriented form.</param>
        /// <param name="powers">The power to which each column variable should be raised.</param>
        /// <returns>The requested multivariate moment.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="columns" />
        /// or <paramref name="powers"/> is <see langword="null"/>, or one of the columns
        /// in <paramref name="columns" /> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The length of <paramref name="powers"/> is not
        /// equal to the length of <paramref name="columns"/>, or lengths of all the component
        /// columns of <paramref name="columns"/> are not equal.</exception>
        public static double RawMoment (this IReadOnlyList<IReadOnlyList<double>> columns, IReadOnlyList<int> powers) {

            if (columns == null) throw new ArgumentNullException(nameof(columns));
            if (powers == null) throw new ArgumentNullException(nameof(powers));
            if (columns.Count != powers.Count) throw new DimensionMismatchException();

            int count = -1;
            foreach (IReadOnlyList<double> column in columns) {
                if (column == null) throw new ArgumentNullException(nameof(columns));
                if (count < 0) {
                    count = column.Count;
                } else {
                    if (count != column.Count) throw new DimensionMismatchException();
                }
            }

            double M = 0.0;
            for (int i = 0; i < count; i++) {
                double t = 1.0;
                for (int j = 0; j < columns.Count; j++) {
                    t *= MoreMath.Pow(columns[j][i], powers[j]);
                }
                M += t;
            }
            M = M / count;

            return (M);

        }

        /// <summary>
        /// Performs a principal component analysis of the multivariate sample.
        /// </summary>
        /// <returns>The result of the principal component analysis.</returns>
        /// <param name="columns">The set of columns to analyze.</param>
        /// <exception cref="ArgumentNullException"><paramref name="columns"/> or one of its members is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">The columns do not all have the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">The number of entries is less than the number of columns.</exception>
        /// <remarks>
        /// <para>For background on principal component analysis, see the remarks at <see cref="Meta.Numerics.Statistics.PrincipalComponentAnalysis"/>.</para>
        /// <para>Note that the <paramref name="columns"/> argument is column-oriented, i.e. it is a list of columns,
        /// not row-oriented, i.e. a list of rows. Either of these would match the method signature, but only
        /// column-oriented inputs will produce correct results.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="columns"/>, or one of the columns in it, is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all of the columns in <paramref name="columns"/> have the same length.</exception>
        /// <exception cref="InsufficientDataException">There are fewer rows than columns.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Principal_component_analysis"/>
        public static PrincipalComponentAnalysis PrincipalComponentAnalysis (this IReadOnlyList<IReadOnlyList<double>> columns) {

            if (columns == null) throw new ArgumentNullException(nameof(columns));

            int dimension = columns.Count;
            int count = -1;
            for (int c = 0; c < dimension; c++) {
                IReadOnlyList<double> column = columns[c];
                if (column == null) throw new ArgumentNullException(nameof(columns));
                if (count < 0) {
                    count = column.Count;
                } else {
                    if (column.Count != count) throw new DimensionMismatchException();
                }
            }

            if (count < dimension) throw new InsufficientDataException();

            // construct a (Count X Dimension) matrix of mean-centered data
            double[] store = MatrixAlgorithms.AllocateStorage(count, dimension);
            int i = 0;
            for (int c = 0; c < dimension; c++) {
                IReadOnlyList<double> column = columns[c];
                double mu = column.Mean();
                for (int r = 0; r < count; r++) {
                    store[i] = column[r] - mu;
                    i++;
                }
            }

            // bidiagonalize it
            double[] a, b;
            MatrixAlgorithms.Bidiagonalize(store, count, dimension, out a, out b);

            // form the U and V matrices
            double[] left = MatrixAlgorithms.AccumulateBidiagonalU(store, count, dimension);
            double[] right = MatrixAlgorithms.AccumulateBidiagonalV(store, count, dimension);

            // find the singular values of the bidiagonal matrix
            MatrixAlgorithms.ExtractSingularValues(a, b, left, right, count, dimension);

            // sort them
            MatrixAlgorithms.SortValues(a, left, right, count, dimension);

            PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis(left, a, right, count, dimension);

            return (pca);

        }

        /// <summary>
        /// Performs a principal component analysis of the columns.
        /// </summary>
        /// <param name="columns">The columns on which to perform the analysis.</param>
        /// <returns>The result of the principal component analysis.</returns>
        public static PrincipalComponentAnalysis PrincipalComponentAnalysis (params IReadOnlyList<double>[] columns) {
            return (PrincipalComponentAnalysis((IReadOnlyList<IReadOnlyList<double>>) columns));
        }

        /// <summary>
        /// Compute k-means clusters of the data.
        /// </summary>
        /// <param name="columns">The data columns on which to perform the analysis.</param>
        /// <param name="m">The number of clusters to compute.</param>
        /// <returns>A description of the identified clusters.</returns>
        /// <remarks>
        /// <para>Note that the <paramref name="columns"/> argument is column-oriented, i.e. it is a list of columns,
        /// not row-oriented, i.e. a list of rows. Either of these would match the method signature, but only
        /// column-oriented inputs will produce correct results.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="columns"/>, or one of the columns in it, is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all of the columns in <paramref name="columns"/> have the same length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">There are fewer data points than the request cluster number <paramref name="m"/>.</exception>
        public static MeansClusteringResult MeansClustering (this IReadOnlyList<IReadOnlyList<double>> columns, int m) {

            if (columns == null) throw new ArgumentNullException(nameof(columns));
 
            int d = columns.Count;
            int n = -1;
            for (int c = 0; c < d; c++) {
                IReadOnlyList<double> column = columns[c];
                if (column == null) throw new ArgumentNullException(nameof(columns));
                if (n < 0) {
                    n = column.Count;
                } else {
                    if (column.Count != n) throw new DimensionMismatchException();
                }
            }
            if ((m <= 0) || (m >= n)) throw new ArgumentOutOfRangeException(nameof(m));


            double[][] centroids = InitializeCentroids(columns, n, m);

            // i = 0 ... n ranges over data rows
            // j = 0 ... d over data columns
            // k = 0 ... m ranges over clusters

            int[] assigments = new int[n];
            int[] counts = new int[m];

            for (int count = 0; count < Global.SeriesMax; count++) {

                bool stable = true;

                // Assignment
                for (int k = 0; k < m; k++) {
                    counts[k] = 0;
                }

                for (int i = 0; i < n; i++) {
                    double minD = Double.MaxValue;
                    int minK = -1;
                    for (int k = 0; k < m; k++) {
                        double D = 0.0;
                        for (int j = 0; j < d; j++) {
                            D += MoreMath.Sqr(centroids[k][j] - columns[j][i]);
                        }
                        if (D < minD) {
                            minD = D;
                            minK = k;
                        }
                    }
                    if (assigments[i] != minK) {
                        assigments[i] = minK;
                        stable = false;
                    }
                    counts[minK]++;
                }

                if (stable) return new MeansClusteringResult(centroids);

                // Update
                for (int k = 0; k < m; k++) {
                    for (int j = 0; j < d; j++) {
                        centroids[k][j] = 0.0;
                    }
                }

                for (int i = 0; i < n; i++) {
                    int k = assigments[i];
                    for (int j = 0; j < d; j++) {
                        centroids[k][j] += columns[j][i];
                    }
                }

                for (int k = 0; k < m; k++) {
                    for (int j = 0; j < d; j++) {
                        centroids[k][j] /= counts[k];
                    }
                }

            }

            throw new NonconvergenceException();
        }

        private static double[][] InitializeCentroids (IReadOnlyList<IReadOnlyList<double>> columns, int n, int m) {

            Debug.Assert(columns != null);
            Debug.Assert((m > 0) && (m < n));

            // Use the first point for the first centroid.
            double[][] centroids = new double[m][];
            for (int k = 0; k < m; k++) centroids[k] = new double[columns.Count];
            for (int j = 0; j < columns.Count; j++) centroids[0][j] = columns[j][0];

            // Choose the remaining centroids
            for (int k = 1; k < m; k++) {

                // For the next centroid, use the point with the largest
                // minimum distance from all the previous centroids.

                int iBest = -1;
                double maxDistance = 0.0;
                for (int i = 0; i < n; i++) {

                    // Determine the smallest distance to an existing centroid
                    double minDistance = Double.PositiveInfinity;
                    for (int k1 = 0; k1 < k; k1++) {
                        double distance = 0.0;
                        for (int j = 0; j < columns.Count; j++) {
                            distance += MoreMath.Sqr(centroids[k1][j] - columns[j][i]);
                        }
                        if (distance < minDistance) minDistance = distance;
                    }

                    // If this distance is larger than any measured before, this
                    // point is the best candidate.
                    if (minDistance > maxDistance) {
                        maxDistance = minDistance;
                        iBest = i;
                    }

                }

                // Use the best point found as the next centroid.
                for (int j = 0; j < columns.Count; j++) {
                    centroids[k][j] = columns[j][iBest];
                }

            }

            // This algorithm produces good candidates, but it's relatively expensive,
            // order m^2 * n * d. This is particularly bad for large m.

            return (centroids);

        }

    }
}
