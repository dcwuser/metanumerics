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
    public static class Multivariate {


        public static MultiLinearRegressionResult MultiLinearRegression(this IReadOnlyList<double> yColumn, IReadOnlyDictionary<string, IReadOnlyList<double>> xColumnDictionary) {
            if (yColumn == null) throw new ArgumentNullException(nameof(yColumn));
            if (xColumnDictionary == null) throw new ArgumentNullException(nameof(xColumnDictionary));

            List<string> xNames = new List<string>();
            List<IReadOnlyList<double>> xColumns = new List<IReadOnlyList<double>>();
            foreach(KeyValuePair<string, IReadOnlyList<double>> xEntry in xColumnDictionary) {
                IReadOnlyList<double> xColumn = xEntry.Value;
                if (xColumn == null) throw new ArgumentNullException("xColumn");
                if (xColumn.Count != yColumn.Count) throw new DimensionMismatchException();
                xNames.Add(xEntry.Key);
                xColumns.Add(xEntry.Value);
            }
            xNames.Add("Intercept");
            xColumns.Add(null);
            return (MultiLinearRegression(yColumn, xColumns, xNames));
        }

        internal static MultiLinearRegressionResult MultiLinearRegression(IReadOnlyList<double> yColumn, IReadOnlyList<IReadOnlyList<double>> xColumns, IReadOnlyList<string> xNames) {
            return (new MultiLinearRegressionResult(yColumn, xColumns, xNames));
            /*
            Debug.Assert(yColumn != null);
            Debug.Assert(xColumns != null);
            Debug.Assert(xColumns.Count > 0);

            int n = yColumn.Count;
            int m = xColumns.Count;
            if (n <= m) throw new InsufficientDataException();

            // Compute the design matrix X.
            bool foundInterceptColumn = false;
            RectangularMatrix X = new RectangularMatrix(n, m);
            for (int c = 0; c < m; c++) {
                IReadOnlyList<double> xColumn = xColumns[c];
                if (xColumn == null) {
                    Debug.Assert(xNames[c] == "Intercept");
                    Debug.Assert(!foundInterceptColumn);
                    for (int r = 0; r < n; r++) {
                        X[r, c] = 1.0;
                    }
                    foundInterceptColumn = true;
                } else {
                    if (xColumn.Count != n) throw new DimensionMismatchException();
                    for (int r = 0; r < n; r++) {
                        X[r, c] = xColumn[r];
                    }
                }
            }
            Debug.Assert(foundInterceptColumn);
            ColumnVector v = new ColumnVector(yColumn);

            // Use X = QR to solve X b = y and compute C.
            ColumnVector b;
            SymmetricMatrix C;
            QRDecomposition.SolveLinearSystem(X, v, out b, out C);

            // For ANOVA, we will need mean and variance of y
            int yn;
            double ym, SST;
            Univariate.ComputeMomentsUpToSecond(yColumn, out yn, out ym, out SST);

            // Compute residuals
            double SSR = 0.0;
            double SSF = 0.0;
            ColumnVector yHat = X * b;
            List<double> residuals = new List<double>(n);
            for (int i = 0; i < n; i++) {
                double z = yColumn[i] - yHat[i];
                residuals.Add(z);
                SSR += z * z;
                SSF += MoreMath.Sqr(yHat[i] - ym);
            }
            double sigma2 = SSR / (n - m);

            // Scale up C by \sigma^2
            // (It sure would be great to be able to overload *=.)
            for (int i = 0; i < m; i++) {
                for (int j = i; j < m; j++) {
                    C[i, j] = C[i, j] * sigma2;
                }
            }

            // Use sums-of-squares to do ANOVA
            AnovaRow fit = new AnovaRow(SSF, m - 1);
            AnovaRow residual = new AnovaRow(SSR, n - m);
            AnovaRow total = new AnovaRow(SST, n - 1);
            OneWayAnovaResult anova = new OneWayAnovaResult(fit, residual, total);

            ParameterCollection parameters = new ParameterCollection(xNames, b, C);
            return (new MultiLinearRegressionResult(parameters, anova, residuals));
            */
        }

        /// <summary>
        /// Performs a multi-variate linear regression.
        /// </summary>
        /// <param name="yColumn">The dependent variable to be predicted by the regression.</param>
        /// <param name="xColumns">The independent variables that serve as inputs to the regression function.</param>
        /// <returns>The regression result.</returns>
        public static MultiLinearRegressionResult MultiLinearRegression(this IReadOnlyList<double> yColumn, params IReadOnlyList<double>[] xColumns) {
            if (xColumns == null) throw new ArgumentNullException(nameof(xColumns));
            if (yColumn == null) throw new ArgumentNullException(nameof(yColumn));
            return (MultiLinearRegression(yColumn, (IReadOnlyList<IReadOnlyList<double>>) xColumns));
        }

        /// <summary>
        /// Performs a multi-variate linear regression.
        /// </summary>
        /// <param name="xColumns">The values of the indepdent variables that serve as inputs to the regression function.</param>
        /// <param name="yColumn">The values of the dependent variable to be predicted by the regression.</param>
        /// <returns>The regression result.</returns>
        public static MultiLinearRegressionResult MultiLinearRegression (IReadOnlyList<double> yColumn, IReadOnlyList<IReadOnlyList<double>> xColumns) {
            if (yColumn == null) throw new ArgumentNullException(nameof(yColumn));
            if (xColumns == null) throw new ArgumentNullException(nameof(xColumns));

            List<string> xNames = new List<string>();
            List<IReadOnlyList<double>> xColumnsCopy = new List<IReadOnlyList<double>>();
            foreach (IReadOnlyList<double> xColumn in xColumns) {
                if (xColumn == null) throw new ArgumentNullException("xColumn");
                if (xColumn.Count != yColumn.Count) throw new DimensionMismatchException();
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

            return (MultiLinearRegression(yColumn, xColumnsCopy, xNames));
            /*
            int count = yColumn.Count;
            int dimension = xColumns.Count;
            foreach (IReadOnlyList<double> xColumn in xColumns) {
                if (xColumn == null) throw new ArgumentNullException(nameof(xColumns));
                if (xColumn.Count != count) throw new DimensionMismatchException();
            }

            if (count <= (dimension + 1)) throw new InsufficientDataException();

            // Compute the design matrix X.
            int n = count;
            int m = dimension + 1;
            RectangularMatrix X = new RectangularMatrix(n, m);
            ColumnVector v = new ColumnVector(n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < dimension; j++) {
                    X[i, j] = xColumns[j][i];
                }
                X[i, dimension] = 1.0;
                v[i] = yColumn[i];
            }

            // Use X = QR to solve X b = y and compute C.
            ColumnVector b;
            SymmetricMatrix C;
            QRDecomposition.SolveLinearSystem(X, v, out b, out C);

            // For ANOVA, we will need mean and variance of y
            int yn;
            double ym, SST;
            Univariate.ComputeMomentsUpToSecond(yColumn, out yn, out ym, out SST);

            // Compute residuals
            double SSR = 0.0;
            double SSF = 0.0;
            ColumnVector yHat = X * b;
            List<double> residuals = new List<double>(count);
            for (int i = 0; i < count; i++) {
                double z = yColumn[i] - yHat[i];
                residuals.Add(z);
                SSR += z * z;
                SSF += MoreMath.Sqr(yHat[i] - ym);
            }
            double sigma2 = SSR / (count - m);


            // Scale up C by \sigma^2
            // (It sure would be great to be able to overload *=.)
            for (int i = 0; i < m; i++) {
                for (int j = i; j < m; j++) {
                    C[i, j] = C[i, j] * sigma2;
                }
            }

            // Use sums-of-squares to do ANOVA
            AnovaRow fit = new AnovaRow(SSF, m - 1);
            AnovaRow residual = new AnovaRow(SSR, n - m);
            AnovaRow total = new AnovaRow(SST, n - 1);
            OneWayAnovaResult anova = new OneWayAnovaResult(fit, residual, total);

            string[] names = new string[m];
            for (int j = 0; j < dimension; j++) {
                names[j] = $"[{j}]";
            }
            names[dimension] = "Intercept";
            ParameterCollection parameters = new ParameterCollection(names, b, C);

            return (new MultiLinearRegressionResult(parameters, anova, residuals));
            */
        }

        /// <summary>
        /// Performs a linear logistic regression analysis.
        /// </summary>
        /// <param name="x">A list of independent variable values.</param>
        /// <param name="y">The boolean dependent variable values.</param>
        /// <returns>A logistic multi-linear model fit. The kth parameter is the slope of the multi-linear model with respect to
        /// the kth column, except for k equal to the <paramref name="outputIndex"/>, for which it is the intercept.</returns>
        /// <remarks>Logistic linear regression is suited to situations where multiple input variables, either continuous or binary indicators, are used to predict
        /// the value of a binary output variable. Like a linear regression, a logistic linear regression tries to find a model that predicts the output variable using
        /// a linear combination of input variables. Unlike a simple linear regression, the model does not assume that this linear
        /// function predicts the output directly; instead it assumes that this function value is then fed into a logit link function, which
        /// maps the real numbers into the interval (0, 1), and interprets the value of this link function as the probability of obtaining success value
        /// for the output variable.</remarks>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> or <paramref name="y"/> is <see langword="null"/>, or one of the columns
        /// in <paramref name="x"/> is <see langword="null"/>.</exception>
        /// <exception cref="DimensionMismatchException">Not all the variable lists have the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">There are not more rows in the sample than columns.</exception>
        /// <exception cref="DivideByZeroException">The curvature matrix is singular, indicating that the data is independent of
        /// one or more of the independent variables, or that two or more variables are linearly dependent.</exception>
        public static MultiLinearLogisticRegressionResult MultiLinearLogisticRegression (IReadOnlyList<IReadOnlyList<double>> x, IReadOnlyList<bool> y) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));

            int inputDimension = x.Count;
            int count = y.Count;
            foreach (IReadOnlyList<double> xColumn in x) {
                if (xColumn == null) throw new ArgumentNullException(nameof(x));
                if (xColumn.Count != count) throw new DimensionMismatchException();
            }
            int parameterDimension = inputDimension + 1;

            if (count <= parameterDimension) throw new InsufficientDataException();

            // Define the log likelihood as a function of the parameter set
            Func<IReadOnlyList<double>, double> logLikelihood = (IReadOnlyList<double> a) => {
                double L = 0.0;
                for (int k = 0; k < count; k++) {
                    double z = a[inputDimension];
                    for (int i = 0; i < inputDimension; i++) {
                        z += a[i] * x[i][k];
                    }
                    double ez = Math.Exp(z);

                    if (y[k]) {
                        L -= MoreMath.LogOnePlus(1.0 / ez);
                    } else {
                        L -= MoreMath.LogOnePlus(ez);
                    }
                }
                return (L);
            };

            // we need  a better starting value
            double[] start = new double[parameterDimension];

            MultiExtremum maximum = MultiFunctionMath.FindLocalMaximum(logLikelihood, start);
            CholeskyDecomposition CD = maximum.HessianMatrix.CholeskyDecomposition();
            if (CD == null) throw new DivideByZeroException();

            string[] names = new string[parameterDimension];
            for (int i = 0; i < inputDimension; i++) {
                names[i] = $"[{i}]";
            }
            names[inputDimension] = "Intercept";
            ParameterCollection parameters = new ParameterCollection(names, maximum.Location, CD.Inverse());
            MultiLinearLogisticRegressionResult result = new MultiLinearLogisticRegressionResult(parameters, maximum.Value);

            return (result);

        }

        /// <summary>
        /// Performs a principal component analysis of the multivariate sample.
        /// </summary>
        /// <returns>The result of the principal component analysis.</returns>
        /// <param name="sample">The set of columns to analyze.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> or one of its members is null.</exception>
        /// <exception cref="DimensionMismatchException">The columns do not all have the same number of entries.</exception>
        /// <exception cref="InsufficientDataException">The number of entries is less than the number of columns.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Principal_component_analysis"/>
        public static PrincipalComponentAnalysis PrincipalComponentAnalysis (IReadOnlyList<IReadOnlyList<double>> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));

            int dimension = sample.Count;
            int count = -1;
            for (int c = 0; c < dimension; c++) {
                IReadOnlyList<double> column = sample[c];
                if (column == null) throw new ArgumentNullException(nameof(column));
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
                IReadOnlyList<double> column = sample[c];
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
        public static PrincipalComponentAnalysis PrincipalComponentAnalysis (params IReadOnlyList<double>[] columns) {
            return (PrincipalComponentAnalysis((IReadOnlyList<IReadOnlyList<double>>) columns));
        }

    }
}
