using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Contains methods for analyzing multivariate samples.
    /// </summary>
    public static class Multivariate {


        public static MultiLinearRegressionResult LinearRegression (IReadOnlyList<IReadOnlyList<double>> x, IReadOnlyList<double> y) {

            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));

            int dimension = x.Count;
            int count = y.Count;
            foreach (IReadOnlyList<double> xColumn in x) {
                if (xColumn == null) throw new ArgumentNullException(nameof(x));
                if (xColumn.Count != count) throw new DimensionMismatchException();
            }

            if (count <= (dimension + 1)) throw new InsufficientDataException();

            // Compute the design matrix X.
            int n = count;
            int m = dimension + 1;
            RectangularMatrix X = new RectangularMatrix(count, dimension + 1);
            ColumnVector v = new ColumnVector(n);
            for (int i = 0; i < count; i++) {
                for (int j = 0; j < dimension; j++) {
                    X[i, j] = x[j][i];
                }
                X[i, dimension] = 1.0;
                v[i] = y[i];
            }

            // Use X = QR to solve X b = y and compute C.
            ColumnVector b;
            SymmetricMatrix C;
            QRDecomposition.SolveLinearSystem(X, v, out b, out C);

            // For ANOVA, we will need mean and variance of y
            int yn;
            double ym, SST;
            Univariate.ComputeMomentsUpToSecond(y, out yn, out ym, out SST);

            // Compute residuals
            double SSR = 0.0;
            double SSF = 0.0;
            ColumnVector yHat = X * b;
            List<double> residuals = new List<double>(count);
            for (int i = 0; i < count; i++) {
                double z = y[i] - yHat[i];
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
        public static MultiLinearLogisticRegressionResult LinearLogisticRegression (IReadOnlyList<IReadOnlyList<double>> x, IReadOnlyList<bool> y) {

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
