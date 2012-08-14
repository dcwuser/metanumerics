using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents the result of a fit procedure.
    /// </summary>
    /// <remarks>
    /// <para>All fit methods in the Meta.Numerics library return their results as an instance of the FitResult class.
    /// This includes methods that fit a sample to a distribution
    /// (e.g. <see cref="Meta.Numerics.Statistics.Distributions.NormalDistribution.FitToSample"/>),
    /// regression methods for bivariate and multivariate data (e.g. <see cref="BivariateSample.LinearLogisticRegression"/>),
    /// and least-squares fits of data with error bars to a model function (e.g. <see cref="UncertainMeasurementSample{T}.FitToFunction"/>).
    /// </para>
    /// <para>A FitResult instance contains not only the parameter values, but also covariances and a goodness-of-fit test.</para> 
    /// <para>The vector of best-fit parameter values can be obtained as an array using the <see cref="Parameters"/> method.
    /// Individual parameter values can be obtained using the <see cref="Parameter"/> method; this method gives not only a best-fit
    /// value but also an uncertainty by returning a <see cref="UncertainValue"/>.</para>
    /// <para>The matrix of covariances can be obtained using the <see cref="CovarianceMatrix"/> method. Covariances between
    /// specific pairs of parameters van be obtained using the <see cref="Covariance"/> method.</para>
    /// <para>The goodness-of-fit test result stored in the <see cref="GoodnessOfFit"/> measures the quality of the fit.
    /// For fits to distributions, it is a Kolmogorov-Smirnov test. For regressions, it is an F-test. For fits to
    /// data with error bars, it is a chi-square test.</para>
    /// <para>Fits are done using the maximum likelyhood method, with results corrected for any small-sample bias.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Maximum_likelihood"/>
    public class FitResult {

        private IList<double> parameters;
        private SymmetricMatrix covarianceMatrix;
        private TestResult test;

        /// <summary>
        /// Gets the number of fit parameters.
        /// </summary>
        public int Dimension {
            get {
                return (parameters.Count);
            }
        }

        /// <summary>
        /// Gets an array of the best fit parameter set.
        /// </summary>
        /// <returns>An independent array of the best fit parameter values.</returns>
        public double[] Parameters () {
            double[] result = new double[Dimension];
            parameters.CopyTo(result, 0);
            return (result);
        }

        /// <summary>
        /// Get an estimate of a fit parameter.
        /// </summary>
        /// <param name="n">The (zero-based) parameter number.</param>
        /// <returns>An estimate of the parameter.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is not within [0,<see cref="Dimension"/>-1].</exception>
        public UncertainValue Parameter (int n) {
            if ((n < 0) || (n >= Dimension)) throw new ArgumentOutOfRangeException("n");
            return (new UncertainValue(parameters[n], Math.Sqrt(covarianceMatrix[n, n])));
        }

        /// <summary>
        /// Gets the covariance of two fit parameters.
        /// </summary>
        /// <param name="n">The (zero-based) number of the fist parameter.</param>
        /// <param name="m">The (zero-based) number of the second parameter.</param>
        /// <returns>The covariance of the two fit parameters.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> or <paramref name="m"/> is not within [0,<see cref="Dimension"/>-1].</exception>
        public double Covariance (int n, int m) {
            if ((n < 0) || (n >= Dimension)) throw new ArgumentOutOfRangeException("n");
            if ((m < 0) || (m >= Dimension)) throw new ArgumentOutOfRangeException("m");
            return (covarianceMatrix[n, m]);
        }

        /// <summary>
        /// Gets the coefficient of correlation between two fit parameters.
        /// </summary>
        /// <param name="n">The (zero-based) number of the first parameter.</param>
        /// <param name="m">The (zero-based) number of the second parameter.</param>
        /// <returns>The correlation coefficient between the two parameters.</returns>
        /// <remarks>
        /// <para>The correlation coefficient between two parameters is a re-scaling of their covariance to a number
        /// between -1 and 1 that indicates the strength of their correlation.</para>
        /// <para>The correlation coefficient is also called the Pearson R coefficient.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> or <paramref name="m"/> is not within [0,<see cref="Dimension"/>-1].</exception>
        public double CorrelationCoefficient (int n, int m) {
            if ((n < 0) || (n >= Dimension)) throw new ArgumentOutOfRangeException("n");
            if ((m < 0) || (m >= Dimension)) throw new ArgumentOutOfRangeException("m");
            return (covarianceMatrix[n, m] / Math.Sqrt(covarianceMatrix[n,n] * covarianceMatrix[m,m]));
        }

        /// <summary>
        /// Gets the covariance matrix containing the variances and covariances for all fit parameters.
        /// </summary>
        /// <returns>The covariance matrix.</returns>
        public SymmetricMatrix CovarianceMatrix () {
            return (covarianceMatrix.Copy());
        }

        /// <summary>
        /// Gets a test of the quality of the fit.
        /// </summary>
        /// <remarks>
        /// <para>Which test is used to evaluate goodness-of-fit depends on the type of fit that was performed.</para>
        /// <para>If no goodness-of-fit test was performed as a part of the fit, this property will be null.</para>
        /// </remarks>
        public TestResult GoodnessOfFit {
            get {
                return (test);
            }
        }

        // one-parameter constructor
        internal FitResult (double p1, double dp1, TestResult test) {
            this.parameters = new ColumnVector(1);
            this.parameters[0] = p1;
            this.covarianceMatrix = new SymmetricMatrix(1);
            this.covarianceMatrix[0, 0] = dp1 * dp1;
            this.test = test;
        }

        // two-parameter constructor
        internal FitResult (double p1, double dp1, double p2, double dp2, double cov, TestResult test) {
            this.parameters = new ColumnVector(2);
            this.parameters[0] = p1;
            this.parameters[1] = p2;
            this.covarianceMatrix = new SymmetricMatrix(2);
            this.covarianceMatrix[0, 0] = dp1 * dp1;
            this.covarianceMatrix[1, 1] = dp2 * dp2;
            this.covarianceMatrix[0, 1] = cov;
            this.test = test;
        }

        // n-parameter constructor
        internal FitResult (IList<double> parameters, SymmetricMatrix covariance, TestResult test) {
            Debug.Assert(parameters != null);
            Debug.Assert(covariance != null);
            Debug.Assert(parameters.Count == covariance.Dimension);
            this.parameters = parameters;
            this.covarianceMatrix = covariance;
            this.test = test;
        }

    }
}

