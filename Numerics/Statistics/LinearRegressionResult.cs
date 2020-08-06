using System;
using System.Collections.Generic;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a linear regression fit.
    /// </summary>
    public sealed class LinearRegressionResult : GeneralLinearRegressionResult {

        internal LinearRegressionResult (IReadOnlyList<double> x, IReadOnlyList<double> y) :
            base () {

            double yMean, xxSum, xySum, yySum;
            Bivariate.ComputeBivariateMomentsUpToTwo(x, y, out n, out xMean, out yMean, out xxSum, out yySum, out xySum);

            b = xySum / xxSum;
            a = yMean - b * xMean;

            residuals = new List<double>(n);
            SSR = 0.0;
            SSF = 0.0;
            for(int i = 0; i < n; i++) {
                double yi = y[i];
                double ypi = a + b * x[i];
                double zi = yi - ypi;
                residuals.Add(zi);
                SSR += zi * zi;
                SSF += MoreMath.Sqr(ypi - yMean);
            }
            SST = yySum;

            xVariance = xxSum / n;
            sigmaSquared = SSR / (n - 2);
            cbb = sigmaSquared / xVariance / n;
            cab = -xMean * cbb;
            caa = (xVariance + xMean * xMean) * cbb;

            rTest = new Lazy<TestResult>(() => {
                double r = xySum / Math.Sqrt(xxSum * yySum);
                TestResult rTest = new TestResult("r", r, new PearsonRDistribution(n), TestType.TwoTailed);
                return (rTest);
            });

        }

        private readonly int n;

        private readonly double a, b;

        private readonly double cbb, cab, caa;

        private readonly double SST, SSF, SSR;

        private readonly double xMean, xVariance, sigmaSquared;

        private readonly List<double> residuals;

        private readonly Lazy<TestResult> rTest;

        /// <summary>
        /// Gets an estimate, with uncertainty, of the intercept of the line.
        /// </summary>
        public override UncertainValue Intercept {
            get {
                return new UncertainValue(a, Math.Sqrt(caa));
            }
        }

        /// <summary>
        /// Gets an estimate, with uncertainty, of the slope of the line.
        /// </summary>
        public UncertainValue Slope {
            get {
                return new UncertainValue(b, Math.Sqrt(cbb));
            }
        }

        /// <summary>
        /// Predicts the Y value at a new X value.
        /// </summary>
        /// <param name="x">The new X value.</param>
        /// <returns>The predicted value of Y and its associated uncertainty.</returns>
        public UncertainValue Predict (double x) {
            double y = a + b * x;
            double yVariance = sigmaSquared * (1.0 + (MoreMath.Sqr(x - xMean) / xVariance + 1.0) / n);
            return new UncertainValue(y, Math.Sqrt(yVariance));
        }

        /// <inheritdoc />
        public override IReadOnlyList<double> Residuals {
            get {
                return residuals;
            }
        }

        /// <inheritdoc />
        public override double SumOfSquaredResiduals {
            get {
                return SSR;
            }
        }

        /// <summary>
        /// Gets the Pearson R test of linear correlation.
        /// </summary>
        public TestResult R {
            get {
                return rTest.Value;
            }
        }

        internal override ParameterCollection CreateParameters () {
            ParameterCollection parameters = new ParameterCollection(
                nameof(Intercept), a, caa, nameof(Slope), b, cbb, cab
            );
            return parameters;
        }

        internal override OneWayAnovaResult CreateAnova () {
            AnovaRow fit = new AnovaRow(SSF, 1);
            AnovaRow residual = new AnovaRow(SSR, n - 2);
            AnovaRow total = new AnovaRow(SST, n - 1);
            OneWayAnovaResult anova = new OneWayAnovaResult(fit, residual, total);
            return anova;
        }

    }

}
