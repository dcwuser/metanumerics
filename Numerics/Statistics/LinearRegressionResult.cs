using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Describes the result of a linear regression.
    /// </summary>
    public sealed class LinearRegressionResult : GeneralLinearRegressionResult {

        internal LinearRegressionResult (
            ParameterCollection parameters,
            TestResult rTest,
            OneWayAnovaResult anova,
            Sample residuals,
            Func<double, UncertainValue> predict
        ) : base(parameters, anova, residuals) {
            this.rTest = rTest;
            this.predict = predict;
        }

        private readonly Func<double, UncertainValue> predict;

        private readonly TestResult rTest;

        /// <summary>
        /// Gets the best fit value of the intercept and its associated uncertainty.
        /// </summary>
        public UncertainValue Intercept {
            get {
                return (this.Parameters[0].Estimate);
            }
        }

        /// <summary>
        /// Gets the best-fit value of the slope and its associated uncertainty.
        /// </summary>
        public UncertainValue Slope {
            get {
                return (this.Parameters[1].Estimate);
            }
        }

        /// <summary>
        /// Predicts the Y value at a new X value.
        /// </summary>
        /// <param name="x">The new X value.</param>
        /// <returns>The predicted value of Y and its associated uncertainty.</returns>
        public UncertainValue Predict (double x) {
            return (predict(x));
        }

        /// <summary>
        /// Gets the Pearson R test of linear correlation.
        /// </summary>
        public TestResult R {
            get {
                return (rTest);
            }
        }

    }


}
