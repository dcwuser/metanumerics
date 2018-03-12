using System;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents the result of a linear logistic regression fit.
    /// </summary>
    public class LinearLogisticRegressionResult : BaseFitResult {

        internal LinearLogisticRegressionResult (ParameterCollection parameters, double logLikelihood) : base(parameters, logLikelihood) {

        }

        /// <summary>
        /// Gets the intercept of the regression model.
        /// </summary>
        public virtual UncertainValue Intercept {
            get {
                return (this.Parameters[0].Estimate);
            }
        }

        /// <summary>
        /// Gets the slope of the regression model.
        /// </summary>
        public virtual UncertainValue Slope {
            get {
                return (this.Parameters[1].Estimate);
            }
        }

        /// <summary>
        /// Computes the outcome probability for the given value.
        /// </summary>
        /// <param name="x">The value of the independent variable.</param>
        /// <returns>The probability of a <see langword="true"/> outcome, as predicted by the model.</returns>
        public virtual double Predict (double x) {
            double t = Intercept.Value + Slope.Value * x;
            return (1.0 / (1.0 + Math.Exp(-t)));
        }

    }

}
