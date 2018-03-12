using System;
using System.Collections.Generic;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents the result of a linear logistic regression fit.
    /// </summary>
    public class MultiLinearLogisticRegressionResult : BaseFitResult {

        internal MultiLinearLogisticRegressionResult (ParameterCollection parameters, double logLikelihood) : base(parameters, logLikelihood) {

        }

        /// <summary>
        /// Gets the intercept of the regression model.
        /// </summary>
        public virtual UncertainValue Intercept {
            get {
                return (this.Parameters[this.Parameters.Count - 1].Estimate);
            }
        }

        /// <summary>
        /// Computes the predicted outcome for the given value.
        /// </summary>
        /// <param name="x">The value of the independent variable.</param>
        /// <returns>The probability of a <see langword="true"/> outcome, as predicted by the model.</returns>
        public virtual double Predict (IReadOnlyList<double> x) {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (x.Count + 1 != Parameters.Count) throw new DimensionMismatchException();
            double t = Parameters[x.Count].Estimate.Value;
            for (int i = 0; i < x.Count; i++) {
                t += Parameters[i].Estimate.Value;
            }
            return (1.0 / (1.0 + Math.Exp(-t)));
        }

    }

}
