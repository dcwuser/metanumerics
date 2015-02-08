using System;

using Meta.Numerics;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents the result of a numerical integration.
    /// </summary>
    public sealed class IntegrationResult : EvaluationResult {

        internal IntegrationResult (UncertainValue estimate, int evaluationCount, EvaluationSettings settings) : base(evaluationCount, settings) {
            this.estimate = estimate;
        }

        private readonly UncertainValue estimate;

        /// <summary>
        /// Gets the estimated value of the integral and its associated error bar.
        /// </summary>
        /// <remarks>
        /// <para>Note that the associated error estimate represents an expected deviation, not
        /// a definitive bound on the deviation.</para>
        /// </remarks>
        public UncertainValue Estimate {
            get {
                return (estimate);
            }
        }

        /// <summary>
        /// Gets the estimated value of the integral.
        /// </summary>
        public double Value {
            get {
                return (estimate.Value);
            }
        }

        /// <summary>
        /// Gets the estimated precision of the evaluated value.
        /// </summary>
        public double Precision {
            get {
                return (estimate.Uncertainty);
            }
        }

    }

}
