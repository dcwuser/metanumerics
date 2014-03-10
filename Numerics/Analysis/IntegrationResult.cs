using System;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Represents the result of a numerical integration.
    /// </summary>
    public sealed class IntegrationResult {

        internal IntegrationResult (UncertainValue estimate, int evaluationCount) {
            this.estimate = estimate;
            this.evaluationCount = evaluationCount;
        }

        private readonly UncertainValue estimate;

        private readonly int evaluationCount;

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
        /// Gets the number of function evaluations that were performed.
        /// </summary>
        public int EvaluationCount {
            get {
                return (evaluationCount);
            }
        }

    }

}
