using System;

using Meta.Numerics;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents the result of a numerical integration.
    /// </summary>
    /// <remarks>
    /// <para>This class is returned by various numerical integration methods, including
    /// <see cref="FunctionMath.Integrate(Func{double, double}, Interval, IntegrationSettings)"/>
    /// and <see cref="MultiFunctionMath.Integrate(Func{System.Collections.Generic.IList{double}, double}, System.Collections.Generic.IList{Interval}, IntegrationSettings)"/>.
    /// In addition to an estimate of the integral and the associated uncertainty, it gives a count of the number of function evaluations
    /// that were required and the <see cref="EvaluationSettings"/> that were used for the integration.</para>
    /// </remarks>
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
        /// Gets the estimated value.
        /// </summary>
        public static implicit operator double (IntegrationResult result) {
            return (result.Value);
        }

    }

}
