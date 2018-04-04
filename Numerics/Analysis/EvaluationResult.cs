using System;
using System.Diagnostics;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents the result of a function analysis.
    /// </summary>
    /// <remarks>
    /// <para>This is the base class for all function analysis result classes.</para>
    /// </remarks>
    public class EvaluationResult {

        internal EvaluationResult (int count) {
            Debug.Assert(count > 0);
            this.count = count;
        }

        private readonly int count;

        /// <summary>
        /// Gets the number of function evaluations performed.
        /// </summary>
        public int EvaluationCount {
            get {
                return (count);
            }
        }

    }
}
