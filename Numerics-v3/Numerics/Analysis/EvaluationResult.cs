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

        internal EvaluationResult (int count, EvaluationSettings settings) {
            Debug.Assert(count > 0);
            Debug.Assert(settings != null);
            this.count = count;
            this.settings = settings;
        }

        private readonly int count;

        private readonly EvaluationSettings settings;

        /// <summary>
        /// Gets the number of function evaluations performed.
        /// </summary>
        public int EvaluationCount {
            get {
                return (count);
            }
        }

        /// <summary>
        /// Gets the settings used for the analysis.
        /// </summary>
        public EvaluationSettings Settings {
            get {
                return (settings);
            }
        }

    }
}
