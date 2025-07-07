using System;

using Meta.Numerics;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains settings used to control the search for roots.
    /// </summary>
    internal class ZeroSettings : EvaluationSettings {

        /// <summary>
        /// A handler that will be called with updates as the root is isolated.
        /// </summary>
        public Action<ZeroResult> Listener { get; set; }

    }


    /// <summary>
    /// Represents the result of a root-finding search.
    /// </summary>
    internal sealed class ZeroResult : EvaluationResult {

        internal ZeroResult (double value, double a, double b, int count) : base(count) {
            this.value = value;
            this.bracket = Interval.FromEndpoints(a, b);
        }

        private readonly double value;

        private readonly Interval bracket;
        
        /// <summary>
        /// The argument at which the function as a root.
        /// </summary>
        public double Value {
            get {
                return value;
            }
        }

        /// <summary>
        /// A bracket containing the root.
        /// </summary>
        public Interval Bracket {
            get {
                return bracket;
            }
        }

    }

}
