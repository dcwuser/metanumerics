using System;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Contains settings controling the evaluation of a function.
    /// </summary>
    public class EvaluationSettings {

        /// <summary>
        /// Initializes a new set of default evaulation settings.
        /// </summary>
        public EvaluationSettings () {
            EvaluationBudget = 5000;
            RelativePrecision = Global.Accuracy;
            AbsolutePrecision = Global.Accuracy;
        }

        /// <summary>
        /// Gets or sets the total number of evaluations allowed.
        /// </summary>
        public int EvaluationBudget { get; set; }

        /// <summary>
        /// Gets or sets targeted relative precision.
        /// </summary>
        public double RelativePrecision { get; set; }

        /// <summary>
        /// Gets or sets the targeted absolute precision.
        /// </summary>
        public double AbsolutePrecision { get; set; }

        /// <summary>
        /// Occurs when an updated evaluation is available.
        /// </summary>
        public event Action<object> Update;

        internal void OnUpdate (object result) {
            if (Update != null) Update(result);
        }

    }


}
