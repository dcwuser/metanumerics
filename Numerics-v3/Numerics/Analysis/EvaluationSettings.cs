using System;

using Meta.Numerics;

namespace Meta.Numerics.Analysis {

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

        internal EvaluationSettings (object isDefault) {
            evaluationBudget = -1;
            relativePrecision = -1.0;
            absolutePrecision = -1.0;
        }
         
        private int evaluationBudget;

        private double relativePrecision;

        private double absolutePrecision;


        /// <summary>
        /// Gets or sets the total number of evaluations allowed.
        /// </summary>
        /// <value>The total number of evaluations allowed, which must be non-negative.</value>
        public int EvaluationBudget {
            get {
                return (evaluationBudget);
            }
            set {
                if (value < 0) throw new ArgumentOutOfRangeException(nameof(value));
                evaluationBudget = value;
            }
        }

        /// <summary>
        /// Gets or sets targeted relative precision.
        /// </summary>
        /// <value>The relative precision to which the result should be evaluated, which must be between 0 and 1.</value>
        public double RelativePrecision {
            get {
                return (relativePrecision);
            }
            set {
                if ((value < 0.0) || (value >= 1.0)) throw new ArgumentOutOfRangeException(nameof(value));
                relativePrecision = value;
            }
        }

        /// <summary>
        /// Gets or sets the targeted absolute precision.
        /// </summary>
        /// <value>The absolute precision to which the result should be evaluated, which must be non-negative.</value>
        public double AbsolutePrecision {
            get {
                return (absolutePrecision);
            }
            set {
                if (value < 0.0) throw new ArgumentOutOfRangeException(nameof(value));
                absolutePrecision = value;
            }
        }

        internal double ComputePrecision (double value) {
            return (absolutePrecision + Math.Abs(value) * relativePrecision);
        }

    }


}
