using System;

using Meta.Numerics;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Contains settings governing the evaluation of a function by a solver.
    /// </summary>
    /// <remarks>
    /// <para>Negative values of <see cref="EvaluationBudget"/>,
    /// <see cref="RelativePrecision"/>, and <see cref="AbsolutePrecision"/>
    /// indicate that the solver method should use its defaults for
    /// that property. You can override the default for a property by
    /// setting it explicitly. If you set values for some properties
    /// but not others, your setting will be applied to the property
    /// you set and the others will use defaults.</para>
    /// <para>When a solver method returns an <see cref="EvaluationResult"/>,
    /// its evaluation setting object will contain the specific
    /// settings used, so you can see which default values were applied.</para>
    /// </remarks>
    public class EvaluationSettings {

        /// <summary>
        /// Initializes a new set of default evaluation settings.
        /// </summary>
        public EvaluationSettings () {
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
