using System;
using System.Diagnostics;

using Meta.Numerics.Matrices;


namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents the result of the integration of a set of couple ordinary differential equations.
    /// </summary>
    public sealed class MultiOdeResult : EvaluationResult {

        internal MultiOdeResult (double x, double[] y, double[] yPrime, int count, MultiOdeSettings settings) : base(count) {
            Debug.Assert(settings != null);
            this.x = x;
            this.y = y;
            this.yPrime = yPrime;
            this.settings = settings;
        }

        private readonly double x;
        private readonly double[] y;
        private readonly double[] yPrime;
        private readonly MultiOdeSettings settings;

        /// <summary>
        /// The value of the independent variable.
        /// </summary>
        /// <value>The value of the ordinate.</value>
        public double X {
            get {
                return (x);
            }
        }

        /// <summary>
        /// The value of the dependent variable.
        /// </summary>
        /// <value>The value of the function.</value>
        public ColumnVector Y {
            get {
                return (new ColumnVector(y, 0, 1, y.Length, true));
            }
        }

        /// <summary>
        /// The value of the derivative.
        /// </summary>
        public ColumnVector YPrime {
            get {
                return (new ColumnVector(yPrime, 0, 1, yPrime.Length, true));
            }
        }

        /// <summary>
        /// Gets the settings used during integration.
        /// </summary>
        public MultiOdeSettings Settings {
            get {
                return (settings);
            }
        }
    }

}
