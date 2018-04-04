using System;
using System.Diagnostics;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents the result of the integration of an ordinary differential equation.
    /// </summary>
    public class OdeResult : EvaluationResult {

        internal OdeResult (double x, double y, double yPrime, int count, OdeSettings settings) : base(count) {
            Debug.Assert(settings != null);
            this.x = x;
            this.y = y;
            this.yPrime = yPrime;
            this.settings = settings;
        }

        private readonly double x, y, yPrime;
        private readonly OdeSettings settings;

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
        public double Y {
            get {
                return (y);
            }
        }

        /// <summary>
        /// The value of the derivative.
        /// </summary>
        public double YPrime {
            get {
                return (yPrime);
            }
        }

        /// <summary>
        /// Gets the settings used during integration.
        /// </summary>
        public OdeSettings Settings {
            get {
                return (settings);
            }
        }
    }

}
