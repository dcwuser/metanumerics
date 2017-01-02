using System;
using System.Collections.Generic;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents the result of the integration of an ordinary differential equation.
    /// </summary>
    public class OdeResult : EvaluationResult {

        internal OdeResult (double x, double y, double yPrime, int count, EvaluationSettings settings) : base(count, settings) {
            this.x = x;
            this.y = y;
            this.yPrime = yPrime;
        }

        private double x, y, yPrime;

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
        /// <value>The value of the abcissa.</value>
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

    }

}
