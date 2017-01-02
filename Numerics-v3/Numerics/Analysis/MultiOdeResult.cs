using System;
using System.Collections.Generic;

using Meta.Numerics.Matrices;


namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents the result of the integration of a set of couple ordinary differential equations.
    /// </summary>
    public class MultiOdeResult : EvaluationResult {

        internal MultiOdeResult (double x, double[] y, double[] yPrime, int count, EvaluationSettings settings) : base(count, settings) {
            this.x = x;
            this.y = y;
            this.yPrime = yPrime;
        }

        private double x;

        private double[] y;

        private double[] yPrime;

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
    }

}
