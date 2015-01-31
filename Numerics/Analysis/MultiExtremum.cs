using System;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents a minimum or maximum of a multidimensional function.
    /// </summary>
    public class MultiExtremum : EvaluationResult {

        internal MultiExtremum (int count, EvaluationSettings settings, double[] point, double value, double precision, double[][] hessian) : base(count, settings) {
            Debug.Assert(point != null);
            this.point = point;
            this.value = value;
            this.precision = precision;
            this.hessian = hessian;
        }

        private readonly double[] point;

        private readonly double value;

        private readonly double precision;

        private double[][] hessian;

        /// <summary>
        /// Gets the dimension of the space on which the function is defined.
        /// </summary>
        public int Dimension {
            get {
                return (point.Length);
            }
        }

        /// <summary>
        /// Gets the location of the extremum.
        /// </summary>
        public ColumnVector Location {
            get {
                return (new ColumnVector(point, 0, 1, point.Length, true));
            }
        }

        /// <summary>
        /// Gets the value of the function at the extremum.
        /// </summary>
        public double Value {
            get {
                return (value);
            }
        }

        /// <summary>
        /// Gets the estimated precision of the function value.
        /// </summary>
        /// <remarks>
        /// <para> </para>
        /// </remarks>
        public double Precision {
            get {
                return (precision);
            }
        }

        /// <summary>
        /// Gets the Hessian matrix at the extremum.
        /// </summary>
        /// <remarks>
        /// <para>The Hessian matrix is the matrix of second partial derivatives. It is symmetric because the order of differentiation does not affect the result.</para>
        /// </remarks>
        public SymmetricMatrix HessianMatrix {
            get {
                return (new SymmetricMatrix(hessian, point.Length, true));
            }
        }

    }

}
