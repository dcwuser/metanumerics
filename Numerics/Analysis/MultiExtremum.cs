using System;
using System.Diagnostics;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents a minimum or maximum of a multidimensional function.
    /// </summary>
    public sealed class MultiExtremum : EvaluationResult {

        internal MultiExtremum (int count, MultiExtremumSettings settings, double[] point, double value, double precision, double[][] hessian) : base(count) {
            Debug.Assert(settings != null);
            Debug.Assert(point != null);
            this.point = point;
            this.value = value;
            this.precision = precision;
            this.hessian = hessian;
            this.settings = settings;
        }

        private readonly double[] point;
        private readonly double value;
        private readonly double precision;
        private readonly double[][] hessian;
        private readonly MultiExtremumSettings settings;

        /// <summary>
        /// Gets the dimension of the space on which the function is defined.
        /// </summary>
        public int Dimension {
            get {
                return (point.Length);
            }
        }

        /// <summary>
        /// Gets the location of the optimum.
        /// </summary>
        /// <value>A read-only vector of the coordinates at which the function reaches its optimum.</value>
        public ColumnVector Location {
            get {
                return (new ColumnVector(point, 0, 1, point.Length, true));
            }
        }

        /// <summary>
        /// Gets the value of the function at the optimum.
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
        /// <para>This should be understood as an estimate of the accuracy of <see cref="Value"/>. It should not
        /// be interpreted as a two-sided error bar, because the true minimum (maximum) value, if different,
        /// can only be smaller (larger).</para>
        /// </remarks>
        public double Precision {
            get {
                return (precision);
            }
        }

        /// <summary>
        /// Gets the Hessian matrix at the optimum.
        /// </summary>
        /// <value>A read-only matrix of approximate second partial derivatives at the optimum,
        /// or null if the algorithm does not provide one.</value>
        /// <remarks>
        /// <para>These values, if present, should be considered very approximate. If you require a more precise
        /// estimate of second derivatives, use numerical differentiation.</para>
        /// </remarks>
        public SymmetricMatrix HessianMatrix {
            get {
                return (new SymmetricMatrix(hessian, point.Length, true));
            }
        }

        /// <summary>
        /// Gets the settings used during optimization.
        /// </summary>
        public MultiExtremumSettings Settings {
            get {
                return (settings);
            }
        }

    }

}
