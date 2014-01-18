using System;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Represents a maximum or minimum of a function of one variable.
    /// </summary>
    public class LineExtremum {

        internal LineExtremum (double x, double f, double f2, int count) {
            this.x = x;
            this.f = f;
            this.f2 = f2;
            this.count = count;
        }

        private readonly double x;
        private readonly double f;
        private readonly double f2;
        private readonly int count;

        /// <summary>
        /// Gets the location (x-value) of the extremum.
        /// </summary>
        /// <remarks>
        /// <para>Note that numerical methods for finding typical a maximum or minimum cannot usually determine
        /// its location to full floating-point precision. Near a quadratic extremum, a change in x of ~&#x3B5; will
        /// change f(x) by ~&#x3B5;<sup>2</sup>. Thus the smallest detectable change in f(x) will
        /// typically correspond to a change in x of order of the square root of full precision. Full
        /// <see cref="System.Double"/> precision being ~16 digits, you should expect the location to
        /// be accurate only to ~8 digits.</para>
        /// </remarks>
        public double Location {
            get {
                return (x);
            }
        }

        /// <summary>
        /// Gets the function value (y-value) at the extremum.
        /// </summary>
        public double Value {
            get {
                return (f);
            }
        }

        /// <summary>
        /// Gets the curvature at the extremum.
        /// </summary>
        /// <remarks>
        /// <para>The curvature is the second derivative of the function at the extremum.</para>
        /// <para>At a typical extremum, where the function has vanishing first derivative, the second derivative will be a number
        /// whose magnitude characterizes the "steepness" with which the function increases as one moves away from the extremum.</para>
        /// <para>At an atypical extremum, for example at an interval boundary or of a non-smooth function, this
        /// value may be meaningless.</para>
        /// <para>Even in the case of a typical extremum, the value of the curvature property will typically be accurate only
        /// to a handfull of digits. If you require a highly accurate determination of the curvature,
        /// you should compute the second derivative of the minimzed function explicitly.</para>
        /// </remarks>
        public double Curvature {
            get {
                return (f2);
            }
        }


        /// <summary>
        /// Gets the number of evaluations of the function that were required to isolate the extremum.
        /// </summary>
        public int EvaluationCount {
            get {
                return (count);
            }
        }

        internal LineExtremum Negate () {
            return (new LineExtremum(x, -f, f2, count));
        }

        /// <summary>
        /// Converts a line extremum to a one-dimensional space extremum.
        /// </summary>
        /// <param name="m">The line extremum.</param>
        /// <returns>The corresponding one-dimensional space extremum.</returns>
        public static implicit operator SpaceExtremum (LineExtremum m) {
            if (m == null) {
                return (null);
            } else {
                double[] s_x = new double[1] { m.x };
                double s_f = m.f;
                SymmetricMatrix s_f2 = new SymmetricMatrix(1);
                s_f2[0, 0] = m.f2;
                return (new SpaceExtremum(s_x, s_f, s_f2));
            }
        }
    }

}
