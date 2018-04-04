using System;
using System.Diagnostics;

namespace Meta.Numerics.Analysis {

    /// <summary>
    /// Represents a maximum or minimum of a function of one variable.
    /// </summary>
    public sealed class Extremum : EvaluationResult {

        internal Extremum (double x, double f, double f2, double a, double b, int count, ExtremumSettings settings) : base(count) {
            Debug.Assert(settings != null);
            this.x = x;
            this.f = f;
            this.f2 = f2;
            this.a = a;
            this.b = b;
            this.settings = settings;
        }

        private readonly double x;
        private readonly double f;
        private readonly double f2;
        private readonly double a, b;
        private readonly ExtremumSettings settings;

        /// <summary>
        /// Gets the location (x-value) of the optimum.
        /// </summary>
        /// <remarks>
        /// <para>Note that numerical methods for finding typical a maximum or minimum cannot usually determine
        /// its location to full floating-point precision. Near a quadratic optimum, a change in x of ~&#x3B5; will
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
        /// Gets the function value (y-value) at the optimum.
        /// </summary>
        public double Value {
            get {
                return (f);
            }
        }

        /// <summary>
        /// Gets a bracket indicating how accurately the location of the optimum was determined.
        /// </summary>
        public Interval Bracket {
            get {
                return (Interval.FromEndpoints(a, b));
            }
        }

        /// <summary>
        /// Gets the curvature at the optimum.
        /// </summary>
        /// <remarks>
        /// <para>The curvature is the second derivative of the function at the optimum.</para>
        /// <para>At a typical optimum, where the function has vanishing first derivative, the second derivative will
        /// be a number whose magnitude characterizes the "steepness" with which the function increases as one moves
        /// away from the optimum.</para>
        /// <para>At an atypical optimum, for example at an interval boundary or of a non-smooth function, this
        /// value may be meaningless.</para>
        /// <para>Even in the case of a typical optimum, the value of the curvature property will often be accurate only
        /// to a few digits. If you require a highly accurate determination of the curvature,
        /// you should use numerical differentiation to find the curvature more accurately.</para>
        /// </remarks>
        public double Curvature {
            get {
                return (f2);
            }
        }

        /// <summary>
        /// Gets the settings that were used during optimization.
        /// </summary>
        public ExtremumSettings Settings {
            get {
                return (settings);
            }
        }


        internal Extremum Negate () {
            return (new Extremum(x, -f, f2, a, b, base.EvaluationCount, settings));
        }

    }

}
