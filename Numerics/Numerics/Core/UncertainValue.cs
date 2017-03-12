
using System;
using System.Globalization;

using Meta.Numerics.Functions;

namespace Meta.Numerics {

    /// <summary>
    /// Represents a value with an associated uncertainty.
    /// </summary>
    public struct UncertainValue : IEquatable<UncertainValue> {

		private readonly double v;
		private readonly double u;
	
        /// <summary>
        /// Gets the best estimate.
        /// </summary>
        public double Value {
			get {
				return(v);
			}
		}
		
        /// <summary>
        /// Gets the uncertainty.
        /// </summary>
		public double Uncertainty {
			get {
				return(u);
			}
		}

        /// <summary>
        /// Gets the relative uncertainty.
        /// </summary>
		public double RelativeUncertainty {
			get {
				return(u / Math.Abs(v));
			}
		}

        /// <summary>
        /// Returns a confidence interval.
        /// </summary>
        /// <param name="P">The required confidence level, which must lie in [0, 1].</param>
        /// <returns>The associated confidence interval.</returns>
        /// <remarks><para>This method assumes that the value is normally distributed with a mean equal to <see cref="Value"/>
        /// and a standard deviation equal to <see cref="Uncertainty"/>.</para></remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="P"/> lies outside [0, 1].</exception>
        public Interval ConfidenceInterval (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            double z = Math.Sqrt(2.0) * AdvancedMath.InverseErf(P);
            return (Interval.FromMidpointAndWidth(v, 2.0 * u * z));
        }


        /// <summary>
        /// Initializes a new uncertain value.
        /// </summary>
        /// <param name="value">The best estimate of the value.</param>
        /// <param name="uncertainty">The uncertainty in the value.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="uncertainty"/> is negative.</exception>
		public UncertainValue (double value, double uncertainty) {
			if (uncertainty < 0.0) throw new ArgumentOutOfRangeException(nameof(uncertainty));
			v = value;
			u = uncertainty;
		}

		// printing
	
        /// <summary>
        /// Creates a string representation of the uncertain value.
        /// </summary>
        /// <returns>A string of the format <i>value</i> &#x00B1; <i>uncertainty</i>.</returns>
		public override string ToString() {
			return String.Format(CultureInfo.CurrentCulture, "{0} \u00B1 {1}", v, u);
        }

        // equality testing

        /// <summary>
        /// Determines whether two uncertain values are equal.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns><see langword="true"/> if the two uncertain values are equal, otherwise <see langword="false"/>.</returns>
        public static bool operator == (UncertainValue v1, UncertainValue v2) {
            // dont' need to check for nulls because UncertainValue is a structure
            return ((v1.Value == v2.Value) && (v1.Uncertainty == v2.Uncertainty));
        }

        /// <summary>
        /// Determines whether two uncertain values are not equal.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns><see langword="true"/> if the two uncertain values are unequal, otherwise <see langword="false"/>.</returns>
        public static bool operator != (UncertainValue v1, UncertainValue v2) {
            return (!(v1 == v2));
        }

        /// <summary>
        /// Determines whether the other uncertain value is the same.
        /// </summary>
        /// <param name="other">The other uncertain value</param>
        /// <returns><see langword="true"/> if <paramref name="other"/> represents the same uncertain value,
        /// otherwise <see langword="false"/>.</returns>
        public bool Equals (UncertainValue other) {
            return (this == other);
        }

        /// <summary>
        /// Determines whether the given object represents the same uncertain value.
        /// </summary>
        /// <param name="obj">The object.</param>
        /// <returns><see langword="true"/> if <paramref name="obj"/> represents the same uncertain value, otherwise <see langword="false"/>.</returns>
        public override bool Equals (object obj) {
            if (obj is UncertainValue) {
                UncertainValue uv = (UncertainValue) obj;
                return (this == uv);
            } else {
                return (false);
            }
        }

        /// <summary>
        /// Computes a hash code for the uncertain value.
        /// </summary>
        /// <returns>A hash code.</returns>
        public override int GetHashCode () {
            return (v.GetHashCode() + 31 * u.GetHashCode());
        }


		// operations with uncertain values

        /// <summary>
        /// Negates an uncertain value.
        /// </summary>
        /// <param name="x">The uncertain value.</param>
        /// <returns>The negative of the uncertain value.</returns>
		public static UncertainValue operator- (UncertainValue x) {
			return( new UncertainValue(-x.Value, x.Uncertainty) );
		}

		// binary operations assume no correlation

        /// <summary>
        /// Adds two uncertain values.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>The sum of the two uncertain values.</returns>
		public static UncertainValue operator+ (UncertainValue v1, UncertainValue v2) {
            return (new UncertainValue(v1.Value + v2.Value, MoreMath.Hypot(v1.Uncertainty, v2.Uncertainty)));
		}

        /// <summary>
        /// Subtracts two uncertain values.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>The difference of the two uncertain values.</returns>
        public static UncertainValue operator- (UncertainValue v1, UncertainValue v2) {
            return(new UncertainValue(v1.Value - v2.Value, MoreMath.Hypot(v1.Uncertainty, v2.Uncertainty)));
		}

        /// <summary>
        /// Multiplies two uncertain values.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>The product of the two uncertain values.</returns>
        public static UncertainValue operator* (UncertainValue v1, UncertainValue v2) {
			double value = v1.Value * v2.Value;
			double uncertainty = Math.Abs(value) * MoreMath.Hypot(v1.RelativeUncertainty, v2.RelativeUncertainty);
            return (new UncertainValue(value, uncertainty));
		}

        /// <summary>
        /// Divides two uncertain values.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>The quotient of the two uncertain values.</returns>
        public static UncertainValue operator/ (UncertainValue v1, UncertainValue v2) {
			double value = v1.Value / v2.Value;
            double uncertainty = Math.Abs(value) * MoreMath.Hypot(v1.RelativeUncertainty, v2.RelativeUncertainty);
			return(new UncertainValue(value, uncertainty));		
		}

        /// <summary>
        /// Adds an uncertain value to an certain value.
        /// </summary>
        /// <param name="v1">The certain value.</param>
        /// <param name="u2">The uncertain value.</param>
        /// <returns>The sum of the two values.</returns>
		public static UncertainValue operator+ (double v1, UncertainValue u2) {
			return( new UncertainValue(v1 + u2.Value, u2.Uncertainty) );
		}

        /// <summary>
        /// Adds a certain value to an uncertain value.
        /// </summary>
        /// <param name="u1">The uncertain value.</param>
        /// <param name="v2">The certain value.</param>
        /// <returns>The sum of the two values.</returns>
        public static UncertainValue operator+ (UncertainValue u1, double v2) {
			return(v2 + u1);
		}

        /// <summary>
        /// Subtracts a certain value from an uncertain value.
        /// </summary>
        /// <param name="u1">The uncertain value.</param>
        /// <param name="v2">The certain value.</param>
        /// <returns>The difference between the two values.</returns>
        public static UncertainValue operator- (UncertainValue u1, double v2) {
            return (u1 + (-v2));
        }

        /// <summary>
        /// Subtracts an uncertain value from a certain value.
        /// </summary>
        /// <param name="v1">The certain value.</param>
        /// <param name="u2">The uncertain vlaue.</param>
        /// <returns>The difference between the two values.</returns>
        public static UncertainValue operator- (double v1, UncertainValue u2) {
            return (new UncertainValue(v1 - u2.Value, u2.Uncertainty));
        }

        /// <summary>
        /// Multiplies a certain value by an uncertain value.
        /// </summary>
        /// <param name="v1">The certain value.</param>
        /// <param name="u2">The uncertain value.</param>
        /// <returns>The product of the two values.</returns>
		public static UncertainValue operator* (double v1, UncertainValue u2) {
			return (new UncertainValue(v1 * u2.Value, Math.Abs(v1) * u2.Uncertainty));
		}

        /// <summary>
        /// Multiplies an uncertain value by a certain value.
        /// </summary>
        /// <param name="u1">The uncertain value.</param>
        /// <param name="v2">The certain value.</param>
        /// <returns>The product of the two values.</returns>
		public static UncertainValue operator* (UncertainValue u1, double v2) {
			return (v2 * u1);
		}

        /// <summary>
        /// Divides an uncertain value by a certain value.
        /// </summary>
        /// <param name="u1">The uncertain value.</param>
        /// <param name="v2">The certain value.</param>
        /// <returns>The quotient of the two values.</returns>
		public static UncertainValue operator/ (UncertainValue u1, double v2) {
			return (new UncertainValue(u1.Value / v2, u1.Uncertainty / Math.Abs(v2)));
		}

	}

}

