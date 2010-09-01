
using System;
using System.Globalization;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a value with an associated uncertainty.
    /// </summary>
	[Serializable]
	public struct UncertainValue {

		private double v;
		private double u;
	
        /// <summary>
        /// Gets the best estimate.
        /// </summary>
        public double Value {
			get {
				return(v);
			}
			internal set {
				v = value;
			}
		}
		
        /// <summary>
        /// Gets the uncertainty.
        /// </summary>
		public double Uncertainty {
			get {
				return(u);
			}
			internal set {
				if (value < 0.0) throw new InvalidOperationException();
				u = value;
			}
		}

        /// <summary>
        /// Gets the relative uncertainty.
        /// </summary>
		public double RelativeUncertainty {
			get {
				return(u/v);
			}
		}

        /// <summary>
        /// Returns a confidence interval.
        /// </summary>
        /// <param name="P">The required confidence level.</param>
        /// <returns>The associated confidence interval.</returns>
        /// <remarks><para>This method assumes </para></remarks>
        public Interval ConfidenceInterval (double P) {
            if ((P <= 0.0) || (P >= 1.0)) throw new ArgumentOutOfRangeException("P");
            double z = Math.Sqrt(2.0) * AdvancedMath.InverseErf(P);
            return (Interval.FromMidpointAndWidth(v, 2.0 * u * z));
        }


        /// <summary>
        /// Initializes a new uncertain value.
        /// </summary>
        /// <param name="value">The best estimate of the value.</param>
        /// <param name="uncertainty">The uncertainty in the value.</param>
		public UncertainValue (double value, double uncertainty) {
			if (uncertainty < 0.0) throw new ArgumentOutOfRangeException("uncertainty");
			v = value;
			u = uncertainty;
		}

		// printing
	
        /// <summary>
        /// Creates a string representation of the uncertain value.
        /// </summary>
        /// <returns>A string of the format <i>value</i> &#x00B1; <i>uncertainty</i>.</returns>
		public override string ToString() {
			return String.Format(CultureInfo.CurrentCulture, "{0} \u00B1 {1}", Value, Uncertainty);
        }

#if SHO
        public string __repr__ () {
            return(ToString());
        }
#endif

        // equality testing

        /// <summary>
        /// Determines whether two uncertain values are equal.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>True if the two uncertain values are equal, otherwise false.</returns>
        public static bool operator == (UncertainValue v1, UncertainValue v2) {
            if (Object.ReferenceEquals(v1, null)) {
                if (Object.ReferenceEquals(v2, null)) {
                    return (true);
                } else {
                    return (false);
                }
            } else {
                if (Object.ReferenceEquals(v2, null)) {
                    return (false);
                } else {
                    return ((v1.Value == v2.Value) && (v1.Uncertainty == v2.Uncertainty));
                }
            } 
        }

        /// <summary>
        /// Determines whether two uncertain values are not equal.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>True if the two uncertain values not are equal, otherwise false.</returns>
        public static bool operator != (UncertainValue v1, UncertainValue v2) {
            return (!(v1 == v2));
        }

        /// <summary>
        /// Determines whether the given object represents the same uncertain value.
        /// </summary>
        /// <param name="obj">The object.</param>
        /// <returns>True if the object represents the same reference point, otherwise false.</returns>
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
            return (v.GetHashCode() ^ u.GetHashCode());
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
			UncertainValue v = new UncertainValue();
			v.Value = v1.Value + v2.Value;
			if (v1.Uncertainty > v2.Uncertainty) {
				v.Uncertainty = v1.Uncertainty * Math.Sqrt( 1.0 + Math.Pow(v2.Uncertainty/v1.Uncertainty,2) );
            } else if (v2.Uncertainty > v1.Uncertainty) {
                v.Uncertainty = v2.Uncertainty * Math.Sqrt(1.0 + Math.Pow(v1.Uncertainty / v2.Uncertainty, 2));
            } else {
                // must handle this as a seperate case to avoid division by zero when v1.Uncertainty == v2.Uncertainty == 0
                v.Uncertainty = Math.Sqrt(2.0) * v1.Uncertainty;
            }
			return(v);
		}

        /// <summary>
        /// Subtracts two uncertain values.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>The difference of the two uncertain values.</returns>
        public static UncertainValue operator - (UncertainValue v1, UncertainValue v2) {
			UncertainValue v = new UncertainValue();
			v.Value = v1.Value - v2.Value;
			if (v1.Uncertainty > v2.Uncertainty) {
				v.Uncertainty = v1.Uncertainty * Math.Sqrt( 1.0 + Math.Pow(v2.Uncertainty/v1.Uncertainty,2) );
            } else if (v2.Uncertainty > v1.Uncertainty) {
                v.Uncertainty = v2.Uncertainty * Math.Sqrt(1.0 + Math.Pow(v1.Uncertainty / v2.Uncertainty, 2));
            } else {
                // must handle this as a seperate case to avoid division by zero when v1.Uncertainty == v2.Uncertainty == 0
                v.Uncertainty = Math.Sqrt(2.0) * v1.Uncertainty;
            }
			return(v);
		}

        /// <summary>
        /// Multiplies two uncertain values.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>The product of the two uncertain values.</returns>
        public static UncertainValue operator * (UncertainValue v1, UncertainValue v2) {
			UncertainValue v = new UncertainValue();
			v.Value = v1.Value * v2.Value;
			v.Uncertainty = Math.Abs(v.Value) * Math.Sqrt( Math.Pow(v1.RelativeUncertainty,2) + Math.Pow(v2.RelativeUncertainty,2) );
			return(v);
		}

        /// <summary>
        /// Divides two uncertain values.
        /// </summary>
        /// <param name="v1">The first uncertain value.</param>
        /// <param name="v2">The second uncertain value.</param>
        /// <returns>The quotient of the two uncertain values.</returns>
        public static UncertainValue operator / (UncertainValue v1, UncertainValue v2) {
			UncertainValue v = new UncertainValue();
			v.Value = v1.Value/v2.Value;
			v.Uncertainty = Math.Abs(v.Value) * Math.Sqrt( Math.Pow(v1.RelativeUncertainty,2) + Math.Pow(v2.RelativeUncertainty,2) );
			return(v);		
		}

        /// <summary>
        /// Adds an uncertain value to an certain value.
        /// </summary>
        /// <param name="v1">The certain value.</param>
        /// <param name="u2">The uncertain value.</param>
        /// <returns>The sum of the two values.</returns>
		public static UncertainValue operator + (double v1, UncertainValue u2) {
			return( new UncertainValue(v1 + u2.Value, u2.Uncertainty) );
		}

        /// <summary>
        /// Adds a certain value to an uncertain value.
        /// </summary>
        /// <param name="u1">The uncertain value.</param>
        /// <param name="v2">The certain value.</param>
        /// <returns>The sum of the two values.</returns>
        public static UncertainValue operator + (UncertainValue u1, double v2) {
			return( v2+u1 );
		}

        /// <summary>
        /// Subtracts a certain value from an uncertain value.
        /// </summary>
        /// <param name="u1">The uncertain value.</param>
        /// <param name="v2">The certain value.</param>
        /// <returns>The difference between the two values.</returns>
        public static UncertainValue operator - (UncertainValue u1, double v2) {
            return (u1 + (-v2));
        }

        /// <summary>
        /// Subtracts an uncertain value from a certain value.
        /// </summary>
        /// <param name="v1">The certain value.</param>
        /// <param name="u2">The uncertain vlaue.</param>
        /// <returns>The difference between the two values.</returns>
        public static UncertainValue operator - (double v1, UncertainValue u2) {
            return (new UncertainValue(v1 - u2.Value, u2.Uncertainty));
        }

        /// <summary>
        /// Multiplies a certain value by an uncertain value.
        /// </summary>
        /// <param name="v1">The certain value.</param>
        /// <param name="u2">The uncertain value.</param>
        /// <returns>The product of the two values.</returns>
		public static UncertainValue operator * (double v1, UncertainValue u2) {
			return( new UncertainValue(v1*u2.Value, Math.Abs(v1)*u2.Uncertainty) );
		}

        /// <summary>
        /// Multiplies an uncertain value by a certain value.
        /// </summary>
        /// <param name="u1">The uncertain value.</param>
        /// <param name="v2">The certain value.</param>
        /// <returns>The product of the two values.</returns>
		public static UncertainValue operator * (UncertainValue u1, double v2) {
			return( v2 * u1 );
		}

        /// <summary>
        /// Divides an uncertain value by a certain value.
        /// </summary>
        /// <param name="u1">The uncertain value.</param>
        /// <param name="v2">The certain value.</param>
        /// <returns>The quotient of the two values.</returns>
		public static UncertainValue operator / (UncertainValue u1, double v2) {
			return( new UncertainValue(u1.Value/v2, u1.Uncertainty/Math.Abs(v2)) );
		}

	}


    /// <summary>
    /// Contains methods for computing basic mathematical functions of uncertain values.
    /// </summary>
    /// <remarks><para>The methods in this static class perform the same basic mathematical operations as the methods of
    /// the <see cref="System.Math"/> class, accounting for the uncertainty in the inputs to produce a corresponding
    /// uncertainty in the output.</para>
    /// <para>As with operations on uncertain values, the methods assume that the uncertainty in input parameters represents the
    /// standard deviation of a distribution of measurements, and produce a value for the uncertainty in the output which
    /// represent a corresponding standard deviation, under the assumption that the standard deviations are small relative to
    /// the best values.</para></remarks>
	public static class UncertainMath {

		// basic functions of uncertain values

        /// <summary>
        /// Computes the square root of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The square root of the argument.</returns>
		public static UncertainValue Sqrt (UncertainValue x) {
			double v = Math.Sqrt(x.Value);
			double u = 0.5 * x.Uncertainty / v;
			return( new UncertainValue(v, u) );
		}

        /// <summary>
        /// Computes the sine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The sine of the argument.</returns>
		public static UncertainValue Sin (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Sin(x.Value);
			y.Uncertainty = Math.Abs( Math.Cos(x.Value) ) * x.Uncertainty;
			return(y);
		}

        /// <summary>
        /// Computes the cosine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The cosine of the argument.</returns>
		public static UncertainValue Cos (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Cos(x.Value);
			y.Uncertainty = Math.Abs ( Math.Sin(x.Value) ) * x.Uncertainty;
			return(y);
		}

        /// <summary>
        /// Computes the tangent of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The tanget of the argument.</returns>
		public static UncertainValue Tan (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Tan(x.Value);
			y.Uncertainty = ( 1.0 + Math.Pow(y.Value, 2) ) * x.Uncertainty;
			return(y);
		}

        /// <summary>
        /// Computes the arcsine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The arcsine of the argument.</returns>
		public static UncertainValue Asin (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Asin(x.Value);
			y.Uncertainty = x.Uncertainty / Math.Abs( Math.Cos(y.Value) );
			return(y);
		}

        /// <summary>
        /// Computes the arccosine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The arccosine of the argument.</returns>
		public static UncertainValue Acos (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Acos(x.Value);
			y.Uncertainty = x.Uncertainty / Math.Abs( Math.Sin(y.Value) );
			return(y);
		}

        /// <summary>
        /// Computes the arctangent of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The arctanget of the argument.</returns>
        public static UncertainValue Atan (UncertainValue x) {
			double xv = x.Value;
			double v = Math.Atan(xv);
			double u = x.Uncertainty/(1+xv*xv);
			return( new UncertainValue(v, u) );
		}

        /// <summary>
        /// Computes the arctangent of the ratio of two uncertain values.
        /// </summary>
        /// <param name="x">The argument of the numerator.</param>
        /// <param name="y">The argument of the denominator.</param>
        /// <returns>The arctangent of the quotient.</returns>
		public static UncertainValue Atan2 (UncertainValue x, UncertainValue y) {
			double xv = x.Value;
			double yv = y.Value;
			double v = Math.Atan2(xv, yv);
			double u1 = x.Uncertainty * yv;
			double u2 = y.Uncertainty * xv;
			double u = Math.Sqrt( u1*u1 + u2*u2) / (xv*xv + yv*yv);
			return( new UncertainValue(v, u) );
		}

        /// <summary>
        /// Computes e to the power of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of e^<sup>x1</sup>.</returns>
        public static UncertainValue Exp (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Exp(x.Value);
			y.Uncertainty = y.Value * x.Uncertainty;
			return(y);
		}

        /// <summary>
        /// Computes the natural logarithm of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of ln(x1).</returns>
        public static UncertainValue Log (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Log(x.Value);
			y.Uncertainty = x.RelativeUncertainty;
			return(y);
		}

        /// <summary>
        /// Computes an uncertain value raised to an arbitrary power.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <param name="p">The power.</param>
        /// <returns>The argument raised to the specified power.</returns>
        public static UncertainValue Pow (UncertainValue x, double p) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Pow(x.Value, p);
			y.Uncertainty = Math.Abs(p * y.Value * x.RelativeUncertainty );
			return(y);
		}

        /// <summary>
        /// Computes the hyperbolic sine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The hyperbolic sine of the argument.</returns>
        public static UncertainValue Sinh (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Sinh(x.Value);
			y.Uncertainty = Math.Cosh(x.Value) * x.Uncertainty;
			return(y);
		}

        /// <summary>
        /// Computes the hyperbolic cosine of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The hyperbolic cosine of the argument.</returns>
		public static UncertainValue Cosh (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Cosh(x.Value);
			y.Uncertainty = Math.Abs(Math.Sinh(x.Value)) * x.Uncertainty;
			return(y);
		}

        /// <summary>
        /// Computes the hyperbolic tangent of an uncertain value.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The hyperbolic tanget of the argument.</returns>
        public static UncertainValue Tanh (UncertainValue x) {
			UncertainValue y = new UncertainValue();
			y.Value = Math.Tanh(x.Value);
			y.Uncertainty = y.Uncertainty / Math.Pow(Math.Cosh(x.Value), 2);
			return(y);
		}

	}

}

