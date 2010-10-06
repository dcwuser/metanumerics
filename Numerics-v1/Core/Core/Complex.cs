using System;
using System.Globalization;

using Meta.Numerics.Functions;

namespace Meta.Numerics {

    /// <summary>Represents a complex number.</summary>
    [Serializable]
    public struct Complex : IEquatable<Complex> {

        private double re;
        private double im;

        /// <summary>
        /// Gets the real part of the complex number.
        /// </summary>
        public double Re {
            get {
                return (re);
            }
            set {
                re = value;
            }
        }

        /// <summary>
        /// Gets the imaginary part of the complex number.
        /// </summary>
        public double Im {
            get {
                return (im);
            }
            set {
                im = value;
            }
        }

        /// <summary>
        /// Gets the complex conjugate of the complex number.
        /// </summary>
		public Complex Conjugate {
			get {
				return( new Complex(re,-im) );
			}
		}

        /// <summary>
        /// Initializes a new complex number.
        /// </summary>
        /// <param name="re">The real part of the complex number.</param>
        /// <param name="im">The imaginary part of the complex number.</param>
		public Complex (double re, double im) {
			this.re = re;
			this.im = im;
		}

		// conversions

        /// <summary>
        /// Converts the complex number to a double-precision real number.
        /// </summary>
        /// <param name="z">The complex number to covert.</param>
        /// <returns>The corresponding double-precision real number.</returns>
        /// <remarks><para>This explicit cast will fail if the complex number has a non-zero imaginary part.
        /// If you just want to obtain the real part of a complex number, use the <see cref="Re" /> property.</para></remarks>
        /// <exception cref="InvalidCastException">z.Im &#x2260; 0</exception>
		public static explicit operator double (Complex z) {
			if (z.Im != 0.0) throw new InvalidCastException("Complex casts to real must have vanishing imaginary parts.");
			return(z.Re);
		}

        /// <summary>
        /// Converts a double-precision real number to a complex number.
        /// </summary>
        /// <param name="x">The double-precision real number to convert.</param>
        /// <returns>The corresponding complex number.</returns>
        /// <remarks><para>The complex number output has a zero imaginary part and real part equal to the input number.</para>
        /// <para>This is an implicit cast; the compiler will apply it automatically whenever a real number is given in a situation
        /// where a complex number is required.</para></remarks>
		public static implicit operator Complex (double x) {
			return( new Complex(x,0.0) );
		}

		// printing

        /// <summary>
        /// Produces a string representation of the complex number.
        /// </summary>
        /// <returns>A string represenation of the complex number.</returns>
		public override string ToString () {
            return (String.Format(CultureInfo.CurrentCulture, "({0},{1})", re, im));
		}


#if SHO
        /// <summary>
        /// Produces the representation of the complex number for the Python interactive console.
        /// </summary>
        /// <returns>A string representation of the complex number.</returns>
        public string __repr__ () {
            return(ToString());
        }
#endif

		// equality operations are right by default

		// static unary operators

        /// <summary>
        /// Negates a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The argument times -1.</returns>
		public static Complex operator- (Complex z) {
			return( new Complex(-z.Re,-z.Im) );
		}

		// equality operators

        /// <summary>
        /// Tests the equality of two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>True if the two complex numbers are equal, otherwise false.</returns>
        public static bool operator == (Complex z1, Complex z2) {
            return ((z1.Re == z2.Re) && (z1.Im == z2.Im));
        }

        /// <summary>
        /// Tests the inequality of two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>False if the two complex numbers are equal, otherwise true.</returns>
        public static bool operator != (Complex z1, Complex z2) {
            return ((z1.Re != z2.Re) || (z1.Im != z2.Im));
        }

        /// <summary>
        /// Determines whether the given object represents the same complex number.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if the object represents the same complex number, otherwise false.</returns>
        public override bool Equals (object obj) {
            return(((Complex) obj) == this);
        }

        /// <summary>
        /// Determines whether the given complex number is the same.
        /// </summary>
        /// <param name="other">The complex number to compare.</param>
        /// <returns>True if the complex number is the same, otherwise false.</returns>
        public bool Equals (Complex other) {
            return (this == other);
        }

        /// <summary>
        /// Returns a hash code for the complex number.
        /// </summary>
        /// <returns>A hash code.</returns>
        public override int GetHashCode () {
            return (Re.GetHashCode() ^ Im.GetHashCode());
        }

        /// <summary>
        /// Adds two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>The sum of the complex numbers.</returns>
		public static Complex operator+ (Complex z1, Complex z2) {
			return( new Complex(z1.Re+z2.Re, z1.Im+z2.Im) );
		}

        /// <summary>
        /// Subtracts the second complex number from the first.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>The difference of the complex numbers.</returns>
		public static Complex operator- (Complex z1, Complex z2) {
			return( new Complex(z1.Re-z2.Re, z1.Im-z2.Im) );
		}

        /// <summary>
        /// Multiplies two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>The product of the two complex numbers.</returns>
		public static Complex operator* (Complex z1, Complex z2) {
			return( new Complex(z1.Re*z2.Re - z1.Im*z2.Im, z1.Re*z2.Im + z1.Im*z2.Re) );
		}

        /// <summary>
        /// Divides two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>The quotient of the two complex numbers.</returns>
		public static Complex operator/ (Complex z1, Complex z2) {
			if (Math.Abs(z2.Re) > Math.Abs(z2.Im)) {
				double x = z2.Im/z2.Re;
				double w = z2.Re + x*z2.Im;
				return( new Complex((z1.Re + x*z1.Im)/w, (z1.Im-x*z1.Re)/w) );
			} else {
				double x = z2.Re/z2.Im;
				double w = x*z2.Re + z2.Im;
				return( new Complex((x*z1.Re+z1.Im)/w, (x*z1.Im-z1.Re)/w) );
			}
		}

		// mixed double-complex binary operations
		
		// these are not strictly necessary, since we have
		// defined an implicit double->Complex converter

        /*
		public static Complex operator+ (Complex z, double a) {
			return( new Complex(z.Re + a, z.Im) );
		}

		public static Complex operator+ (double a, Complex z) {
			return( z + a );
		}

		public static Complex operator- (Complex z, double a) {
			return( new Complex(z.Re -a, z.Im) );
		}

		public static Complex operator- (double a, Complex z) {
			return( new Complex(a - z.Re, -z.Im) );
		}

		public static Complex operator* (double a, Complex z) {
			return( new Complex(a*z.Re, a*z.Im) );
		}

		public static Complex operator* (Complex z, double a) {
			return( a * z );
		}

		public static Complex operator/ (double a, Complex z) {
			if (Math.Abs(z.Re)>Math.Abs(z.Im)) {
				double x = z.Im/z.Re;
				double w = z.Re + x*z.Im;
				return( new Complex(a/w, -a*x/w) );
			} else {
				double x = z.Re/z.Im;
				double w = x*z.Re + z.Im;
				return( new Complex(a*x/w, -a/w) );
			}
		}

		public static Complex operator/ (Complex z, double a) {
			return( new Complex(z.Re/a, z.Im/a) );
		}
        */

	}

}
