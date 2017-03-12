using System;
using System.Globalization;

using Meta.Numerics.Functions;

namespace Meta.Numerics {

    // In Debug mode, math operations (e.g. z1 + z2) are faster if we directly access the backing fields instead of using the
    // property accessors. I suspect that in Release mode, this difference is optimized away, but i haven't checked.

    /// <summary>Represents a complex number.</summary>
    /// <remarks>
    /// <para>Version 4.0 of the .NET Framework introduced a Complex structure equivalent to this one. To maintain compatibility
    /// with earlier versions of the .NET Framework, Meta.Numerics maintains its own Complex structure.</para>
    /// </remarks>
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
        }

        /// <summary>
        /// Gets the imaginary part of the complex number.
        /// </summary>
        public double Im {
            get {
                return (im);
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


        /// <summary>
        /// Converts a System.Numerics.Complex number to a Meta.Numerics.Complex number.
        /// </summary>
        /// <param name="value">The System.Numerics.Complex number.</param>
        /// <returns>The Meta.Numerics.Complex number.</returns>
        /// <remarks>
        /// <para>The <see cref="System.Numerics.Complex"/> data type has been available in the .NET Framework since version 4.0.
        /// Because Meta.Numerics offered its own <see cref="Complex"/> data type before the .NET Framework, and since even in
        /// the latest version of the .NET Framework, its Complex type has some notable deficiencies (e.g., compute
        /// <see cref="System.Numerics.Complex.Sqrt(System.Numerics.Complex)"/> of -1.0 and note that it does
        /// not equal <see cref="System.Numerics.Complex.ImaginaryOne"/>), our own Complex type persists. Eventually, we
        /// expect the deficiencies of <see cref="System.Numerics.Complex"/> to be corected. Until that time,
        /// to ease interoperation, we provide an implicit cast that converts the .NET Framework Complex type into
        /// the Meta.Numerics Complex type.</para>
        /// </remarks>
        //public static implicit operator Complex (System.Numerics.Complex value) {
        //    return (new Complex(value.Real, value.Imaginary));
        //}

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
        [CLSCompliant(false)]
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
			return( new Complex(-z.re, -z.im) );
		}

		// equality operators

        private static bool Equals (Complex z1, Complex z2) {
            return ((z1.re == z2.re) && (z1.im == z2.im));
        }

        /// <summary>
        /// Tests the equality of two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>True if the two complex numbers are equal, otherwise false.</returns>
        public static bool operator == (Complex z1, Complex z2) {
            return (Equals(z1, z2));
        }

        /// <summary>
        /// Tests the inequality of two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>False if the two complex numbers are equal, otherwise true.</returns>
        public static bool operator != (Complex z1, Complex z2) {
            return (!Equals(z1, z2));
       }

        /// <summary>
        /// Determines whether the given object represents the same complex number.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if the object represents the same complex number, otherwise false.</returns>
        public override bool Equals (object obj) {
            if (obj is Complex) {
                return (Equals(this, (Complex) obj));
            } else {
                return (false);
            }
        }

        /// <summary>
        /// Determines whether the given complex number is the same.
        /// </summary>
        /// <param name="other">The complex number to compare.</param>
        /// <returns>True if the complex number is the same, otherwise false.</returns>
        public bool Equals (Complex other) {
            return (Equals(this, other));
        }

        /// <summary>
        /// Returns a hash code for the complex number.
        /// </summary>
        /// <returns>A hash code.</returns>
        public override int GetHashCode () {
            return (Re.GetHashCode() ^ (Im.GetHashCode() * 31));
        }

        /// <summary>
        /// Adds two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>The sum of the complex numbers.</returns>
		public static Complex operator + (Complex z1, Complex z2) {
            return (new Complex(z1.re + z2.re, z1.im + z2.im));
		}

        /// <summary>
        /// Subtracts the second complex number from the first.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>The difference of the complex numbers.</returns>
		public static Complex operator - (Complex z1, Complex z2) {
            return (new Complex(z1.re - z2.re, z1.im - z2.im));
		}

        /// <summary>
        /// Multiplies two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>The product of the two complex numbers.</returns>
		public static Complex operator * (Complex z1, Complex z2) {
            return (new Complex(z1.re * z2.re - z1.im * z2.im, z1.re * z2.im + z1.im * z2.re));
		}

        /// <summary>
        /// Divides two complex numbers.
        /// </summary>
        /// <param name="z1">The first complex number.</param>
        /// <param name="z2">The second complex number.</param>
        /// <returns>The quotient of the two complex numbers.</returns>
		public static Complex operator / (Complex z1, Complex z2) {
            return (Divide(z1, z2));
		}

        private static Complex Divide (Complex z1, Complex z2) {

            // In math class we are taught
            //   \frac{a + i b}{c + i d} = \frac{(a + i b)(c - i d)}{(c + i d)(c - i d)} =
            //   \frac{(a c + b d) + i (b c  - a d)}{c^2 + d^2} 
            // This is a terrible formula to use numerically, because c^2 + d^2 can easily overflow.

            // Smith's 1962 aglorithm factors out the larger of c or d, e.g. if c is larger
            //   \frac{c (a + b d / c) + i c (b - a d / c)}{c (c + d d / c)} =
            //   \frac{(a + b r) + i (b - a r)}{c + d r}
            // where r = d / c. This is where most implementations of complex division stop.

            // In "A Robust Complex Division in Scilab" (2012) (https://arxiv.org/pdf/1210.4539v2.pdf),
            // Baudin & Smith note this algorithm also has problems when r underflows, which
            // can happen easily since it has been chosen to be small, and suggest in that
            // case changing the order of operations in the numerator terms.
            //   \frac{(a + d (b / c)) + i (b - d (a / c))}{c}
            // They give several examples where this improves the result and find that it
            // reduces by rate of problematic results over their space of test cases by
            // nearly two orders of magnitude over Smith's original method. Aside from a single
            // extra test, its common path is no less performant than Smith's original method.

            double re, im;
            if (Math.Abs(z2.im) <= Math.Abs(z2.re)) {
                Divide_Internal(z1.re, z1.im, z2.re, z2.im, out re, out im);
            } else {
                // The d > c version of the formula is obtained by a <-> b, c <-> d, f <-> -f.
                Divide_Internal(z1.im, z1.re, z2.im, z2.re, out re, out im);
                im = -im;
            }

            return (new Complex(re, im));
             
        }

        private static void Divide_Internal (double a, double b, double c, double d, out double e, out double f) {
            double r = d / c;
            double t = 1.0 / (c + d * r);
            const double normalLimit = 1.0 / Double.MaxValue;
            if (Math.Abs(r) > normalLimit) {
                e = (a + b * r) * t;
                f = (b - a * r) * t;
            } else {
                e = (a + d * (b / c)) * t;
                f = (b - d * (a / c)) * t;
            }
        }

        // Mixed double-complex binary operations

        // These are not strictly necessary, since we have defined an implicit double -> Complex cast,
        // but they are presumably faster than doing a cast and then an operation with zero im parts.

        /// <summary>
        /// Computes the sum of a complex and a real number.
        /// </summary>
        /// <param name="z">The complex number.</param>
        /// <param name="a">The real number.</param>
        /// <returns>The sum z + a.</returns>
		public static Complex operator+ (Complex z, double a) {
            // This is 1 flop instead of a cast plus 2 flops.
			return (new Complex(z.Re + a, z.Im));
		}

        /// <summary>
        /// Computes the sum of a real and a complex number.
        /// </summary>
        /// <param name="a">The real number.</param>
        /// <param name="z">The complex number.</param>
        /// <returns>The sum a + z.</returns>
		public static Complex operator+ (double a, Complex z) {
            // This is 1 flop instead a cast plus 2 flops.
            return (new Complex(a + z.Re, z.Im));
		}

        /*
		public static Complex operator- (Complex z, double a) {
			return( new Complex(z.Re -a, z.Im) );
		}

		public static Complex operator- (double a, Complex z) {
			return( new Complex(a - z.Re, -z.Im) );
		}
        */

        // Mixed complex/double arithmetic

        /// <summary>
        /// Multiplies a complex number by a real number.
        /// </summary>
        /// <param name="a">The real number.</param>
        /// <param name="z">The complex number.</param>
        /// <returns>The product az.</returns>
        public static Complex operator * (double a, Complex z) {
            // This is 2 flops instead of a cast plus 6 flops.
			return( new Complex(a * z.re, a * z.im) );
		}

        /// <summary>
        /// Multiplies a real number by a complex number.
        /// </summary>
        /// <param name="z">The complex number.</param>
        /// <param name="a">The real number.</param>
        /// <returns>The product za.</returns>
		public static Complex operator * (Complex z, double a) {
			return( a * z );
		}

        /*
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
        */

        /// <summary>
        /// Divides a complex number by a real number.
        /// </summary>
        /// <param name="z">The complex dividend.</param>
        /// <param name="a">The real divisor.</param>
        /// <returns>The quotient z / a.</returns>
		public static Complex operator / (Complex z, double a) {
            // This is 2 flops instead of a cast plus 9 (!) flops.
			return( new Complex(z.Re / a, z.Im / a) );
		}

        /// <summary>
        /// Gets the complex value of zero.
        /// </summary>
        public static readonly Complex Zero = new Complex(0.0, 0.0);

        /// <summary>
        /// Gets the complex value of one.
        /// </summary>
        public static readonly Complex One = new Complex(1.0, 0.0);

        /// <summary>
        /// Gets the square root of negative one.
        /// </summary>
        public static readonly Complex I = new Complex(0.0, 1.0);

        /// <summary>
        /// Determines if the given complex number is not-a-number.
        /// </summary>
        /// <param name="z">The complex number.</param>
        /// <returns><see langword="true"/> if the argument is not-a-number,
        /// otherwise <see langword="false"/>.</returns>
        /// <remarks>
        /// <para>The <see cref="Double.NaN"/> value is used to signal the result of a failed or impossible calculation,
        /// such as dividing zero by zero. A <see cref="Complex"/> value is considered to be not-a-number if
        /// either its real or imaginary part is not-a-number. This method tests for such an occurance.</para>
        /// </remarks>
        public static bool IsNaN (Complex z) {
            return (Double.IsNaN(z.Re) || Double.IsNaN(z.Im));
        }

	}

}
