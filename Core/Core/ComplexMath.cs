using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics {

    /// <summary>
    /// Provides simple functions of complex arguments. 
    /// </summary>
    public static class ComplexMath {

        // static pure imaginary

        /// <summary>
        /// Gets the unit imaginary number I.
        /// </summary>
        public static Complex I {
            get {
                return (new Complex(0.0, 1.0));
            }
        }

        // basic functions of complex arguments

        /// <summary>
        /// Computes the absolute value of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of |z|.</returns>
        /// <remarks>
        /// <para>The absolute value of a complex number is the distance of the number from the origin
        /// in the complex plane. This is a compatible generalization of the definition of the absolute
        /// value of a real number.</para>
        /// </remarks>
        /// <seealso cref="Math.Abs(Double)"/>
        public static double Abs (Complex z) {
            return (MoreMath.Hypot(z.Re, z.Im));
        }

        /// <summary>
        /// Computes the phase of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of arg(z).</returns>
        /// <remarks>
        /// <para>The phase of a complex number is the angle between the line joining it to the origin and the real axis of the complex plane.</para>
        /// <para>The phase of complex numbers in the upper complex plane lies between 0 and &#x3C0;. The phase of complex numbers
        /// in the lower complex plane lies between 0 and -&#x3C0;. The phase of a real number is zero.</para>
        /// </remarks>
        public static double Arg (Complex z) {
            // returns 0 to PI in the upper complex plane (Im>=0),
            // 0 to -PI in the lower complex plane (Im<0)
            return (Math.Atan2(z.Im, z.Re));
        }

        /// <summary>
        /// Computes e raised to the power of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of e<sup>z</sup>.</returns>
        public static Complex Exp (Complex z) {
            double m = Math.Exp(z.Re);
            return (new Complex(m * AdvancedMath.Cos(z.Im, 0.0), m * AdvancedMath.Sin(z.Im, 0.0)));
        }

        /// <summary>
        /// Computes the natrual logarithm of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of ln(z).</returns>
        public static Complex Log (Complex z) {
            return (new Complex(Math.Log(Abs(z)), Arg(z)));
        }

        /// <summary>
        /// Computes the square root of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The square root of the argument.</returns>
        public static Complex Sqrt (Complex z) {
            if (z.Im == 0) {
                return (Math.Sqrt(z.Re));
            } else {
                return (Pow(z, 0.5));
            }
        }

        /// <summary>
        /// Computes the sine of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of sin(z).</returns>
        public static Complex Sin (Complex z) {
            double p = Math.Exp(z.Im);
            double q = 1 / p;
            double sinh = (p - q) / 2.0;
            double cosh = (p + q) / 2.0;
            return (new Complex(Math.Sin(z.Re) * cosh, Math.Cos(z.Re) * sinh));
        }

        /// <summary>
        /// Computes the hyperbolic sine of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of sinh(z).</returns>
        public static Complex Sinh (Complex z) {
            // sinh(z) = -i sin(iz)
            Complex sin = Sin(new Complex(-z.Im, z.Re));
            return (new Complex(sin.Im, -sin.Re));
        }

        /// <summary>
        /// Computes the cosine of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of cos(z).</returns>
        public static Complex Cos (Complex z) {
            double p = Math.Exp(z.Im);
            double q = 1 / p;
            double sinh = (p - q) / 2.0;
            double cosh = (p + q) / 2.0;
            return (new Complex(Math.Cos(z.Re) * cosh, -Math.Sin(z.Re) * sinh));
        }

        /// <summary>
        /// Computes the hyperbolic cosine of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of cosh(z).</returns>
        public static Complex Cosh (Complex z) {
            // cosh(z) = cos(iz)
            return (Cos(new Complex(-z.Im, z.Re)));
        }

        /// <summary>
        /// Computes the tangent of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of tan(z).</returns>
        public static Complex Tan (Complex z) {
            // tan z = [sin(2x) + I sinh(2y)]/[cos(2x) + I cosh(2y)]
            double x2 = 2.0 * z.Re;
            double y2 = 2.0 * z.Im;
            double p = Math.Exp(y2);
            double q = 1 / p;
            double cosh = (p + q) / 2.0;
            if (Math.Abs(z.Im) < 4.0) {
                double sinh = (p - q) / 2.0;
                double D = Math.Cos(x2) + cosh;
                return (new Complex(Math.Sin(x2) / D, sinh / D));
            } else {
                // when Im(z) gets too large, sinh and cosh individually blow up
                // but ratio is still ~1, so rearrage to use tanh instead
                double F = (1.0 + Math.Cos(x2) / cosh);
                return (new Complex(Math.Sin(x2) / cosh / F, Math.Tanh(y2) / F));
            }
        }

        /// <summary>
        /// Computes the hyperbolic tangent of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of tanh(z).</returns>
        public static Complex Tanh (Complex z) {
            return (Sinh(z) / Cosh(z));
        }

        // pure complex binary operations

        /// <summary>
        /// Raises a complex number to an arbitrary real power.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <param name="p">The power.</param>
        /// <returns>The value of z<sup>p</sup>.</returns>
        public static Complex Pow (Complex z, double p) {
            double m = Math.Pow(Abs(z), p);
            double t = Arg(z) * p;
            return (new Complex(m * Math.Cos(t), m * Math.Sin(t)));
        }

        /// <summary>
        /// Raises a real number to an arbitrary complex power.
        /// </summary>
        /// <param name="x">The real base, which must be non-negative.</param>
        /// <param name="z">The complex exponent.</param>
        /// <returns>The value of x<sup>z</sup>.</returns>
        public static Complex Pow (double x, Complex z) {
            if (x < 0.0) throw new ArgumentOutOfRangeException("x");
            if (z == 0.0) return (1.0);
            if (x == 0.0) return (0.0);
            double m = Math.Pow(x, z.Re);
            double t = Math.Log(x) * z.Im;
            return (new Complex(m * Math.Cos(t), m * Math.Sin(t)));
        }

        /// <summary>
        /// Raises a complex number to an integer power.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <param name="n">The power.</param>
        /// <returns>The value of z<sup>n</sup>.</returns>
        public static Complex Pow (Complex z, int n) {

            // this is a straight-up copy of MoreMath.Pow with x -> z, double -> Complex

            if (n < 0) return (1.0 / Pow(z, -n));

            switch (n) {
                case 0:
                    // we follow convention that 0^0 = 1
                    return (1.0);
                case 1:
                    return (z);
                case 2:
                    // 1 multiply
                    return (z * z);
                case 3:
                    // 2 multiplies
                    return (z * z * z);
                case 4: {
                        // 2 multiplies
                        Complex z2 = z * z;
                        return (z2 * z2);
                    }
                case 5: {
                        // 3 multiplies
                        Complex z2 = z * z;
                        return (z2 * z2 * z);
                    }
                case 6: {
                        // 3 multiplies
                        Complex z2 = z * z;
                        return (z2 * z2 * z2);
                    }
                case 7: {
                        // 4 multiplies
                        Complex z3 = z * z * z;
                        return (z3 * z3 * z);
                    }
                case 8: {
                        // 3 multiplies
                        Complex z2 = z * z;
                        Complex z4 = z2 * z2;
                        return (z4 * z4);
                    }
                case 9: {
                        // 4 multiplies
                        Complex z3 = z * z * z;
                        return (z3 * z3 * z3);
                    }
                case 10: {
                        // 4 multiplies
                        Complex z2 = z * z;
                        Complex z4 = z2 * z2;
                        return (z4 * z4 * z2);
                    }
                default:
                    return (ComplexMath.Pow(z, (double)n));
            }

        }


    }
}
