using System;
using System.Diagnostics;
using System.Runtime.CompilerServices;

using Meta.Numerics.Functions;

namespace Meta.Numerics {

    /// <summary>
    /// Contains methods that compute basic functions of complex arguments. 
    /// </summary>
    public static class ComplexMath {

        // static pure imaginary

        /// <summary>
        /// Gets the unit imaginary number I.
        /// </summary>
        /// <value>The unit imaginary number.</value>
        /// <seealso href="https://en.wikipedia.org/wiki/Imaginary_unit"/>
        public static Complex I {
            get {
                return (Complex.I);
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
        /// in the lower complex plane lies between 0 and -&#x3C0;. The phase of a positive real number is zero.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Argument_(complex_analysis)"/>
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
            return (new Complex(m * MoreMath.Cos(z.Im), m * MoreMath.Sin(z.Im)));
        }

        /// <summary>
        /// Computes the natrual logarithm of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of ln(z).</returns>
        /// <remarks>
        /// <para>The image below shows the complex log function near the origin, using domain coloring.</para>
        /// <img src="../images/ComplexLogPlot.png" />
        /// <para>You can see the zero at (0, 1) and the branch cut extending along the negative real axis from the pole at the origin.</para>
        /// </remarks>
        public static Complex Log (Complex z) {
            return (new Complex(Math.Log(Abs(z)), Arg(z)));
        }


        /// <summary>
        /// Computes the square of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of z<sup>2</sup>.</returns>
        /// <remarks>
        /// <para>Unlike <see cref="MoreMath.Sqr(double)"/>, there is slightly more to this method than shorthand for z * z.
        /// In terms of real and imaginary parts z = x + i y, the product z * z = (x * x - y * y) + i(x * y + x * y),
        /// which not only requires 6 flops to evaluate, but also computes the real part as an expression that can easily overflow
        /// or suffer from significant cancelation error. By instead computing z<sup>2</sup> via as (x - y) * (x + y) + i 2 * x * y,
        /// this method not only requires fewer flops but is also less subject to overflow and cancelation error. You should
        /// therefore generally favor the computation of z<sup>2</sup> using this method over its computation as z * z.</para> 
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Complex Sqr (Complex z) {
            // This form has one less flop than z * z, and, more importantly, evaluates
            // (x - y) * (x + y) instead of x * x - y * y; the latter would be more
            // subject to cancelation errors and overflow or underflow 
            return (new Numerics.Complex((z.Re - z.Im) * (z.Re + z.Im), 2.0 * z.Re * z.Im));
        }

        /// <summary>
        /// Computes the square root of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The square root of the argument.</returns>
        /// <remarks>
        /// <para>The image below shows the complex square root function near the origin, using domain coloring.</para>
        /// <img src="../images/ComplexSqrtPlot.png" />
        /// <para>You can see the branch cut extending along the negative real axis from the zero at the origin.</para>
        /// </remarks>
        public static Complex Sqrt (Complex z) {

            // Handle the degenerate case quickly.
            if (z.Im == 0.0) {
                if (z.Re < 0.0) {
                    return (new Complex(0.0, Math.Sqrt(-z.Re)));
                } else {
                    return (new Complex(Math.Sqrt(z.Re), 0.0));
                }
            }
            // This also eliminates need to worry about Im(z) = 0 in subsequent formulas.

            // One way to compute Sqrt(z) is just to call Pow(z, 0.5), which coverts to polar coordinates
            // (sqrt + atan), halves the phase, and reconverts to cartesian coordinates (cos + sin).
            // Not only is this more expensive than necessary, it also fails to preserve certain expected
            // symmetries, such as that the square root of a pure negative is a pure imaginary, and that the
            // square root of a pure imaginary has exactly equal real and imaginary parts. This all goes
            // back to the fact that Math.PI is not stored with infinite precision, so taking half of Math.PI
            // does not land us on an argument with sine exactly equal to zero.

            // To find a fast and symmetry-respecting formula for complex square root,
            // note x + i y = \sqrt{a + i b} implies x^2 + 2 i x y - y^2 = a + i b,
            // so x^2 - y^2 = a and 2 x y = b. Cross-substitute and use the quadratic formula to obtain
            //   x = \sqrt{\frac{\sqrt{a^2 + b^2} + a}{2}}  y = \pm \sqrt{\frac{\sqrt{a^2 + b^2} - a}{2}}
            // There is just one complication: depending on the sign on a, either x or y suffers from
            // cancelation when |b| << |a|. We can get aroud this by noting that our formulas imply
            // x^2 y^2 = b^2 / 4, so |x| |y| = |b| / 2. So after computing the one that doesn't suffer
            // from cancelation, we can compute the other with just a division. This is basically just
            // the right way to evaluate the quadratic formula without cancelation.

            // All this reduces our total cost to two sqrts and a few flops, and it respects the desired
            // symmetries. Much better than atan + cos + sin!

            // The signs are a matter of choice of branch cut, which is traditionally taken so x > 0 and sign(y) = sign(b).

            double a = z.Re;
            double b = z.Im;

            // If the components are too large, Hypot(a, b) will overflow, even though the subsequent sqrt would
            // make the result representable. To avoid this, we re-scale (by exact powers of 2 for accuracy)
            // when we encounter very large components to avoid intermediate infinities.
            bool rescale = false;
            if ((Math.Abs(a) >= sqrtRescaleThreshold) || (Math.Abs(b) >= sqrtRescaleThreshold)) {
                if (Double.IsInfinity(b) && !Double.IsNaN(a)) {
                    return (new Complex(Double.PositiveInfinity, b));
                }
                a *= 0.25;
                b *= 0.25;
                rescale = true;
            }
            // 1. Note that, if one component is very large and the other is very small, this re-scale could cause
            // the very small component to underflow. This is not a problem, though, because the neglected term
            // is ~ (very small) / (very large), which would underflow anyway.
            // 2. Note that our test for re-scaling requires two abs and two floating point comparisons.
            // For the typical non-overflowng case, it might be faster to just test for IsInfinity after Hypot,
            // then re-scale and re-do Hypot if overflow is detected. That version makes the code more
            // complicated and simple profiling experiments show no discernible perf improvement, so
            // stick with this version for now.

            double x, y;
            if (a >= 0.0) {
                x = Math.Sqrt((MoreMath.Hypot(a, b) + a) / 2.0);
                y = b / (2.0 * x);
            } else {
                y = Math.Sqrt((MoreMath.Hypot(a, b) - a) / 2.0);
                if (b < 0.0) y = -y;
                x = b / (2.0 * y);
            }
            // Note that, because we have re-scaled when any components are very large,
            // (2.0 * x) and (2.0 * y) are guaranteed not to overflow.

            if (rescale) {
                x *= 2.0;
                y *= 2.0;
            }


            return (new Complex(x, y));

        }

        // This is the largest value x for which (Hypot(x,x) + x) does not overflow.
        private static readonly double sqrtRescaleThreshold = Double.MaxValue / (1.0 + Math.Sqrt(2.0));

        /// <summary>
        /// Computes the sine of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of sin(z).</returns>
        public static Complex Sin (Complex z) {
            // Sine of a complex mixes sin, cos, sinh, and cosh. By computing sinh and cosh ourselves
            // from exp, we prevent having to evaluate exp twice inside seperate calls to sinh and cosh.
            double p = Math.Exp(z.Im);
            double q = 1.0 / p;
            double sinh = 0.5 * (p - q);
            double cosh = 0.5 * (p + q);
            return (new Complex(MoreMath.Sin(z.Re) * cosh, MoreMath.Cos(z.Re) * sinh));
            // For large z.Im, it's pretty easy for cosh and sinh to overflow. There is a very tiny space
            // of arguments for which z.Im causes cosh and sinh barely overflow, but the sin and cos of z.Re
            // are small enough to bring the result back into the representable range. We don't
            // handle this, so we return inf for those values.
        }

        internal static Complex SinPi (Complex z) {
            double y = Math.PI * z.Im;
            double p = Math.Exp(y);
            double q = 1.0 / p;
            double sinh = 0.5 * (p - q);
            double cosh = 0.5 * (p + q);
            return (new Complex(MoreMath.SinPi(z.Re) * cosh, MoreMath.CosPi(z.Re) * sinh));
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
            double q = 1.0 / p;
            double sinh = (p - q) / 2.0;
            double cosh = (p + q) / 2.0;
            return (new Complex(MoreMath.Cos(z.Re) * cosh, -MoreMath.Sin(z.Re) * sinh));
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
            // tan x = sin z / cos z, but to avoid unnecessary repeated trig computations, use
            //   tan z = \frac{\sin(2x) + i \sinh(2y)}{\cos(2x) + \cosh(2y)}
            // from A&S 4.3.57, and compute trig functions here.

            // Even this can only be used for |y| not-too-big, beacuse sinh and cosh (correctly) overflow
            // for quite moderate y, even though their ratio does not. In that case, divide
            // through by cosh to get:
            //   tan z = \frac{\frac{\sin(2x)}{\cosh(2y)} + i \tanh(2y)}{1 + \frac{\cos(2x)}{\cosh(2y)}}
            // which correctly computes the (tiny) real part and the (normal-sized) imaginary part.

            double x2 = 2.0 * z.Re;
            double y2 = 2.0 * z.Im;
            double p = Math.Exp(y2);
            double q = 1.0 / p;
            double cosh = (p + q) / 2.0;
            if (Math.Abs(z.Im) < 4.0) {
                double sinh = (p - q) / 2.0;
                double D = MoreMath.Cos(x2) + cosh;
                return (new Complex(MoreMath.Sin(x2) / D, sinh / D));
            } else {
                // when Im(z) gets too large, sinh and cosh individually blow up
                // but ratio is still ~1, so rearrage to use tanh instead
                double D = (1.0 + Math.Cos(x2) / cosh);
                return (new Complex(MoreMath.Sin(x2) / cosh / D, Math.Tanh(y2) / D));
            }
        }

        /// <summary>
        /// Computes the hyperbolic tangent of a complex number.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of tanh(z).</returns>
        public static Complex Tanh (Complex z) {
            // tanh(z) = -i tan(iz)
            Complex tan = Tan(new Complex(-z.Im, z.Re));
            return (new Complex(tan.Im, -tan.Re));
        }

        // Using sin(w) = \frac{e^{iw} - e^{-iw}}{2i}, write z = sin(w) and e^{iw} = v and solve for z.
        //    v^2 - 2 i v z - 1 = 0
        //    v = i z \pm \sqrt{1-z^2}
        //    w = -i \ln (i z \pm \sqrt{1-z^2})
        // Basically, we just need to compute i z \pm \sqrt{1-z^2}. The log determines its magnitude and
        // phase, multiplying by i switches their places.

        /// <summary>
        /// Computes the inverse sine (arcsine) of a complex number.
        /// </summary>
        /// <param name="z">The number.</param>
        /// <returns>The value of arcsin(z).</returns>
        public static Complex Asin (Complex z) {

            Asin_Internal(Math.Abs(z.Re), Math.Abs(z.Im), out double b, out double bPrime, out double q);

            double p;
            if (bPrime < 0.0) {
                p = Math.Asin(b);
            } else {
                p = Math.Atan(bPrime);
            }
            if (z.Re < 0.0) p = -p;

            if (z.Im < 0.0) q = -q;

            return (new Complex(p, q));

        }

        /// <summary>
        /// Computes the inverse cosine (arccosine) of a complex number.
        /// </summary>
        /// <param name="z">The number.</param>
        /// <returns>The value of arccos(z).</returns>
        public static Complex Acos (Complex z) {

            Asin_Internal(Math.Abs(z.Re), Math.Abs(z.Im), out double b, out double bPrime, out double q);

            double p;
            if (bPrime < 0.0) {
                p = Math.Acos(b);
            } else {
                p = Math.Atan(1.0 / bPrime);
            }
            if (z.Re < 0.0) p = Math.PI - p;

            if (z.Im > 0.0) q = -q;

            return (new Complex(p, q));

        }

        private static void Asin_Internal (double x, double y, out double b, out double bPrime, out double v) {

            // This method for the inverse complex sine (and cosine) is described in
            // Hull and Tang, "Implementing the Complex Arcsine and Arccosine Functions Using Exception Handling",
            // ACM Transactions on Mathematical Software (1997)
            // (https://www.researchgate.net/profile/Ping_Tang3/publication/220493330_Implementing_the_Complex_Arcsine_and_Arccosine_Functions_Using_Exception_Handling/links/55b244b208ae9289a085245d.pdf)

            // First, the basics: start with sin(w) = \frac{e^{iw} - e^{-iw}}{2i} = z. Here z is the input
            // and w is the output. To solve for w, define t = e^{i w} and multiply through by t to
            // get the quadratic equation t^2 - 2 i z t - 1 = 0. The solution is t = i z + \sqrt{1 - z^2}, so
            //   w = arcsin(z) = - i \ln ( i z + \sqrt{1 - z^2} )
            // Decompose z = x + i y, multiply out i z + \sqrt{1 - z^2}, use \ln s = |s| + i arg(s), and do a
            // bunch of algebra to get the components of w = arcsin(z) = u + i v
            //   u = arcsin(\beta)  v = sign(y) \ln (\alpha + \sqrt{\alpha^2 - 1})
            // where
            //   \alpha = \frac{\rho + \sigma}{2}  \beta = \frac{\rho - \sigma}{2}
            //   \rho = \sqrt{(x + 1)^2 + y^2}  \sigma = \sqrt{(x - 1)^2 + y^2}
            // This appears in DLMF section 4.23. (http://dlmf.nist.gov/4.23), along with the analogous
            //   arccos(w) = arccos(\beta) - i sign(y) \ln (\alpha + \sqrt{\alpha^2 - 1})
            // So \alpha and \beta together give us arcsin(w) and arccos(w).

            // As written, \alpha is not susceptable to cancelation errors, but \beta is. To avoid cancelation, note
            //   \beta = \frac{\rho^2 - \sigma^2}{2(\rho + \sigma)} = \frac{2 x}{\rho + \sigma} = \frac{x}{\alpha}
            // which is not subject to cancelation. Note \alpha >= 1 and |\beta| <= 1.

            // For \alpha ~ 1, the argument of the log is near unity, so we compute (\alpha - 1) instead,
            // and write the argument as 1 + (\alpha - 1) + \sqrt{(\alpha - 1)(\alpha + 1)}.
            // For \beta ~ 1, arccos does not accurately resolve small angles, so we compute the tangent of the angle
            // instead. Hull and Tang derive formulas for (\alpha - 1) and \beta' = \tan(u) that do not suffer
            // cancelation for these cases.

            // For simplicity, we assume all positive inputs and return all positive outputs. The caller should
            // assign signs appropriate to the desired cut conventions. We return v directly since its magnitude
            // is the same for both arcsin and arccos. Instead of u, we usually return \beta and sometimes \beta'.
            // If \beta' is not computed, it is set to -1; if it is computed, it should be used instead of \beta
            // to determine u. Compute u = \arcsin(\beta) or u = \arctan(\beta') for arcsin, u = \arccos(\beta)
            // or \arctan(1/\beta') for arccos.

            Debug.Assert((x >= 0.0) || Double.IsNaN(x));
            Debug.Assert((y >= 0.0) || Double.IsNaN(y));

            // For x or y large enough to overflow \alpha^2, we can simplify our formulas and avoid overflow.
            if ((x > Global.SqrtMax / 2.0) || (y > Global.SqrtMax / 2.0)) {

                b = -1.0;
                bPrime = x / y;

                double small, big;
                if (x < y) {
                    small = x;
                    big = y;
                } else {
                    small = y;
                    big = x;
                }
                v = Global.LogTwo + Math.Log(big) + 0.5 * MoreMath.LogOnePlus(MoreMath.Sqr(small / big));

            } else {

                double r = MoreMath.Hypot((x + 1.0), y);
                double s = MoreMath.Hypot((x - 1.0), y);

                double a = (r + s) / 2.0;
                b = x / a;
                //Debug.Assert(a >= 1.0);
                //Debug.Assert(Math.Abs(b) <= 1.0);

                if (b > 0.75) {
                    //double amx;
                    if (x <= 1.0) {
                        double amx = (y * y / (r + (x + 1.0)) + (s + (1.0 - x))) / 2.0;
                        bPrime = x / Math.Sqrt((a + x) * amx);
                    } else {
                        // In this case, amx ~ y^2. Since we take the square root of amx, we should
                        // pull y out from under the square root so we don't loose its contribution
                        // when y^2 underflows.
                        double t = (1.0 / (r + (x + 1.0)) + 1.0 / (s + (x - 1.0))) / 2.0;
                        bPrime = x / Math.Sqrt((a + x) * t) / y;
                    }
                    //Debug.Assert(amx >= 0.0);
                    //bPrime = x / Math.Sqrt((a + x) * amx);
                } else {
                    bPrime = -1.0;
                }

                if (a < 1.5) {
                    if (x < 1.0) {
                        // This is another case where our expression is proportional to y^2 and
                        // we take its square root, so again we pull out a factor of y from
                        // under the square root. Without this trick, our result has zero imaginary
                        // part for y smaller than about 1.0E-155, while the imaginary part
                        // should be about y.
                        double t = (1.0 / (r + (x + 1.0)) + 1.0 / (s + (1.0 - x))) / 2.0;
                        double am1 = y * y * t;
                        v = MoreMath.LogOnePlus(am1 + y * Math.Sqrt(t * (a + 1.0)));
                    } else {
                        double am1 = (y * y / (r + (x + 1.0)) + (s + (x - 1.0))) / 2.0;
                        v = MoreMath.LogOnePlus(am1 + Math.Sqrt(am1 * (a + 1.0)));
                    }
                    //Debug.Assert(am1 >= 0.0);
                } else {
                    // Because of the test above, we can be sure that a * a will not overflow.
                    v = Math.Log(a + Math.Sqrt((a - 1.0) * (a + 1.0)));
                }

            }

        }

        // pure complex binary operations

        /// <summary>
        /// Raises a complex number to an arbitrary real power.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <param name="r">The power.</param>
        /// <returns>The value of z<sup>r</sup>.</returns>
        public static Complex Pow (Complex z, double r) {
            double m = Math.Pow(Abs(z), r);
            double t = Arg(z) * r;
            return (new Complex(m * MoreMath.Cos(t), m * MoreMath.Sin(t)));
        }

        /// <summary>
        /// Raises a real number to an arbitrary complex power.
        /// </summary>
        /// <param name="x">The real base, which must be non-negative.</param>
        /// <param name="z">The complex exponent.</param>
        /// <returns>The value of x<sup>z</sup>.</returns>
        public static Complex Pow (double x, Complex z) {
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if (z == Complex.Zero) return (Complex.One);
            if (x == 0.0) return (Complex.Zero);
            double m = Math.Pow(x, z.Re);
            double t = Math.Log(x) * z.Im;
            return (new Complex(m * MoreMath.Cos(t), m * MoreMath.Sin(t)));
        }

        /// <summary>
        /// Raises a complex number to a complex power.
        /// </summary>
        /// <param name="z">The base.</param>
        /// <param name="r">The exponent.</param>
        /// <returns>The value of z<sup>r</sup>.</returns>
        /// <seealso href="http://mathworld.wolfram.com/ComplexExponentiation.html"/>
        public static Complex Pow (Complex z, Complex r) {
            if (r == Complex.Zero) return (Complex.One);
            if (z == Complex.Zero) return (Complex.Zero);
            if (z == Complex.One) return (Complex.One);
            if (r == Complex.One) return (z);
            double m = Abs(z);
            double t = Arg(z);
            double n = Math.Pow(m, r.Re) * Math.Exp(-r.Im * t);
            double u = r.Re * t + r.Im * Math.Log(m);
            return (new Complex(n * MoreMath.Cos(u), n * MoreMath.Sin(u)));
        }

        /// <summary>
        /// Raises a complex number to an integer power.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <param name="n">The power.</param>
        /// <returns>The value of z<sup>n</sup>.</returns>
        public static Complex Pow (Complex z, int n) {

            if (n < 0) return (1.0 / Pow(z, -n));

            switch (n) {
                case 0:
                    // We follow convention that 0^0 = 1
                    return (1.0);
                case 1:
                    return (z);
                case 2:
                    // 1 multiply
                    return (Sqr(z));
                case 3:
                    // 2 multiplies
                    return (Sqr(z) * z);
                case 4: {
                        // 2 multiplies
                        Complex z2 = Sqr(z);
                        return (Sqr(z2));
                    }
                case 5: {
                        // 3 multiplies
                        Complex z2 = Sqr(z);
                        return (Sqr(z2) * z);
                    }
                case 6: {
                        // 3 multiplies
                        Complex z2 = Sqr(z);
                        return (Sqr(z2) * z2);
                    }
                case 7: {
                        // 4 multiplies
                        Complex z3 = Sqr(z) * z;
                        return (Sqr(z3) * z);
                    }
                case 8: {
                        // 3 multiplies
                        Complex z2 = Sqr(z);
                        Complex z4 = Sqr(z2);
                        return (Sqr(z4));
                    }
                case 9: {
                        // 4 multiplies
                        Complex z3 = Sqr(z) * z;
                        return (Sqr(z3) * z3);
                    }
                case 10: {
                        // 4 multiplies
                        Complex z2 = Sqr(z);
                        Complex z4 = Sqr(z2);
                        return (Sqr(z4) * z2);
                    }
                case 12: {
                        // 4 multiplies
                        Complex z3 = Sqr(z) * z;
                        Complex z6 = Sqr(z3);
                        return (Sqr(z6));
                    }
                case 16: {
                        // 4 multiplies
                        Complex z2 = Sqr(z);
                        Complex z4 = Sqr(z2);
                        Complex z8 = Sqr(z4);
                        return (Sqr(z8));
                    }
                // that's all the cases do-able in 4 or fewer complex multiplies
                default:
                    return (ComplexMath.Pow(z, (double)n));
            }

        }


    }
}
