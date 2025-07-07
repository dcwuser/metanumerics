using System;
using System.Diagnostics;

namespace Meta.Numerics.Functions {

	public static partial class AdvancedMath {

        // one-argument functions

        /// <summary>
        /// Computes the natural logarithm of the Gamma function.
        /// </summary>
        /// <param name="x">The argument, which must be positive.</param>
        /// <returns>The log Gamma function ln(&#x393;(x)).</returns>
        /// <remarks>
        /// <para>Mathematically, this this just the natural logarithm of the value of &#x393;(x), the
        /// function that is computed by <see cref="Gamma(double)"/>.</para>
        /// <para>Because &#x393;(x) grows rapidly for increasing positive x, the function <see cref="Gamma(double)"/>
        /// overflows even for moderatlely large arguments. This function provides accurate values of ln(&#x393;(x)) even
        /// for values for which &#x393;(x) would overflow. Additionally, &#x393;(x) ~ 1 near x ~ 1 and x ~ 2. This function provides fully accurate values
        /// of ln(&#x393;(x)) even for those values of x for which &#x393;(x) is indistinguishable from one or only distinguished by the last few digits.</para>
        /// <para>If you need to compute a product or quotient of several Gamma functions, computing it by exponentiating
        /// a sum or difference of evaluations of this function will allow computation of the ratio even in cases
        /// for which an individual Gamma function overflows. However, if the value you want is either a <see cref="Pochhammer(double, double)"/>
        /// function (the ratio of two Gamma functions) or a <see cref="Beta(double, double)"/> function, it will be more efficient and more accurate
        /// to use those specific functions.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="Gamma(double)" />
        public static double LogGamma (double x) {
            if (x < 2.5) {
                if (x < 0.0) {
                    throw new ArgumentOutOfRangeException(nameof(x));
                } else if (x < 0.5) {
                    // For small arguments, use reflection to move us to the series. The log term dominates.
                    return Math.Log(Math.PI / MoreMath.SinPi(x)) - GammaSeries.LogGammaOnePlus(-x);
                } else if (x < 1.5) {
                    // Use the series expansion near 1.
                    return GammaSeries.LogGammaOnePlus(x - 1.0);
                } else {
                    // The series expansion can be adapted near 2, too.
                    return GammaSeries.LogGammaTwoPlus(x - 2.0);
                }
            } else {
                if (x < 16.0) {
                    // In between, we still use Lanczos.
                    return Lanczos.LogGamma(x);
                } else if (x < Double.PositiveInfinity) {
                    // For large arguments, the asymptotic series is even faster than the Lanczos approximation.
                    return Stirling.LogGamma(x);
                } else if (x == Double.PositiveInfinity) {
                    // Precisely at infinity x * Math.Log(x) - x => NaN, so special-case it.
                    return Double.PositiveInfinity;
                } else {
                    Debug.Assert(Double.IsNaN(x));
                    return x;
                }
            }
		}

        /// <summary>
        /// Computes the Gamma function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of &#x393;(x).</returns>
        /// <remarks>
        /// <para>The Gamma function is a generalization of the factorial (see <see cref="AdvancedIntegerMath.Factorial"/>) to arbitrary real values.</para>
        /// <img src="../images/GammaIntegral.png" />
        /// <para>For positive integer arguments, this integral evaluates to &#x393;(n+1)=n!, but it can also be evaluated for non-integer z.</para>
        /// <para>Like the factorial, &#x393;(x) grows rapidly with increasing x; &#x393;(x) overflows <see cref="System.Double" /> 
        /// for all x larger than ~171. For arguments in this range, you may find it useful to work with the <see cref="LogGamma" /> method, which
        /// returns accurate values for ln(&#x393;(x)) even in the range for which &#x393;(x) overflows.</para>
        /// <para>To evaluate the Gamma function for a complex argument, use <see cref="AdvancedComplexMath.Gamma" />.</para>
        /// <h2>Domain, Range, and Accuracy</h2>
        /// <para>The function is defined for all x. It has poles at all negative integers and at zero; the method returns <see cref="Double.NaN"/> for these arguments. For positive
        /// arguments, the value of the function increases rapidly with increasing argument. For values of x greater than about 170, the value of the function exceeds
        /// <see cref="Double.MaxValue"/>; for these arguments the method returns <see cref="Double.PositiveInfinity"/>. The method is accurate to full precision over its entire
        /// domain.</para>
        /// </remarks>
        /// <seealso cref="AdvancedIntegerMath.Factorial" />
        /// <seealso cref="LogGamma" />
        /// <seealso cref="AdvancedComplexMath.Gamma" />
        /// <seealso href="http://en.wikipedia.org/wiki/Gamma_function" />
        /// <seealso href="http://mathworld.wolfram.com/GammaFunction.html" />
        /// <seealso href="http://dlmf.nist.gov/5">DLMF on the Gamma Function</seealso>
        public static double Gamma (double x) {
            // Break up decision tree into two levels to try to have a more uniform number of branch tests.
            if (x < 2.5) {
                if (x < -0.5) {
                    // Use \Gamma(x) \Gamma(1-x) = \frac{\pi}{\sin(\pi x)} and -x \Gamma(-x) = \Gamma(1-x) to move negative arguments to positive arguments.
                    return -Math.PI / (MoreMath.SinPi(x) * x * Gamma(-x));
                } else if (x < 0.5) {
                    // For |x| < 1/2, specialize the reflection formula to the series around x = 1 to avoid loss of accuracy in computing 1 - x.
                    return Math.PI / MoreMath.SinPi(x) / GammaSeries.GammaOnePlus(-x);
                } else if (x < 1.5) {
                    // Unlike the Lanczos approximation, the series guarantees \Gamma(1) = \Gamma(2) = 1 exactly
                    return GammaSeries.GammaTwoPlus(x - 1.0) / x;
                    // Don't call GammaOnePlus(x - 1) here because it computes (x - 1) + 1, which just reverses the subtraction with lost accuracy
                } else {
                    return GammaSeries.GammaTwoPlus(x - 2.0);
                }
            } else {
                if (x < 16.0) {
                    return Lanczos.Gamma(x);
                } else if (x < 172.0) {
                    return Stirling.Gamma(x);
                } else if (x <= Double.PositiveInfinity) {
                    // For x >~ 172, Gamma(x) overflows.
                    return Double.PositiveInfinity;
                } else {
                    Debug.Assert(Double.IsNaN(x));
                    return x;
                }
            }
		}

        /// <summary>
        /// Computes the digamma function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of &#x3C8;(x).</returns>
        /// <remarks>
        /// <para>The psi function, also called the digamma function, is the logrithmic derivative of the &#x393; function.</para>
        /// <img src="../images/DiGamma.png" />
        /// <para>Because it is defined as a <i>logarithmic</i> derivative, the digamma function does not overflow <see cref="System.Double"/>
        /// even for arguments for which <see cref="Gamma(double)"/> does.</para>
        /// <para>To evaluate the psi function for complex arguments, use <see cref="AdvancedComplexMath.Psi" />.</para>
        /// </remarks>
        /// <seealso cref="Gamma(double)"/>
        /// <seealso cref="AdvancedComplexMath.Psi"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Digamma_function" />
        /// <seealso href="http://mathworld.wolfram.com/DigammaFunction.html" />
		public static double Psi (double x) {
            if (x < 2.5) {
                if (x < -0.5) {
                    return Psi(1.0 - x) - Math.PI / MoreMath.TanPi(x);
                } else if (x < 0.5) {
                    // Specializing to the series here avoids loosing accuracy by computing 1-x.
                    return GammaSeries.PsiOnePlus(-x) - Math.PI / MoreMath.TanPi(x);
                } else if (x < 1.5) {
                    return GammaSeries.PsiTwoPlus(x - 1.0) - 1.0 / x;
                } else {
                    return GammaSeries.PsiTwoPlus(x - 2.0);
                }
            } else {
                if (x < 16.0) {
                    return Lanczos.Psi(x);
                } else if (x <= Double.PositiveInfinity) {
                    // For large arguments, the Stirling asymptotic expansion is faster than the Lanzcos approximation
                    return Stirling.Psi(x);
                } else {
                    Debug.Assert(Double.IsNaN(x));
                    return x;
                }
            }
		}

        /// <summary>
        /// Computes the polygamma function.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value of &#968;<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The polygamma function gives higher logarithmic derivatives of the Gamma function.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Polygamma_function"/>
        /// <seealso href="http://mathworld.wolfram.com/PolygammaFunction.html"/>
        /// <seealso href="http://dlmf.nist.gov/5.15">DLMF on the Polygamma Function</seealso>
        public static double Psi (int n, double x) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(x));

            // for n=0, use normal digamma algorithm
            if (n == 0) return Psi(x);

            // for small x, use the reflection formula
            if (x <= 0.0) {
                if (x == Math.Ceiling(x)) {
                    return (Double.NaN);
                } else {
                    // use the reflection formula
                    // this requires that we compute the nth derivative of cot(pi x)
                    // i was able to do this, but my algorithm is O(n^2), so for large n and not-so-negative x this is probably sub-optimal
                    double y = Math.PI * x;
                    double d = -EvaluateCotDerivative(ComputeCotDerivative(n), y) * MoreMath.Pow(Math.PI, n + 1);
                    if (n % 2 == 0) {
                        return (d + Psi(n, 1.0 - x));
                    } else {
                        return (d - Psi(n, 1.0 - x));
                    }
                }
            }

            // compute the minimum x required for our asymptotic series to converge
            // my original approximation was to say that, for 1 << k << n, the factorials approach e^(2k),
            // so to achieve convergence to 10^(-16) at the 8'th term we should need x ~ 10 e ~ 27
            // unfortunately this estimate is both crude and an underestimate; for now i use an emperical relationship instead
            double xm = 16.0 + 2.0 * n;

            // by repeatedly using \psi(n,z) = \psi(n,z+1) - (-1)^n n! / z^(n+1),
            // increase x until it is large enough to use the asymptotic series
            // keep track of the accumulated shifts in the variable s
            double s = 0.0;
            while (x < xm) {
                s += 1.0 / MoreMath.Pow(x, n + 1);
                x += 1.0;
            }

            // now that x is big enough, use the asymptotic series
            // \psi(n,z) = - (-1)^n (n-1)! / z^n [ 1 + n / 2 z + \sum_{k=1}^{\infty} (2k+n-1)! / (2k)! / (n-1)! * B_{2k} / z^{2k} ]
            double t = 1.0 + n / (2.0 * x);
            double x2 = x * x;
            double t1 = n * (n + 1) / 2.0 / x2;
            for (int i = 1; i < AdvancedIntegerMath.Bernoulli.Length; i++) {
                double t_old = t;
                t += AdvancedIntegerMath.Bernoulli[i] * t1;
                if (t == t_old) {
                    double g = AdvancedIntegerMath.Factorial(n - 1) * (t / MoreMath.Pow(x, n) + n * s);
                    if (n % 2 == 0) g = -g;
                    return (g);
                }
                int i2 = 2 * i;
                t1 *= 1.0 * (n + i2) * (n + i2 + 1) / (i2 + 2) / (i2 + 1) / x2;
            }
            throw new NonconvergenceException();

            // for small x, the s-part strongly dominates the t-part, so it isn't actually necessary for us to determine t
            // very accurately; in the future, we should modify this code to allow the t-series to converge when its overall
            // contribution no longer matters, rather than requiring t to converge to full precision

        }

        // The nth derivative of cot(x) can be obtained by noting that cot(x) = cos(x) / sin(x)
        // Define r_n = (c/s)^n. Then by explicit differentiation D_x r_n = - n (r_{n-1} + r_{n+1})
        // Since the result is expressed in terms of r's, differentiation can be repeated using the same formula
        // The next two methods simply implement this machinery. Is there a way we can turn this into a recursion formula so it is O(n) instead of O(n^2)?

        private static double[] ComputeCotDerivative (int n) {

            // make an array to hold coefficients of 1, r, r^2, ..., r_{n+1}
            double[] p = new double[n + 2];

            // start with one power of r
            p[1] = 1.0;

            // differentiate n times
            for (int i = 1; i <= n; i++) {

                // only even or odd powers of r appear at any given order; this fact allows us to use just one array:
                // the entries of one parity are the source (and are set to zero after being used), the of the other parity are the target
                for (int j = i; j >= 0; j -= 2) {
                    // add -j times our coefficient to the coefficient above
                    p[j + 1] += -j * p[j];
                    // same for the coefficient below; we need not add since no one else has addressed it yet (since we are moving down)
                    if (j > 0) p[j - 1] = -j * p[j];
                    // we are done with this coefficient; make it zero for the next time
                    p[j] = 0.0;
                }

            }

            // return the resulting array coefficients
            return (p);
        }

        private static double EvaluateCotDerivative (double[] p, double x) {

            // compute cot(x), which is what our polynomial is in powers of
            double r = 1.0 / Math.Tan(x);

            // compute its square, which we will multiply by to move ahead by powers of two as we evaluate
            double r2 = r * r;

            // compute the first term and set the index of the next term
            int i; double rp;
            if (p.Length % 2 == 0) {
                i = 1;
                rp = r;
            } else {
                i = 0;
                rp = 1.0;
            }

            double f = 0.0;
            while (i < p.Length) {
                // add the current term
                f += p[i] * rp;
                // prepare for next term
                i += 2;
                rp *= r2;
            }

            return (f);

        }

        // This function computes x^n / n! or x^{\nu} / \Gamma(\nu + 1), which can easily become
        // Infinity/Infinity=NaN for large n if computed naively.

        // Poisson probability is \mu^k e^{-\mu} / k!. This can easily result in Infinity/Infinity = NaN
        // if computed naively even for values of \mu and k that produce a normal-sized result.

        internal static double PoissonProbability (double mu, int k) {
            Debug.Assert(mu > 0.0);
            Debug.Assert(k >= 0);
            double P;
            if (k < 170 && mu < 707.0 && k * Math.Log(mu) < 707.0) {
                // If neither k! or e^{\mu} nor \mu^k will overflow, it is fastest and most
                // accurate to directly apply the Poisson formula.
                P = MoreMath.Pow(mu, k) * Math.Exp(-mu) / AdvancedIntegerMath.Factorial(k);
            } else if (k < 20) {
                // If our overflow conditions are not met even though k < 20
                // we must have \mu > 707, and we can be sure P < ~10^{-270}.
                // Since k << \mu, the cancellation in log space is not delicate,
                // so we accept it.
                Debug.Assert(mu >= 707.0);
                P = Math.Exp(k * Math.Log(mu) - AdvancedIntegerMath.LogFactorial(k) - mu);
                Debug.Assert(P < 2.0E-270);
            } else {
                // For large k not too far from \mu, the balance between factors can become
                // quite delicate, producing a normal-sized result even for factors
                // that would massively overflow or underflow. In this region,
                // we apply a specialization of the Stirling approximation which
                // deals with most of this cancellation analytically.
                Debug.Assert(k >= 20);
                P = Stirling.PoissonProbability(mu, k);
            }
            Debug.Assert(0.0 <= P);
            Debug.Assert(P <= 1.0);
            return P;
        }

        // Smame logic but for non-integer k, so x^a e^{-x} / \Gamma(a + 1).

        // You might think that for small a, we should adapt the Gamma series logic to produce this directly,
        // but we don't really get anything from that: we eliminate one Pow and one Exp but add two Logs, and
        // since it converts multiplication to addtion/subtraction with cancellation, accuracy suffers.

        internal static double PoissonProbability(double x, double a) {
            Debug.Assert(x >= 0.0);
            Debug.Assert(a >= 0.0);
            if (x == 0.0) {
                return a == 0.0 ? 1.0 : 0.0;
            } else if (Double.IsInfinity(x)) {
                // We need to handle infinity explicitly b/c otherwise it produces infinity * 0.0 or infinity - infinity
                return 0.0;
            } else if (a < 170.0 && x < 707.0 && a * Math.Log(x) < 707.0) {
                return Math.Pow(x, a) * Math.Exp(-x) / AdvancedMath.Gamma(a + 1.0);
            } else if (a < 20.0) {
                return Math.Exp(a * Math.Log(x) - AdvancedMath.LogGamma(a + 1.0) - x);
            } else {
                return Stirling.PoissonProbability(x, a);
            }
        }

        internal static double ExpTimesPower (double x, double a) {

            double e = (x - a) / a;
            double d = e - MoreMath.LogOnePlus(e);
            double f = Math.Pow(a / Math.E, a) * Math.Exp(-a * d);

            return Math.Exp(-x) * Math.Pow(x, a);
        }



        internal static double PowerOverFactorial (double x, int n) {
            if (n < 16) {
                // There is a small range (x, n) for which the power overflows,
                // but the factorial would bring the true result back into the representable
                // range. Since 15! ~ 1.3E12, the result overflows for values of x^n about
                // 12 orders of magnitude lower than it should. Would be good to fix this.
                return MoreMath.Pow(x, n) / AdvancedIntegerMath.Factorial(n);
            } else {
                return Stirling.PowerOverFactorial(x, n);
            }
        }

        internal static double PowerOverFactorial (double x, double nu) {
            if (nu < 16.0) {
                //return Lanczos.PowerOverFactorial(x, nu);
                return Math.Pow(x, nu) / AdvancedMath.Gamma(nu + 1.0);
            } else {
                return Stirling.PowerOverFactorial(x, nu);
            }
        }

        // x^n / (2n + 1)!!

        internal static double PowerOverDoubleFactorial (double x, int n) {
            if (n < 16) {
                return MoreMath.Pow(x, n) / AdvancedIntegerMath.DoubleFactorial(2 * n + 1);
            } else {
                // Would be good to create a dedicated method for this to avoid \sqrt{\pi} cancelation
                return(Math.Sqrt(Math.PI / 2.0 / x) * PowerOverFactorial(0.5 * x, n + 0.5));
            }
        }


        /*
        public static double LogPochhammer (double x, double y) {

            double z = x + y;
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if (z < 0.0) throw new ArgumentOutOfRangeException(nameof(y));

            if ((x < 2.5) && (z < 2.5)) {
                // Double Series
            } else if ((x > 16.0) && (z > 16.0)) {
                // Double Asymptotic
            } else if (x < 0.5) {
                // Avoid Lanczos with small x
            } else {
                // Lanczos
            }

        }
        */
        /*
        /// <summary>
        /// Computes the Pochammer symbol (x)<sub>y</sub>.
        /// </summary>
        /// <param name="x">The first argument.</param>
        /// <param name="y">The second argument.</param>
        /// <returns>The value of (x)<sub>y</sub>.</returns>
        /// <remarks>
        /// <para>The Pochhammer symbol is defined as a ratio of Gamma functions.</para>
        /// <para>For positive integer y, this is equal to a rising factorial product.</para>
        /// <para>Note that while the Pochhammer symbol notation is very common, combinatorialists sometimes
        /// use the same notation for a falling factorial.</para>
        /// <para>If you need to compute a ratio of Gamma functions, in most cases you should compute it by calling this
        /// function instead of taking the quotient of two calls to the <see cref="Gamma(double)"/> function,
        /// because there is a large fraction of the parameter space where Gamma(x+y) and Gamma(x) overflow,
        /// but their quotient does not, and is correctly computed by this method.</para>
        /// <para>This function accepts both positive and negative values for both x and y.</para>
        /// </remarks>
        /// <seealso href="http://mathworld.wolfram.com/PochhammerSymbol.html"/>
        public static double Pochhammer1 (double x, double y) {

            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if (y < 0.0) throw new ArgumentOutOfRangeException(nameof(y));

            // Handle simplest cases quickly
            if (y == 0.0) {
                return 1.0;
            } else if (y == 1.0) {
                return x;
            }

            // Should we also handle more small integer y explicitly?
            if (0.0 < y && y < 16.0) {
                double yRound = Math.Round(y);
                if (y == yRound) {
                    double p = x;
                    for (int i = ((int)yRound) - 1; i > 0; i--) {
                        x += 1.0;
                        p *= x;
                    }
                    return p;
                }
            }

            // The key case of x and y positive is pretty simple. Both the Lanczos and Stirling forms of \Gamma
            // admit forms that preserve accuracy in \Gamma(x+y) / \Gamma(x) - 1 as y -> 0. It turns
            // out, however, to be ridiculously complicated to handle all the corner cases that appear in
            // other regions of the x-y plane.

            double z = x + y;

            if (y < 0.5) { 

            }

            if (x < 0.25) {

                bool xNonPositiveInteger = (Math.Round(x) == x);
                bool zNonPositiveInteger = (z <= 0.0) && (Math.Round(z) == z);

                if (xNonPositiveInteger) {
                    if (zNonPositiveInteger) {
                        long m;
                        double p;
                        if (z >= x) {
                            m = (long) (z - x);
                            p = Pochhammer(-z + 1, m);
                        } else {
                            m = (long) (x - z);
                            p = 1.0 / Pochhammer(-x + 1, m);
                        }
                        if (m % 2 != 0) p = -p;
                        return (p);
                    } else {
                        return (0.0);
                    }
                }

                if (y < 0.0) {
                    // x negative, y negative; we can transform to make x and y positive
                    return (1.0 / Pochhammer(1.0 - x, -y) * MoreMath.SinPi(x) / MoreMath.SinPi(z));
                } else if (z <= 0.25) {
                    // x negative, y positive, but not positive enough, so z still negative
                    // we can transform to make x and y positive
                    return (Pochhammer(1.0 - z, y) * MoreMath.SinPi(x) / MoreMath.SinPi(z));
                } else {
                    // x negative, y positive enough to make z positive
                    return (-x * Gamma(-x) * Gamma(z) * MoreMath.SinPi(x) / Math.PI);
                }

            } else {

                if (z >= 0.25) {
                    // This is the standard case: x positive, y may be negative, but not negative enough to make z negative
                    //double P;
                    if ((x > 16.0) && (z > 16.0)) {
                        return Math.Exp(Stirling.LogPochhammer(x, y));
                        //P = Stirling.ReducedLogPochhammer(x, y);
                    //} else if ((x < 2.5) && (z < 2.5)) {
                    //    return AdvancedMath.Gamma(z) / AdvancedMath.Gamma(x);
                    } else {
                        //P = Lanczos.ReducedLogPochhammer(x, y);
                        return Lanczos.Pochammer(x, y);
                        //return Math.Exp(Lanczos.LogPochammer(x, y));
                    }
                    //return (Math.Exp(y * P));
                } else {
                    // if y is very negative, z will also be negative, so we can transform one gamma
                    if (Math.Round(z) == z) {
                        return (Double.PositiveInfinity);
                    } else {
                        return (Math.PI / (Gamma(x) * Gamma(1.0 - z) * MoreMath.SinPi(z)));
                    }
                }

            }

        }
        */

        // two-argument functions



        // "Efficient and accurate algorithms for the computation and inversion of the incomplete gamma function ratios"
        // https://arxiv.org/abs/1306.1754
        // Division of 1st quadrant, \alpha(x), Q series, another aysmptotic fomrulation

        // "Computation of the incomplete gamma function for negative values of the argument" 2016
        // https://arxiv.org/abs/1608.04152



        

        /// <summary>
        /// Computes the upper incomplete Gamma function.
        /// </summary>
        /// <param name="a">The shape parameter, which must be positive.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of &#x393;(a,x).</returns>
        /// <remarks>
        /// <para>The incomplete Gamma function is defined by the same integrand as the Gamma function (<see cref="Gamma(double)"/>),
        /// but the integral is not taken over the full positive real axis.</para>
        /// <img src="../images/UpperIncompleteGammaIntegral.png" />
        /// <para>Like the &#x393; function itself, this function gets large very quickly. For many
        /// purposes, you will prefer to use the regularized incomplete gamma functions <see cref="LeftRegularizedGamma"/> and
        /// <see cref="RightRegularizedGamma"/>, which are accurately computed even in regions for which this function overflows.</para>
        /// <para>To make clearer which incomplete gamma function is meant, this function is be retired. Call <see cref="UpperIncompleteGamma(double, double)"/> instead.</para>
        /// </remarks>
        /// <seealso cref="Gamma(double)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Incomplete_Gamma_function"/>
        public static double Gamma (double a, double x) {
            return AdvancedMath.UpperIncompleteGamma(a, x);
            /*
            if (a < 0.0) throw new ArgumentOutOfRangeException(nameof(a));
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            if ((a > 128.0) && (Math.Abs(x - a) < 0.25 * a)) {
                Gamma_Temme(a, x, out _, out double Q);
                return Q * AdvancedMath.Gamma(a);
            } else {
                if (a >= IncompleteGammaAlpha(x)) {
                    // P is smaller, so compute P
                    double Q = 1.0 - GammaP_Series(a, x);
                    return Q * AdvancedMath.Gamma(a);
                } else {
                    // Q is smaller, so compute Q
                    if (x > 1.5) {
                        Debug.Assert(x >= a);
                        return GammaQ_ContinuedFraction(a, x) * AdvancedMath.ExpTimesPower(x, a);
                    } else {
                        double Q = GammaQ_Series(a, x);
                        return Q * AdvancedMath.Gamma(a);
                    }
                }
            }
            */
            //return RightRegularizedGamma(a, x) * Gamma(a);
        }



        // three-argument functions

        /// <summary>
        /// Computes the incomplete Beta function.
        /// </summary>
        /// <param name="a">The left shape parameter, which must be non-negative.</param>
        /// <param name="b">The right shape parameter, which must be non-negative.</param>
        /// <param name="x">The integral endpoint, which must lie in [0,1].</param>
        /// <returns>The value of B<sub>x</sub>(a, b).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> or <paramref name="b"/> is negative, or
        /// <paramref name="x"/> lies outside [0, 1].</exception>
        /// <remarks>
        /// <para>This function is defined by the same integral that defines the Beta function (<see cref="Beta(double, double)"/>, but
        /// with the integral taken from 0 to x instead of from 0 to 1.</para>
        /// <para>Note that in most mathematical literature, <paramref name="x"/> is regarded as the first argument. We, however,
        /// follow the programming convention that optional arguments follow required ones, so <paramref name="x"/> is the last
        /// argument of the method.</para>
        /// <para>If you require the regularized function, use <see cref="LeftRegularizedBeta(double, double, double)"/>
        /// to obtain it directly.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function"/>
        /// <seealso href="http://mathworld.wolfram.com/IncompleteBetaFunction.html"/>
        public static double Beta (double a, double b, double x) {
            if (a < 0.0) throw new ArgumentOutOfRangeException(nameof(a));
            if (b < 0.0) throw new ArgumentOutOfRangeException(nameof(b));
            if ((x < 0.0) || (x > 1.0)) throw new ArgumentOutOfRangeException(nameof(x));
			if (x == 0.0) return(0.0);
            double xtp = (a + 1.0) / (a + b + 2.0);
            if (x > xtp) {
                return (Beta(a, b) - Beta(b, a, 1.0 - x));
            } else {
                return (Math.Pow(x, a) * Math.Pow(1.0 - x, b) * RegularizedBeta_ContinuedFraction(a, b, x));
            }
		}

        internal static double RegularizedBeta_ContinuedFraction (double a, double b, double x) {

            // Use the continued fraction (DLMF 8.17.22,  http://dlmf.nist.gov/8.17#v)
            //   B(a, b, x) = x^a (1-x)^b BCF(a, b, x)
            //   BCF(a, b, x) = \frac{1}{a} \left[ \frac{1}{1+} \frac{d_1}{1+} \frac{d_2}{1+} \cdots \right]
            // where
            //   d_{2m} = \frac{x m (b-m)}{(a + 2m - 1)(a + k)}
            //   d_{2m+1} = -\frac{x (a + m) (a + b + m)}{(a + m)(a + 2m + 1)}
            // which is good for x < (a + 1) / (a + b + 2). For larger x, compute compliment and subtract.

            // Evaluate via Steed's method.
            // We can simplify slightly due to this being a Stieltjes CF (denominator always 1).
            // Do k = 0 and k = 1 explicitly, since simplifications are possible for those terms.

            double ab = a + b;
            double p = -x * ab / (a + 1.0);
            double D = 1.0 / (1.0 + p);
            double Df = -p * D;
            double f = 1.0 + Df;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double f_old = f;
                int m = k / 2;
                p = x / ((a + (k - 1)) * (a + k));
                if ((k % 2) == 0) {
                    p *= m * (b - m);
                } else {
                    p *= -(a + m) * (ab + m);
                }
                D = 1.0 / (1.0 + p * D);
                Df = (D - 1.0) * Df;
                f += Df;
                 if (f == f_old) {
                    return (f / a);
                }
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes the regularized incomplete Beta function.
        /// </summary>
        /// <param name="a">The left shape parameter, which must be non-negative.</param>
        /// <param name="b">The right shape parameter, which must be non-negative.</param>
        /// <param name="x">The integral endpoint, which must lie in [0,1].</param>
        /// <returns>The value of I<sub>x</sub>(a, b) = B<sub>x</sub>(a, b) / B(a, b).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> or <paramref name="b"/> is negative or zero,
        /// or <paramref name="x"/> lies outside [0, 1].</exception>
        public static double LeftRegularizedBeta (double a, double b, double x) {
            if (a <= 0.0) throw new ArgumentOutOfRangeException(nameof(a));
            if (b <= 0.0) throw new ArgumentOutOfRangeException(nameof(b));
            if ((x < 0.0) || (x > 1.0)) throw new ArgumentOutOfRangeException(nameof(x));
            double xtp = (a + 1.0) / (a + b + 2.0);
            if (x > xtp) {
                return (1.0 - LeftRegularizedBeta(b, a, 1.0 - x));
            } else {
                return (RegularizedBeta_ContinuedFraction(a, b, x) * PowOverBeta(a, b, x));
            }
        }

        private static double PowOverBeta (double a, double b, double x) {
            if ((a > 16.0) && (b > 16.0)) {
                return (Stirling.PowOverBeta(a, b, x));
            } else {
                return (Math.Pow(x, a) * Math.Pow(1.0 - x, b) / Beta(a, b));
            }
        }

	}

    /// <summary>
    /// Contains methods that compute advanced functions of complex arguments.
    /// </summary>
    public static partial class AdvancedComplexMath {

        /// <summary>
        /// Computes the complex Gamma function.
        /// </summary>
        /// <param name="z">The complex argument.</param>
        /// <returns>The complex value of &#x393;(z).</returns>
        /// <remarks>
        /// <para>The image below shows the complex &#x393; function near the origin using domain coloring.</para>
        /// <img src="../images/ComplexGammaPlot.png" />
        /// </remarks>
        /// <seealso cref="AdvancedMath.Gamma(double)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Gamma_function" />
        /// <seealso href="http://mathworld.wolfram.com/GammaFunction.html" />
        public static Complex Gamma (Complex z) {
            if (z.Re < 0.5) {
                // 1-z form
                return (Math.PI / Gamma(1.0 - z) / ComplexMath.SinPi(z));
            } else if (ComplexMath.Abs(z) < 16.0) {
                return (Lanczos.Gamma(z));
            } else {
                // Add flag to do z-reduction
                return (ComplexMath.Exp(LogGamma_Stirling(z)));
            }
        }

        /// <summary>
        /// Compute the complex log Gamma function.
        /// </summary>
        /// <param name="z">The complex argument.</param>
        /// <returns>The principal complex value y for which exp(y) = &#x393;(z).</returns>
        /// <seealso cref="AdvancedMath.LogGamma" />
        /// <seealso href="http://mathworld.wolfram.com/LogGammaFunction.html"/>
        public static Complex LogGamma (Complex z) {
            if (z.Im == 0.0 && z.Re < 0.0) {
                // Handle the pure negative case explicitly.
                double re = Math.Log(Math.PI / Math.Abs(MoreMath.SinPi(z.Re))) - AdvancedMath.LogGamma(1.0 - z.Re);
                double im = Math.PI * Math.Floor(z.Re);
                return new Complex(re, im);
            } else if (z.Re > 16.0 || Math.Abs(z.Im) > 16.0) {
                // According to https://dlmf.nist.gov/5.11, the Stirling asymptoic series is valid everywhere
                // except on the negative real axis. So at first I tried to use it for |z.Im| > 0, |z| > 16. But in practice,
                // it exhibits false convergence close to the negative real axis, e.g. z = -16 + i. So I have
                // moved to requiring |z| large and reasonably far from the negative real axis.
                return (Stirling.LogGamma(z));
            } else if (z.Re >= 0.125) {
                return (Lanczos.LogGamma(z));
            } else {
                // For the remaining z < 0, we need to use the reflection formula.
                // For large z.Im, SinPi(z) \propto e^{\pi |z.Im|} overflows even though its log does not.
                // It's possible to do some algebra to get around that problem, but it's not necessary
                // because for z.Im that big we would have used the Stirling series.
                Complex f = ComplexMath.Log(Math.PI / ComplexMath.SinPi(z));
                Complex g = Lanczos.LogGamma(1.0 - z);
                // The reflection formula doesn't stay on the principal branch, so we need to add a multiple of 2 \pi i
                // to fix it up. See Hare, "Computing the Principal Branch of Log Gamma" for how to do this.
                // https://pdfs.semanticscholar.org/1c9d/8865836a312836500126cb47c3cbbed3043e.pdf
                Complex h = new Complex(0.0, 2.0 * Math.PI * Math.Floor(0.5 * (z.Re + 0.5)));
                if (z.Im < 0.0) h = -h;
                return (f - g + h);
            }
        }

        private static Complex LogGamma_Stirling (Complex z) {

            // work in the upper complex plane; i think this isn't actually necessary
            //if (z.Im < 0.0) return (LogGamma_Stirling(z.Conjugate).Conjugate);

            Complex f = (z - 0.5) * ComplexMath.Log(z) - z + 0.5 * Math.Log(Global.TwoPI);

            // reduce f.Im modulo 2*PI
            // result is cyclic in f.Im modulo 2*PI, but if f.Im starts off too big, the corrections
            // applied below will be lost because they are being added to a big number
            //f = new Complex(f.Re, AdvancedMath.Reduce(f.Im, 0.0));
            // we should still do an optional reduction by 2 pi so that we get phase of Gamma right even
            // for large z.

            Complex zz = ComplexMath.Sqr(z);
            Complex zp = z;
            for (int i = 1; i < AdvancedIntegerMath.Bernoulli.Length; i++) {
                Complex f_old = f;
                f += AdvancedIntegerMath.Bernoulli[i] / ((2 * i) * (2 * i - 1)) / zp;
                if (f == f_old) return (f);
                zp *= zz;
            }
            throw new NonconvergenceException();
        }

        /// <summary>
        /// Computes the complex digamma (&#x3C8;) function.
        /// </summary>
        /// <param name="z">The complex argument.</param>
        /// <returns>The value of &#x3C8;(z).</returns>
        /// <remarks>
        /// <para>The image below shows the complex &#x3C8; function near the origin using domain coloring.</para>
        /// <img src="../images/ComplexPsiPlot.png" />
        /// </remarks>
        /// <seealso cref="AdvancedMath.Psi(double)" />
        /// <seealso href="https://en.wikipedia.org/wiki/Digamma_function"/>
        public static Complex Psi (Complex z) {
            if (z.Re < 0.5) {
                // reduce z.Re in order to handle large real values!
                return (Psi(1.0 - z) - Math.PI / ComplexMath.Tan(Math.PI * z));
            } else {
                // add Stirling for large z
                return (Lanczos.Psi(z));
            }
        }

    }

}
