
using System;

namespace Meta.Numerics.Functions {

	public static partial class AdvancedMath {

		// one-argument functions

        /// <summary>
        /// Computes the natural logrithm of the Gamma function.
        /// </summary>
        /// <param name="x">The argument, which must be positive.</param>
        /// <returns>The log Gamma function ln(&#x393;(x)).</returns>
        /// <remarks>
        /// <para>Because &#x393;(x) grows rapidly for increasing positive x, it is often necessary to
        /// work with its logarithm in order to avoid overflow. This function returns accurate
        /// values of ln(&#x393;(x)) even for values of x which would cause &#x393;(x) to overflow.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative or zero.</exception>
        /// <seealso cref="Gamma(double)" />
		public static double LogGamma (double x) {
            if (x <= 0.0) {
                throw new ArgumentOutOfRangeException("x");
            } else if (x > 16.0) {
                // for large arguments, the asymptotic series is even faster than the Lanczos approximation
                return (LogGamma_Stirling(x));
            } else {
                return (LanczosLogGamma(x));
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
        /// <para>Because &#x393;(x) grows beyond the largest value that can be represented by a <see cref="System.Double" /> at quite
        /// moderate values of x, you may find it useful to work with the <see cref="LogGamma" /> method, which returns ln(&#x393;(x)).</para>
        /// <para>To evaluate the Gamma function for a complex argument, use <see cref="AdvancedComplexMath.Gamma" />.</para>
        /// </remarks>
        /// <seealso cref="AdvancedIntegerMath.Factorial" />
        /// <seealso cref="LogGamma" />
        /// <seealso cref="AdvancedComplexMath.Gamma" />
        /// <seealso href="http://en.wikipedia.org/wiki/Gamma_function" />
        /// <seealso href="http://mathworld.wolfram.com/GammaFunction.html" />
        public static double Gamma (double x) {
			if (x <= 0.0) {
                if (x == Math.Ceiling(x)) {
                    // poles at zero and negative integers
                    return (Double.NaN);
                } else {
                    return (Math.PI / Gamma(-x) / (-x) / AdvancedMath.Sin(0.0, x / 2.0));
                }
			}
			return( Math.Exp(LogGamma(x)) );
		}

        /// <summary>
        /// Computes the digamma function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of &#x3C8;(x).</returns>
        /// <remarks>
        /// <para>The psi function, also called the digamma function, is the logrithmic derivative of the &#x393; function.</para>
        /// <img src="../images/DiGamma.png" />
        /// <para>To evaluate the Psi function for complex arguments, use <see cref="AdvancedComplexMath.Psi" />.</para>
        /// </remarks>
        /// <seealso cref="Gamma(double)"/>
        /// <seealso cref="AdvancedComplexMath.Psi"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Digamma_function" />
        /// <seealso href="http://mathworld.wolfram.com/DigammaFunction.html" />
		public static double Psi (double x) {
            if (x <= 0.0) {
                if (x == Math.Ceiling(x)) {
                    // there are poles at zero and negative integers
                    return (Double.NaN);
                } else {
                    // use the reflection formula to change to a positive x
                    return (Psi(1.0 - x) - Math.PI / Math.Tan(Math.PI * x));
                }
            } else if (x > 16.0) {
                // for large arguments, the Stirling asymptotic expansion is faster than the Lanzcos approximation
                return (Psi_Stirling(x));
            } else {
                return (LanczosPsi(x));
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
        public static double Psi (int n, double x) {

            if (n < 0) throw new ArgumentOutOfRangeException("n");

            // for n=0, use normal digamma algorithm
            if (n == 0) return (Psi(x));

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
                    // add -j times our coeffcient to the coefficient above
                    p[j + 1] += -j * p[j];
                    // same for the coefficient below; we need not add since no one else has addressed it yet (since we are moving down)
                    if (j > 0) p[j - 1] = -j * p[j];
                    // we are done with this coefficeint; make it zero for the next time
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

        // two-argument functions

        /// <summary>
        /// Computes the Beta function.
        /// </summary>
        /// <param name="a">The first parameter.</param>
        /// <param name="b">The second parameter.</param>
        /// <returns>The beta function B(a,b).</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Beta_function"/>
        public static double Beta (double a, double b) {
			return( Math.Exp( LogGamma(a) + LogGamma(b) - LogGamma(a+b) ) );
		}

        /// <summary>
        /// Computes the normalized lower (left) incomplete Gamma function.
        /// </summary>
        /// <param name="a">The shape parameter, which must be positive.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of &#x3B3;(a,x)/&#x393;(x).</returns>
        /// <remarks><para>The incomplete Gamma function is obtained by carrying out the Gamma function integration from zero to some
        /// finite value x, instead of to infinity. The function is normalized by dividing by the complete integral, so the
        /// function ranges from 0 to 1 as x ranges from 0 to infinity.</para>
        /// <para>For large values of x, this function becomes 1 within floating point precision. To determine its deviation from 1
        /// in this region, use the complementary function <see cref="RightRegularizedGamma"/>.</para>
        /// <para>For a=&#x3BD;/2 and x=&#x3C7;<sup>2</sup>/2, this function is the CDF of the &#x3C7;<sup>2</sup> distribution with &#x3BD; degrees of freedom.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is negative or zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="RightRegularizedGamma" />
        public static double LeftRegularizedGamma (double a, double x) {
			if (a <= 0) throw new ArgumentOutOfRangeException("a");
			if (x < 0) throw new ArgumentOutOfRangeException("x");
            if ((a > 128.0) && (Math.Abs(x - a) < 0.25 * a)) {
                double P, Q;
                Gamma_Temme(a, x, out P, out Q);
                return (P);
            } else if (x<(a+1.0)) {
				return( GammaP_Series(a, x) );
			} else {
				return( 1.0 - GammaQ_ContinuedFraction(a, x) );
			}
		}

        /// <summary>
        /// Computes the normalized upper (right) incomplete Gamma function.
        /// </summary>
        /// <param name="a">The shape paraemter, which must be positive.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of &#x393;(a,x)/&#x393;(x).</returns>
        /// <remarks>
        /// <para>This function is the complement of the left incomplete Gamma function <see cref="LeftRegularizedGamma"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is negative or zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="LeftRegularizedGamma"/>
		public static double RightRegularizedGamma (double a, double x) {
			if (a <= 0) throw new ArgumentOutOfRangeException("a");
			if (x < 0) throw new ArgumentOutOfRangeException("x");
            if ((a > 128.0) && (Math.Abs(x - a) < 0.25 * a)) {
                double P, Q;
                Gamma_Temme(a, x, out P, out Q);
                return (Q);
            } else if (x < (a+1.0)) {
				return( 1.0 - GammaP_Series(a, x) );
			} else {
				return( GammaQ_ContinuedFraction(a, x) );
			}
		}

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
        /// <para>Like the &#x393; function itself, this function gets large very quickly. For most
        /// purposes, you will prefer to use the regularized incomplete gamma functions <see cref="LeftRegularizedGamma"/> and
        /// <see cref="RightRegularizedGamma"/>.</para>
        /// </remarks>
        /// <seealso cref="Gamma(double)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Incomplete_Gamma_function"/>
        public static double Gamma (double a, double x) {
            return (RightRegularizedGamma(a, x) * Gamma(a));
        }



		// three-argument functions
        
        /// <summary>
        /// Computes the incomplete Beta function.
        /// </summary>
        /// <param name="a">The left shape parameter, which must be non-negative.</param>
        /// <param name="b">The right shape paraemter, which must be non-negative.</param>
        /// <param name="x">The integral endpoint, which must lie in [0,1].</param>
        /// <returns>The value of B<sub>x</sub>(a, b).</returns>
		public static double Beta (double a, double b, double x) {
            if (a < 0.0) throw new ArgumentOutOfRangeException("a");
            if (b < 0.0) throw new ArgumentOutOfRangeException("b");
            if ((x < 0.0) || (x > 1.0)) throw new ArgumentOutOfRangeException("x");
			if (x == 0.0) return(0.0);
            double xtp = (a + 1.0) / (a + b + 2.0);
			if (x > xtp) {
				return(Beta(a,b) - Beta(b, a, 1.0-x));
			}
			// evaluate via continued fraction via Steed's method
			double aa = 1.0;			// a_1
			double bb = 1.0;			// b_1
			double D = 1.0;			    // D_1 = 1 / b_1
			double Df = aa / bb;		// Df_1 = f_1 - f_0 = a_1 / b_1
			double f = 0.0 + Df;		// f_1 = f_0 + Df_1 = b_0 + Df_1
			for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
				int m = k/2;
				if ((k % 2) == 0) {
					aa = x * m * (b - m) / ((a+k-1)*(a+k));
				} else {
					aa = -x * (a + m) * (a + b + m) / ((a+k-1)*(a+k));
				}
				D = 1.0 / (bb + aa * D);
				Df = (bb * D - 1.0) * Df;
				f += Df;
                if (f == f_old) {
                    return (Math.Pow(x, a) * Math.Pow(1.0 - x, b) / a * f);
                }
			}
            throw new NonconvergenceException();
		}

        /// <summary>
        /// Computes the regularized incomplete Beta function.
        /// </summary>
        /// <param name="a">The left shape parameter, which must be non-negative.</param>
        /// <param name="b">The right shape paraemter, which must be non-negative.</param>
        /// <param name="x">The integral endpoint, which must lie in [0,1].</param>
        /// <returns>The value of I<sub>x</sub>(a, b) = B<sub>x</sub>(a, b) / B(a, b).</returns>
        public static double LeftRegularizedBeta (double a, double b, double x) {
            return (Beta(a, b, x) / Beta(a, b));
        }

		// helper functions

        // Stirling's asymptotic series, which does well for large x
        //   \log \Gamma(z) = (z - 1/2) \log z - z + \log (2 \pi) / 2 +
        //     \sum_{k=1}^{\infty} \frac{B_{2k}}{2k (2k-1) z^{2k-1}}
        // this converges in ~8 terms for z > ~10

        private static double LogGamma_Stirling (double x) {

            double f = (x - 0.5) * Math.Log(x) - x + Math.Log(Global.TwoPI) / 2.0;

            double xx = x * x;
            double xp = x;
            for (int i = 1; i < AdvancedIntegerMath.Bernoulli.Length; i++) {
                double f_old = f;
                f += AdvancedIntegerMath.Bernoulli[i] / (2 * i) / (2 * i - 1) / xp;
                if (f == f_old) return(f);
                xp *= xx;
            }

            throw new NonconvergenceException();
        }

        // by differentiating the Stirling asymptotic series for log Gamma, we obtain an analogous asymptotic series for psi
        //   \psi(x) = \log z - 1 / (2 z) - \sum_{k=1}^{\infty} \frac{B_{2k}{2k z^{2k}}

        private static double Psi_Stirling (double x) {

            double f = Math.Log(x) - 1.0 / (2.0 * x);

            double x2 = x * x;
            double xp = x2;
            for (int i = 1; i < AdvancedIntegerMath.Bernoulli.Length; i++) {
                double f_old = f;
                f -= AdvancedIntegerMath.Bernoulli[i] / (2 * i) / xp;
                if (f == f_old) return (f);
                xp *= x2;
            }
            throw new NonconvergenceException();

        }

        private static double LanczosLogGamma (double x) {

            // compute the Lanczos series
            double s = LanczosD[0];
            for (int i = 1; i < LanczosD.Length; i++) {
                s += LanczosD[i] / (x + i);
            }
            s = 2.0 / Global.SqrtPI * s / x;

            // compute the leading terms
            double xx = x + 0.5;
            double t = xx * Math.Log(xx + LanczosR) - x;

            return (t + Math.Log(s));

        }

        private static double LanczosPsi (double x) {

            // compute the Lanczos series
            double s0 = LanczosD[0];
            double s1 = 0.0;
            for (int i = 1; i < LanczosD.Length; i++) {
                double xi = x + i;
                double st = LanczosD[i] / xi;
                s0 += st;
                s1 += st / xi;
            }

            // compute the leading terms
            double xx = x + LanczosR + 0.5;
            double t = Math.Log(xx) - LanczosR / xx - 1.0 / x;

            return (t - s1/s0);

        }

        // These values are from Table 8.5 of Glendon Ralph Pugh's 2004 PhD thesis,
        // available at http://web.viu.ca/pughg/phdThesis/phdThesis.pdf
        // He claims that they "guarantee 16 digit floating point accuracy in the right-half plane"

        internal static readonly double[] LanczosD = new double[] {
             2.48574089138753565546e-5,
             1.05142378581721974210,
            -3.45687097222016235469,
             4.51227709466894823700,
            -2.98285225323576655721,
             1.05639711577126713077,
            -1.95428773191645869583e-1,
             1.70970543404441224307e-2,
            -5.71926117404305781283e-4,
             4.63399473359905636708e-6,
            -2.71994908488607703910e-9
        };

        internal const double LanczosR = 10.900511;

        /*
        private static double LanczosLogGamma (double[] f, double g, double c, double x) {
			double p = x + 0.5;
			double q = g + p;
			double s = f[0];
			for (int i=1; i<f.Length; i++) {
				s += f[i] / (x + i);
			}
			return( p * Math.Log(q) - q + Math.Log(c * s / x) );
		}

		private static double LanczosPsi (double[] f, double g, double x) {
			double p = x + 0.5;
			double q = g + p;
			double s = f[0];
			double t = 0.0;
			for (int i=1; i<f.Length; i++) {
				double xi = x + i;
				double ds = f[i] / xi;
				s += ds;
				t += ds/xi;
			}
			return( Math.Log(q) - g/q - t/s - 1.0/x );
		}
        */

		// Compute GammaP(a,x) for x < a+1
		private static double GammaP_Series (double a, double x) {
            if (x == 0.0) return (0.0);
			double ap = a;
			double ds = Math.Exp( a * Math.Log(x) - x - LogGamma(a + 1.0) );
			double s = ds;
            for (int i=0; i<Global.SeriesMax; i++) {
				ap += 1.0;
				ds *= (x / ap);
                double s_old = s;
				s += ds;
                if (s == s_old) {
                    return (s);
                }
			}
            throw new NonconvergenceException();
		}

		// Compute GammaQ(a,x) for x > a+1
		private static double GammaQ_ContinuedFraction (double a, double x) {
            if (Double.IsPositiveInfinity(x)) return (0.0);
			double aa = 1.0;			// a_1
			double bb = x - a + 1.0;	// b_1
			double D = 1.0/bb;		    // D_1 = b_0/b_1
			double Df = aa/bb;		    // Df_1 = f_1 - f_0
			double f = 0.0 + Df;		// f_1 = f_0 + Df_1 = b_0 + Df_1
            // entering this loop with bb infinite (as caused e.g. by infinite x) will cause a
            // NonconvergenceException instead of the expected convergence to zero
			for (int k=1; k<Global.SeriesMax; k++) {
				double f_old = f;
				aa = -k * (k-a);
				bb += 2.0;
				D = 1.0 / (bb + aa * D);
				Df = (bb * D - 1.0) * Df;
				f += Df;
                if (f == f_old) {
                    return (Math.Exp(a * Math.Log(x) - x - LogGamma(a)) * f);
                }
			}
			throw new NonconvergenceException();
		}

        // For large a, for x ~ a, the convergence of both the series and the continued fraction for the incomplete gamma
        // function is very slow; the problem is that the kth term goes like x/(a+k), and for large a, adding k
        // makes little difference until k ~ a.

        // In this region, NR uses ~15-point Gaussian quadrature of the peaked integrand, which should be about as good
        // as one iteration of the 15-point Gauss-Kronrod integrator used by our adaptive integrator. But when I tried
        // our adaptive integrator, it wanted to subdivide and repeat, requiring hundreds of evaluations to achieve full
        // accuracy. This leads me to believe that the NR algorithm is unlikely to achieve full accuracy.

        // So in this region we use instead a rather strange expansion due to Temme. It looks simple enough at first:
        //   P = erfc(-z) / 2 -  R     Q = erfc(z) / 2 + R
        // where the erfc term is (nearly) the Normal(a,sqrt(a)) approximation and R is a correction term.

        // The first odditity is that z is not quite (x-a)/sqrt(2a). Instead
        //   z^2 / a = eta^2 / 2 = (x-a)/a - log(x/a) = e - log(1+e) = exp(u) - 1 - u
        // where e = (x-a)/a and u = x/a. Note for x ~ a, this makes z have nearly the expected value, but with O(e)~O(u) corrections.

        // R can be expressed as a double power series in 1/a and u (or e).
        //   R = exp(-z^2) / sqrt(2 pi a) \sum_{ij} D_{ij} u^{j} / a_{i}
        // To obtain the coefficients, expand inverse powers of eta in powers of e
        //   C_0 = -1/eta = singular -1/3 + 1/12 e - 23/540 e^2 + ...
        //   C_1 = 1/eta^3 = singular - 1/540 - 1/288 e + ...
        //   C_2 = -3/eta^5 = singular + 25/6048 + ...
        //   C_k = (-1)^(k+1) (2k-1)!! / eta^(2k+1)
        // Discard the singular terms. The remaining terms give the coefficients for R. This weird prescription is the
        // second oddity.

        // Since the coefficients decrease faster in terms of u than in terms of e, we convert to a power series in u after
        // dropping the singular terms. Note this is not the same as dropping the singular terms in a direct expansion in u.

        // We record enough terms to obtain full accuracy when a > 100 and |e| < 0.25.

        private static readonly double[][] TemmeD = new double[][] {
            new double[] { - 1.0 / 3.0, 1.0 / 12.0, - 1.0 / 1080.0, - 19.0 / 12960.0, 1.0 / 181440.0, 47.0 / 1360800.0,
                1.0 / 32659200.0, - 221.0 / 261273600.0, - 281.0 / 155196518400.0,  857.0 / 40739086080.0, 1553.0 / 40351094784000.0 },
            new double[] { -1.0 / 540.0, - 1.0 / 288.0, 25.0 / 12096.0, - 223.0 / 1088640.0, - 89.0 / 1088640.0,
                757.0 / 52254720.0,  445331.0 / 155196518400.0, - 1482119.0 / 2172751257600.0, - 7921307.0 / 84737299046400.0 },
            new double[] { 25.0 / 6048.0, - 139.0 / 51840.0, 101.0 / 311040.0, 1379.0 / 7464960.0, - 384239.0 / 7390310400.0,
                - 1007803.0 / 155196518400.0, 88738171.0 / 24210656870400.0, 48997651.0 / 484213137408000.0 },
            new double[] { 101.0 / 155520.0, 571.0 / 2488320.0, - 3184811.0 / 7390310400.0, 36532751.0 / 310393036800.0,
                10084279.0 / 504388684800.0, - 82273493.0 / 5977939968000.0 },
            new double[] { - 3184811.0 / 3695155200.0, 163879.0 / 209018880.0, - 2745493.0 / 16303472640.0,
                - 232938227.0 / 2934625075200.0, 256276123.0 / 5869250150400.0 },
            new double[] { - 2745493.0 / 8151736320.0, - 5246819.0 / 75246796800.0, 119937661.0 / 451480780800.0,
                - 294828209.0 / 2708884684800.0 },
            new double[] { 119937661.0 / 225740390400.0, - 534703531.0 / 902961561600.0 },
            new double[] { 8325705316049.0 / 24176795811840000.0 , 4483131259.0 / 86684309913600.0 }
        };

        private static void Gamma_Temme (double a, double x, out double P, out double Q) {

            double u = Math.Log(x / a);
            
            // compute argument of error function, which is almost (x-a)/sqrt(a)

            double dz = 1.0;
            double z = dz;
            for (int i = 3; true; i++) {
                if (i > Global.SeriesMax) throw new NonconvergenceException();
                double z_old = z;
                dz *= u / i;
                z += dz;
                if (z == z_old) break;
            }
            z = u * Math.Sqrt(a * z / 2.0);

            // the first approximation is just the almost-Gaussian one

            if (z > 0) {
                Q = AdvancedMath.Erfc(z) / 2.0;
                P = 1.0 - Q;
            } else {
                P = AdvancedMath.Erfc(-z) / 2.0;
                Q = 1.0 - P;
            }

            // compute Temme's correction to the Gaussian approximation

            double R0 = Math.Exp(-z*z) / Math.Sqrt(Global.TwoPI * a);

            double S0 = 0.0;
            double ai = 1.0;
            for (int i=0; i < TemmeD.Length; i++) {
                double dS = 0.0;
                double uj = 1.0;
                for (int j = 0; j < TemmeD[i].Length; j++) {
                    dS += TemmeD[i][j] * uj;
                    uj *= u;
                }
                S0 += dS / ai;
                ai *= a;
            }

            double R = R0 * S0;
            Q = Q + R;
            P = P - R;

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
                return (Math.PI / Gamma(1.0 - z) / ComplexMath.Sin(Math.PI * z));
                // -z form
                //return (-Math.PI / Gamma(-z) / z / ComplexMath.Sin(Math.PI * z));
            }
            return (ComplexMath.Exp(LogGamma(z)));
        }

        /// <summary>
        /// Compute the complex log Gamma function.
        /// </summary>
        /// <param name="z">The complex argument, which must have a non-negative real part.</param>
        /// <returns>The complex value ln(&#x393;(z)).</returns>
        /// <exception cref="ArgumentOutOfRangeException">The real part of <paramref name="z"/> is negative.</exception>
        /// <seealso cref="AdvancedMath.LogGamma" />
        public static Complex LogGamma (Complex z) {
            if (z.Re < 0.0) throw new ArgumentOutOfRangeException("z");
            if (ComplexMath.Abs(z) > 15.0) {
                return (LogGamma_Stirling(z));
            } else {
                return (LanczosLogGamma(z));
            }
        }


        private static Complex LogGamma_Stirling (Complex z) {

            // work in the upper complex plane; i think this isn't actually necessary
            if (z.Im < 0.0) return (LogGamma_Stirling(z.Conjugate).Conjugate);

            Complex f = (z - 0.5) * ComplexMath.Log(z) - z + Math.Log(Global.TwoPI) / 2.0;

            // reduce f.Im modulo 2*PI
            // result is cyclic in f.Im modulo 2*PI, but if f.Im starts off too big, the corrections
            // applied below will be lost because they are being added to a big number
            f = new Complex(f.Re, AdvancedMath.Reduce(f.Im, 0.0));

            Complex zz = z * z;
            Complex zp = z;
            for (int i = 1; i < AdvancedIntegerMath.Bernoulli.Length; i++) {
                Complex f_old = f;
                f += AdvancedIntegerMath.Bernoulli[i] / (2 * i) / (2 * i - 1) / zp;
                if (f == f_old) return (f);
                zp *= zz;
            }
            throw new NonconvergenceException();
        }

        private static Complex LanczosLogGamma (Complex z) {

            // compute the Lanczos series
            Complex s = AdvancedMath.LanczosD[0];
            for (int i = 1; i < AdvancedMath.LanczosD.Length; i++) {
                s += AdvancedMath.LanczosD[i] / (z + i);
            }
            s = (2.0 / Global.SqrtPI) * (s / z);

            // compute the leading terms
            Complex zz = z + 0.5;
            Complex t = zz * ComplexMath.Log(zz + AdvancedMath.LanczosR) - z;

            return (t + ComplexMath.Log(s));

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
        public static Complex Psi (Complex z) {
            if (z.Re < 0.5) {
                // reduce z.Re in order to handle large real values!
                return (LanczosPsi(1.0 - z) - Math.PI / ComplexMath.Tan(Math.PI * z));
            } else {
                return (LanczosPsi(z));
            }
        }

        private static Complex LanczosPsi (Complex z) {

            // compute the Lanczos series
            Complex s0 = AdvancedMath.LanczosD[0];
            Complex s1 = 0.0;
            for (int i = 1; i < AdvancedMath.LanczosD.Length; i++) {
                Complex zi = z + i;
                Complex st = AdvancedMath.LanczosD[i] / zi;
                s0 += st;
                s1 += st / zi;
            }

            // compute the leading terms
            Complex zz = z + AdvancedMath.LanczosR + 0.5;
            Complex t = ComplexMath.Log(zz) - AdvancedMath.LanczosR / zz - 1.0 / z;

            return (t - s1 / s0);

        }

        /*
        private static Complex LanczosLogGamma (double[] f, double g, double c, Complex z) {
            Complex p = z + 0.5;
            Complex q = g + p;
            Complex s = f[0];
            for (int i = 1; i < f.Length; i++) {
                s += f[i] / (z + i);
            }
            return (p * ComplexMath.Log(q) - q + ComplexMath.Log(c * s / z));
        }
        */


    }
	
}
