
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
        /// <para>Because the Gamma function grows rapidly for increasing real, positive arguments, it is often necessary to
        /// work with its logrithm in order to avoid overflow.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="Gamma(Double)" />
		public static double LogGamma (double x) {
			if (x <= 0.0) throw new ArgumentOutOfRangeException("x");
            if (x > 15.0) {
                return (StirlingLogGamma(x));
            } else {
                return (LanczosLogGamma(x));
                //return (LanczosLogGamma(Lanczos.Lanczos9_C, Lanczos.Lanczos9_G, Lanczos.SqrtTwoPi, x));
                //return (LanczosLogGamma(Lanczos.Lanczos15_F, Lanczos.Lanczos15_G, Lanczos.SqrtTwoPi, x));
                //return (LanczosLogGamma(Lanczos.LanczosF, Lanczos.LanczosG, x));
            }
		}

        /// <summary>
        /// Computes the Gamma function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The Gamma function &#x393;(x).</returns>
        /// <remarks>
        /// <para>The Gamma function is a generalization of the factorial (see <see cref="AdvancedIntegerMath.Factorial"/>) to arbitrary real values. If we define &#x393;(x) =
        /// <sub>0</sub>&#x222B;<sup>&#x221E;</sup>dt t<sup>x-1</sup> e<sup>-t</sup>, then for integer values &#x393;(n+1)=n!, but the
        /// the integral can also be evaluated for non-integer x.</para>
        /// <para>Because &#x393;(x) grows beyond the largest value that can be represented by a <see cref="System.Double" /> at quite
        /// moderate values of x, you may find it useful to work with the <see cref="LogGamma" /> method, which returns ln(&#x393;(x)).</para>
        /// </remarks>
        /// <seealso cref="AdvancedIntegerMath.Factorial" />
        /// <seealso cref="LogGamma" />
		public static double Gamma (double x) {
			if (x < 0.0) {
                return(Math.PI / Gamma(-x) / (-x) / AdvancedMath.Sin(0.0, x/2.0));
                //double y = 1.0 - x;
                //double py = Math.PI * y;
                //return( py / Math.Sin(py) / Gamma(1.0+y) );
			}
			return( Math.Exp(LogGamma(x)) );
		}

        /// <summary>
        /// Computes the Psi function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value of &#x3C8;(x).</returns>
        /// <remarks>The Psi function is the derivative of the LogGamma function.</remarks>
        /// <seealso cref="LogGamma"/>
		public static double Psi (double x) {
            if (x < 0.0) {
                return (Psi(1.0 - x) - Math.PI / Math.Tan(Math.PI * x));
            }
//			if (x < 1.0) {
//				double y = 1.0 - x;
//				return( Psi(y) + Math.PI / Math.Tan(Math.PI * y) );
//			}
            return (LanczosPsi(x));
			//return( LanczosPsi(Lanczos.LanczosF, Lanczos.LanczosG, x) );
		}

		// need special handling of negative arguments


        // two-argument functions

        /// <summary>
        /// Computes the Beta function.
        /// </summary>
        /// <param name="a">The first parameter.</param>
        /// <param name="b">The second parameter.</param>
        /// <returns>The beta function B(a,b).</returns>
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
        /// in this region, use the complementary function <see cref="RightGamma"/>.</para>
        /// <para>For a=&#x3BD;/2 and x=&#x3C7;<sup>2</sup>/2, this function is the CDF of the &#x3C7;<sup>2</sup> distribution with &#x3BD; degrees of freedom.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is negative.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="RightGamma" />
        public static double LeftGamma (double a, double x) {
			if (a <= 0) throw new ArgumentOutOfRangeException("a");
			if (x < 0) throw new ArgumentOutOfRangeException("x");
			if (x<(a+1)) {
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
        /// <remarks>This function is the complement of the left incomplete Gamma function <see cref="LeftGamma"/>. </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is negative.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="LeftGamma"/>
		public static double RightGamma (double a, double x) {
			if (a <= 0) throw new ArgumentOutOfRangeException("a");
			if (x < 0) throw new ArgumentOutOfRangeException("x");
			if (x < (a+1)) {
				return( 1.0 - GammaP_Series(a, x) );
			} else {
				return( GammaQ_ContinuedFraction(a, x) );
			}
		}

        /// <summary>
        /// Computes the incomplete Gamma function.
        /// </summary>
        /// <param name="a">The shape parameter, which must be positive.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of &#x393;(a,x).</returns>
        /// <remarks><para>Like the &#x393; function itself, this function gets large very quickly. For most
        /// purposes, you will prefer to use the normalized incomplete gamma functions <see cref="LeftGamma"/> and
        /// <see cref="RightGamma"/>.</para></remarks>
        public static double Gamma (double a, double x) {
            return (RightGamma(a, x) * Gamma(a));
        }



		// three-argument functions
        
        /// <summary>
        /// Computes the incomplete Beta function.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <returns></returns>
		public static double Beta (double a, double b, double x) {
			if (a<0.0) throw new ArgumentOutOfRangeException("a");
			if (b<0.0) throw new ArgumentOutOfRangeException("b");
			if ((x<0.0) || (x>1.0)) throw new ArgumentOutOfRangeException("x");
			if (x == 0.0) return(0.0);
			double xtp = (a+1.0)/(a+b+2.0);
			if (x > xtp) {
				return(Beta(a,b) - Beta(b, a, 1.0-x));
			}
			// evaluate via continued fraction via Steed's method
			double aa = 1.0;			// a_1
			double bb = 1.0;			// b_1
			double D = 1.0;			// D_1 = b_0/b_1
			double Df = aa/bb;		// Df_1 = f_1 - f_0
			double f = 0.0 + Df;		// f_1 = f_0 + Df_1 = b_0 + Df_1
			int k = 0;
			do {
				k++;
				int m = k/2;
				if ((k % 2) == 0) {
					aa = x * m * (b - m) / ((a+k-1)*(a+k));
				} else {
					aa = -x * (a + m) * (a + b + m) / ((a+k-1)*(a+k));
				}
				D = 1.0 / (bb + aa * D);
				Df = (bb * D - 1.0) * Df;
				f += Df;
			} while ((f+Df) != f);
			return( Math.Pow(x, a) * Math.Pow(1.0-x, b) / a * f );
			
		}
        

		// helper functions

        // Sterling's asymptotic series, which does well for large x

        private static double StirlingLogGamma (double x) {

            double f = (x - 0.5) * Math.Log(x) - x + 0.5 * Math.Log(2.0 * Math.PI);

            double xx = x*x;
            double fx = 1.0/x;
            for (int i = 0; i < bernoulli.Length; i++) {
                double df = fx * bernoulli[i] / (2*i+1) / (2*i+2);
                double f_old = f;
                f += df;
                if (f == f_old) return(f);
                fx = fx / xx;
            }

            throw new NonconvergenceException();
        }

        private static readonly double[] bernoulli =
            new double[] { 1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66, -691.0 / 2730, 7.0 / 6, -3617.0 / 510 };

        /*
        private static double[] ComputeBernoulliNumbers (int n) {
            double[] b = new double[n];
            b[0] = 1.0;
            for (int k = 1; k < n; k++) {
                double s = 0;
                for (int m = 0; m < k; m++) {
                    s += AdvancedIntegerMath.BinomialCoefficient(k + 1, m) * b[m];
                }
                b[k] = -s / AdvancedIntegerMath.BinomialCoefficient(k + 1, k);
            }
            return (b);
        }
        */

        /*
        private static double LanczosLogGamma (double[] f, double g, double x) {
            return (LanczosLogGamma(f, g, 1.0, x));
        }
        */

        private static double LanczosLogGamma (double x) {

            // compute the Lanczos series
            double s = LanczosD[0];
            for (int i = 1; i < LanczosD.Length; i++) {
                s += LanczosD[i] / (x + i);
            }
            s = 2.0 / Math.Sqrt(Math.PI) * s / x;

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

        internal static readonly double LanczosR = 10.900511;

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
			double ds = Math.Exp( a * Math.Log(x) - x - LogGamma(a+1) );
			double s = ds;
            for (int i=0; i<Global.SeriesMax; i++) {
				ap += 1.0;
				ds *= (x / ap);
                double s_old = s;
				s += ds;
                if (s == s_old) return (s);
			}
            throw new NonconvergenceException();
		}

		// Compute GammaQ(a,x) for x > a+1
		private static double GammaQ_ContinuedFraction (double a, double x) {
			double aa = 1.0;			// a_1
			double bb = x - a + 1.0;	// b_1
			double D = 1.0/bb;		// D_1 = b_0/b_1
			double Df = aa/bb;		// Df_1 = f_1 - f_0
			double f = 0.0 + Df;		// f_1 = f_0 + Df_1 = b_0 + Df_1
			for (int k=1; k<Global.SeriesMax; k++) {
				double f_old = f;
				aa = -k * (k-a);
				bb += 2.0;
				D = 1.0 / (bb + aa * D);
				Df = (bb * D - 1.0) * Df;
				f += Df;
				if (f == f_old) return( Math.Exp(a * Math.Log(x) - x - LogGamma(a)) * f );
			}
			throw new NonconvergenceException();
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
        /// <returns>The complex value &#x393;(z).</returns>
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
        /// <param name="z">The complex argument.</param>
        /// <returns>The complex value ln(&#x393;(z)).</returns>
        /// <exception cref="ArgumentOutOfRangeException">The real part of <paramref name="z"/> is negative.</exception>
        public static Complex LogGamma (Complex z) {
            if (z.Re < 0.0) throw new ArgumentOutOfRangeException("z");
            if (ComplexMath.Abs(z) > 15.0) {
                return (StirlingLogGamma(z));
            } else {
                return (LanczosLogGamma(z));
                //return (LanczosLogGamma(Lanczos.Lanczos9_C, Lanczos.Lanczos9_G, Lanczos.SqrtTwoPi, z));
                //return (LanczosLogGamma(Lanczos.Lanczos15_F, Lanczos.Lanczos15_G, Lanczos.SqrtTwoPi, z));
                //return (LanczosLogGamma(Lanczos.LanczosF, Lanczos.LanczosG, z));
            }
        }


        private static Complex StirlingLogGamma (Complex z) {

            //Console.WriteLine("Stirling");

            // work in the upper complex plane
            if (z.Im < 0.0) return (StirlingLogGamma(z.Conjugate).Conjugate);

            Complex f = (z - 0.5) * ComplexMath.Log(z) - z + 0.5 * Math.Log(2.0 * Math.PI);
            //Console.WriteLine("f={0}", f);

            // reduce f.Im modulo 2*PI
            // result is cyclic in f.Im modulo 2*PI, but if f.Im starts off too big, the corrections
            // applied below will be lost because they are being added to a big number
            f = new Complex(f.Re, AdvancedMath.Reduce(f.Im, 0.0));
            //Console.WriteLine("f={0}", f);

            Complex zz = z * z;
            Complex fz = 1.0 / z;
            for (int i = 0; i < bernoulli.Length; i++) {
                Complex df = fz * bernoulli[i] / (2 * i + 1) / (2 * i + 2);
                Complex f_old = f;
                f += df;
                //Console.WriteLine("f={0}", f);
                if (f == f_old) return (f);
                fz = fz / zz;
            }

            throw new NonconvergenceException();
        }

        private static readonly double[] bernoulli =
            new double[] { 1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66, -691.0 / 2730, 7.0 / 6, -3617.0 / 510 };


        /*
        private static Complex LanczosLogGamma (double[] f, double g, Complex z) {
            return (LanczosLogGamma(f, g, 1.0, z));
        }
        */


        private static Complex LanczosLogGamma (Complex z) {

            // compute the Lanczos series
            Complex s = AdvancedMath.LanczosD[0];
            for (int i = 1; i < AdvancedMath.LanczosD.Length; i++) {
                s += AdvancedMath.LanczosD[i] / (z + i);
            }
            s = (2.0 / Math.Sqrt(Math.PI)) * (s / z);

            // compute the leading terms
            Complex zz = z + 0.5;
            Complex t = zz * ComplexMath.Log(zz + AdvancedMath.LanczosR) - z;

            return (t + ComplexMath.Log(s));

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
