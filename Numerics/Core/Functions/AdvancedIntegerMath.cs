
using System;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Contains methods that compute advanced functions of integer arguments.
    /// </summary>
	public static class AdvancedIntegerMath {

		private static long[]  factorialTable = CreateFactorialTable(20);

		private static long[] CreateFactorialTable (int n) {
			long[] factorials = new long[n];
			factorials[0] = 1;
			for (int k=1; k<n; k++) {
				factorials[k] = k * factorials[k-1];
			}
			return(factorials);
		}

        /// <summary>
        /// Computes the factorial of an integer.
        /// </summary>
        /// <param name="n">The argument, which must be non-negative.</param>
        /// <returns>The factorial <paramref name="n"/>!.</returns>
        /// <remarks>
        /// <para>The factorial of an integer n is the product of all integers from 1 to n.</para>
        /// <para>Because n! becomes too large to be representable as a double-precision floating point number for quite
        /// moderate values of n, you may find it convenient to use the <see cref="LogFactorial"/> in order to avoid
        /// overflow when computing expression in which large factorials will cancel with other large factors.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso cref="LogFactorial"/>
		public static double Factorial (int n) {
			if (n<0) throw new ArgumentOutOfRangeException("n");
			if (n<factorialTable.Length) {
				return( (double) factorialTable[n]);
			} else {
				return( Math.Round( AdvancedMath.Gamma(n+1) ) );
			}
		}

        /// <summary>
        /// Computes the logrithm of the factorial of an integer.
        /// </summary>
        /// <param name="n">The argument, which must be non-negative.</param>
        /// <returns>The log factorial ln(<paramref name="n"/>!).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso cref="Factorial"/>
		public static double LogFactorial (int n) {
			if (n<0) throw new ArgumentOutOfRangeException("n");
			if (n<factorialTable.Length) {
				return( Math.Log( (double) factorialTable[n] ) );
			} else {
				return( AdvancedMath.LogGamma(n+1) );
			}
		}

        /// <summary>
        /// Computes binomial coefficients.
        /// </summary>
        /// <param name="n">The upper argument, which must be non-negative. (The order of the polynomial.)</param>
        /// <param name="m">The lower argument, which must be non-negative and less than or equal to <paramref name="n"/>. (The order of the term.)</param>
        /// <returns>The binomial coefficent C(<paramref name="n"/>,<paramref name="m"/>),
        /// also denoted "<paramref name="n"/> choose <paramref name="m"/>".</returns>
        /// <remarks>
        /// <para>The binomial coefficient C(n,m) is the coefficient of x<sup><paramref name="m"/></sup> in the expansion
        /// of (1+x)<sup><paramref name="n"/></sup>.</para>
        /// <para>C(n,m) can also be given a combinatoric intrepretation as the total number of distinct subsets of m items in a set of n items.</para>
        /// <para>Pascal's triangle is a classic representation of binomial coefficients.</para>
        /// <table style="text-align: center;">
        /// <tr><td colspan="3"></td><td colspan="2">C(0,0)</td></tr>
        /// <tr><td colspan="2"></td><td colspan="2">C(1,0)</td><td colspan="2">C(1,1)</td></tr>
        /// <tr><td colspan="1"></td><td colspan="2">C(2,0)</td><td colspan="2">C(2,1)</td><td colspan="2">C(2,2)</td></tr>
        /// <tr><td colspan="2">C(3,0)</td><td colspan="2">C(3,1)</td><td colspan="2">C(3,2)</td><td colspan="2">C(3,3)</td></tr>
        /// </table>
        /// <para>The relation of an element in Pascal's triangle to its two parent elements is C(n+1,m) = C(n,m-1) + C(n,m).
        /// There are many other relationships among binomial coefficients. Among the most computationally useful is
        /// B(n,m+1) = (n-m)/(m+1) B(n,m), which can be used to generate all the binomial coefficients in a row of Pascal's
        /// triangle (i.e., all the coefficients for a given order polynomial) starting from an outer values B(n,0) = 1 =B(n,n).
        /// If you need a series of binomial coefficients, using a recurrion will be more computationally efficient than
        /// calling this method for each one.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative, or <paramref name="m"/> lies outside [0,<paramref name="n"/>].</exception>
		public static long BinomialCoefficient (int n, int m) {
			if (n<0) throw new ArgumentOutOfRangeException("n");
			if (m<0) throw new ArgumentOutOfRangeException("m");
			if (m>n) throw new ArgumentOutOfRangeException("m");
			if ((n<factorialTable.Length) && (m<factorialTable.Length)) {
				return( factorialTable[n] / factorialTable[m] / factorialTable[n-m] );
			} else {
                return (BinomialCoefficientLoop(n, m));
				//return( (long) Math.Round( Math.Exp(LogFactorial(n) - LogFactorial(m) - LogFactorial(n-m)) ) );
			}
		}

        // this is an O(m) algorithm for computing (n m)

        private static long BinomialCoefficientLoop (int n, int m) {

            int km;
            if (m < n / 2) {
                km = n - m;
            } else {
                km = m;
            }

            long t = 1;
            for (int i = n; i > km; i--) {
                t = t * i / (n - i + 1);
            }
            return (t);
        }

        /*
        A suggestion from http://delab.csd.auth.gr/papers/SBI02m.pdf for computing binomial coefficient
        (n k) in O(min(k,n-k)

        FUNCTION Comb6(n,k:INTEGER): INTEGER;
        VAR t: INTEGER;
        BEGIN
            t:=1
            IF k<n-k
                THEN FOR i:=n DOWNTO n-k+1 DO
                    t:=t*i/(n-i+1)
                ELSE FOR i:=n DOWNTO k+1 DO
                    t:=t*i/(n-i+1);
            Comb6:=t
        END;
        
        */


        private static long DoubleFactorial_Multiply (int n) {
            long f = 1;
            for (int k = n; k > 1; k = k - 2) {
                f = f * k;
            }
            return (f);
        }

        private static double LogDoubleFactorial_Gamma (int n) {
            if (n % 2 == 0) {
                // m = n/2, n!! = 2^m Gamma(m+1)
                int m = n / 2;
                return (m * Math.Log(2.0) + AdvancedMath.LogGamma(m + 1.0));
            } else {
                // m = (n+1)/2, n!! = 2^m Gamma(m+1/2) / Sqrt(PI)
                int m = (n + 1) / 2;
                return (m * Math.Log(2.0) + AdvancedMath.LogGamma(m + 0.5) - 0.5 * Math.Log(Math.PI));
            }
        }

        internal static double DoubleFactorial (int n) {
            if (n < 30) {
                return ((double) DoubleFactorial_Multiply(n));
            } else {
                return (Math.Round(Math.Exp(LogDoubleFactorial_Gamma(n))));
            }
        }

        internal static double LogDoubleFactorial (int n) {
            if (n < 30) {
                return (Math.Log((double) DoubleFactorial_Multiply(n)));
            } else {
                return (LogDoubleFactorial_Gamma(n));
            }
        }

        /*
		public static long StirlingS1 (int n, int m) {
			throw new NotImplementedException();
		}

		public static long StirlingS2 (int n, int m) {
			throw new NotImplementedException();
		}

		public static double BernoulliNumberB (int n) {
			throw new NotImplementedException();
		}

		public static long EulerNumberE (int n) {
			throw new NotImplementedException();
		}
        */

		// 2-integer functions

		// Euclid's algorithm; Knuth gives a binary algorithm that may be faster
        /// <summary>
        /// Computes the greatest common factor of two integers.
        /// </summary>
        /// <param name="u">The first integer.</param>
        /// <param name="v">The second integer.</param>
        /// <returns>The greatest common factor of <paramref name="u"/> and <paramref name="v"/>.</returns>
        /// <seealso cref="LCM" />
		public static long GCF (long u, long v) {
			while (v != 0) {
				long t = u % v;
				u = v;
				v = t;
			}
			return(u);
		}
		// Knuth gives example GCD(40902,24140) = 34

        /// <summary>
        /// Computes the least common multiple of two integers.
        /// </summary>
        /// <param name="u">The first integer.</param>
        /// <param name="v">The second integer.</param>
        /// <returns>The least common multiple of <paramref name="u"/> and <paramref name="v"/>.</returns>
        /// <seealso cref="GCF"/>
		public static long LCM (long u, long v) {
			return( u / GCF(u,v) * v );
		}

		// primality testing
		// 3J 6J 9J symbols

		// Bernoulli and Euler polynomials

        // this doesn't work due to integer overflow; we need an arbitrary-precision integer
        // structure to make it work

        /*
        public static uint PowMod (uint b, uint e, uint m) {

            uint r = 1;

            while (e > 0) {
                if ((e & 1) == 1) {
                    r = checked((r * b) % m);
                }
                e = e >> 1;
                b = checked((b * b) % m);
            }

            return (r);

        }
        */

        // see http://en.wikipedia.org/wiki/Modular_exponentiation

	}
	
}
