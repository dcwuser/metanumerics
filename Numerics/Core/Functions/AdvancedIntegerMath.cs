using System;
using System.Collections;
using System.Collections.Generic;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Contains methods that compute advanced functions of integer arguments.
    /// </summary>
	public static class AdvancedIntegerMath {

		private static readonly long[]  factorialTable = CreateFactorialTable(20);

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
        /// <returns>The factorial n!.</returns>
        /// <remarks>
        /// <para>The factorial of an integer n is the product of all integers from 1 to n. For example, 4! = 4 * 3 * 2 * 1 = 24.</para>
        /// <para>n! also has a combinatorial intrepretation as the number of permutations of n objects. For example, a set of 3
        /// objects (abc) has 3! = 6 permutations: (abc), (bac), (cba), (acb), (cab), (bca).</para>
        /// <para>Because n! grows extremely quickly with increasing n, we return the result as a double, even though
        /// the value is always an integer. (13! would overlow an int. 21! would overflow a long. 171! overflows even a double.)</para>
        /// <para>In order to deal with factorials of larger runbers, you can use the <see cref="LogFactorial"/> method, which
        /// returns accurate values of ln(n!) even for values of n for which n! would overflow a double.</para>
        /// <para>The factorial is generalized to non-integer arguments by the &#x393; function (<see cref="AdvancedMath.Gamma(double)"/>).</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso cref="LogFactorial"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Factorial"/>
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
        /// <returns>The natrual logarithm of the factorial: ln(n!).</returns>
        /// <remarks>
        /// <para>This function provides accurate values of ln(n!) even for values of n which would cause n! to overflow.</para>
        /// </remarks>
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
        /// Computes a binomial coefficient.
        /// </summary>
        /// <param name="n">The upper argument, which must be non-negative.</param>
        /// <param name="m">The lower argument, which must be non-negative and less than or equal to <paramref name="n"/>.</param>
        /// <returns>The binomial coefficent C(<paramref name="n"/>,<paramref name="m"/>),
        /// also denoted "<paramref name="n"/> choose <paramref name="m"/>".</returns>
        /// <remarks>
        /// <para>The binomial coefficient C(n,m) is the coefficient of x<sup>m</sup> in the expansion of (1+x)<sup>n</sup>.</para>
        /// <img src="../images/BinomialExpansion.png" />
        /// <para>C(n,m) can also be given a combinatoric intrepretation as the total number of ways to pick m items from a set of
        /// n distinct items.</para>
        /// <para>For example C(4,2) = 6. This can be seen by expanding (1+x)<sup>4</sup> = 
        /// 1 + 4 x + 6 x<sup>2</sup> + 4 x<sup>3</sup> + x<sup>4</sup> and noting that the coefficient of the x<sup>2</sup> term is 6.
        /// It can also be seen by considering the four-member set (abcd) and noting that there are 6 possible two-member subests:
        /// (ab), (ac), (ad), (bc), (bd), (cd).</para>
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
        /// triangle (i.e., all the coefficients for a given order polynomial) starting from an outer values B(n,0) = 1 = B(n,n).
        /// If you need a series of binomial coefficients, using a recursion will be more computationally efficient than
        /// calling this method for each one. The <see cref="BinomialCoefficients"/> method provides a fast enumeration of
        /// all the binomial coefficients in a given row of Pascal's triangle.</para>
        /// <para>Binomial coefficients are always integers, but we return a the result as double because the value can exceed
        /// the capacity of an int or even a long for even quite moderate values of n and m.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative, or <paramref name="m"/> lies outside [0,<paramref name="n"/>].</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Binomial_coefficient"/>
		public static double BinomialCoefficient (int n, int m) {
			if (n < 0) throw new ArgumentOutOfRangeException("n");
			if ((m < 0) || (m > n)) throw new ArgumentOutOfRangeException("m");

            // make lower argument small
            if (m > n / 2) m = n - m;
            // this increases the likelyhood that we will be able to do factorial table lookup and increases the accuracy
            // of the upper bound (e n / m)^k

			if ((n<factorialTable.Length) && (m<factorialTable.Length)) {
                // do table look-up, if possible
				return( factorialTable[n] / factorialTable[m] / factorialTable[n-m] );
			} else if ((n < 60) || (m < 4) || (m * (1.0 + Math.Log(n / m)) < 0.75 * LogLongMax)) {
                // use integer arithmetic recursion, if we will not overflow a long
                // i tested n specifically: n=61 does not overflow but n=62 does
                long B = 1;
                for (int i = 1; i <= m; i++) {
                    B = B * (n - i + 1) / i;
                    // we must do multiplication before division; it would be nice to do division first, in order to keep quantities
                    // as small as possible, but the previous coefficient is not guaranteed to be exactly divisible by i
                }
                return (B);
             } else {
                // our result will be very big, so we probably won't loose accuracy in subtracting big factorial logs
                return (Math.Round(Math.Exp(LogFactorial(n) - LogFactorial(m) - LogFactorial(n - m))));
			}
		}

        private static readonly double LogLongMax = Math.Log(Int64.MaxValue);

#if PAST

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
#endif

        /// <summary>
        /// Evaluates a row of binomial coefficients.
        /// </summary>
        /// <param name="n">The upper argument.</param>
        /// <returns>An enumeration of the binomial coefficients in the nth row of Pascal's triangle.</returns>
        public static IEnumerable<double> BinomialCoefficients (int n) {

            if (n < 0) throw new ArgumentOutOfRangeException("n");

            if (n < 60) {
                // use integer arithmetic if we won't overflow a long
                long B = 1;
                yield return (B);
                for (int k = 1; k <= n; k++) {
                    B = B * (n - k + 1) / k;
                    yield return (B);
                }
            } else {
                // otherwise use a double
                double B = 1.0;
                yield return(B);
                for (int k = 1; k <= n; k++) {
                    B = Math.Round(B / k * (n - k + 1));
                    yield return (B);
                }
                // at first i was worried that the falling B's in the second half of the row would not be reproduced accurately
                // by floating point arithmetic,  but that doesn't appear to be the case; if it does become a problem, just
                // store the first half in an array and move down it in the second half
            }

        }

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
                return (m * Global.LogTwo + AdvancedMath.LogGamma(m + 1.0));
            } else {
                // m = (n+1)/2, n!! = 2^m Gamma(m+1/2) / Sqrt(PI)
                int m = (n + 1) / 2;
                return (m * Global.LogTwo + AdvancedMath.LogGamma(m + 0.5) - 0.5 * Math.Log(Math.PI));
            }
        }

        /// <summary>
        /// Computes the double factorial of the given number.
        /// </summary>
        /// <param name="n">The argument, which must be positive.</param>
        /// <returns>The double factorial n!!.</returns>
        /// <remarks>
        /// <para>The double factorial of an integer is the product all integers of the same parity, up to an including the integer.
        /// Thus 5! = 5 * 3 * 1 = 15 and 6! = 6 * 4 * 2 = 48.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso href="http://mathworld.wolfram.com/DoubleFactorial.html"/>
        public static double DoubleFactorial (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (n < 32) {
                return ((double) DoubleFactorial_Multiply(n));
            } else {
                return (Math.Round(Math.Exp(LogDoubleFactorial_Gamma(n))));
            }
        }

        /// <summary>
        /// Computes the natural logarithm of the double factorial of the given number.
        /// </summary>
        /// <param name="n">The argument.</param>
        /// <returns>The value of ln(n!!).</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso cref="DoubleFactorial"/>
        public static double LogDoubleFactorial (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n"); 
            if (n < 32) {
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

		// Bernoulli and Euler polynomials
        
        // Stirling numbers


        /// <summary>
        /// Computes a power of an integer in modular arithmetic.
        /// </summary>
        /// <param name="b">The base, which must be positive.</param>
        /// <param name="e">The exponent, which must be positive.</param>
        /// <param name="m">The modulus, which must be positive.</param>
        /// <returns>The value of b<sup>e</sup> mod m.</returns>
        /// <remarks>
        /// <para>Modular exponentiation is used in many number-theory applications, including
        /// primality testing, prime factorization, and cryptography.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="b"/>,  <paramref name="e"/>, or <paramref name="m"/> is not positive.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Modular_exponentiation"/>
        public static int PowMod (int b, int e, int m) {
            if (b < 1) throw new ArgumentOutOfRangeException("b");
            if (e < 1) throw new ArgumentOutOfRangeException("e");
            if (m < 1) throw new ArgumentOutOfRangeException("m");

            // use long internally
            // since the "worst" we ever do before modding is to square, and since a long should
            // hold twice as many digits as an int, this algorithm should not overflow
            long bb = Convert.ToInt64(b);
            long mm = Convert.ToInt64(m);
            long rr = 1;

            while (e > 0) {
                if ((e & 1) == 1) {
                    rr = checked((rr * bb) % mm);
                }
                e = e >> 1;
                bb = checked((bb * bb) % mm);
            }

            return (Convert.ToInt32(rr));

        }

        /// <summary>
        /// Enumerates all partitions of the given integer
        /// </summary>
        /// <param name="n">The integer to partition, which must be positive.</param>
        /// <returns>An enumeration of all partitions of the given integer.</returns>
        /// <remarks>
        /// <para>Integer partitions are ways to write an integer as a sum of smaller integers. For example, the integer 4 has 5 partitions: 4,
        /// 3 + 1, 2 + 2, 2 + 1 + 1, and 1 + 1 + 1 + 1.</para>
        /// <para>Integer partitions appear in combinatoric problems and solutions to problems that may be mapped into combinatoric problems.
        /// For example, the terms which appear in <a href="http://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula">Faà di Bruno's formula</a>
        /// correspond to integer partitions.</para>
        /// <para>The number of partitions grows very rapidly with n. Since enumerating through partitions does not require us to count them,
        /// no overflows will occur even for large values of <paramref name="n"/>. However, completing the enumeration of
        /// such a large number of paritions will take a long time, even though our algorithm produces each partition very quickly. For
        /// example, there are about two hundred million partitions of the integer 100.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is not positive.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Integer_partition"/>
        public static IEnumerable<int[]> Partitions (int n) {

            if (n < 1) throw new ArgumentOutOfRangeException("n");

            // initialize the state
            int[] a = new int[n + 1];
            int k = 1;
            a[0] = 0;
            a[1] = n;

            while (k != 0) {

                // advance to the next partition
                int y = a[k] - 1;
                k--;
                int x = a[k] + 1;
                while (x <= y) {
                    a[k] = x;
                    y -= x;
                    k++;
                }
                a[k] = x + y;

                // return a copy so our state array can't be disturbed
                int[] p = new int[k + 1];
                Array.Copy(a, p, k + 1);
                yield return (p);
            }

        }

#if FUTURE

        private static readonly int[] primes = new int[] { 2, 3, 5, 7, 11, 13 };

        public static List<int> FactorByTrial (ref int n) {

            List<int> factors = new List<int>();

            foreach (int p in primes) {

                if (n % p == 0) {
                    n = n / p;
                    factors.Add(p);
                    if (n < p) return (factors);
                }

            }

            return (factors);

        }

        private static int pollardE = 360360;
        // lcm(1...5) = 60
        // lcm(1...10) = 2520
        // lcm(1...15) = 360360
        // lcm(1...20) = 232792560

        public static List<int> FactorByPollard1 (ref int n) {

            List<int> factors = new List<int>();

            while (true) {
                int g = PowMod(2, pollardE, n) - 1;
                g = (int) GCF(g, n);
                if ((g > 1) && (g < n)) {
                    factors.Add(g);
                    n = n / g;
                } else {
                    return (factors);
                }
            }

        }

#endif

	}

}
