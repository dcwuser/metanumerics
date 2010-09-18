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
        /// <para>Because n! grows extremely quickly with increasing n, we return the result as a double, even though
        /// the value is always an integer. (13! would overlow an int. 21! would overflow a long. 171! overflows even a double.)</para>
        /// <para>In order to deal with factorials of larger runbers, you can use the <see cref="LogFactorial"/> method, which
        /// returns accurate values of ln(n!) even for values of n for which n! itself overflows a double.</para>
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
        /// Computes binomial coefficients.
        /// </summary>
        /// <param name="n">The upper argument, which must be non-negative. (The order of the polynomial.)</param>
        /// <param name="m">The lower argument, which must be non-negative and less than or equal to <paramref name="n"/>. (The order of the term.)</param>
        /// <returns>The binomial coefficent C(<paramref name="n"/>,<paramref name="m"/>),
        /// also denoted "<paramref name="n"/> choose <paramref name="m"/>".</returns>
        /// <remarks>
        /// <para>The binomial coefficient C(n,m) is the coefficient of x<sup>m</sup> in the expansion
        /// of (1+x)<sup>n</sup>.</para>
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
        /// <seealso href="http://en.wikipedia.org/wiki/Binomial_coefficient"/>
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
        /// <returns></returns>
        /// <remarks>
        /// <para>Integer partitions are ways to write an integer as a sum of smaller integers. For example, the integer 4 has 5 partitions: 4,
        /// 3 + 1, 2 + 2, 2 + 1 + 1, and 1 + 1 + 1 + 1.</para>
        /// <para>Integer partitions appear in combinatoric problems and solutions to problems that may be mapped into combinatoric problems.
        /// For example, the terms which appear in <a href="http://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno%27s_formula">Faà di Bruno's formula</a>
        /// correspond to integer partitions.</para>
        /// <para>The number of partitions grows very rapidly with n. Already for n = 122, the number of partitions is
        /// 2,291,320,912, larger than that the capacity of an int. Since enumerating through partitions does not require us to count them,
        /// this method will work for this and larger values of <paramref name="n"/>. Just keep in mind that completing the enumeration of
        /// such a large number of paritions is liable to take a long time, even though our algorithm produces each partition very quickly.</para>
        /// </remarks>
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

        public static void Permute (int n) {

            int[] p = new int[n + 1];
            p[n] = 1;

            // print the permutation
            for (int j = n; j >= 1; j--) {
                for (int i = 1; i <= p[j]; i++) {
                    Console.Write("+{0}", j);
                }
            }
            Console.WriteLine();

            int count = 1;

            while (p[1] != n) {

                // compute next permutation

                // find lowest number over one
                int k = 2;
                while (p[k] == 0) {
                    k++;
                }

                // reduce it by one
                p[k]--;

                // determine how much we have to sum to
                int s = n;
                for (int j = k; j <= n; j++) {
                    s -= j * p[j];
                }

                // produce that sum as quickly as we can with numbers less than or equal to k
                for (int j = k - 1; j >= 1; j--) {
                    p[j] = 0;
                    while (s >= j) {
                        p[j]++;
                        s -= j;
                    }
                }

                // print the permutation
                for (int j = n; j >= 1; j--) {
                    for (int i = 1; i <= p[j]; i++) {
                        Console.Write("+{0}", j);
                    }
                }
                Console.WriteLine();

                count++;

            }

            Console.WriteLine(count);

        }

        /*
        def ruleAsc(n):
    a = [0 for i in range(n + 1)]
    k = 1
    a[0] = 0
    a[1] = n
    while k != 0:
        x = a[k - 1] + 1
        y = a[k] - 1
        k -= 1
        while x <= y:
            a[k] = x
            y -= x
            k += 1
        a[k] = x + y
        yield a[:k + 1]
        */

        public static void P2 (int n) {

            int[] a = new int[n + 2];
            //int[] a = new int[n + 1];

            int k = 2;
            //int k = 1;
            a[1] = 0;
            //a[0] = 0;
            a[2] = n;
            //a[1] = n;

            while (k != 1) {
            //while (k != 0) {
                int y = a[k] - 1;
                k--;
                int x = a[k] + 1;
                while (x <= y) {
                    a[k] = x;
                    y -= x;
                    k++;
                }
                a[k] = x + y;

                // print the permutation
                Console.WriteLine("k={0}",k);
                for (int j = 0; j < a.Length; j++) {
                    Console.Write("{0} ", a[j]);
                }
                Console.WriteLine();

            }


        }
#endif

	}

#if FUTURE
    public class IntegerPartitionEnumerator : IEnumerator<int[]> {

        public IntegerPartitionEnumerator (int n) {
            this.n = n;
            a = new int[n + 1];
            Reset();
        }

        // the integer being partioned
        int n;

        // store for the current partition (only the first k elements are relevent)
        int[] a;

        // the number of terms in the current partition
        int k;

        public int[] Current {
            get {
                // if (a[0] == 0) throw new InvalidOperationException();
                int[] p = new int[k+1];
                Array.Copy(a, p, k + 1);
                //Array.Copy(a, 1, p, 0, k+1);
                return (p);
            }
        }

        object IEnumerator.Current {
            get {
                return (Current);
            }
        }

        public void Reset () {
            k = 1;
            a[0] = 0;
            a[1] = n;
        }

        public bool MoveNext () {

            if (k == 0) {
                return (false);
            } else {
                int y = a[k] - 1;
                k--;
                int x = a[k] + 1;
                while (x <= y) {
                    a[k] = x;
                    y -= x;
                    k++;
                }
                a[k] = x + y;

                /*
                // print the permutation
                Console.WriteLine("k={0}", k);
                for (int j = 0; j < a.Length; j++) {
                    Console.Write("{0} ", a[j]);
                }
                Console.WriteLine();
                */

                return (true);
            }

        }

        void IDisposable.Dispose () {
        }

    }
#endif
	
}
