using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

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
        /// <returns>The value of n!.</returns>
        /// <remarks>
        /// <para>The factorial of an integer n is the product of all integers from 1 to n. For example, 4! = 4 * 3 * 2 * 1 = 24.</para>
        /// <para>n! also has a combinatorial intrepretation as the number of permutations of n objects. For example, a set of 3
        /// objects (abc) has 3! = 6 permutations: (abc), (bac), (cba), (acb), (cab), (bca).</para>
        /// <para>Because n! grows extremely quickly with increasing n, we return the result as a double, even though
        /// the value is always an integer. (13! would overlow an int, 21! would overflow a long, 171! overflows even a double.)</para>
        /// <para>In order to deal with factorials of larger numbers, you can use the <see cref="LogFactorial"/> method, which
        /// returns accurate values of ln(n!) even for values of n for which n! would overflow a double.</para>
        /// <para>The factorial is generalized to non-integer arguments by the &#x393; function (<see cref="AdvancedMath.Gamma(double)"/>).</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso cref="LogFactorial"/>
        /// <seealso cref="AdvancedMath.Gamma(double)"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Factorial"/>
		public static double Factorial (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n < factorialTable.Length) {
				return( (double) factorialTable[n]);
			} else {
				return( Math.Round( AdvancedMath.Gamma(n+1) ) );
			}
		}

        /// <summary>
        /// Computes the logrithm of the factorial of an integer.
        /// </summary>
        /// <param name="n">The argument, which must be non-negative.</param>
        /// <returns>The value of ln(n!).</returns>
        /// <remarks>
        /// <para>This function provides accurate values of ln(n!) even for values of n which would cause n! to overflow.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso cref="Factorial"/>
		public static double LogFactorial (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n < factorialTable.Length) {
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
        /// <para>Binomial coefficients are always integers, but we return the result as double because the value can exceed
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

        /// <summary>
        /// Enumerates the binomial coefficients of a given order.
        /// </summary>
        /// <param name="n">The upper argument, which must be non-negative.</param>
        /// <returns>An enumeration of the binomial coefficients in the nth row of Pascal's triangle.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
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
                return (m * Global.LogTwo + AdvancedMath.LogGamma(m + 0.5) - Math.Log(Math.PI) / 2.0);
            }
        }

        /// <summary>
        /// Computes the double factorial of the given integer.
        /// </summary>
        /// <param name="n">The argument, which must be positive.</param>
        /// <returns>The value of n!!.</returns>
        /// <remarks>
        /// <para>The double factorial of an integer is the product all integers of the same parity, up to and including the integer.
        /// Thus 5! = 5 * 3 * 1 = 15 and 6! = 6 * 4 * 2 = 48.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso href="http://mathworld.wolfram.com/DoubleFactorial.html"/>
        public static double DoubleFactorial (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n < 32) {
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
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n < 32) {
                return (Math.Log((double) DoubleFactorial_Multiply(n)));
            } else {
                return (LogDoubleFactorial_Gamma(n));
            }
        }


        /// <summary>
        /// Computes the given harmonic number.
        /// </summary>
        /// <param name="n">The index of the harmonic number to compute, which must be non-negative.</param>
        /// <returns>The harmonic number H<sub>n</sub>.</returns>
        /// <remarks>
        /// <para>H<sub>n</sub> is the nth partial sum of the harmonic series.</para>
        /// <para>Since the harmonic series diverges, H<sub>n</sub> grows without bound as n increases, but
        /// it does so extremely slowly, approximately as log(n).</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Harmonic_series_(mathematics)"/>
        public static double HarmonicNumber (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n < 32) {
                // for small values, just add up the harmonic series
                double H = 0.0;
                for (int i = 1; i <= n; i++) {
                    H += 1.0 / i;
                }
                return (H);
            } else {
                // for large values, use the digamma function
                // this will route to the the Stirling asymptotic expansion
                return (AdvancedMath.Psi(n+1) + AdvancedMath.EulerGamma);
            }

        }

        /// <summary>
        /// Computes a Fibonacci number.
        /// </summary>
        /// <param name="n">The index of the Fibonacci number to compute, which must be non-negative.</param>
        /// <returns>The nth Fibonacci number F<sub>n</sub>.</returns>
        /// <see href="http://en.wikipedia.org/wiki/Fibonacci_number"/>
        public static double FibonacciNumber (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else {
                // Begin with relation F_n = \frac{\phi^n - (-\phi)^{-n}}{\sqrt{5}}, which is exact. For large n, second
                // term becomes negligable, so we certainly expect F_n = \left[ \frac{\phi^n}{\sqrt{5}} \right] for large enough n.
                // What's surprising is that n = 0 is already large enough, giving us a really fast way to get the nth Fibonacci
                // number (at least to floating point precision) without computing any lower Fibonacci numbers, much less all of them.
                return (Math.Round(MoreMath.Pow(AdvancedMath.GoldenRatio, n) / Math.Sqrt(5.0)));
            }
        }

        /// <summary>
        /// Computes the given Bernoulli number.
        /// </summary>
        /// <param name="n">The index of the Bernoulli number to compute, which must be non-negative.</param>
        /// <returns>The Bernoulli number B<sub>n</sub>.</returns>
        /// <remarks>
        /// <para>B<sub>n</sub> vanishes for all odd n except n=1. For n about 260 or larger, B<sub>n</sub> overflows a double.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Bernoulli_number"/>
        /// <seealso href="http://mathworld.wolfram.com/BernoulliNumber.html"/>
        public static double BernoulliNumber (int n) {

            if (n < 0) throw new ArgumentOutOfRangeException("n");

            // B_1 is the only odd Bernoulli number.
            if (n == 1) return (-1.0 / 2.0);

            // For all other odd arguments, return zero.
            if (n % 2 != 0) return (0.0);

            // If the argument is small enough, look up the answer in our stored array.
            int m = n / 2;
            if (m < Bernoulli.Length) return (Bernoulli[m]);

            // Otherwise, use the relationship with the Riemann zeta function to get the Bernoulli number.
            // Since this is only done for large n, it would probably be faster to just sum the zeta series explicitly here.
            double B = 2.0 * AdvancedMath.RiemannZeta(n) / AdvancedMath.PowOverGammaPlusOne(Global.TwoPI, n);
            if (m % 2 == 0) B = -B;
            return (B);

        }

        // the even Bernoulli numbers B_2n = Bernoulli[n]
        // the only nonvanishing odd Bernoulli number is B_1 = -1/2, which must be handled seperately if you
        // use these numbers in any series expansion

        internal static readonly double[] Bernoulli = new double[] {
            1.0, 1.0 / 6.0, -1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0, 5.0 / 66.0, -691.0 / 2730.0, 7.0 / 6.0,
            -3617.0 / 510.0, 43867.0 / 798.0, -174611.0 / 330.0, 854513.0 / 138.0, -236364091.0 / 2730.0, 8553103.0 / 6.0, -23749461029.0 / 870.0, 8615841276005.0 / 14322.0
        };

        // the expansions in which they appear are asymptotic; the numbers grow rapidly after ~B_16

        /// <summary>
        /// Computes a Stirling number of the first kind.
        /// </summary>
        /// <param name="n">The upper argument, which must be non-negative.</param>
        /// <param name="k">The lower argument, which must lie between 0 and n.</param>
        /// <returns>The value of the unsigned Stirling number of the second kind.</returns>
        public static double StirlingNumber1 (int n, int k) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
            if ((k < 0) || (k > n)) throw new ArgumentOutOfRangeException(nameof(k));

            if (k == n) {
                return (1.0);
            } else if (k == 0) {
                return (0.0);
            } else if (k == 1) {
                return (AdvancedIntegerMath.Factorial(n - 1));
            } else if (k == (n - 1)) {
                return (AdvancedIntegerMath.BinomialCoefficient(n, 2));
            } else {
                double[] s = Stirling1_Recursive(n, k);
                return (s[k]);
            }

        }

        /// <summary>
        /// Computes a row of Sterling numbers of the first kind.
        /// </summary>
        /// <param name="n">The upper argument, which must be non-negative.</param>
        /// <returns>An array with n+1 elements. The element with (zero-based) index k contains the
        /// unsigned Sterling number of the first kind with upper argument n and lower argument k.</returns>
        public static double[] StirlingNumbers1 (int n) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));

            if (n == 0) {
                return (new double[] { 1.0 });
            } else {
                return (Stirling1_Recursive(n, n));
            }

        }

        private static double[] Stirling1_Recursive (int n, int kMax) {

            Debug.Assert(n > 0);
            Debug.Assert(kMax > 0);
            Debug.Assert(kMax <= n);

            double[] s = new double[kMax + 1];

            // Seed j = 1 row
            s[0] = 0.0;
            s[1] = 1.0;

            // Compute higher rows
            for (int j = 2; j <= n; j++) {

                int jkMax;
                if (j <= kMax) {
                    s[j] = 1.0;
                    jkMax = j - 1;
                } else {
                    jkMax = kMax;
                }

                for (int k = jkMax; k > 0; k--) {
                    s[k] = (j - 1) * s[k] + s[k - 1];
                }

            }

            return (s);

        }

        /// <summary>
        /// Computes a Stirling number of the second kind.
        /// </summary>
        /// <param name="n">The upper argument, which must be non-negative.</param>
        /// <param name="k">The lower argument, which must lie between 0 and n.</param>
        /// <returns>The value of the Stirling number of the second kind.</returns>
        public static double StirlingNumber2 (int n, int k) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
            if ((k < 0) || (k > n)) throw new ArgumentOutOfRangeException(nameof(k));

            if ((k == 1) || (k == n)) {
                return (1.0);
            } else if (k == 0) {
                return (0.0);
                // The exceptional value 1 for n = k = 0 will already have been returned by the previous case. 
            } else if (k == 2) {
                return (Math.Round(MoreMath.Pow(2.0, n - 1) - 1.0));
            } else if (k == (n - 1)) {
                return (AdvancedIntegerMath.BinomialCoefficient(n, 2));
            } else {
                double[] s = Stirling2_Recursive(n, k);
                return (s[k]);
            }

            // There is a formula for Stirling numbers
            //   { n \brace k } = \frac{1}{k!} \sum{j=0}^{k} (-1)^j { k \choose j} (k - j)^n
            // which would be faster than recursion, but it has large cancelations between
            // terms. We could try to use it when all values are less than 2^52, for which
            // double arithmetic is exact for integers. For k!, that means k < 18. For 
            // largest term in sum for all k, that means n < 14.

        }

        /// <summary>
        /// Computes a row of Stirling numbers of the second kind.
        /// </summary>
        /// <param name="n">The upper argument, which must be non-negative.</param>
        /// <returns>An array with n+1 elements. The element with (zero-based) index k contains the
        /// Stirling number of the second kind with upper argument n and lower argument k.</returns>
        public static double[] StirlingNumbers2 (int n) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));

            if (n == 0) {
                return (new double[] { 1.0 });
            } else {
                return (Stirling2_Recursive(n, n));
            }

        }

        private static double[] Stirling2_Recursive (int n, int kMax) {

            Debug.Assert(n > 0);
            Debug.Assert(kMax > 0);
            Debug.Assert(kMax <= n);

            double[] s = new double[kMax + 1];

            // Seed j = 1 row
            s[0] = 0.0;
            s[1] = 1.0;

            // Compute higher rows
            for (int j = 2; j <= n; j++) {

                int jkMax;
                if (j <= kMax) {
                    s[j] = 1.0;
                    jkMax = j - 1;
                } else {
                    jkMax = kMax;
                }

                for (int k = jkMax; k > 1; k--) {
                    s[k] = k * s[k] + s[k - 1];
                }

            }

            return (s);

        }

        // There are two well-known algorithms for computing Bell numbers. One is Dobinski's formula
        //   B_n = \frac{1}{e} \sum_{k=0}^{\infty} \frac{k^n}{k!}
        // which is basically just the relationship between the Bell numbers and the Poisson moments.
        // The other is the triangle relationship. It computes all Bell numbers up to the one sought,
        // but is O(n^2) in time and O(n) in memory and is slower when a single Bell number is sought.

        // B_n excees the capacity of a double for n >= 219.

        // Very approximately, the terms in Dobinski's formula increase up to k ~ n / W(n) ~ n / log n
        // before they start to decrease, where W(n) is the Lambert W function. It might be worthwhile
        // to start at this value of k and move outwards.

        // As implemented, this returns Infinity too soon, at n=165, because the power in the numerator
        // overflows. We need to deal with that by calling PowOverGamma instead.

        /// <summary>
        /// Computes the given Bell number.
        /// </summary>
        /// <param name="n">The index of the Bell number to compute, which must be non-negative.</param>
        /// <returns>The Bell number B<sub>n</sub>.</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Bell_number"/>
        public static double BellNumber (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");

            // Dobinski's formula only good for n > 0
            if (n == 0) return (1.0);

            double s = 1.0;
            double q = 1.0;
            for (int k = 2; k < Global.SeriesMax; k++) {
                double s_old = s;
                q *= k;
                s += MoreMath.Pow(k, n) / q;
                if (s == s_old) {
                    return (Math.Round(s / Math.E));
                }
            }
            
            throw new NonconvergenceException();
        }


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

        internal static int GCD (int u, int v) {
            while (v != 0) {
                int t = u % v; u = v; v = t;
            }
            return (u);
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
            if (b < 0) throw new ArgumentOutOfRangeException("b");
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

        /// <summary>
        /// Determines whether the given integer is prime.
        /// </summary>
        /// <param name="n">The integer, which must be positive.</param>
        /// <returns>True if the integer is prime, otherwise false.</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Prime_number"/>
        /// <seealso href="http://mathworld.wolfram.com/PrimeNumber.html"/>
        public static bool IsPrime (int n) {
            if (n <= 0) {
                throw new ArgumentOutOfRangeException("n");
            } if (n < 8) {
                return ((n == 2) || (n == 3) || (n == 5) || (n == 7));
            } else {

                // even numbers (larger than 2) aren't prime
                if (n % 2 == 0) return (false);

                // we will use the Miller-Rabin probabilistic prime testing; to do so we must write m = n - 1 = 2^s * d
                int m = n - 1;
                int d = m;
                int s = 0;
                while (d % 2 == 0) {
                    s++;
                    d = d / 2;
                }

                // according to articles quoted at http://primes.utm.edu/prove/prove2_3.html
                // witnesses 2 and 3 are sufficient for all n < 1,373,653 
                // witnesses 2, 3, and 5 are sufficient for all n < 25,326,001
                // witnesses 2, 3, 5, and 7 are sufficient for all n < 118,670,087,467 except n = 3,215,031,751
                // witnesses 2, 3, 5, 7, and 11 are sufficient for all n < 2,152,302,898,747
                // witnesses 31 and 73 are sufficient for all n < 9,080,191
                // witnesses 2, 7, and 61 are sufficient for all n < 4,759,123,141

                if (n < 1373653) {
                    return (IsProbablyPrime(n, m, s, d, 2) && IsProbablyPrime(n, m, s, d, 3));
                } else {
                    return (IsProbablyPrime(n, m, s, d, 2) && IsProbablyPrime(n, m, s, d, 7) && IsProbablyPrime(n, m, s, d, 61));
                }

            }
        }

        // The Miller-Rabin probable prime test http://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test

        // The odd integer n = m + 1 = 2^s * d + 1
        // The witness is w, which must be 2 <= w <= n - 2

        // should loop be 0...(s-1) or 1...(s-1)?

        private static bool IsProbablyPrime (int n, int m, int s, int d, int w) {
            int x = PowMod(w, d, n);
            if ((x == 1) || (x == m)) return (true);
            for (int i = 0; i < s; i++) {
                x = PowMod(x, 2, n);
                if (x == 1) return (false);
                if (x == m) return (true);
            }
            return (false);
        }

        // Prime factorization. Leave this internal for now, until it is cleaned up.

        // As currently implemented, it does not actually guarantee full prime factorization! Pollard's rho
        // method can yield non-prime factors, and this appears to occur for about 0.25% of all integers under 1,000,000.
        // For example, "factors" of 1681 = 41 * 41, 6751 = 43 * 157, and 9167 = 89 * 103 are claimed.
        // These composite "factors" are, however, still co-prime to the other factors, so the almost-factorization
        // will still work for reduction of Fourier transforms, which is how we are currently using it.

        internal static List<Element> Factor (int n) {
            if (n < 1) throw new ArgumentOutOfRangeException("n");

            List<Element> factors = new List<Element>();

            if (n > 1) FactorByTrialDivision(factors, ref n);

            if (n > 1) FactorByPollardsRhoMethod(factors, ref n);

            if (n > 1) factors.Add(new Element(n, 1));

            return(factors);
        }

        internal static readonly int[] SmallPrimes = new int[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };

        // Trial division is the simplest prime factorization method. It consists of attempted to divide by known primes.
        // It is a good way to eliminate known small prime factors before proceeding on to bigger and more difficult prime factors.

        private static void FactorByTrialDivision (List<Element> factors, ref int n) {

            foreach (int p in SmallPrimes) {

                int m = 0;
                while (n % p == 0) {
                    n = n / p;
                    m++;
                }
                if (m > 0) factors.Add(new Element(p, m));

                if (n == 1) return;

            }

        }

        private static void FactorByPollardsRhoMethod (List<Element> factors, ref int n) {

            int x = 5; int y = 2; int k = 1; int l = 1;

            for (int c = 0; c < Global.SeriesMax; c++) {
                //while (true) {
                int g = AdvancedIntegerMath.GCD(Math.Abs(y - x), n);
                if (g == n) {
                    // the factor n will repeat itself indefinitely; either n is prime or the method has failed
                    return;
                } else if (g == 1) {
                    k--;
                    if (k == 0) {
                        y = x;
                        l = 2 * l;
                        k = l;
                    }
                    // take x <- (x^2 + 1) mod n
                    x = AdvancedIntegerMath.PowMod(x, 2, n) + 1;
                    if (x == n) x = 0;
                } else {
                    // g is a factor of n; in all likelyhood, it is prime, although this isn't guaranteed
                    // for our current approximate-factoring purposes, we will assume it is prime
                    // it is at least co-prime to all other recognized factors
                    int m = 0;
                    while (n % g == 0) {
                        n = n / g;
                        x = x % n;
                        y = y % n;
                        m++;
                    }
                    factors.Add(new Element(g, m));
                }
            }

        }

	}


    internal struct Element {

        internal Element (int value, int multiplicity) {
            this.value = value;
            this.multiplicity = multiplicity;
        }

        private readonly int value, multiplicity;

        public int Value {
            get {
                return (value);
            }
        }

        public int Multiplicity {
            get {
                return (multiplicity);
            }
        }

    }


}
