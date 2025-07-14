
using System;
using System.Diagnostics;
using System.Collections.Generic;

namespace Meta.Numerics.Functions {

    /// <summary>
    /// Contains methods that compute the values of orthogonal polynomials.
    /// </summary>
    /// <remarks>
    /// <para>Orthogonal polynomials are complete families of polynomials that are orthogonal on a given interval with
    /// a given integration weight. Because of this property, any function on the interval can be expanded
    /// in the polynomials in a unique way.</para>
    /// <para>Because their fluctuations across the interval are driven by cancellations among multiple terms,
    /// polynomials are notoriously difficult to evaluate with guaranteed numerical accuracy over all
    /// points in their domain. As order increases, our methods loose some accuracy. Each method description
    /// gives an example of the accuracy that can be expected for different argument regiemes.</para>
    /// <para>By the way, the Wikipedia article on the classical orthogonal polynomials is particularly good.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Classical_orthogonal_polynomials"/>
    /// <seealso href="https://mathworld.wolfram.com/OrthogonalPolynomials.html"/>
	public static class OrthogonalPolynomials {


        /// <summary>
        /// Computes the value of a (physicists') Hermite polynomial.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value H<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The Hermite polynomials are orthogonal on the interval (-&#8734;,+&#8734;) with the
        /// weight e<sup>-x<sup>2</sup></sup>.</para>
        /// <img src="../images/HermiteHOrthonormality.png" />
        /// <para>They appear in the solution of the one-dimensional, quantum mechanical, harmonic oscillator.</para>
        /// <para>Statisticans' Hermite polynomials (see <see cref="HermiteHe"/>) are related to physicists' Hermite
        /// polynomials via H<sub>n</sub>(x) = 2<sup>n</sup>H<sub>n</sub>(x &#x221A;2). Statisticians' Hermite polynomials
        /// do not grow as quickly as physicists', and may therefore be preferable for large values of <paramref name="n"/>
        /// and <paramref name="x"/>, for which the return value of this method can overflow.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso cref="HermiteHe"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Hermite_polynomials"/>
        /// <seealso href="http://mathworld.wolfram.com/HermitePolynomial.html" />
		public static double HermiteH (int n, double x) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            } else if (n == 0) {
                return 1.0;
            } else {
                // Use the recurrence H_{n+1} = 2x H_n - 2n H_{n-1}.
                // The recurrence is unstable, but H is the dominant solution.
                double H0 = 1.0;
                double H1 = 2.0 * x;
                for (int k = 1; k < n; k++) {
                    double H2 = 2.0 * (x * H1 - k * H0);
                    H0 = H1;
                    H1 = H2;
                }
                return H1;
            }
		}

        /// <summary>
        /// Computes the value of a (statisticians') Hermite polynomial.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="x">The argument.</param>
        /// <returns>The value He<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>The (statisticians') Hermite polynomials are orthogonal on the interval (-&#8734;,+&#8734;) with a
        /// weight function equal to the standard normal probability distribution.</para>
        /// <img src="../images/HermiteHeOrthonormality.png" />
        /// <para>Their ortho-normality relation makes them a useful basis for expressing perturbations
        /// around a normal distribution.</para>
        /// <para>Physicists' Hermite polynomials (<see cref="HermiteH"/>) are related to statisticians' Hermite
        /// polynomials via H<sub>n</sub>(x) = 2<sup>n</sup>H<sub>n</sub>(x &#x221A;2).</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        /// <seealso cref="HermiteH"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Hermite_polynomials" />
        public static double HermiteHe (int n, double x) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            } else if (n == 0) {
                return 1.0;
            } else {
                double H0 = 1.0;
                double H1 = x;
                for (int k = 1; k < n; k++) {
                    double H2 = x * H1 - k * H0;
                    H0 = H1;
                    H1 = H2;
                }
                return H1;
            }
        }

		// orthogonal on [0,Infinity] with weight e^{-x}
        /// <summary>
        /// Computes the value of a Laguerre polynomial.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value L<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>Laguerre functions are orthogonal on the interval [0,+&#8734;) with the weight e<sup>-x</sup>.</para>
        /// <img src="../images/LaguerreLOrthonormality.png" />
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> or <paramref name="x"/> is negative.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Laguerre_polynomial" />
        /// <seealso href="http://mathworld.wolfram.com/LaguerrePolynomial.html" />
        /// <seealso cref="LaguerreL(int,double,double)"/>
		public static double LaguerreL (int n, double x) {
			if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
			if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

			if (n == 0) return 1.0;

            // Recurrence (n+1)L_{n+1} = (2n+1-x)L_{n} - nL_{n-1}
            double L0 = 1.0;
			double L1 = 1.0 - x;

			for (int k = 1; k < n; k++) {
				double L2 = ( (2 * k + 1 - x) * L1 - k * L0 ) / (k + 1);
				L0 = L1;
				L1 = L2;
			}

			return L1;
		}



        /// <summary>
        /// Computes the value of an associated Laguerre polynomial.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="alpha">The associated order, which must be greater than -1.</param>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value L<sub>n</sub><sup>a</sup>(x).</returns>
        /// <remarks>
        /// <para>The associated Laguerre polynomials are orthogonal on the interval [0,+&#8734;) with the weight
        /// x<sup>a</sup> e<sup>-x</sup>.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> of <paramref name="x"/> is negative,
        /// or <paramref name="alpha"/> is less than or equal to -1.</exception>
        /// <seealso href="https://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html" />
        public static double LaguerreL (int n, double alpha, double x) {
			if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
			if (alpha <= -1) throw new ArgumentOutOfRangeException(nameof(alpha)); 
			if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));

            double L0 = 0.0;
            double L1 = 1.0;
            for (int k = 0; k < n; k++) {
                double L2 = ((2 * k + 1 + alpha - x) * L1 - (k + alpha) * L0) / (k + 1);
                L0 = L1;
                L1 = L2;
            }
            return L1;

		}

        /// <summary>
        /// Computes the value of a Legendre polynomial.
        /// </summary>
        /// <param name="ell">The order, which must be non-negative.</param>
        /// <param name="x">The argument, which must lie on the closed interval between -1 and +1.</param>
        /// <returns>The value of P<sub>l</sub>(x).</returns>
        /// <remarks>
        /// <para>Legendre polynomials are orthogonal on the interval [-1,1].</para>
        /// <img src="../images/LegendrePOrthonormality.png" />
        /// <para>The values returned by this are fully accurate (14-16 decimal digits) over
        /// the full range of arguments for orders up to one million.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> lies outside [-1,+1].</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Legendre_polynomials"/>
        /// <seealso href="http://mathworld.wolfram.com/LegendrePolynomial.html"/>
        public static double LegendreP (int ell, double x) {

			if (Math.Abs(x) > 1.0) throw new ArgumentOutOfRangeException(nameof(x));

            if (ell < 0) {
                return LegendreP(-(ell + 1), x);
            } else if (ell == 0) {
                return 1.0;
            } else {
                // use recurrence (n+1) P_{n+1} = (2n+1) x P_{n} - n P_{n-1}
                double P0 = 1.0;
                double P1 = x;
                for (int n = 1; n < ell; n++) {
                    double P2 = ((2 * n + 1) * x * P1 - n * P0) / (n + 1);
                    P0 = P1;
                    P1 = P2;
                }
                return P1;
            }
		}


        // Legendre polynomials normalized for their use in the spherical harmonics

        /// <summary>
        /// Computes the value of an associated Legendre polynomial.
        /// </summary>
        /// <param name="ell">The order, which must be non-negative.</param>
        /// <param name="m">The associated order, which must lie between -&#x2113; and +&#x2113; inclusive.</param>
        /// <param name="x">The argument, which must lie on the closed interval between -1 and +1.</param>
        /// <returns>The value of P<sub>l,m</sub>(x).</returns>
        /// <remarks>
        /// <para>Associated Legendre polynomials appear in the definition of the <see cref="AdvancedMath.SphericalHarmonic"/> functions.</para>
        /// <para>For values of l and m over about 150, values of this polynomial can exceed the capacity of double-wide floating point numbers.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="ell"/> is negative, <paramref name="m"/> lies outside {-&#x2113; .. +&#x2113;}, or
        /// <paramref name="x"/> lies outside [-1, +1].</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Associated_Legendre_polynomials"/>
        public static double LegendreP (int ell, int m, double x) {

            if (ell < 0) throw new ArgumentOutOfRangeException(nameof(ell));
            if (Math.Abs(m) > ell) throw new ArgumentOutOfRangeException(nameof(m));
            if (Math.Abs(x) > 1.0) throw new ArgumentOutOfRangeException(nameof(x));

            //if (l < 0) {
            //    return (LegendreP(-l + 1, m, x));
            //}

            double f;
            if (ell < 10) {
                // For low enough orders, we can can get the factorial quickly from a table look-up and without danger of overflow.
                f = Math.Sqrt(AdvancedIntegerMath.Factorial(ell + m) / AdvancedIntegerMath.Factorial(ell - m));
            } else {
                // For higher orders, use Pochhammer evaluation to avoid overflow. (Note that if the Pochhammer symbol is > ~1.0E154, this still overflows.)
                //f = Math.Sqrt(AdvancedMath.Pochhammer(l - m + 1, 2 * m));
                f = Math.Exp((AdvancedIntegerMath.LogFactorial(ell + m) - AdvancedIntegerMath.LogFactorial(ell - m)) / 2.0);
            }

            if (m < 0) {
                m = -m;
                if (m % 2 != 0) f = -f;
            }

            return f * LegendrePe(ell, m, x);

        }


        // renormalized associated Legendre polynomials Pe{l,m} = sqrt( (l-m)! / (l+m)! ) P{l,m}
        // unlike the un-renormalized P{l,m}, the renormalized Pe{l,m} do not get too big

        // this is not quite the same renormalization used by NR; it omits a factor sqrt((2l+1)/4Pi)
        // by omitting this factor, we avoid some unnecessary factors and divisions by 4Pi

        // the l-recursion (l-m) P{l,m} = x (2l-1) P{l-1,m} - (l+m-1) P{l-2,m} becomes
        // sqrt((l-m)(l+m)) P{l,m} = x (2l-1) P{l-1,m} - sqrt((l-1-m)(l-1+m)) P{l-2,m}
        // this is stable for increasing l

        // the initial value P{m,m} = (-1)^m (2m-1)!! (1-x^2)^(m/2) becomes
        // Pe{m,m} = (-1)^m (2m-1)!! sqrt( (1-x^2)^m / (2m)! ) = (-1)^m sqrt( prod_{k=1}^{m} (2k-1) (1-x^2) / (2k) )

        internal static double LegendrePe (int l, int m, double x) {

            Debug.Assert(l >= 0);
            Debug.Assert(0 <= m); Debug.Assert(m <= l);
            Debug.Assert(Math.Abs(x) <= 1.0);

            double xx = (1.0 + x) * (1.0 - x);
            // determine P{m,m}
            double P0 = 1.0;
            for (int k = 1; k <= m; k++) {
                P0 *= (1.0 - 1.0 / (2 * k)) * xx;
            }
            P0 = Math.Sqrt(P0);
            if (m % 2 != 0) P0 = -P0;
            if (l == m) return P0;
            // determine P{m+1,m}
            double s0 = Math.Sqrt(2*m+1);
            double P1 = x * s0 * P0;
            // iterate up to P{l,m}
            for (int k = m + 2; k <= l; k++) {
                double s2 = Math.Sqrt((k - m) * (k + m));
                double P2 = (x * (2 * k - 1) * P1 - s0 * P0) / s2;
                // prepare for next iteration
                s0 = s2;
                P0 = P1;
                P1 = P2;
            }
            return P1;

        }

        /// <summary>
        /// Computes the value of a Cebyshev polynomial of the frist kind.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="x">The argument, which must lie in the closed interval between -1 and +1.</param>
        /// <returns>The value of T<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>Chebyshev polynomials of the first kind are orthogonal on the interval [-1,1] with the weight (1-x<sup>2</sup>)<sup>-1/2</sup>.</para>
        /// <img src="../images/ChebyshevOrthonormality.png" />
        /// <para>Values returned are fully accurate (14-16 decimal digits) over the full range of argument for orders up to thousands.
        /// By orders up to a million, ~11 decimal digits remain accurate.</para>
        /// <para>For high orders, close to the endpoints, the oscillation of Chebyshev polynomials becomes extremely rapid. Particularly
        /// in these regions, keep in mind that the nearest representable floating point number might be far enough from the decimal number
        /// argument you want that the Chebyshev value shifts significantly. For example, you might claim that our value for T(1000000,0.9999)
        /// is accurate only to ~9 digits. But that is because 0.9999 is actually stored as 4503149267407759 X 2^(-52) = 0.99990000000000001101...,
        /// and the value returned is indeed accurate for that argument to ~11 digits.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative, or <paramref name="x"/> lies outside [-1,+1].</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Chebyshev_polynomials"/>
        /// <seealso href="http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html"/>
        public static double ChebyshevT (int n, double x) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
            if (Math.Abs(x) > 1.0) throw new ArgumentOutOfRangeException(nameof(x));

			if (n == 0) return 1.0;
            if (x == 1.0) return 1.0;
            if (x == -1.0) return (n % 2 == 0) ? 1.0 : -1.0;

            // very close to the endpoints, the recurrence looses accuracy for high n
            // use a series expansion there instead
            //if ((n > 10) && ((1.0 - Math.Abs(x) * n * n) < 1.0)) return ChebyshevT_Series1(n,x);
            //if ((n > 10) && (Math.Abs(n * x) < 0.1)) return (ChebyshevT_Series0(n, x));

            // Use the recurrence T_{n+1} = 2xT_{n} - T_{n-1}.
            // This is the same recurrence as for U_n, just with different starting values.
            // Amazingly, the recurence appears to work well for both: neither dominates over most
            // of parameter space.

            double two_x = 2.0 * x;

			double T0 = 1.0;
			double T1 = x;

			for (int k = 1; k < n; k++) {
				double T2 = two_x * T1 - T0;
				T0 = T1;
				T1 = T2;
			}

			return T1;
		}
		// Use: approximation of functions with minimum error

        /// <summary>
        /// Computes the value of a Cebyshev polynomial of the second kind.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="x">The argument, which must lie in the closed interval between -1 and +1.</param>
        /// <returns>The value of U<sub>n</sub>(x).</returns>
        /// <remarks>
        /// <para>Chebyshev polynomials of the second kind are orthogonal on the interval [-1,1] with the weight (1-x<sup>2</sup>)<sup>1/2</sup>.</para>
        /// <para>Values returned are fully accurate (14-16 decimal digits) over the full range of argument for orders up to thousands.
        /// For orders up to a million, ~11 decimal digits remain accurate.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative, or <paramref name="x"/> lies outside [-1,+1].</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Chebyshev_polynomials"/>
        /// <seealso href="https://mathworld.wolfram.com/ChebyshevPolynomialoftheSecondKind.html"/>
        public static double ChebyshevU(int n, double x) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
            if (Math.Abs(x) > 1.0) throw new ArgumentOutOfRangeException(nameof(x));

            if (n == 0) return 1.0;

            double two_x = 2.0 * x;

            double T0 = 1.0;
            double T1 = two_x;

            for (int k = 1; k < n; k++) {
                double T2 = two_x * T1 - T0;
                T0 = T1;
                T1 = T2;
            }

            return T1;
        }

        // orthogonal on [-1,1] with weight (1-x^2)^{1/2}
        // use recurrence U_{n+1} = 2xU_{n} - U_{n-1}
        /*
		public static double ChebyshevU (int n, double x) {
			if (x<-1.0 || x>1.0) throw new ArgumentOutOfRangeException("x",x,"The argument of a Chebyshew polynomial must be between -1 and 1.");
			if (n<0) throw new ArgumentOutOfRangeException("n", n, "The order of a Chebyshev polynomial must be non-negative.");
			if (n==0) return(1.0);
			double U0 = 1.0;
			double U1 = 2.0 * x;
			for (int k=1; k<n; k++) {
				double U2 = 2*x*U1 - U0;
				U0 = U1;
				U1 = U2;
			}
			return(U1);
		}
        */

        /// <summary>
        /// Computes the value of a Gegenbauer polynomial.
        /// </summary>
        /// <param name="n">The order, which must be non-negative.</param>
        /// <param name="alpha">The weight parameter, which must be non-negative.</param>
        /// <param name="x">The argument, which must lie in the closed interval [-1,+1].</param>
        /// <returns>The vaue of C<sup>(alpha)</sup><sub>n</sub>(x).</returns>
        /// <remarks><para>The Gegenbauer polynomials, also called the ultraspherical polynomials, are orthogonal on the interval [-1,+1] with the
        /// weight (1 - x^2)^(alpha - 1/2).</para>
        /// <para>The Legendre (<see cref="LegendreP(int, double)"/>) and Chebyshev (<see cref="ChebyshevT(int, double)"/> polynomials are
        /// thus special cases of the Gegenbauer polynomials.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative, or <paramref name="alpha"/> is negative, or
        /// <paramref name="x"/> lies outside [-1, +1].</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Gegenbauer_polynomials" />
        public static double GegenbauerC (int n, double alpha, double x) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
            if (alpha < 0.0) throw new ArgumentOutOfRangeException(nameof(alpha));
            if (Math.Abs(x) > 1.0) throw new ArgumentOutOfRangeException(nameof(x));

            if (n == 0) return 1.0;

            double two_alpha = 2.0 * alpha;
            double two_x = 2.0 * x;

            double C0 = 1.0;
            double C1 = two_alpha * x;

            for (int k = 1; k < n; k++) {
                double C2 = (two_x * (k + alpha) * C1 - (k - 1 + two_alpha) * C0) / (k + 1);
                C0 = C1;
                C1 = C2;
            }

            return C1;
        }

        // associated Legendre, Laguerre

        //  R00
        //      R11
        //  R20     R22
        //      R31     R33
        //  R40     R42     R44
        //      R51     R53     R55

        // 

        /// <summary>
        /// Computes the value of a Zernike polynomial.
        /// </summary>
        /// <param name="n">The order paramter, which must be non-negative.</param>
        /// <param name="m">The index parameter, which must lie between 0 and n inclusive.</param>
        /// <param name="rho">The argument, which must lie between 0 and 1.</param>
        /// <returns>The value of R<sub>n</sub><sup>m</sup>(&#x3C1;).</returns>
        /// <remarks>
        /// <para>Zernike polynomials are orthononal on the interval [0,1] with the weight &#x3C1;.</para>
        /// <para>They are often used in optics to characterize the imperfections in a lens. In
        /// this context, the amplitude of each is associated with a name given in the following table.</para>
        /// <table>
        ///     <tr><th>n</th><th>m</th><th>name</th></tr>
        ///     <tr><td>1</td><td>1</td><td>tilt</td></tr>
        ///     <tr><td>2</td><td>0</td><td>defocus</td></tr>
        ///     <tr><td>2</td><td>2</td><td>astigmatism</td></tr>
        ///     <tr><td>3</td><td>1</td><td>coma</td></tr>
        ///     <tr><td>3</td><td>3</td><td>trefoil</td></tr>
        /// </table>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative, or <paramref name="m"/>
        /// lies outside {0 .. n}, or <paramref name="rho"/> lies outside [0, 1].</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Zernike_polynomials"/>
        /// <seealso href="https://mathworld.wolfram.com/ZernikePolynomial.html"/>
        public static double ZernikeR (int n, int m, double rho) {

            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
            if ((m < 0) || (m > n)) throw new ArgumentOutOfRangeException(nameof(m));
            if ((rho < 0.0) || (rho > 1.0)) throw new ArgumentOutOfRangeException(nameof(rho));

            // n and m have the same parity
            if ((n - m) % 2 != 0) return 0.0;

            // R00
            if (n == 0) return 1.0; 

            // R^{m}_m
            double r2 = Math.Pow(rho, m);
            if (n == m) return r2;

            // R^{m+1}_{m+1}
            int k = m;
            double r1 = r2 * rho;

            while (true) {

                k += 2;

                // *
                //  \
                //   * recurrence involving two lesser m's
                //  /
                // *
                // 2n R^{m+1}_{n-1} = (n+m) R^{m}_{n-2} + (n-m) R^{m}_{n}

                double r0 = ((2 * k) * rho * r1 - (k + m) * r2) / (k - m);

                if (k == n) return r0;

                //   *
                //  /
                // * recurrence involving two greater m's
                //  \
                //   *
                // 

                double rp = (2 * (k + 1) * rho * r0 - (k - m) * r1) / (k + m + 2);

                r2 = r0;
                r1 = rp;

            }

        }

    }

}
