using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Functions {

    public partial class AdvancedMath {

        /// <summary>
        /// Computes the Riemann zeta function.
        /// </summary>
        /// <param name="x">The argument.</param>
        /// <returns>The value &#x3B6;(s).</returns>
        /// <remarks>
        /// <para>The Riemann &#x3B6; function can be defined as the sum of the <paramref name="x"/>th inverse power of the natural numbers.</para>
        /// <img src="../images/ZetaSeries.png" />
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Riemann_zeta_function"/>
        public static double RiemannZeta (double x) {
            if (x < 0.0) {
                // for negative numbers, use the reflection formula
                double t = 1.0 - x;
                return (2.0 * Math.Pow(Global.TwoPI, -t) * MoreMath.Cos(Global.HalfPI * t) * AdvancedMath.Gamma(t) * RiemannZeta(t));
            } else {
                double xm1 = x - 1.0;
                if (Math.Abs(xm1) < 0.25) {
                    // near the singularity, use the Stjielts expansion
                    return (RiemannZeta_LaurentSeries(xm1));
                } else {
                    // call Dirichlet function, which converges faster
                    return (DirichletEta(x) / (1.0 - Math.Pow(2.0, 1.0 - x)));
                }
            }
        }

        /// <summary>
        /// Computes the Dirichlet eta function.
        /// </summary>
        /// <param name="x">The argument, which must be non-negative.</param>
        /// <returns>The value of &#x3B7;(s).</returns>
        /// <remarks>
        /// <para>The Dirichlet eta function is the sum of the <paramref name="x"/>th inverse power of the natural numbers,
        /// with alternating signs.</para>
        /// <img src="../images/DirichletEtaSeries.png" />
        /// <para>Because these are just the terms of the Riemann zeta function (<see cref="RiemannZeta"/>) with
        /// alternating signs, it is also called the alternating zeta function.</para>
        /// <para>It can be related to the Riemann &#x3B6; function.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative.</exception>
        /// <seealso cref="RiemannZeta"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Dirichlet_eta_function"/>
        public static double DirichletEta (double x) {
            if (x < 0.0) throw new ArgumentOutOfRangeException(nameof(x));
            return (DirichletEta_Borwein(x));
        }

        // Borwein's amazing method for computing eta is detailed at http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html.
        // The original paper is http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf.

        // An amazing, fast, fixed-length sequence approximation to eta that is good to 1/8^(DirichletEta_Coefficients.Length) for all s > 0.

        private static double DirichletEta_Borwein (double x) {
            double sum1 = 0.0;
            bool sign = true;
            for (int k = 0; k < DirichletEta_BorweinCoefficients.Length; k++) {
                double term = Math.Pow(k + 1, -x);
                if (!sign) term = -term;
                sum1 += term;
                sign = !sign;
            }
            double sum2 = 0.0;
            for (int k = 0; k < DirichletEta_BorweinCoefficients.Length; k++) {
                double term = DirichletEta_BorweinCoefficients[k] * Math.Pow(k + DirichletEta_BorweinCoefficients.Length + 1, -x);
                if (!sign) term = -term;
                sum2 += term;
                sign = !sign;
            }
            return (sum1 + sum2);
        }

        private static readonly double[] DirichletEta_BorweinCoefficients = ComputeBorweinEtaCoefficients(16);

        // The Borwein coefficients e_k = \sum_{j=k}^{n} \binom{n, j}, i.e.
        //   e_n = \binom{n, n}
        //   e_{n-1} = \binom{n, n} + \binom{n, n-1 }
        // etc.

		private static double[] ComputeBorweinEtaCoefficients (int n) {
            double norm = MoreMath.Pow(2.0, n);
			double[] e = new double[n];
			double sum = 0.0;

            IEnumerator<double> binomials = AdvancedIntegerMath.BinomialCoefficients(n).GetEnumerator();
            for (int k = n - 1; k >= 0; k--) {
                binomials.MoveNext();
                sum += binomials.Current;
                e[k] = sum / norm;
            }

			return(e);
		}

        // The Laurent expansion of zeta near s = 1 is
        //   \zeta(s) = 1/(s-1) + \sum_{k=0}^{\infty} (-1)^k \gamma_k (s-1)^k / k!
        // where \gamma_k are Stieltjes constants
        // Note that the argument of this function is x = s - 1, not s

        private static double RiemannZeta_LaurentSeries (double x) {
            double df = 1.0;
            double f = 1.0 / x + StieltjesConstants[0];
            for (int i = 1; i < StieltjesConstants.Length; i++) {
                double f_old = f;
                df = - df * x / i;
                f += StieltjesConstants[i] * df;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        // Here are the first 16 Stieltjes constants, from http://pi.lacim.uqam.ca/piDATA/stieltjesgamma.txt
        // Since the last term in the Laurent expansion of zeta goes like \gamma_n (s-1)^n / n!, this should
        // be enough to allow us to use the expansion up to (s-1) ~ 1

        internal static readonly double[] StieltjesConstants = new double[] {
            0.577215664901532860607,
           -0.072815845483676724861,
           -0.009690363192872318484,
            0.002053834420303345866,
            0.002325370065467300058,
            0.000793323817301062702,
           -0.000238769345430199610,
           -0.000527289567057751046,
           -0.000352123353803039510,
           -0.000034394774418088048,
            0.000205332814909064795,
            0.000270184439543903527,
            0.000167272912105140193,
           -0.000027463806603760159,
           -0.000209209262059299946,
           -0.000283468655320241447
        };

    }

    public static partial class AdvancedComplexMath {

        /// <summary>
        /// Computes the Riemann zeta function for complex values.
        /// </summary>
        /// <param name="z">The argument.</param>
        /// <returns>The value of &#x3B6;(z).</returns>
        /// <remarks>
        /// <para>As the imaginary part of the argument increases, the computation of the zeta function becomes slower and more difficult.
        /// The computation time is approximately proportional to the imaginary part of z. The result also slowly looses accuracy for arguments with
        /// very large imaginary parts; for arguments with z.Im of order 10^d, approximately the last d digits of the result are suspect.</para>
        /// <para>The image below shows the complex &#x393; function near the origin using domain coloring. You can see the first non-trivial
        /// zeros at (1/2, &#177;14.13...) as well as the trivial zeros along the negative real axis.</para>
        /// <img src="../images/ComplexRiemannZetaPlot.png" />
        /// </remarks>
        public static Complex RiemannZeta (Complex z) {

            // Use conjugation and reflection symmetry to move to the first quadrant.
            if (z.Im < 0.0) return (RiemannZeta(z.Conjugate).Conjugate);
            if (z.Re < 0.0) {
                Complex zp = Complex.One - z;
                return (2.0 * ComplexMath.Pow(Global.TwoPI, -zp) * ComplexMath.Cos(Global.HalfPI * zp) * AdvancedComplexMath.Gamma(zp) * RiemannZeta(zp)); 
            }

            // Close to pole, use Laurent series.
            Complex zm1 = z - Complex.One;
            if (ComplexMath.Abs(zm1) < 0.50) {
                return (RiemannZeta_LaurentSeries(zm1));
            }

            // Fall back to Euler-Maclaurin summation.
            int n = RiemannZeta_EulerMaclaurin_N(z.Re, z.Im);
            return (RiemannZeta_EulerMaclaurin(z, n));

        }

        private static Complex RiemannZeta_LaurentSeries (Complex z) {
            Complex df = 1.0;
            Complex f = 1.0 / z + AdvancedMath.StieltjesConstants[0];
            for (int i = 1; i < AdvancedMath.StieltjesConstants.Length; i++) {
                Complex f_old = f;
                df = -df * z / i;
                f += AdvancedMath.StieltjesConstants[i] * df;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        // An important decision when using Euler-Maclaurin summation is how far to sum the series explicitly before
        // applying the Euler-Maclaurin approximation to the remaining terms. If you have an arbitrary number of
        // pre-computed Bernoulli numbers, this question is quite subtle, since there are trade-offs between the number
        // of explicit and Bernoulli terms. For us the question is simpler because we only have a few (32 as of this comment)
        // pre-computed Bernoulli numbers and computing higher ones on the fly is quite expensive. Therefore we want to
        // be sure we converge within a fixed number of Bernoulli terms (say 16 to be safe). We approximate the error
        // by the last Bernoulli term
        //    \frac{B_j (z)_j}{j! n^(z+j)}
        // Using \frac{B_j}{j!} \frac{(2\pi)^j}{2} \rightarrow 1 as j increases and n^z ~ n^s and (z)_j ~ (s + t + j)^j,
        // this becomes
        //    \frac{2}{n^s} \left( \frac{s + t + j}{2 \pi n} \right)^j
        // setting this to \epsilon and solving for n gives
        //    n^{s + j} ~ \frac{2}{\epsilon} \left( \frac{s + t + j}{2 \pi} \right)^j
        //    (s + j) \ln n ~ \ln(2/\epsilon) + j \ln \left ( \frac{s + t + j}{2 \pi} \right)
        // We apply this to obtain n. At a couple logs and exp it's costly to just obtain an integer, but every
        // term is already a log and a couple of trig functions, so it's cost is still trivial compared to the
        // Euler-Maclaurin algorithm that follows.

        private static int RiemannZeta_EulerMaclaurin_N (double s, double t) {
            Debug.Assert(s >= 0.0); Debug.Assert(t >= 0.0);

            // target convergence half-way through our tabulated Bernoulii numbers
            const int j = 16;

            double n = Math.Exp((baseN + j * Math.Log((s + t + j) / Global.TwoPI)) / (s + j));
            if (n > Int32.MaxValue) {
                throw new NonconvergenceException();
            } else {
                // add 4 to ensure a lower bound
                return (((int) Math.Ceiling(n)) + 4);
            }
        }

        private static readonly double baseN = Math.Log(2.0 / Global.Accuracy);

        // Euler-Maclaurin summation approximates a sum by an integral.
        //   \sum_{k=a}^{b} f(k) = \int_a^b \! f(x) \, dx  - B_1 \left[ f(b) + f(a) \right] + 
        //     \sum_{j=0}^{\infty} \frac{B_{2j}}{(2j)!} \left[ f^{(2j-1)}(b) - f^{(2j-1)}(a) \right]

        private static Complex RiemannZeta_EulerMaclaurin (Complex z, int n) {

            Complex f = 0.0;

            // Sum terms up to (n-1)^{-z}
            for (int k = 1; k < n; k++) {
                f += ComplexMath.Pow(k, -z);
            }

            // Integral term
            Complex f_old = f;
            Complex zm1 = z - Complex.One;
            f += ComplexMath.Pow(n, -zm1) / zm1;
            if (f == f_old) return (f);

            // B_1 term
            f_old = f;
            Complex df = ComplexMath.Pow(n, -z) / 2.0;
            f += df;
            if (f == f_old) return (f);

            // B_2 term
            df *= z / n;
            f += AdvancedIntegerMath.BernoulliNumber(2) * df;
            if (f == f_old) return (f);

            // Higher Bernoulli terms
            for (int k = 4; k < 2 * AdvancedIntegerMath.Bernoulli.Length; k += 2) {
                f_old = f;
                df *= (z + (k - 2)) / n * (z + (k - 3)) / n / k / (k - 1);
                f += AdvancedIntegerMath.Bernoulli[k / 2] * df;
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();

        }

    }

}
