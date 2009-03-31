using System;

namespace Meta.Numerics.Functions {

    public partial class AdvancedMath {

        /// <summary>
        /// Compute the Riemann zeta function.
        /// </summary>
        /// <param name="s">The argument.</param>
        /// <returns>The value &#x3B6;(s).</returns>
        /// <remarks>
        /// <para>The Riemann &#x3B6; function can be defined as the sum of the <paramref name="s"/>th inverse power of the natural numbers.</para>
        /// </remarks>
        public static double RiemannZeta (double s) {
            if (s < 0.0) {
                // for negative numbers, use the reflection formula
                double t = 1.0 - s;
                double z = 2.0 * Math.Pow(2.0 * Math.PI, -t) * Math.Cos(t * Math.PI / 2.0) * AdvancedMath.Gamma(t) * RiemannZeta(t);
                return (z);
            } else {
                if (Math.Abs(s - 1.0) < 0.25) {
                    // near the sigularity, use the Stjielts expansion
                    return (RiemannZeta_Series(s - 1.0));
                } else {
                    // call Dirichlet function, which converges faster
                    return (DirichletEta(s) / (1.0 - Math.Pow(2.0, 1.0 - s)));
                }
            }
        }

        /// <summary>
        /// Computes the Dirichlet eta function.
        /// </summary>
        /// <param name="s">The argument, which must be non-negative.</param>
        /// <returns>The value of &#x3B7;(s).</returns>
        /// <remarks>
        /// <para>The Dirichlet eta function is the sum of the <paramref name="s"/>th inverse power of the natural numbers,
        /// with alternating signs. It can be related to the Riemann &#x3B6; function.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="s"/> is negative.</exception>
        /// <seealso cref="RiemannZeta"/>
        public static double DirichletEta (double s) {
            if (s < 0.0) throw new ArgumentOutOfRangeException("s");
            return (DirichletEta_Sequence(s));
        }

        // an amazing, fast, fixed-length sequence approximation to eta that is good to 1/8^(DirichletEta_Coefficients.Length) for all s > 0
        private static double DirichletEta_Sequence (double s) {
            double sum1 = 0.0;
            bool sign = true;
            for (int k = 0; k < DirichletEta_Coefficients.Length; k++) {
                double term = Math.Pow(k + 1, -s);
                if (!sign) term = -term;
                sum1 += term;
                sign = !sign;
            }
            double sum2 = 0.0;
            for (int k = 0; k < DirichletEta_Coefficients.Length; k++) {
                double term = DirichletEta_Coefficients[k] * Math.Pow(k + DirichletEta_Coefficients.Length + 1, -s);
                if (!sign) term = -term;
                sum2 += term;
                sign = !sign;
            }
            return (sum1 + sum2 / Math.Pow(2.0, DirichletEta_Coefficients.Length));
        }

        private static double[] DirichletEta_Coefficients = ReimannCoefficients(15);

		private static double[] ReimannCoefficients (int n) {
			double[] e = new double[n];
			double sum = 0.0;
			for (int k = n; k>0; k--) {
				sum += (double) AdvancedIntegerMath.BinomialCoefficient(n, k);
				e[k-1] = sum;
			}
			return(e);
		}

        // an expansion near s=1 using Stieltjes constants; here x = s-1
        private static double RiemannZeta_Series (double x) {
            double dz = 1.0;
            double z = 1.0 / x + StieltjesConstants[0];
            for (int i = 1; i < StieltjesConstants.Length; i++) {
                double z_old = z;
                dz = - dz * x / i;
                z = z_old + StieltjesConstants[i] * dz;
                if (z == z_old) return (z);
            }
            throw new NonconvergenceException();
        }

        private static double[] StieltjesConstants = new double[] {
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
           -0.000027463806603760159
        };

    }

}
