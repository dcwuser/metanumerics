using System;


namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        /// <summary>
        /// Computes the Lambert W function.
        /// </summary>
        /// <param name="x">The argument, which must be greater than or equal to -1/e.</param>
        /// <returns>The value W(x).</returns>
        /// <remarks>
        /// <para>The Lambert W function solves the transcendental equation W e<sup>W</sup> = x.
        /// The function appears in a number of contexts, including the solution of differential
        /// equations and the enumeration of trees.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is less than -1/e.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Lambert_W_function" />
        /// <seealso href="https://mathworld.wolfram.com/LambertW-Function.html"/>
        /// <seealso href="https://dlmf.nist.gov/4.13"/>
        public static double LambertW (double x) {

            if (x < -EI) throw new ArgumentOutOfRangeException(nameof(x));

            // use an initial approximation
            double W;
            if (x < -EI / 2.0) {
                // We need to be careful not to move on to Halley's method too close to the branch point,
                // because Halley's method has (W-1) in the denominator and therefore blows up close to
                // the branch point. Fortunately, the branch point series converges quickly close to the
                // branch point, and if it converges, we don't need to move on to Halley's method.
                if (Lambert_SeriesBranch(x, out W)) return (W);
            } else if (x < EI) {
                W = Lambert_SeriesZero(x);
            } else if (x > Math.E) {
                W = Lambert_SeriesLarge(x);
            } else {
                // It would be nice to have something better in the transition region, but
                // this converges after just a few iterations, so it's hard to justify doing much more.
                W = 0.5;
            }

            // Use Halley's method (Newton's method plus second derivative) to hone in on solution.
            W = Lambert_Halley(x, W);

            return (W);

        }

        private const double EI = 1.0 / Math.E;

        private static double Lambert_Halley (double x, double W) {

            for (int i = 0; i < Global.SeriesMax; i++) {
                double e = Math.Exp(W);
                double f = e * W - x;
                double dW = f / ((W + 1.0) * e - ((W + 2.0) / (W + 1.0)) * f / 2.0);
                W = W - dW; 
                if (Math.Abs(dW) <= Global.Accuracy * Math.Abs(W)) return (W);
            }

            throw new NonconvergenceException();

        }

        // series useful near 0
        private static double Lambert_SeriesZero (double x) {
            double W = x * (1.0 - x + (3.0 / 2.0) * x * x - (8.0 / 3.0) * x * x * x + (125.0 / 24.0) * x * x * x * x);
            return (W);
        }


        // This is a series development around the branch point (i.e. x = -1/e) that is documented in DLMF.
        // The series is W(-e^{-(1 + t^2 / 2)}) = -\sum_{k=0}^{\infty} c_k (-t)^k
        // The coefficients are c_0 = 1, c_1 = 1, c_k = \frac{1}{k+1} \left( c_{k-1} \sum_{j=2}^{k-1} j c_j c_{k+1-j} \right)

        private static bool Lambert_SeriesBranch (double x, out double W) {

            double t = Math.Sqrt(-2.0 * (1.0 + Math.Log(-x)));
            W = -1.0 + t;
            double tk = t;
            for (int k = 2; k < Lambert_SeriesBranch_Coefficients.Length; k++) {
                double W_old = W;
                tk *= -t;
                W += Lambert_SeriesBranch_Coefficients[k] * tk;
                if (W == W_old) {
                    return (true);
                }
            }
            return (false);

        }

        private static readonly double[] Lambert_SeriesBranch_Coefficients = ComputeLambertBranchCoefficients(8);

        private static double[] ComputeLambertBranchCoefficients (int n) {
            double[] c = new double[n];
            c[0] = 1.0;
            c[1] = 1.0;
            for (int k = 2; k < c.Length; k++) {
                double t = 0.0;
                for (int j = 2; j < k; j++) {
                    t += j * c[j] * c[k + 1 - j];
                }
                c[k] = (c[k - 1] - t) / (k + 1);
            }
            return (c);
        }

        // series useful for large x
        private static double Lambert_SeriesLarge (double x) {
            double L1 = Math.Log(x);
            double L2 = Math.Log(L1);
            double W = L1 - L2 + L2 / L1 + L2 * (L2 - 2.0) / (2.0 * L1 * L1) + L2 * (2.0 * L2 * L2 - 9.0 * L2 + 6.0) / (6.0 * L1 * L1 * L1);
            return (W);
        }


    }

}
