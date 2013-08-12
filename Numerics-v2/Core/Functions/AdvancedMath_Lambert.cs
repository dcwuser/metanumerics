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
        /// <seealso href="http://en.wikipedia.org/wiki/Lambert_W_function" />
        /// <seealso href="http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf"/>
        public static double LambertW (double x) {

            if (x < -EI) throw new ArgumentOutOfRangeException("x");

            // use an initial approximation
            double W;
            if (x < -EI / 2.0) {
                W = Lambert_SeriesSmall(x);
                // don't try Halley extremely close to the edge; it will blow up due to the large derivatives
                if ((x + EI) < 1.0E-6) return (W);
            } else if (x < EI) {
                W = Lambert_SeriesZero(x);
            } else if (x > Math.E) {
                W = Lambert_SeriesLarge(x);
            } else {
                // it would be nice to have something better in the transition region, but
                // this converges after just a few iterations, so it's hard to justify doing much more
                W = 0.5;
            }

            // use Halley's method (Newton's method plus second derivative) to hone in on solution
            W = Lambert_Halley(x, W);

            return (W);

        }

        private const double EI = 1.0 / Math.E;

        private static double Lambert_Halley (double x, double w0) {

            for (int i = 0; i < Global.SeriesMax; i++) {
                double e = Math.Exp(w0);
                double f = e * w0 - x;
                double dw = f / ((w0 + 1.0) * e - ((w0 + 2.0) / (w0 + 1.0)) * f / 2.0);
                double w1 = w0 - dw;
                if (w1 == w0) {
                    return (w1);
                }
                w0 = w1;
            }

            throw new NonconvergenceException();

        }


        // series useful near 0
        private static double Lambert_SeriesZero (double x) {
            double W = x - x * x + (3.0 / 2.0) * x * x * x - (8.0 / 3.0) * x * x * x * x;
            return (W);
        }

        // series useful near -1/e
        private static double Lambert_SeriesSmall (double x) {
            double p = Math.Sqrt(2.0 * (Math.E * x + 1.0));
            double W = -1.0 + p - p * p / 3.0 + (11.0 / 72.0) * p * p * p - (43.0 / 540.0) * p * p * p * p;
            return (W);
        }

        // series useful for large x
        private static double Lambert_SeriesLarge (double x) {
            double L1 = Math.Log(x);
            double L2 = Math.Log(L1);
            double W = L1 - L2 + L2 / L1 + L2 * (L2 - 2.0) / L1 / L1 / 2.0;
            return (W);
        }


    }

}
