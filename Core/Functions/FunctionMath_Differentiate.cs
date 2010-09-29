using System;
using System.Collections.Generic;
using Meta.Numerics.Statistics;

namespace Meta.Numerics.Functions {

    public static partial class FunctionMath {

#if FUTURE

        /// <summary>
        /// Estimates the derivative of the given function at the given point.
        /// </summary>
        /// <param name="f">The function to differentiate.</param>
        /// <param name="x">The point at which to evaluate the derivative.</param>
        /// <returns>The estimate, with uncertainty, of the value of the derivative.</returns>
        public static UncertainValue Differentiate (Function<double, double> f, double x) {

            if (f == null) throw new ArgumentNullException("f");

            double step = Math.Pow(2, -4) * Math.Abs(x) + Math.Pow(2, -8);
            double factor = Math.Pow(2.0, 1.0 / 2.0);
            UncertainValue v = Differentiate(f, x, step, factor);
            return (v);

        }

        private static UncertainValue Differentiate (Function<double, double> f, double x, double h, double s) {

            // choose an initial step size
            //double h = Math.Abs(x) / 16.0 + Math.Pow(2.0, -16);
            //double t = x + h;
            //h = t - x;

            // keep track of the best value
            UncertainValue best = new UncertainValue(0.0, Double.MaxValue);

            // create a tableau
            int max = 12;
            double[][] D = new double[max][];

            // fill out the tableau
            for (int j = 0; j < max; j++) {
                // create new row in the tableau
                D[j] = new double[j + 1];

                //Console.WriteLine(j);

                // add our next evaluation
                double fp = f(x + h);
                double fm = f(x - h);
                D[j][0] = (fp - fm) / (2.0 * h);
                //Console.WriteLine(D[j][0]);

                // fill out the row by extrapolation
                double e = 1.0;
                for (int k = 1; k <= j; k++) {
                    e = e / (s*s);
                    D[j][k] = (D[j][k - 1] - e * D[j - 1][k - 1]) / (1 - e);
                }

                //Console.WriteLine(D[j][j]);
                if (j > 0) {
                    double u = Math.Max(
                        Math.Abs(D[j][j] - D[j - 1][j - 1]),
                        Math.Abs(D[j][j] - D[j][j-1])
                     );

                    // check for 
                    if (u < best.Uncertainty) {
                        best = new UncertainValue(D[j][j], u);
                    } else if (u > 2.0 * best.Uncertainty) {
                        return (best);
                    }
                }

                // reduce the step-size for the next iteration
                h = h / s;
            }

            return (best);
            //throw new NonconvergenceException();

        }

#endif

    }

}