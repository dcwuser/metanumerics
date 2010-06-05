using System;
using System.Collections.Generic;
using Meta.Numerics.Statistics;

namespace Meta.Numerics.Functions {

    public static partial class FunctionMath {

        public static double Differentiate (Function<double, double> f, double x) {

            // choose an initial step size
            double h = Math.Abs(x) / 16.0 + Math.Pow(2.0, -16);
            double t = x + h;
            h = t - x;

            // create a tableau
            int max = 12;
            double[][] D = new double[max][];

            // fill out the tableau
            for (int j = 0; j < max; j++) {
                // create new row in the tableau
                D[j] = new double[j + 1];

                Console.WriteLine(j);

                // add our next evaluation
                double fp = f(x + h);
                double fm = f(x - h);
                D[j][0] = (fp - fm) / (2.0 * h);
                Console.WriteLine(D[j][0]);

                // fill out the row by extrapolation
                double e = 1.0;
                for (int k = 1; k <= j; k++) {
                    e = e / 2.0;
                    //e = e / 4.0;
                    D[j][k] = (D[j][k - 1] - e * D[j - 1][k - 1]) / (1 - e);
                }

                Console.WriteLine(D[j][j]);

                // check for convergence
                if ((j > 0) && (Math.Abs(D[j][j] - D[j-1][j-1]) <= Math.Pow(2.0,-42) * Math.Abs(D[j][j]))) {
                    return (D[j][j]);
                }

                // halve the step-size for the next iteration
                h = h / Math.Sqrt(2.0);
                //h = h / 2.0;
            }

            throw new NonconvergenceException();

        }

    }

}