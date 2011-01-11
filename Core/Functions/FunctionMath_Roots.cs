using System;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Functions {

    public static partial class FunctionMath {

        /// <summary>
        /// Isolates a root in the vicinity of a given point.
        /// </summary>
        /// <param name="f">The function whoose zero is sought.</param>
        /// <param name="x">A ordinate believed to be near the sought zero.</param>
        /// <returns>An ordinate at which the function has a zero.</returns>
        public static double FindZero (Func<double, double> f, double x) {

            if (f == null) throw new ArgumentNullException("f");

            // take a step
            double fx = f(x);
            double dx = 0.1 * Math.Abs(x) + 0.01;

            // keep stepping until we change sign
            for (int n = 0; n < 100; n++) {

                double y = x + dx;
                double fy = f(y);

                // if we the function changed sign, we have bracketed a root
                if (Math.Sign(fy) != Math.Sign(fx)) {
                    return (FindZero(f, x, fx, y, fy));
                }

                // the the function got bigger, reverse direction
                if (Math.Abs(fy) > Math.Abs(fx)) dx = -dx;

                // prepare for the next loop
                dx = 2.0 * dx;
                x = y;
                fx = fy;

            }

            throw new NonconvergenceException();

        }

        /// <summary>
        /// Isolates a root within a given interval.
        /// </summary>
        /// <param name="f">The function whoose zero is sought.</param>
        /// <param name="bracket">An interval bracketing the root.</param>
        /// <returns>An ordinate within the bracket at which the function has a zero.</returns>
        /// <exception cref="InvalidOperationException">The function does not change sign across the given interval.</exception>
        public static double FindZero (Func<double, double> f, Interval bracket) {

            if (f == null) throw new ArgumentNullException("f");

            double x1 = bracket.LeftEndpoint;
            double x2 = bracket.RightEndpoint;
            double f1 = f(x1);
            double f2 = f(x2);

            // make sure the bracket points really do bracket a root
            if (Math.Sign(f1) == Math.Sign(f2)) throw new InvalidOperationException();

            return (FindZero(f, x1, f1, x2, f2));
        }

        private static double FindZero (Func<double, double> f, double x1, double f1, double x2, double f2) {

            double x0 = Double.MaxValue;

            for (int n = 0; n < Global.SeriesMax; n++) {

                //Debug.WriteLine(String.Format("Bracket ({0},{1}) ({2},{3})", x1, f1, x2, f2));
                //Console.WriteLine("Bracket ({0},{1}) ({2},{3})", x1, f1, x2, f2);

                // evaluate at the bracket mid-point
                double x3 = (x1 + x2) / 2.0;
                double f3 = f(x3);
                //Debug.WriteLine(String.Format("Midpoint f({0}) = {1}", x3, f3));
                if (f3 == 0.0) return (x3);

                // find the Ridder point...
                double x4;
                double q = Math.Sqrt(f3 * f3 - f1 * f2);
                if (q == 0.0) return (x3);
                double dx = (x3 - x1) * f3 / q;
                if (f1 > f2) {
                    x4 = x3 + dx;
                } else {
                    x4 = x3 - dx;
                }

                // ...and evaluate there
                double f4 = f(x4);
                //Debug.WriteLine(String.Format("Ridder f({0}) = {1}", x4, f4));
                if (f4 == 0.0) return (x4);

                if (Math.Abs(x4 - x0) <= RelativePrecision * Math.Abs(x0)) return (x4);
                if (Math.Abs(x4 - x0) <= AbsolutePrecision) return (x4);

                // update the bracket
                if (Math.Sign(f3) == Math.Sign(f4)) {
                    // x3 and x4 are on the same side of the root
                    // replace the same-signed bracket point with x4
                    if (Math.Sign(f4) == Math.Sign(f1)) {
                        x1 = x4;
                        f1 = f4;
                    } else {
                        x2 = x4;
                        f2 = f4;
                    }
                } else {
                    // x3 and x4 are on oppisite sides of the root
                    // use them to replace the bracket
                    x1 = x3;
                    f1 = f3;
                    x2 = x4;
                    f2 = f4;
                }

                // add end criteria to see whether we have moved from the last iteration
                x0 = x4;

            }

            throw new NonconvergenceException();

        }

        // multidimensional

#if FUTURE
        public static void FindZero (Function<double[], double[]> f, double[] x0) {

            if (f == null) throw new ArgumentNullException("f");
            if (x0 == null) throw new ArgumentNullException("x0");

            int d = x0.Length;
            SquareMatrix B = new SquareMatrix(d);
            ColumnVector x = new ColumnVector(x0);

            for (int n = 0; n < Global.SeriesMax; n++) {

                // determine
                ColumnVector F1 = new ColumnVector(f(x.ToArray()));

                // determine Newton step
                SquareLUDecomposition LU = B.LUDecomposition();
                ColumnVector dx = -LU.Solve(F1);

                // take newton step
                x = x + dx;

                // update B

            }

            throw new NonconvergenceException();

        }
#endif

    }

}
