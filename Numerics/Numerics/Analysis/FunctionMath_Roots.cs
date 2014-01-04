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
            if (fx == 0.0) return (x);
            double dx = (Math.Abs(x) + 1.0 / 16.0) / 16.0;
            //double dx = 0.1 * Math.Abs(x) + 0.01;

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

            // x1 and x2 are bracket points, with function values f1 and f2

            // x0 is the current best estimate of the root; it will be updated in the first loop
            double x0 = Double.MaxValue;

            for (int n = 0; n < Global.SeriesMax; n++) {

                //Debug.WriteLine(String.Format("Bracket ({0},{1}) ({2},{3})", x1, f1, x2, f2));

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

                // stopping criteria
                if (f4 == 0.0) return (x4);
                if (Math.Abs(x4 - x0) <= Global.Accuracy * Math.Abs(x0)) return (x4);
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

        /// <summary>
        /// Finds a vector argument which makes a vector function zero.
        /// </summary>
        /// <param name="f">The vector function.</param>
        /// <param name="x0">The vector argument.</param>
        /// <returns>The vector argument which makes all components of the vector function zero.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> or <paramref name="x0"/> is null.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="f"/> is not equal to the
        /// dimension of <paramref name="x0"/>.</exception>
        public static double[] FindZero (Func<double[], double[]> f, double[] x0) {

            if (f == null) throw new ArgumentNullException("f");
            if (x0 == null) throw new ArgumentNullException("x0");

            // we will use Broyden's method, which is a generalization of the secant method to the multi-dimensional problem
            // just as the secant method is essentially Newton's method with a crude, numerical value for the slope,
            // Broyden's method is a multi-dimensional Newton's method with a crude, numerical value for the Jacobian matrix

            int d = x0.Length;

            // we should re-engineer this code to work directly on vector/matrix storage rather than go back and forth
            // between vectors and arrays, but for the moment this works; implementing Blas2 rank-1 update would help

            // starting values
            ColumnVector x = new ColumnVector(x0);
            //double[] x = new double[d]; Blas1.dCopy(x0, 0, 1, x, 0, 1, d);
            ColumnVector F = new ColumnVector(f(x.ToArray()));
            //double[] F = f(x);
            if (F.Dimension != d) throw new DimensionMismatchException();
            SquareMatrix B = ApproximateJacobian(f, x0);
            double g = MoreMath.Pow2(F.Norm());
            //double g = Blas1.dDot(F, 0, 1, F, 0, 1, d);

            for (int n = 0; n < Global.SeriesMax; n++) {

                // determine the Newton step
                LUDecomposition LU = B.LUDecomposition();
                ColumnVector dx = -LU.Solve(F);

                //for (int i = 0; i < d; i++) {
                //    Console.WriteLine("F[{0}]={1} x[{0}]={2} dx[{0}]={3}", i, F[i], x[i], dx[i]);
                //}

                // determine how far we will move along the Newton step
                ColumnVector x1 = x + dx;
                ColumnVector F1 = new ColumnVector(f(x1.ToArray()));
                double g1 = MoreMath.Pow2(F1.Norm());

                // check whether the Newton step decreases the function vector magnitude
                // NR suggest that it's necessary to ensure that it decrease by a certain amount, but I have yet to see that make a diference
                double gm = g - 0.0;
                if (g1 > gm) {
                    // the Newton step did not reduce the function vector magnitude, so we won't step that far
                    // determine how far along the descent direction we will step by parabolic interpolation
                    double z = g / (g + g1);
                    //Console.WriteLine("z={0}", z);
                    // take at least a small step in the descent direction
                    if (z < (1.0 / 16.0)) z = 1.0 / 16.0;
                    dx = z * dx;
                    x1 = x + dx;
                    F1 = new ColumnVector(f(x1.ToArray()));
                    g1 = MoreMath.Pow2(F1.Norm());
                    // NR suggest that this be repeated, but I have yet to see that make a difference
                }

                // take the step
                x = x + dx;

                // check for convergence
                if ((F1.InfinityNorm() < Global.Accuracy) || (dx.InfinityNorm() < Global.Accuracy)) {
                //if ((g1 < Global.Accuracy) || (dx.Norm() < Global.Accuracy)) {
                    return (x.ToArray());
                }

                // update B
                ColumnVector dF = F1 - F;
                RectangularMatrix dB = (dF - B * dx) * (dx / (dx.Transpose() * dx)).Transpose();
                //RectangularMatrix dB = F1 * ( dx / MoreMath.Pow(dx.Norm(), 2) ).Transpose();
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < d; j++) {
                            B[i,j] += dB[i, j];
                    }
                }
                
                // prepare for the next iteration
                F = F1;
                g = g1;

            }

            throw new NonconvergenceException();

        }

        // we just need a very approximate Jacobian to start off the Broyden's root-finder,
        // so we don't employ all the right tricks to maximize accuracy here

        private static SquareMatrix ApproximateJacobian (Func<double[], double[]> f, double[] x0) {

            int d = x0.Length;

            double[] f0 = f(x0);

            double[] xp = new double[d];

            SquareMatrix B = new SquareMatrix(d);
            for (int j = 0; j < d; j++) {
                Array.Copy(x0, xp, d);
                double dx = (Math.Abs(x0[j]) + 1.0 / 64.0) / 64.0;
                xp[j] += dx;
                double[] fp = f(xp);
                for (int i = 0; i < d; i++) {
                    B[i, j] = (fp[i] - f0[i]) / dx;
                }
            }

            return (B);

        }

    }

}
