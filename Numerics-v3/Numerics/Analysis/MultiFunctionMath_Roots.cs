using System;
using System.Collections.Generic;

using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    public static partial class MultiFunctionMath {

        /// <summary>
        /// Finds a vector argument which makes a vector function zero.
        /// </summary>
        /// <param name="f">The vector function.</param>
        /// <param name="x0">The vector argument.</param>
        /// <returns>The vector argument which makes all components of the vector function zero.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> or <paramref name="x0"/> is null.</exception>
        /// <exception cref="DimensionMismatchException">The dimension of <paramref name="f"/> is not equal to the
        /// dimension of <paramref name="x0"/>.</exception>
        public static ColumnVector FindZero (Func<IList<double>, IList<double>> f, IList<double> x0) {

            if (f == null) throw new ArgumentNullException("f");
            if (x0 == null) throw new ArgumentNullException("x0");

            // we will use Broyden's method, which is a generalization of the secant method to the multi-dimensional problem
            // just as the secant method is essentially Newton's method with a crude, numerical value for the slope,
            // Broyden's method is a multi-dimensional Newton's method with a crude, numerical value for the Jacobian matrix

            int d = x0.Count;

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
                ColumnVector F1 = new ColumnVector(f(x1));
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
                    return (x);
                }

                // update B
                ColumnVector dF = F1 - F;
                RectangularMatrix dB = (dF - B * dx) * (dx / (dx.Transpose() * dx)).Transpose();
                //RectangularMatrix dB = F1 * ( dx / MoreMath.Pow(dx.Norm(), 2) ).Transpose();
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < d; j++) {
                        B[i, j] += dB[i, j];
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

        private static SquareMatrix ApproximateJacobian (Func<IList<double>, IList<double>> f, IList<double> x0) {

            int d = x0.Count;

            IList<double> f0 = f(x0);

            double[] xp = new double[d];

            SquareMatrix B = new SquareMatrix(d);
            for (int j = 0; j < d; j++) {
                x0.CopyTo(xp, 0);
                double dx = (Math.Abs(x0[j]) + 1.0 / 64.0) / 64.0;
                xp[j] += dx;
                IList<double> fp = f(xp);
                for (int i = 0; i < d; i++) {
                    B[i, j] = (fp[i] - f0[i]) / dx;
                }
            }

            return (B);

        }


    }

}
