using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    public static partial class FunctionMath {

        // Multidimensional minimization

        // 1. Iterate on axes
        // 2. Powell's method with reset
        // 3. Powell's heuristic method
        // 4. Powell's method with Brent modification

        /// <summary>
        /// Minimizes a function on a multi-dimensional space in the vicinity of a given point. 
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="x">The starting point for the search.</param>
        /// <returns>The minimum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> or <paramref name="x"/> is null</exception>
        /// <exception cref="NonconvergenceException">The minimum was not found to the required precision within the budgeted number of function evaluations.</exception>
        internal static SpaceExtremum FindMinimum (Func<double[], double> f, double[] x) {

            return (FindMinimum(f, x, new EvaluationSettings() {
                RelativePrecision = Global.Accuracy,
                AbsolutePrecision = Global.Accuracy * Global.Accuracy,
                EvaluationBudget = 1024 * x.Length
            }));

        }

        /// <summary>
        /// Minimizes a function on a multi-dimensional space in the vicinity of a given point, subject to the given settings. 
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="x">The starting point for the search.</param>
        /// <param name="settings">The evaluation settings.</param>
        /// <returns>The minimum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/>, <paramref name="x"/>, or <paramref name="settings"/> is null.</exception>
        /// <exception cref="NonconvergenceException">The minimum was not found to the required precision within the budgeted number of function evaluations.</exception>
        internal static SpaceExtremum FindMinimum (Func<double[], double> f, double[] x, EvaluationSettings settings) {

            int d = x.Length;

            // put the function into a Functor we will use for line searches
            LineSearchFunctor fo = new LineSearchFunctor(f);

            // keep track of the (conjugate) minimization directions
            double[][] Q = new double[d][];
            for (int i = 0; i < d; i++) {
                Q[i] = new double[d];
                for (int j = 0; j < d; j++) {
                    Q[i][j] = 0.0;
                }
                // pick a step size in each direction that represents a fraction of the input value
                Q[i][i] = (1.0 / 16.0) * (Math.Abs(x[i]) + (1.0 / 16.0));
            }
            // keep track of the curvature in these directions
            double[] r = new double[d];

            // keep track of the function value
            double y = f(x);

            // keep track of the total number of evaluations
            //int count = 0;

            bool skip = false;

            // interate until convergence
            while (fo.EvaluationCount < settings.EvaluationBudget) {

                // remember our starting position
                double[] x0 = new double[d];
                Array.Copy(x, x0, d);
                double y0 = y;

                // keep track of direction of largest decrease
                double max_dy = 0.0;
                int max_i = 0;

                // now minimize in each direction
                for (int i = 0; i < d; i++) {

                    // if we placed the net direction in the first slot last time,
                    // we are already at the minimum along that direction, so don't
                    // minimize along it
                    if (skip) {
                        skip = false;
                        continue;
                    }
                    //if ((n > 0) && (i == 0) && (max_i == 0)) continue;

                    //Console.WriteLine("i = {0}", i);
                    //WriteVector(Q[i]);

                    // form a line function
                    //LineFunction f1 = new LineFunction(f, x, Q[i]);
                    fo.Origin = x;
                    fo.Direction = Q[i];

                    // minimize it
                    Extremum m = FindMinimum(fo, 0.0, y, 1.0, new ExtremumSettings() { EvaluationBudget = settings.EvaluationBudget, AbsolutePrecision = 0.0, RelativePrecision = 0.0 });
                    //LineExtremum m = FindMinimum(new Func<double,double>(f1.Evaluate), 0.0, y, 1.0);

                    // add to the evaluation count
                    //count += f1.Count;

                    // update the current position
                    x = fo.ComputeLocation(m.Location);
                    //x = f1.Position(m.Location);
                    //WriteVector(x);
                    r[i] = m.Curvature;

                    // keep track of how much the function dropped, and
                    // if this is the direction of largest decrease
                    double dy = y - m.Value;
                    //Console.WriteLine("dy = {0}", dy);
                    if (dy > max_dy) {
                        max_dy = dy;
                        max_i = i;
                    }
                    y = m.Value;


                }

                //Console.WriteLine("max_i = {0}, max_dy = {1}", max_i, max_dy);
                //Console.WriteLine("y0 = {0}, y = {1}", y0, y);

                // figure out the net direction we have moved
                double[] dx = new double[d];
                for (int i = 0; i < d; i++) {
                    dx[i] = x[i] - x0[i];
                }
                //Console.WriteLine("Finish:");
                //WriteVector(x);
                //Console.WriteLine("Net direction:");
                //WriteVector(dx);

                // check termination criteria
                // we do this before minimizing in the net direction because if dx=0 it loops forever
                double Dy = Math.Abs(y0 - y);
                if ((Dy < settings.AbsolutePrecision) || (2.0 * Dy < (Math.Abs(y) + Math.Abs(y0)) * settings.RelativePrecision)) {
                    SymmetricMatrix A = ComputeCurvature(f, x);
                    return (new SpaceExtremum(x, y, A));
                }

                // attempt a minimization in the net direction
                fo.Origin = x;
                fo.Direction = dx;
                //LineFunction f2 = new LineFunction(f, x, dx);
                //LineExtremum mm = FindMinimum(new Func<double,double>(f2.Evaluate), 0.0, y, 1.0);
                Extremum mm = FindMinimum(fo, 0.0, y, 1.0, new ExtremumSettings() { EvaluationBudget = settings.EvaluationBudget, RelativePrecision = 0.0, AbsolutePrecision = 0.0 });
                //count += f2.Count;
                //x = f2.Position(mm.Location);
                x = fo.ComputeLocation(mm.Location);
                y = mm.Value;

                // rotate this direction into the direction set
                /*
                for (int i = 0; i < (d - 1); i++) {
                    Q[i] = Q[i + 1];
                    r[i] = r[i + 1];
                }
                Q[d - 1] = dx;
                r[d - 1] = mm.Curvature;
                */
                // this is the basic Powell procedure, and it leads to linear dependence

                // replace the direction of largest decrease with the net direction
                Q[max_i] = dx;
                r[max_i] = mm.Curvature;
                if (max_i == 0) skip = true;
                // this is powell's modification to avoid linear dependence

                // reset

            }

            throw new NonconvergenceException();

        }

        // numerical approximation of Hessian
        // requires 3 evaluations for diagonals, there are d of those
        // requires 4 evaluations for off-diagonals; there are d(d-1)/2 of those
        // total of 2d^2 + d evaluations required

        private static SymmetricMatrix ComputeCurvature (Func<double[], double> f, double[] x) {

            int d = x.Length;

            double e = Math.Pow(2.0, -15.0);
            double[] dx = new double[d];
            for (int i = 0; i < d; i++) {
                double h = e * (Math.Abs(x[i]) + 1.0);
                // ensure that step is exactly representable
                double xh = x[i] + h;
                h = xh - x[i];
                // record d
                dx[i] = h;
            }

            SymmetricMatrix H = new SymmetricMatrix(d);
            for (int i = 0; i < d; i++) {
                double[] xp = (double[]) x.Clone(); xp[i] += dx[i];
                double[] xm = (double[]) x.Clone(); xm[i] -= dx[i];
                double f0 = f(x);
                double fp = f(xp);
                double fm = f(xm);
                H[i, i] = (fm - 2.0 * f0 + fp) / (dx[i]*dx[i]);
                for (int j = 0; j < i; j++) {
                    double[] xpp = (double[]) x.Clone(); xpp[i] += dx[i]; xpp[j] += dx[j];
                    double[] xpm = (double[]) x.Clone(); xpm[i] += dx[i]; xpm[j] -= dx[j];
                    double[] xmm = (double[]) x.Clone(); xmm[i] -= dx[i]; xmm[j] -= dx[j];
                    double[] xmp = (double[]) x.Clone(); xmp[i] -= dx[i]; xmp[j] += dx[j];
                    double fpp = f(xpp);
                    double fpm = f(xpm);
                    double fmm = f(xmm);
                    double fmp = f(xmp);
                    H[i, j] = (fpp - fpm - fmp + fmm) / dx[i] / dx[j] / 4.0; ;
                }
            }

            return (H);
        }

#if PAST
        /// <summary>
        /// Maximizes a function on a multi-dimensional space in the vicinity of a given point, subject to the given settings. 
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="x">The starting point for the search.</param>
        /// <param name="settings">The evaluation settings.</param>
        /// <returns>The maximum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/>, <paramref name="x"/>, or <paramref name="settings"/> is null.</exception>
        /// <exception cref="NonconvergenceException">The maximum was not found to the required precision within the budgeted number of function evaluations.</exception>
        internal static SpaceExtremum FindMaximum (Func<double[], double> f, double[] x, EvaluationSettings settings) {
            return (FindMinimum((double[] p) => -f(p), x, settings));
        }
#endif

    }


    internal class LineSearchFunctor : Functor {

        public LineSearchFunctor (Func<double[], double> f) : base(null) {
                this.Function = LineFunction;
                this.SpaceFunction = f;
        }

        public Func<double[], double> SpaceFunction { get; private set; }

        public double[] Origin { get; set; }

        public double[] Direction { get; set; }

        public double[] ComputeLocation (double s) {
            double[] x = new double[Origin.Length];
            for (int i = 0; i < x.Length; i++) {
                x[i] = Origin[i] + Direction[i] * s;
            }
            return (x);
        }

        private double LineFunction (double s) {
            return (SpaceFunction(ComputeLocation(s)));
        }


    }

    /*
    internal class LineFunction {

        public LineFunction (Func<double[], double> f, double[] x, double[] dx) {

            if (x.Length != dx.Length) throw new DimensionMismatchException();

            this.f = f;
            this.x = x;
            this.dx = dx;
            this.d = x.Length;
            this.Count = 0;
        }

        private Func<double[], double> f;

        private double[] x;

        private double[] dx;

        private int d;

        public int Count { get; private set; }

        public double[] Position (double s) {
            double[] y = new double[d];
            for (int i = 0; i < d; i++) {
                y[i] = x[i] + dx[i] * s;
            }
            return (y);
        }

        public double Evaluate (double s) {
            Count++;
            return (f(Position(s)));
        }

    }
    */

    // multi-dimensional minimum

    /// <summary>
    /// Represents a maximum or minimum of a function on a multi-dimensional space.
    /// </summary>
    internal class SpaceExtremum {

        internal SpaceExtremum (double[] x, double f, SymmetricMatrix f2) {
            this.x = x;
            this.f = f;
            this.f2 = f2;
        }

        private double[] x;
        private double f;
        private SymmetricMatrix f2;

        /// <summary>
        /// Gets the location of the extremum.
        /// </summary>
        /// <returns>The coordinates of the extremum.</returns>
        public double[] Location () {
            return ((double[]) x.Clone());
        }

        /// <summary>
        /// Gets the value of the function at the extremum.
        /// </summary>
        public double Value {
            get {
                return (f);
            }
        }

        /// <summary>
        /// Gets the curvature matrix at the extremum.
        /// </summary>
        /// <returns>The curvature matrix.</returns>
        public SymmetricMatrix Curvature () {
            return (f2.Copy());
        }

        /// <summary>
        /// Gets the dimension of the space on which the function is defined. 
        /// </summary>
        public int Dimension {
            get {
                return (x.Length);
            }
        }

    }

}
