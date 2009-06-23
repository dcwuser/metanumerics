using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Functions {

    public static partial class FunctionMath {

        // One-dimensional minimization

        /// <summary>
        /// Minimizes a function in the vicinity of a given point.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="x">A point suspected to be near the minimum; the search begins at this point.</param>
        /// <returns>The minimum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null.</exception>
        /// <remarks>
        /// <para>Since the search algorithm begins by evaluating <paramref name="f"/> near <paramref name="x"/>,
        /// it can fail if <paramref name="x"/> is near a singularity or other point at which the evaluation
        /// of <paramref name="f"/> could fail. If you can reliably bracket a minimum, the other
        /// overload of this method is safer and, if your bracket is any good, slightly faster.</para>
        /// </remarks>
        public static LineExtremum FindMinimum (Function<double, double> f, double x) {
            if (f == null) throw new ArgumentNullException("f");
            double fx = f(x);
            double d = 0.1 * Math.Abs(x) + 0.01;
            return (FindMinimum(f, x, fx, d));
        }

        // finds a minimum given an initial step size

        private static LineExtremum FindMinimum (Function<double,double> f, double x, double fx, double d) {

            // take a step
            // add a while loop to make sure we have moved enough that fy != fx
            double y = x + d;
            double fy = f(y);

            // if the step was uphill, reverse
            if (fy > fx) {
                d = -d;

                double t;

                t = y;
                y = x;
                x = t;

                t = fy;
                fy = fx;
                fx = t;
            }

            // keep stepping downhill...
            while (true) {

                double z = y + d;
                double fz = f(z);

                if (fz > fy) {

                    // ...until we go uphill

                    // that means we have bracketed a minimum
                    double a, b;
                    if (x < z) {
                        a = x;
                        b = z;
                    } else {
                        a = z;
                        b = x;
                    }

                    // organize our three points in order from lowest to highest values
                    double u, v, w, fu, fv, fw;
                    u = y;
                    fu = fy;
                    if (fx < fz) {
                        v = x;
                        fv = fx;
                        w = z;
                        fw = fz;
                    } else {
                        v = z;
                        fv = fz;
                        w = x;
                        fw = fx;
                    }

                    // hand the brackets and the three lowest points to Brent's algorithm
                    LineExtremum min = FindMinimum(f, a, b, u, fu, v, fv, w, fw);
                    return (min);

                }

                // ...still going downhill, so shuffle our points
                x = y;
                fx = fy;

                y = z;
                fy = fz;

                // and double the step size 
                d = 2.0 * d;

            }

        }

        /// <summary>
        /// Minimizes a function on a given interval.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="r">The interval.</param>
        /// <returns>The minimum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null.</exception>
        public static LineExtremum FindMinimum (Function<double, double> f, Interval r) {

            if (f == null) throw new ArgumentNullException("f");

            // sample three points within the interval

            double x1 = (3.0 * r.LeftEndpoint + r.RightEndpoint) / 4.0;
            double x2 = (r.LeftEndpoint + r.RightEndpoint) / 2.0;
            double x3 = (r.LeftEndpoint + 3.0 * r.RightEndpoint) / 4.0;

            double f1 = f(x1);
            double f2 = f(x2);
            double f3 = f(x3);

            // order them

            double u, v, w, fu, fv, fw, a, b;
            if ((f1 < f2) && (f1 < f3)) {
                // x1 is minimum
                u = x1;
                fu = f1;
                if (f2 < f3) {
                    v = x2;
                    fv = f2;
                    w = x3;
                    fw = f3;
                } else {
                    v = x3;
                    fv = f3;
                    w = x2;
                    fw = f2;
                }
                a = r.LeftEndpoint;
                b = x3;
            } else if ((f2 < f1) && (f2 < f3)) {
                // x2 is minimum
                u = x2;
                fu = f2;
                if (f1 < f3) {
                    v = x1;
                    fv = f1;
                    w = x3;
                    fw = f3;
                } else {
                    v = x3;
                    fv = f3;
                    w = x1;
                    fw = f1;
                }
                a = x1;
                b = x3;
            } else {
                // x3 is minimum
                u = x3;
                fu = f3;
                if (f1 < f2) {
                    v = x1;
                    fv = f1;
                    w = f2;
                    fw = f2;
                } else {
                    v = x2;
                    fv = f2;
                    w = x1;
                    fw = f1;
                }
                a = x1;
                b = r.RightEndpoint;
            }

            // pass them to Brent's algorithm

            LineExtremum min = FindMinimum(f, a, b, u, fu, v, fv, w, fw);
            return(min);

        }

        // Brent's algorithm: use 3-point parabolic interpolation,
        // switching to golden section if interval does not shrink fast enough
        // see Richard Brent, "Algorithms for Minimization Without Derivatives"

        private static LineExtremum FindMinimum (Function<double,double> f, double a, double b, double u, double fu, double v, double fv, double w, double fw) {

            if (f == null) throw new ArgumentNullException("f");
            
            // keep track of previous step sizes
            double dd = 0.5*(b-a);
            double ddd = (b-a);

            for (int n = 0; n < Global.SeriesMax; n++) {

                //Console.WriteLine("(u,fu)=({0},{1}) (v,fv)=({2},{3}) (w,fw)=({4},{5})", u, fu, v, fv, w, fw);

                // do a parameteric fit to the parabola:
                //   f = f0 + 0.5 * f2 (x - x0)^2
                // through the points (u,fu) (v,fv) (w,fw)

                // the solution is:
                //   f2 = (x1-x2)(f3-f1) + (x3-x1)(f2-f1) / (x1-x2)(x3-x1)(x3-x2)
                //   x0 = x1 - 0.5 * [ (f3-f1)(x1-x2)^2 - (f2-f1)(x3-x1)^2 ] / [ (f3-f1)(x1-x2) + (f2-f1)(x3-x1) ]

                double vu = u - v;
                double uw = w - u;
                double t1 = (fw - fu) * vu;
                double t2 = (fv - fu) * uw;
                double p = t1 * vu - t2 * uw;
                double q = 2.0 * (t1 + t2);

                double d = 0.0;
                double x = 0.0;
                
                // if denominator vanishes, prefer golden section
                if (q != 0.0) {

                    d = - p / q;

                    //double cu = q / (vu * uw * (w - v));
                    //Console.WriteLine("cu = {0}", cu);

                    // if step is not getting smaller, prefer golden section
                    if (Math.Abs(d) < 0.5 * Math.Abs(ddd)) {

                        x = u + d;

                        // if step takes us out of bounds, prefer golden section
                        if ((x <= a) || (x >= b)) d = 0.0;

                    } else {

                        d = 0.0;
                    }
                }

                if (d == 0.0) {

                    //Console.WriteLine("GS");

                    // the parabolic fit didn't work out, so use golden section instead

                    double au = u - a;
                    double ub = b - u;

                    if (au > ub) {
                        d = - au * Gold;
                    } else {
                        d = ub * Gold;
                    }

                    x = u + d;

                }

                //Console.WriteLine("d={0}", d);

                // evaluate the function
                double fx = f(x);

                //Console.WriteLine("(x,fx)=({0},{1})", x, fx);

                // check for convergence
                // check for change in f
                // to do: add check change in x given an input precision target
                // note: the = part of <= is important, because the RHS can be zero
                // note: by adding |f(x)| and |f(u)|, we eliminate the danger of a zero RHS
                double df = fx - fu;
                if (2.0*Math.Abs(df) <= RelativePrecision * (Math.Abs(fu) + Math.Abs(fx))) {

                    // compute curvature
                    double f2 = q / (vu * uw * (w - v));
                    //Console.WriteLine("  f''={0}",f2);

                    LineExtremum minimum = new LineExtremum(x, fx, f2);
                    // minimum.EvaluationCount = n;
                    return (minimum);
                }

                if (fx < fu) {

                    // we have a new lowest point

                    // adjust the bracket
                    if (x < u) {
                        b = u;
                    } else {
                        a = u;
                    }
                    // adjust our set of three lowest points
                    w = v;
                    fw = fv;
                    v = u;
                    fv = fu;
                    u = x;
                    fu = fx;

                } else {

                    // our candidate point was higher

                    // adjust the bracket
                    if (x < u) {
                        a = x;
                    } else {
                        b = x;
                    }

                    // adjust our set of three lowest points
                    if (fx < fv) {

                        w = v;
                        fw = fv;
                        v = x;
                        fv = fx;

                    } else {

                        if (fx < fw) {

                            w = x;
                            fw = fx;
                        }
                    }

                }

                // remember our step size
                ddd = dd;
                dd = d;

            }

            throw new NonconvergenceException();

        }

        // golden ratio ~0.3819660;
        private static readonly double Gold = (3.0 - Math.Sqrt(5.0)) / 2.0;

        //private static readonly double AbsolutePrecision = Math.Pow(2.0, -512);


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
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null.</exception>
        /// <exception cref="ArgumentNullException"><paramref name="x"/> is null.</exception>
        public static SpaceExtremum FindMinimum (Function<double[], double> f, double[] x) {

            if (f == null) throw new ArgumentNullException("f");
            if (x == null) throw new ArgumentNullException("x");

            int d = x.Length;

            // keep track of the (conjugate) minimization directions
            double[][] Q = new double[d][];
            for (int i = 0; i < d; i++) {
                Q[i] = new double[d];
                for (int j = 0; j < d; j++) {
                    Q[i][j] = 0.0;
                }
                // pick a step size in each direction that represents a fraction of the input value
                Q[i][i] = 0.075 * (Math.Abs(x[i]) + 0.0001);
            }
            // keep track of the curvature in these directions
            double[] r = new double[d];

            // keep track of the function value
            double y = f(x);

            bool skip = false;

            // interate until convergence
            for (int n = 0; n < 33; n++ ) {

                //Console.WriteLine("N = {0}", n);

                // remember our starting position
                double[] x0 = new double[d];
                Array.Copy(x, x0, d);
                //Console.WriteLine("Start:");
                //WriteVector(x0);
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
                    LineFunction f1 = new LineFunction(f, x, Q[i]);

                    // minimize it
                    LineExtremum m = FindMinimum(new Function<double,double>(f1.Evaluate), 0.0, y, 1.0);

                    // update the current position
                    x = f1.Position(m.Location);
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
                if (2.0 * Math.Abs(y0 - y) <= (Math.Abs(y) + Math.Abs(y0)) * RelativePrecision) {
                    //Console.WriteLine("Terminated");
                    SymmetricMatrix A = ComputeCurvature(f, x);
                    //Console.WriteLine("A=");
                    //WriteMatrix(A);
                    /*
                    SymmetricMatrix H = ComputeHessian(Q, r);
                    Console.WriteLine("H=");
                    WriteMatrix(H);
                    CholeskyDecomposition CD = H.CholeskyDecomposition();
                    if (CD == null) Console.WriteLine("not positive definite");
                    SymmetricMatrix HI = CD.Inverse();
                    Console.WriteLine("HI=");
                    WriteMatrix(HI);
                    */
                    return (new SpaceExtremum(x, y, A));
                }

                //Console.WriteLine("i = P");

                // attempt a minimization in the net direction
                LineFunction f2 = new LineFunction(f, x, dx);
                LineExtremum mm = FindMinimum(new Function<double,double>(f2.Evaluate), 0.0, y, 1.0);
                x = f2.Position(mm.Location);
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

        private static SymmetricMatrix ComputeCurvature (Function<double[], double> f, double[] x) {

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

        private static void WriteVector (double[] x) {
            for (int i = 0; i < x.Length; i++) {
                Console.WriteLine("  {0}", x[i]);
            }
        }

        private static void WriteMatrix (IMatrix M) {
            int im = M.RowCount;
            int jm = M.ColumnCount;
            for (int i = 0; i < im; i++) {
                for (int j = 0; j < jm; j++) {
                    Console.Write("{0} ", M[i, j]);
                }
                Console.WriteLine();
            }
        }

    }

    internal class LineFunction {

        public LineFunction (Function<double[],double> f, double[] x, double[] dx) {

            if (x.Length != dx.Length) throw new DimensionMismatchException();

            this.f = f;
            this.x = x;
            this.dx = dx;
            this.d = x.Length;
        }

        private Function<double[],double> f;

        private double[] x;

        private double[] dx;

        private int d;

        public double[] Position (double s) {
            double[] y = new double[d];
            for (int i = 0; i < d; i++) {
                y[i] = x[i] + dx[i] * s;
            }
            return (y);
        }

        public double Evaluate (double s) {

            return (f(Position(s)));

        }

    }

    // one-dimensional minimum

    /// <summary>
    /// Represents a maximum or minimum of a function of one variable.
    /// </summary>
    public class LineExtremum {

        internal LineExtremum (double x, double f, double f2) {
            this.x = x;
            this.f = f;
            this.f2 = f2;
        }

        private double x;
        private double f;
        private double f2;

        /// <summary>
        /// Gets the location (x-value) of the extremum.
        /// </summary>
        /// <remarks>
        /// <para>Note that numerical methods for finding typical a maximum or minimum cannot determine
        /// its location to full precision. Near a quadratic extremum, a change in x of ~&#x3B5; will
        /// change f(x) by ~&#x3B5;<sup>2</sup>. Thus the smallest detectable change in f(x) will
        /// typically correspond to a change in x of order of the square root of full precision. Full
        /// <see cref="System.Double"/> precision being ~16 digits, you should expect the location to
        /// be accurate only to ~8 digits.</para>
        /// </remarks>
        public double Location {
            get {
                return (x);
            }
        }

        /// <summary>
        /// Gets the function value (y-value) at the extremum.
        /// </summary>
        public double Value {
            get {
                return (f);
            }
        }

        /// <summary>
        /// Gets the curvature at the extremum.
        /// </summary>
        /// <remarks>
        /// <para>The curvature is the second derivative of the function at the minimum.</para>
        /// <para>At a typical minimum, where the function has vanishing first derivative, the second derivative will be a positive number
        /// whose magnitude characterizes the "steepness" with which the function increases as one moves away from the minimum.</para>
        /// <para>At an atypical minimum, for example an interval boundary or a non-smooth function, this
        /// value may be meaningless.</para>
        /// <para>Even in the case of a typical minimum, the value of the curvature property will typically be accurate only
        /// to a handfull of digits. If you require a highly accurate determination of the curvature,
        /// you should compute the second derivative of the minimzed function explicitly.</para>
        /// </remarks>
        public double Curvature {
            get {
                return (f2);
            }
        }

        /// <summary>
        /// Converts a line extremum to a one-dimensional space extremum.
        /// </summary>
        /// <param name="m">The line extremum.</param>
        /// <returns>The corresponding one-dimensional space extremum.</returns>
        public static implicit operator SpaceExtremum (LineExtremum m) {
            double[] s_x = new double[1] { m.x };
            double s_f = m.f;
            SymmetricMatrix s_f2 = new SymmetricMatrix(1);
            s_f2[0, 0] = m.f2;
            return (new SpaceExtremum(s_x, s_f, s_f2));
        }

        // keep track of iterations

        /*
        private int count;

        public int EvaluationCount {
            get {
                return (count);
            }
            internal set {
                count = value;
            }
        }
        */

    }

    // multi-dimensional minimum

    /// <summary>
    /// Represents a maximum or minimum of a function on a multi-dimensional space.
    /// </summary>
    public class SpaceExtremum {

        internal SpaceExtremum (double[] x, double f, SymmetricMatrix f2) {
            this.x = x;
            this.f= f;
            this.f2 = f2;
        }

        private double[] x;
        private double f;
        private SymmetricMatrix f2;

        /// <summary>
        /// Gets the location of the extremum.
        /// </summary>
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
        public SymmetricMatrix Curvature () {
            return (f2.Clone());
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

    /*
    public class Minimum : FunctionPoint<double> {

        internal Minimum (double x, double y, double d2xdy2) : base(x, y) {
            this.f2 = d2xdy2;
        }

        private double f2;

        /// <summary>
        /// The curvatre (second derivative) at the minimum.
        /// </summary>
        /// <remarks>
        /// <para>The curvature is the second derivative at the minimum.</para>
        /// <para>A atypical minima (e.g. an interval end-point or for a non-smooth function),
        /// this value may not be meaningful. Even for a typical minimum, the value may be
        /// accurate only to a few significant digits.</para>
        /// </remarks>
        public double Curvature {
            get {
                return (f2);
            }
        }

    }
    */ 

}
