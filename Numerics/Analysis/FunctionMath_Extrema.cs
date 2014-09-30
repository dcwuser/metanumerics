using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    // The basic public 1D extrema seach methods are
    //   Min(f, x), Min(f, x, settings), Min(f, interval), Min(f, interval, settings)
    // and the corresponding methods for maxima. The flow of control from these to private methods is as follows.
    // For the point methods:
    //   Min(f, x) -> Min(f, x, settings) -> Min(functor, x, settings) <- Max(f, x, settings) <- Max(f, x)
    // then
    //   Min(functor, start, settings) -> Min(functor, x, fx, step, settings) -> Min(functor, a, b, u, fu, v, fv, w, fw, settings)
    // For the interval methods:
    //   Min(f, interval) -> Min(f, interval, settings) -> Min(functor, a, b, settings) <- Max(f, interval, settings) <- Max(f, interval)
    // then
    //   Min(functor, a, b, settings) -> Min(functor, a, b, u, fu, v, fv, w, fw, settings)
    // So in the end, they all flow into the Brent algorithm.

    public static partial class FunctionMath {

        // One-dimensional minimization

        // An evaluation budget of 32 is sufficient for all our test cases except for |x|, which requires 82 (!) evaluations to converge. Parabolic fitting just does a very poor job
        // for this function (at all scales, since it is scale invariant).

        private static readonly EvaluationSettings DefaultExtremaSettings = new EvaluationSettings() { EvaluationBudget = 128, RelativePrecision = 0.0, AbsolutePrecision = 0.0 };

        /// <summary>
        /// Maximizes a function in the vicinity of a given point.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="x">A point suspected to be near the maximum. The search begins at this point.</param>
        /// <returns>The maximum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null.</exception>
        /// <exception cref="NonconvergenceException">More than the maximum allowed number of function evaluations occured without a maximum being determined.</exception>
        public static Extremum FindMaximum (Func<double, double> f, double x) {
            return (FindMaximum(f, x, DefaultExtremaSettings));
        }


        /// <summary>
        /// Maximizes a function in the vicinity of a given point, subject to the given evaluation settings.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="x">A point suspected to be near the maximum. The search begins at this point.</param>
        /// <param name="settings">The settings to use when searching for the maximum.</param>
        /// <returns>The maximum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null or <paramref name="settings"/> is null.</exception>
        /// <exception cref="NonconvergenceException">More than the maximum allowed number of function evaluations occured without a maximum being determined to the prescribed precision.</exception>
        /// <remarks>
        /// <para>When you supply <paramref name="settings"/>, note that the supplied <see cref="EvaluationSettings.RelativePrecision"/> and <see cref="EvaluationSettings.AbsolutePrecision"/>
        /// values refer to argument (i.e. x) values, not function (i.e. f) values. Note also that, for typical functions, the best attainable relative precision is of the order of the
        /// square root of machine precision (about 10<sup>-7</sup>), i.e. half the number of digits in a <see cref="Double"/>. This is because to identify an extremum we need to resolve changes
        /// in the function value, and near an extremum  &#x3B4;f &#x223C; (&#x3B4;x)<sup>2</sup>, so changes in the function value &#x3B4;f &#x223C; &#x3B5; correspond to changes in the
        /// argument value &#x3B4;x &#x223C; &#x221A;&#x3B5;. If you supply zero values for both precision settings, the method will adaptively approximate the best attainable precision for
        /// the supplied function and locate the extremum to that resolution. This is our suggested practice unless you know that you require a less precise determination.</para>
        /// </remarks>
        public static Extremum FindMaximum (Func<double, double> f, double x, EvaluationSettings settings) {
            if (f == null) throw new ArgumentNullException("f");
            if (settings == null) throw new ArgumentNullException("settings");
            return (FindMinimum(new Functor(f, true), x, settings).Negate());
        }

        /// <summary>
        /// Minimizes a function in the vicinity of a given point.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="x">A point suspected to be near the minimum. The search begins at this point.</param>
        /// <returns>The minimum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null.</exception>
        /// <exception cref="NonconvergenceException">More than the maximum allowed number of function evaluations occured without a minimum being determined.</exception>
        public static Extremum FindMinimum (Func<double, double> f, double x) {
            return (FindMinimum(f, x, DefaultExtremaSettings));
        }

        /// <summary>
        /// Minimizes a function in the vicinity of a given point subject to the given evaluation settings.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="x">A point suspected to be near the minimum. The search begins at this point.</param>
        /// <param name="settings">The settings to use when searching for the minimum.</param>
        /// <returns>The minimum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null or <paramref name="settings"/> is null.</exception>
        /// <exception cref="NonconvergenceException">More than the maximum allowed number of function evaluations occured without a minimum being determined to the prescribed precision.</exception>
        /// <remarks>
        /// <para>When you supply <paramref name="settings"/>, note that the supplied <see cref="EvaluationSettings.RelativePrecision"/> and <see cref="EvaluationSettings.AbsolutePrecision"/>
        /// values refer to argument (i.e. x) values, not function (i.e. f) values. Note also that, for typical functions, the best attainable relative precision is of the order of the
        /// square root of machine precision (about 10<sup>-7</sup>), i.e. half the number of digits in a <see cref="Double"/>. This is because to identify an extremum we need to resolve changes
        /// in the function value, and near an extremum  &#x3B4;f &#x223C; (&#x3B4;x)<sup>2</sup>, so changes in the function value &#x3B4;f &#x223C; &#x3B5; correspond to changes in the
        /// argument value &#x3B4;x &#x223C; &#x221A;&#x3B5;. If you supply zero values for both precision settings, the method will adaptively approximate the best attainable precision for
        /// the supplied function and locate the extremum to that resolution. This is our suggested practice unless you know that you require a less precise determination.</para>
        /// <para>Since the search algorithm begins by evaluating <paramref name="f"/> at points near <paramref name="x"/>, it can fail if <paramref name="x"/> is near a singularity
        /// or other point at which the evaluation of <paramref name="f"/> could fail. If you can reliably bracket an extremum, the <see cref="FindMinimum(Func{Double,Double},Interval,EvaluationSettings)"/>
        /// overload of this method is safer and, if your bracket is any good, usually slightly faster.</para>
        /// </remarks>
        public static Extremum FindMinimum (Func<double, double> f, double x, EvaluationSettings settings) {
            
            if (f == null) throw new ArgumentNullException("f");
            if (settings == null) throw new ArgumentNullException("settings");

            return (FindMinimum(new Functor(f), x, settings));
        }

        private static Extremum FindMinimum (Functor f, double x, EvaluationSettings settings) {

            Debug.Assert(f != null); Debug.Assert(settings != null);

            // To call the bracketing function, we need an initial value and an initial step.
            double fx = f.Evaluate(x);
            // Pick a relatively small initial value to try to avoid running into any nearly singularities.
            double d = (Math.Abs(x) + 1.0 / 32.0) / 32.0;

            return(FindMinimum(f, x, fx, d, settings));

        }

        private static Extremum FindMinimum (Functor f, double x, double fx, double d, EvaluationSettings settings) {

            // This function brackets a minimum by starting from x and taking increasing steps downhill until it moves uphill again.

            // We write this function assuming f(x) has already been evaluated because when it is called
            // to do a line fit for Powell's multi-dimensional minimization routine, that is the case and
            // we don't want to do a superfluous evaluation.

            Debug.Assert(f != null); Debug.Assert(d > 0.0); Debug.Assert(settings != null);

            // evaluate at x + d
            double y = x + d;
            double fy = f.Evaluate(y);

            // if we stepped uphill, reverse direction of steps and exchange x & y
            if (fy > fx) {
                Global.Swap(ref x, ref y); Global.Swap(ref fx, ref fy);
                d = -d;
            }

            // we now know f(x) >= f(y) and we are stepping downhill
            // continue stepping until we step uphill
            double z, fz;
            while (true) {

                if (f.EvaluationCount >= settings.EvaluationBudget) throw new NonconvergenceException();

                z = y + d;
                fz = f.Evaluate(z);

                Debug.WriteLine(String.Format("f({0})={1} f({2})={3} f({4})={5} d={6}", x, fx, y, fy, z, fz, d));

                if (fz > fy) break;

                // increase the step size each time
                d = AdvancedMath.GoldenRatio * d;

                // x <- y <- z
                x = y; fx = fy; y = z; fy = fz;


            }

            // we x and z now bracket a local minimum, with y the lowest point evaluated so far
            double a = Math.Min(x, z); double b = Math.Max(x, z);
            if (fz < fx) { Global.Swap(ref x, ref z); Global.Swap(ref fx, ref fz);}

            return (FindMinimum(f, a, b, y, fy, x, fx, z, fz, settings));

        }

        /// <summary>
        /// Maximizes a function on the given interval.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="r">The interval.</param>
        /// <returns>The maximum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null.</exception>
        /// <exception cref="NonconvergenceException">More than the maximum allowed number of function evaluations occured without a minimum being determined.</exception>
        public static Extremum FindMaximum (Func<double, double> f, Interval r) {
            return (FindMaximum(f, r, DefaultExtremaSettings));
        }

        /// <summary>
        /// Maximizes a function on the given interval, subject to the given evaluation settings.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="r">The interval.</param>
        /// <param name="settings">The settings used when searching for the maximum.</param>
        /// <returns>The maximum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null or <paramref name="settings"/> is null.</exception>
        /// <exception cref="NonconvergenceException">More than the maximum allowed number of function evaluations occured without a maximum being determined to the prescribed precision.</exception>
        /// <remarks>
        /// <para>When you supply <paramref name="settings"/>, note that the supplied <see cref="EvaluationSettings.RelativePrecision"/> and <see cref="EvaluationSettings.AbsolutePrecision"/>
        /// values refer to argument (i.e. x) values, not function (i.e. f) values. Note also that, for typical functions, the best attainable relative precision is of the order of the
        /// square root of machine precision (about 10<sup>-7</sup>), i.e. half the number of digits in a <see cref="Double"/>. This is because to identify an extremum we need to resolve changes
        /// in the function value, and near an extremum  &#x3B4;f &#x223C; (&#x3B4;x)<sup>2</sup>, so changes in the function value &#x3B4;f &#x223C; &#x3B5; correspond to changes in the
        /// argument value &#x3B4;x &#x223C; &#x221A;&#x3B5;. If you supply zero values for both precision settings, the method will adaptively approximate the best attainable precision for
        /// the supplied function and locate the extremum to that resolution. This is our suggested practice unless you know that you require a less precise determination.</para>
        /// </remarks>
        public static Extremum FindMaximum (Func<double, double> f, Interval r, EvaluationSettings settings) {
            if (f == null) throw new ArgumentNullException("f");
            if (settings == null) throw new ArgumentNullException("settings");
            return (FindMinimum(new Functor(f, true), r.LeftEndpoint, r.RightEndpoint, settings).Negate());
        }

        /// <summary>
        /// Minimizes a function on the given interval.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="r">The interval.</param>
        /// <returns>The minimum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null.</exception>
        /// <exception cref="NonconvergenceException">More than the maximum allowed number of function evaluations occured without a minimum being determined.</exception>
        public static Extremum FindMinimum (Func<double, double> f, Interval r) {
            return (FindMinimum(f, r, DefaultExtremaSettings));
        }

        /// <summary>
        /// Minimizes a function on the given interval, subject to the given evaluation settings.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <param name="r">The interval.</param>
        /// <param name="settings">The settings used when searching for the minimum.</param>
        /// <returns>The minimum.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is null or <paramref name="settings"/> is null.</exception>
        /// <exception cref="NonconvergenceException">More than the maximum allowed number of function evaluations occured without a minimum being determined to the prescribed precision.</exception>
        /// <remarks>
        /// <para>When you supply <paramref name="settings"/>, note that the supplied <see cref="EvaluationSettings.RelativePrecision"/> and <see cref="EvaluationSettings.AbsolutePrecision"/>
        /// values refer to argument (i.e. x) values, not function (i.e. f) values. Note also that, for typical functions, the best attainable relative precision is of the order of the
        /// square root of machine precision (about 10<sup>-7</sup>), i.e. half the number of digits in a <see cref="Double"/>. This is because to identify an extremum we need to resolve changes
        /// in the function value, and near an extremum  &#x3B4;f &#x223C; (&#x3B4;x)<sup>2</sup>, so changes in the function value &#x3B4;f &#x223C; &#x3B5; correspond to changes in the
        /// argument value &#x3B4;x &#x223C; &#x221A;&#x3B5;. If you supply zero values for both precision settings, the method will adaptively approximate the best attainable precision for
        /// the supplied function and locate the extremum to that resolution. This is our suggested practice unless you know that you require a less precise determination.</para>
        /// </remarks>
        public static Extremum FindMinimum (Func<double, double> f, Interval r, EvaluationSettings settings) {
            if (f == null) throw new ArgumentNullException("f");
            if (settings == null) throw new ArgumentNullException("settings");
            return (FindMinimum(new Functor(f), r.LeftEndpoint, r.RightEndpoint, settings));
        }

        private static Extremum FindMinimum (
            Functor f,
            double a, double b,
            EvaluationSettings settings
        ) {

            // evaluate three points within the bracket
            double u = (3.0 * a + b) / 4.0;
            double v = (a + b) / 2.0;
            double w = (a + 3.0 * b) / 4.0;

            double fu = f.Evaluate(u); double fv = f.Evaluate(v); double fw = f.Evaluate(w);

            Debug.WriteLine(String.Format("f({0})={1}  f({2})={3}  f({4})={5}", u, fu, v, fv, w, fw));

            // move in the bracket boundaries, if possible
            if (fv < fu) { a = u; if (fw < fv) a = v; }
            if (fv < fw) { b = w; if (fu < fv) b = v; }

            Debug.WriteLine(String.Format("[{0} {1}]", a, b));

            // sort u, v, w by fu, fv, fw values
            // these three comparisons are the most efficient three-item sort
            if (fv < fu) { Global.Swap(ref v, ref u); Global.Swap(ref fv, ref fu); }
            if (fw < fu) { Global.Swap(ref w, ref u); Global.Swap(ref fw, ref fu); }
            if (fw < fv) { Global.Swap(ref w, ref v); Global.Swap(ref fw, ref fv); }

            // pass to Brent's algorithm to shrink bracket
            return (FindMinimum(f, a, b, u, fu, v, fv, w, fw, settings));

        }

        // Brent's algorithm: use 3-point parabolic interpolation,
        // switching to golden section if interval does not shrink fast enough
        // see Richard Brent, "Algorithms for Minimization Without Derivatives"

        // The bracket is [a, b] and the three lowest points are (u,fu), (v,fv), (w, fw)
        // Note that a or b may be u, v, or w.

        private static Extremum FindMinimum (
            Functor f,
            double a, double b,
            double u, double fu, double v, double fv, double w, double fw,
            EvaluationSettings settings
        ) {

            double tol = 0.0; double fpp = Double.NaN;

            while (f.EvaluationCount < settings.EvaluationBudget) {

                Debug.WriteLine(String.Format("n={0} tol={1}", f.EvaluationCount, tol));
                Debug.WriteLine(String.Format("[{0}  f({1})={2}  f({3})={4}  f({5})={6}  {7}]", a, u, fu, v, fv, w, fw, b));

                // While a < u < b is guaranteed, a < v, w < b is not guaranteed, since the bracket can sometimes be made tight enough to exclude v or w.
                // For example, if u < v < w, then we can set b = v, placing w outside the bracket.

                Debug.Assert(a < b);
                Debug.Assert((a <= u) && (u <= b));
                Debug.Assert((fu <= fv) && (fv <= fw));

                // Expected final situation is a<tol><tol>u<tol><tol>b, leaving no point left to evaluate that is not within tol of an existing point.

                if ((b - a) <= 4.0 * tol) return (new Extremum(u, fu, fpp, f.EvaluationCount));

                double x; ParabolicFit(u, fu, v, fv, w, fw, out x, out fpp);
                Debug.WriteLine(String.Format("parabolic x={0} f''={1}", x, fpp));

                if (Double.IsNaN(fpp) || (fpp <= 0.0) || (x < a) || (x > b)) {

                    // the parabolic fit didn't work out, so do a golden section reduction instead

                    // to get the most reduction of the bracket, pick the larger of au and ub
                    // for self-similarity, pick a point inside it that divides it into two segments in the golden section ratio,
                    // i.e. 0.3820 = \frac{1}{\phi + 1} and 0.6180 = \frac{\phi}{\phi+1}
                    // put the smaller segment closer to u so that x is closer to u, the best minimum so far

                    double au = u - a;
                    double ub = b - u;

                    if (au > ub) {
                        x = u - au / (AdvancedMath.GoldenRatio + 1.0);
                    } else {
                        x = u + ub / (AdvancedMath.GoldenRatio + 1.0);
                    }

                    Debug.WriteLine(String.Format("golden section x={0}", x));

                }

                // ensure we don't evaluate within tolerance of an existing point
                if (Math.Abs(x - u) < tol) { Debug.WriteLine(String.Format("shift from u (x={0})", x)); x = (x > u) ? u + tol : u - tol; }
                if ((x - a) < tol) { Debug.WriteLine(String.Format("shift from a (x={0})", x)); x = a + tol; }
                if ((b - x) < tol) { Debug.WriteLine(String.Format("shift from b (x={0})", x)); x = b - tol; }

                // evaluate the function at the new point x
                double fx = f.Evaluate(x);
                Debug.WriteLine(String.Format("f({0}) = {1}", x, fx));
                Debug.WriteLine(String.Format("delta={0}", fu - fx));

                // update a, b and u, v, w based on new point x

                if (fx < fu) {

                    // the new point is lower than all the others; this is success 

                    // u now becomes a bracket point
                    if (u < x) {
                        a = u;
                    } else {
                        b = u;
                    }

                    // x -> u -> v -> w
                    w = v; fw = fv;
                    v = u; fv = fu;
                    u = x; fu = fx;

                } else {

                    // x now becomes a bracket point
                    if (x < u) {
                        a = x;
                    } else {
                        b = x;
                    }

                    if (fx < fv) {

                        // the new point is higher than u, but still lower than v and w
                        // this isn't what we expected, but we have lower points that before

                        // x -> v -> w
                        w = v; fw = fv;
                        v = x; fv = fx;

                    } else if (fx < fw) {

                        // x -> w
                        w = x; fw = fx;

                    } else {

                        // the new point is higher than all our other points; this is the worst case

                        // we might still want to replace w with x because
                        // (i) otherwise a parabolic fit will reproduce the same x and
                        // (ii) w is quite likely far outside the new bracket and not telling us much about the behavior near u
                        // w = x; fw = fx;
                        // but tests with Rosenbrock function indicate this increases evaluation count

                        Debug.WriteLine("bad point");

                    }

                }

                // if the user has specified a tollerance, use it
                if ((settings.RelativePrecision > 0.0 || settings.AbsolutePrecision > 0.0)) {
                    tol = Math.Max(Math.Abs(u) * settings.RelativePrecision, settings.AbsolutePrecision);
                } else {
                    // otherwise, try to get the tollerance from the curvature
                    if (fpp > 0.0) {
                        tol = Math.Sqrt(2.0 * Global.Accuracy * (Math.Abs(fu) + Global.Accuracy) / fpp);
                    } else {
                        // but if we don't have a useable curvature either, wing it
                        if (tol == 0.0) tol = Math.Sqrt(Global.Accuracy);
                    }
                }

            }

            throw new NonconvergenceException();

        }

        private static void ParabolicFit (
           double x1, double f1, double x2, double f2, double x3, double f3,
           out double x0, out double fpp
       ) {

            // We want to find the parabola
            //   f = f0 + (f'' / 2) (x - x0)^2
            // that passes throught the points (x1,f1), (x2,f2), (x3,f3)

            // The solution, after much algebra is:
            //   x0 = x1 - \frac{1}{2} \frac{(f3-f1)(x1-x2)^2 - (f2-f1)(x3-x1)^2}{(f3-f1)(x1-x2) + (f2-f1)(x3-x1)}
            //   f0 = ?
            //   f'' = 2 \frac{(x1-x2)(f3-f1) + (x3-x1)(f2-f1)}{(x1-x2)(x3-x1)(x3-x2)}

            // compute the differences that appear in our expressions
            double x12 = x1 - x2; double f21 = f2 - f1;
            double x31 = x3 - x1; double f31 = f3 - f1;
            double x32 = x3 - x2;

            // compute the numerator and denominator in the expression for x0
            double t1 = f31 * x12;
            double t2 = f21 * x31;
            double p = t1 * x12 - t2 * x31;
            double q = 2.0 * (t1 + t2);

            if (q == 0.0) {
                // the denominator vanishes when the points are colinear; there is no parabolic fit
                x0 = Double.NaN;
                fpp = Double.NaN;
            } else {
                x0 = x1 - p / q;
                fpp = q / (x12 * x31 * x32);
            }

        }

    }
        

    // This class is used to wrap a function, storing some associated state such as the evaluation count.
    // It isn't truely a functor in the C++ sense, since .NET doesn't allow () to be overloaded, but
    // it is a functor in the sense that it is a class used to represent a function.

    internal class Functor {

        public Functor (Func<double, double> f) : this(f, false) {}

        public Functor (Func<double, double> f, bool negate) {
            Function = f;
            Negate = negate;
        }

        public Func<double, double> Function { get; protected set; }

        public int EvaluationCount { get; protected set; }

        public bool Negate { get; set; }

        public virtual double Evaluate (double x) {

            double f = Function(x);
            if (Negate) f = -f;
            EvaluationCount++;
            return (f);

        }

    }

}
