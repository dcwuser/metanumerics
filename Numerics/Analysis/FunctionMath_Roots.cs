using System;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Analysis {

    public static partial class FunctionMath {

        /// <summary>
        /// Isolates a root in the vicinity of a given point.
        /// </summary>
        /// <param name="f">The function whoose zero is sought.</param>
        /// <param name="x">A ordinate believed to be near the sought zero.</param>
        /// <returns>An ordinate at which the function has a zero.</returns>
        public static double FindZero (Func<double, double> f, double x) {

            if (f is null) throw new ArgumentNullException(nameof(f));

            ZeroSettings settings = new ZeroSettings() { RelativePrecision = Global.Accuracy, AbsolutePrecision = Global.Accuracy * Global.Accuracy, EvaluationBudget = 200 };
            Functor ftr = new Functor(f) { EvaluationBudget = settings.EvaluationBudget };

            // take a step
            double fx = ftr.Evaluate(x);
            if (fx == 0.0) return x;
            double dx = (Math.Abs(x) + 1.0 / 16.0) / 16.0;

            // keep stepping until we change sign
            for (int n = 0; n < settings.EvaluationBudget; n++) {

                double y = x + dx;
                double fy = ftr.Evaluate(y);
                if (fy == 0.0) return y;

                // if we the function changed sign, we have bracketed a root
                if (Math.Sign(fy) != Math.Sign(fx)) {
                    ZeroResult zr = FindZero_Ridder(ftr, x, fx, y, fy, settings);
                    return zr.Value;
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

            if (f is null) throw new ArgumentNullException(nameof(f));

            ZeroSettings settings = new ZeroSettings() { RelativePrecision = Global.Accuracy, AbsolutePrecision = Global.Accuracy * Global.Accuracy, EvaluationBudget = 200 };
            Functor ftr = new Functor(f) { EvaluationBudget = settings.EvaluationBudget };

            double a = bracket.LeftEndpoint;
            double fa = ftr.Evaluate(a);
            if (fa == 0.0) return a;

            double b = bracket.RightEndpoint;
            double fb = ftr.Evaluate(b);
            if (fb == 0.0) return b;

            // Make sure given interval really brackets the root.
            if (Math.Sign(fa) == Math.Sign(fb)) throw new InvalidOperationException();

            ZeroResult zr = FindZero_Ridder(ftr, a, fa, b, fb, settings);
            //ZeroResult zr = FindZero_Brent(ftr, a, fa, b, fb, settings);
            return zr.Value;

        }

        internal static double FindZero (Func<double, double> f, double a, double b, double x) {

            Debug.Assert(f is object);

            ZeroSettings settings = new ZeroSettings() { RelativePrecision = Global.Accuracy, AbsolutePrecision = Global.Accuracy * Global.Accuracy, EvaluationBudget = 200 };
            Functor ftr = new Functor(f) { EvaluationBudget = settings.EvaluationBudget };

            double fa = ftr.Evaluate(a);
            if (fa == 0.0) return a;

            double fb = ftr.Evaluate(b);
            if (fb == 0.0) return b;

            Debug.Assert(Math.Sign(fa) != Math.Sign(fb));

            Debug.Assert((a < x && x < b) || (b < x && x < a));

            //ZeroResult zr = FindZero_Brent(ftr, a, fa, b, fb, x, settings);
            ZeroResult zr = FindZero_Ridder(ftr, a, fa, b, fb, settings);
            return zr.Value;

        }

        private static ZeroResult FindZero_Ridder (Functor f, double a, double fa, double b, double fb, ZeroSettings settings) {

            Debug.Assert(f is object);
            Debug.Assert(settings is object);
            Debug.Assert(Math.Sign(fa) != Math.Sign(fb));

            for (int i = 0; i < settings.EvaluationBudget; i += 2) {

                double m = 0.5 * (a + b);
                double tol = settings.ComputePrecision(m);
                if (Math.Abs(a - b) < tol) return new ZeroResult(m, a, b, f.EvaluationCount);

                double fm = f.Evaluate(m);
                if (fm == 0.0) return new ZeroResult(m, m, m, f.EvaluationCount);

                // Since fm * fm > 0 and fa * fb < 0, fm * fm - fa * fb should be positive and its root non-zero
                // But it is possible that all are so small that their squares are zero, in which case return.
                double s = Math.Sqrt(fm * fm - fa * fb);
                if (s == 0.0) return new ZeroResult(m, a, b, f.EvaluationCount);
                double x = m + (m - a) * ((fa > fb) ? fm : -fm) / s; // a <-> b gives same result, with a couple ulp jiggle.
                //Debug.Assert((a <= x && x <= b) || (b <= x && x <= a));

                double fx = f.Evaluate(x);
                if (fx == 0.0) return new ZeroResult(x, x, x, f.EvaluationCount);

                int sign = Math.Sign(fx);
                if (sign == Math.Sign(fa)) {
                    a = x;
                    fa = fx;
                    if (sign != Math.Sign(fm)) {
                        b = m;
                        fb = fm;
                    }
                } else {
                    Debug.Assert(sign == Math.Sign(fb));
                    b = x;
                    fb = fx;
                    if (sign != Math.Sign(fm)) {
                        a = m;
                        fa = fm;
                    }
                }

            }

            /*
            // x1 and x2 are bracket points, with function values f1 and f2

            // x0 is the current best estimate of the root; it will be updated in the first loop
            double x0 = Double.MaxValue;

            for (int n = 0; n < Global.SeriesMax; n++) {

                // evaluate at the bracket mid-point
                double x3 = 0.5 * (a + b);
                double f3 = f.Evaluate(x3);
                if (f3 == 0.0) return new ZeroResult(x3, x3, x3, f.EvaluationCount);

                // find the Ridder point...
                double x4;
                double q = Math.Sqrt(f3 * f3 - fa * fb);
                if (q == 0.0) return new ZeroResult(x3, x3, x3, f.EvaluationCount); // Is this right?
                double dx = (x3 - a) * f3 / q;
                if (fa > fb) {
                    x4 = x3 + dx;
                } else {
                    x4 = x3 - dx;
                }

                // ...and evaluate there
                double f4 = f.Evaluate(x4);
                if (f4 == 0.0) return new ZeroResult(x4, x4, x4, f.EvaluationCount);

                // stopping criteria
                double tol = settings.ComputePrecision(x4);
                if (Math.Abs(x4 - x0) <= tol) return new ZeroResult(x0, a, b, f.EvaluationCount);
                //if (Math.Abs(x4 - x0) <= Global.Accuracy * Math.Abs(x0)) return (x4);
                //if (Math.Abs(x4 - x0) <= AbsolutePrecision) return (x4);

                // update the bracket
                if (Math.Sign(f3) == Math.Sign(f4)) {
                    // x3 and x4 are on the same side of the root
                    // replace the same-signed bracket point with x4
                    if (Math.Sign(f4) == Math.Sign(fa)) {
                        a = x4;
                        fa = f4;
                    } else {
                        b = x4;
                        fb = f4;
                    }
                } else {
                    // x3 and x4 are on oppisite sides of the root
                    // use them to replace the bracket
                    a = x3;
                    fa = f3;
                    b = x4;
                    fb = f4;
                }

                // add end criteria to see whether we have moved from the last iteration
                x0 = x4;

            }
            */
            throw new NonconvergenceException();

        }


        // This basically the Brent algorithm, but modified slightly based on evaluation counts of our test functions.
        // Richard Brent, "Algorithms for Minimization Without Derivatives", 1973, Chapter 4
        // Brent, R. P., "An Algorithm With Guaranteed Convergence for Finding a Zero of a Function",
        // The Computer Journal, 14 (1971) 422 (https://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)

        // Up until v5 we used Ridder's algorithm, which is simpler but, according to our test cases,
        // not as efficient: 366 vs 286 function evaluations for all test functions.

        // One way we differ from the original Brent algorithm is that it throws away the third point
        // and falls back to secant interpolation using two points rather easily. We don't.
        // An implementation I found on the web (https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.cpp)
        // took 339 evaluations for our test functions, so our changes do appear to work.

        private static ZeroResult FindZero_Brent (Functor f, double a, double fa, double b, double fb, double? z, ZeroSettings settings) {

            Debug.Assert(f is object);
            Debug.Assert(Math.Sign(fa) != Math.Sign(fb));
            Debug.Assert(settings is object);

            // At start, without a 3rd point, we can't use inverse quadratic interpolation.
            // We could use bisection or secant. Secant would likely be better if the bracket
            // were already very good, but for the wide, human-guessed brackets of my test cases,
            // it was worse, so choose bisection.
            double x = z.HasValue ? z.Value : 0.5 * (a + b);
            double fx = f.Evaluate(x);
            if (fx == 0.0) return new ZeroResult(x, a, b, f.EvaluationCount);

            // We will track the last and next-to-last points evaluated.
            double fy = Math.Abs(fa) > Math.Abs(fb) ? fa : fb;

            // We will track the last and next-to-last steps.
            double d = a - b;
            double d1 = d;

            // We will track 3 points. A bracket [a, b] and a 3rd point c outside the bracket.
            double c, fc;

            for (int i = 0; i < settings.EvaluationBudget; i++) {

                // Choose new points a, b, and c according to the sign of fx.
                // We always discard the old c, which will lie outside both the old and the new bracket.
                if (Math.Sign(fx) == Math.Sign(fa)) {
                    c = a;
                    fc = fa;
                    a = x;
                    fa = fx;
                } else {
                    Debug.Assert(Math.Sign(fx) == Math.Sign(fb));
                    c = b;
                    fc = fb;
                    b = x;
                    fb = fx;
                }

                // The bracket is always [a, b].
                Debug.Assert(Math.Sign(fa) != Math.Sign(fb));

                // Third point c always lies outside the bracket (otherwise it would be part of bracket).
                // (Can we get rid of equals? I don't think we can force yet.)
                Debug.Assert(c <= Math.Min(a, b) || c >= Math.Max(a, b));

                // Make a the side of the bracket with a function value closest to zero.
                if (Math.Abs(fb) < Math.Abs(fa)) {
                    double t = a; a = b; b = t;
                    t = fa; fa = fb; fb = t;
                }

                // a should now have the smaller function value.
                Debug.Assert(Math.Abs(fa) <= Math.Abs(fb));

                // Note that Brent and most implementations use [b,c] for the bracket and b for the
                // point with smallest value. This unnatural notation drove me crazy, so I changed
                // it, but my naming won't agree with other implementors.

                // Check for convergence.
                double tol = settings.ComputePrecision(a);
                double ell = b - a;
                if (Math.Abs(ell) < tol) return new ZeroResult(a, a, b, f.EvaluationCount);

                // Report progress if requested.
                if (settings.Listener is object) settings.Listener(new ZeroResult(a, a, b, f.EvaluationCount));

                // We will now compute the next step.
                // Shift previous steps to the appropriate registers.
                double d2 = d1;
                d1 = d;
                d = 0.0;

                // We hope to use inverse quadratic interpolation to guess a step that will take us very close to the root.
                // In order to use the interpolate, we require a few conditions:
                //   The function values are decreasing.
                //   The interpolate is within the bracket.
                //   The interpolate is not unexpectedly far toward the wrong side of the bracket.
                //   The interpolate steps are decreasing by at least a factor 2 per two steps.
                //   (The interpolate step is not too small. I tried with a 1/4 tol minimum, but many functions needed more steps.)

                // Brent forces an occasional bisection with the criteria |d2| > tol, which is eventually violated if
                // convergence is good. I am somewhat skeptical of this idea: why force a bisection if interpolation is working well?
                // But it makes no difference for any of my test functions so for now I leave it in. One point in my favor: if
                // you change it to |d1| > tol, so it is triggered more often, it makes us do worse for most test functions
                // (although slightly better for y = x^9).

                if (Math.Abs(fx) < Math.Abs(fy) && Math.Abs(d2) > tol) {

                    double r = fa / fb;
                    double s = fa / fc;
                    double t = fc / fb;

                    double p = s * (t * (r - t) * ell - (1.0 - r) * (a - c));
                    double q = (r - 1.0) * (s - 1.0) * (t - 1.0);

                    // It is important that these conditions not allow equality, because that allows p = q = 0 => d = NaN.
                    if (Math.Sign(p) == Math.Sign(ell * q) && 2.0 * Math.Abs(p) < 3.0 * Math.Abs(ell * q) && Math.Abs(p) < 0.5 * Math.Abs(d2 * q)) {
                        d = p / q;
                    }
                }

                // If interpolation didn't work, fall back to bisection.
                // Testing a + d == a instead of d == 0 is a way to catch non-zero d's so small that they do not
                // move the evaluation point. Before this change, I found cases where we would repeatedly cycle between
                // interpolations that would not move the point and bisections that would, ultimately doing
                // twice as many evaluations once we hit that phase of the search.
                if (a + d == a) {
                    d = 0.5 * ell;
                }

                x = a + d;

                fy = fx;
                fx = f.Evaluate(x);
                if (fx == 0.0) return new ZeroResult(x, x, x, f.EvaluationCount);

            }

            throw new NonconvergenceException();
        }


    }

}
