using System;
using System.Collections.Generic;
using System.Diagnostics;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using System.Runtime.Remoting.Messaging;
using System.Reflection.Emit;

namespace Test {

    public class RootResult {

        public RootResult (double r, int n) {
            this.Root = r;
            this.EvaluationCount = n;
        }

        public double Root { get; private set; }

        public int EvaluationCount { get; private set; }

    }

    public readonly struct ValueAndDerivative {

        public ValueAndDerivative (double value, double derivative) {
            this.Value = value;
            this.Derivative = derivative;
        }

        public double Value { get; }

        public double Derivative { get; }

    }

    public class TestFunctor {

        public TestFunctor(Func<double, double> f, bool negate) {
            this.f = f;
            this.negate = negate;
        }

        private readonly Func<double, double> f;
        private readonly bool negate;

        public double Evaluate (double x) {
            double y = f(x);
            return negate ? -y : y;
        }
    }

    [TestClass]
    public class FutureRootTest {

        [TestMethod]
        public  void AiryTiming () {

            Console.WriteLine(AdvancedMath.AiryAiZero(1));

            double za = 0.0;
            Stopwatch sa = Stopwatch.StartNew();
            for(int k = 1; k < 1000; k++) {
                za += AdvancedMath.AiryAiZero(k);
            }
            sa.Stop();
            Console.WriteLine(sa.ElapsedMilliseconds);

            double zb = 0.0;
            Stopwatch sb = Stopwatch.StartNew();
            for (int k = 1; k < 1000; k++) {
                zb += AdvancedMath.AiryBiZero(k);
            }
            sb.Stop();
            Console.WriteLine(sb.ElapsedMilliseconds);


        }

        [TestMethod]
        public void FunctorVsLambda () {

            Func<double, double> f = x => x * x;
            TestFunctor t = new TestFunctor(f, true);
            Func<double, double> u = x => -f(x);

            double yt = 0.0;
            Stopwatch st = Stopwatch.StartNew();
            for (double x = 1.0E-1; x < 10.0; x += 0.00001) {
                yt += t.Evaluate(x);
            }
            st.Stop();
            Console.WriteLine(st.ElapsedMilliseconds);

            double yu = 0.0;
            Stopwatch su = Stopwatch.StartNew();
            for (double x = 1.0E-1; x < 10.0; x += 0.00001) {
                yu += u(x);
            }
            su.Stop();
            Console.WriteLine(su.ElapsedMilliseconds);

            Console.WriteLine($"{yt} {yu}");

        }

        private static double tol = 1.0E-14;

        public RootResult Newton (Func<double,ValueAndDerivative> f, double a, double b) {

            ValueAndDerivative fa = f(a);
            if (fa.Value == 0.0) return new RootResult(a, 1);

            ValueAndDerivative fb = f(b);
            if (fb.Value == 0.0) return new RootResult(b, 2);

            if (Math.Sign(fa.Value) == Math.Sign(fb.Value)) throw new InvalidOperationException();

            for (int i = 0; i < 200; i++) {

                // Make a the smaller value
                if (Math.Abs(fb.Value) < Math.Abs(fa.Value)) {
                    double t = a; a = b; b = t;
                    ValueAndDerivative ft = fa; fa = fb; fb = ft;
                }

                double ell = b - a;
                if (Math.Abs(ell) < tol) return new RootResult(a, 2 + i);

                double d = 0.0;
                if (Math.Sign(fa.Value) == Math.Sign(-ell * fa.Derivative) && Math.Abs(fa.Value) < Math.Abs(ell * fa.Derivative)) {
                    d = -fa.Value / fa.Derivative;
                    if (Math.Abs(d) < 0.25 * tol) d = 0.25 * tol * Math.Sign(ell);
                }

                if (a + d == a) {
                    d = 0.5 * ell;
                }

                double x = a + d;
                ValueAndDerivative fx = f(x);
                if (fx.Value == 0.0) return new RootResult(b, 3 + i);

                if (Math.Sign(fx.Value) == Math.Sign(fa.Value)) {
                    a = x;
                    fa = fx;
                } else {
                    Debug.Assert(Math.Sign(fx.Value) == Math.Sign(fb.Value));
                    b = x;
                    fb = fx;
                }

            }

            throw new NonconvergenceException();

        }

        private IReadOnlyList<Func<Func<double, ValueAndDerivative>, double, double, RootResult>> derivativeSolvers = new List<Func<Func<double, ValueAndDerivative>, double, double, RootResult>> {
        };

        private IReadOnlyList<(string, Func<double, ValueAndDerivative>, double, double)> derivativeTestCases = new List<(string, Func<double, ValueAndDerivative>, double, double)> {
            ("Cosine", (double x) =>new ValueAndDerivative(Math.Cos(x), -Math.Sin(x)), 0.0, 2.0),
            ("LogGamma", (double x) => new ValueAndDerivative(AdvancedMath.LogGamma(x), AdvancedMath.Psi(x)), 0.1, 1.6),
            ("Psi", (double x) => new ValueAndDerivative(AdvancedMath.Psi(x), AdvancedMath.Psi(1, x)), 0.1, 3.0),
            ("Ei", (double x) => new ValueAndDerivative(AdvancedMath.IntegralEi(x), Math.Exp(x) / x ), 0.01, 1.0),
            ("J0", (double x) => new ValueAndDerivative(AdvancedMath.BesselJ(0, x), -AdvancedMath.BesselJ(1, x)), 0.4, 4.0),
            //("W-", (double w) => w * Math.Exp(w) + 0.25, -1.0 / Math.E, 2.0), // This has very large slope close to endpoint.
            //("W+", (double w) => w * Math.Exp(w) - 1.0, -1.0 / Math.E, 2.0),
            //("P1", (double x) => x * (3.0 + x * (-4.0 + x * 2.0)), -1.0, +1.0), // This polynomial is known to cause naive false position to misbehave.
            ("Ci", (double x) => new ValueAndDerivative(AdvancedMath.IntegralCi(x), Math.Cos(x) / x), 0.03, 3.0),
            ("Flat", (double x) => (x == 0.0) ? new ValueAndDerivative(0.0, 0.0) : new ValueAndDerivative(x * Math.Exp(-1.0 / Math.Abs(x)), Math.Exp(-1.0 / Math.Abs(x) * (1.0 + Math.Abs(x)) / Math.Abs(x))), -0.6, 1.5), // Ridiculously flat function similiar to example by Brent
            ("x9", (double x) => new ValueAndDerivative(MoreMath.Pow(x, 9), 9.0 * MoreMath.Pow(x, 8)), -1.0, 4.0), // Test function used by Brent
            //("Step", (double x) => (x < 0.0) ? -1.0 : +1.0, -3.0, +2.0),
            //("Well", (double u) => Math.Sqrt(1.0 - u * u) - u * Math.Tan(u), 0.0, 1.0), // QM Square Well Bound State
            //("Near", (double x) => (x - 0.1) * (x * x + 0.01), -1.5, +1.0) // Cubic with almost double root near real root
        };

        [TestMethod]
        public void DerivativeRootFinderTest() {

            Func<Func<double, ValueAndDerivative>, double, double, RootResult> derivativeSolver = Newton;

            //foreach (Func<Func<double, double>, double, double, RootResult> solver in solvers) {

                Console.WriteLine(derivativeSolver.Method.Name);

                int count = 0;
                foreach ((string name, Func<double, ValueAndDerivative> f, double a, double b) in derivativeTestCases) {
                    try {
                        RootResult r = derivativeSolver(f, a, b);
                        Console.WriteLine($"  {name} {r.Root} {r.EvaluationCount}");
                        count += r.EvaluationCount;
                    } catch (NonconvergenceException) {
                        Console.WriteLine($"  {name} Noncovergence");
                    }
                }
                Console.WriteLine(count);

            //}
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

        public static RootResult Brent (Func<double, double> f, double a, double b) {

            double fa = f(a);
            if (fa == 0.0) return new RootResult(a, 1);

            double fb = f(b);
            if (fb == 0.0) return new RootResult(b, 2);

            if (Math.Sign(fa) == Math.Sign(fb)) throw new InvalidOperationException();

            // At start, without a 3rd point, we can't use inverse quadratic interpolation.
            // We could use bisection or secant. Secant would likely be better if the bracket
            // were already very good, but for the wide, human-guessed brackets of my test cases,
            // it was worse, so choose bisection.
            double x = 0.5 * (a + b);
            double fx = f(x);
            if (fx == 0.0) return new RootResult(x, 3);

            // We will track the last and next-to-last points evaluated.
            double fy = Math.Abs(fa) > Math.Abs(fb) ? fa : fb;

            // We will track the last and next-to-last steps.
            double d = a - b;
            double d1 = d;

            // We will track 3 points. A bracket [a, b] and a 3rd point c outside the bracket.
            double c, fc;

            for (int i = 0; i < 200; i++) {

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
                double ell = b - a;
                if (Math.Abs(ell) < tol) return new RootResult(a, 3 + i);

                // We will not compute the next step.
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
                if (a + d == a ) {
                    d = 0.5 * ell;
                }

                x = a + d;

                fy = fx;
                fx = f(x);
                if (fx == 0.0) return new RootResult(x, 4 + i);

            }

            throw new NonconvergenceException();
        }

        public static RootResult BrentCPP (Func<double, double> f, double a, double b) {

            // https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.cpp

            double fa = f(a);
            double fb = f(b);

            double c = a;
            double fc = fa;
            double e = b - a;
            double d = e;

            for (int i = 0; i < 200 ; i++ )
            {
                if (Math.Abs(fc) < Math.Abs(fb)) {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }

                double m = 0.5 * (c - b);

                if (Math.Abs(m) <= tol || fb == 0.0) {
                    return new RootResult(b, 2 + i);
                }

                if (Math.Abs(e) < tol || Math.Abs(fa) <= Math.Abs(fb)) {
                    e = m;
                    d = e;
                } else {
                    double s = fb / fa;

                    double p, q;
                    if (a == c) {
                        p = 2.0 * m * s;
                        q = 1.0 - s;
                    } else {
                        q = fa / fc;
                        double r = fb / fc;
                        p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                    }

                    if (0.0 < p) {
                        q = -q;
                    } else {
                        p = -p;
                    }

                    s = e;
                    e = d;

                    if (2.0 * p < 3.0 * m * q - Math.Abs(tol * q) &&
                      p < Math.Abs(0.5 * s * q)) {
                        d = p / q;
                    } else {
                        e = m;
                        d = e;
                    }
                }
                a = b;
                fa = fb;

                if (tol < Math.Abs(d)) {
                    b = b + d;
                } else if (0.0 < m) {
                    b = b + tol;
                } else {
                    b = b - tol;
                }

                fb = f(b);

                if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)) {
                    c = a;
                    fc = fa;
                    e = b - a;
                    d = e;
                }
            }

            throw new NonconvergenceException();

        }

        public static RootResult BrentWiki (Func<double, double> f, double a, double b) {

            double fa = f(a);
            double fb = f(b);

            if (fa * fb > 0.0) throw new InvalidOperationException();
           
            double c = a;
            double fc = fa;

            bool mflag = true;
            double d = 0.0;

            for (int i = 0; i < 200; i++) {

                if (Math.Abs(fa) < Math.Abs(fb)) {
                    double t = a; a = b; b = t;
                    t = fa; fa = fb; fb = t;
                }

                if (fb == 0.0 || Math.Abs(b - a) <= tol) return new RootResult(b, 2 + i);

                double s;
                if ((fa != fc) && (fb != fc)) {
                    s = a * fb * fc / (fa - fb) / (fa - fc) +
                        b * fa * fc / (fb - fa) / (fb - fc) +
                        c * fa * fb / (fc - fa) / (fc - fb);
                } else {
                    s = b - fa * (b - a) / (fb - fa);
                }

                Interval r = Interval.FromEndpoints(b, (3.0 * a + b) / 4.0);
                if (!r.Contains(s) || 
                    (mflag && Math.Abs(s-b) >= Math.Abs(b-c) / 2.0) ||
                    (!mflag && Math.Abs(s - b) >= Math.Abs(c - d) / 2.0) ||
                    (mflag && Math.Abs(b-c) < tol) ||
                    (!mflag && Math.Abs(c-d) < tol)) {
                    s = (a + b) / 2.0;
                    mflag = true;
                } else {
                    mflag = false;
                }

                double fs = f(s);

                d = c;
                c = b;

                if (fa * fs < 0.0) {
                    b = s;
                } else {
                    a = s;
                }

            }

            throw new NonconvergenceException();

        }

        public static RootResult BrentNR (Func<double, double> f, double a, double b) {

            double fa = f(a);
            double fb = f(b);

            if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) throw new InvalidOperationException();

            double c = b;
            double fc = fb;

            double d = 0.0;
            double e = 0.0;

            for (int i = 0; i < 100; i++) {

                if ((fb > 0.0) && (fc > 0.0) || (fb < 0.0) && (fc < 0.0)) {
                    c = a;
                    fc = fa;
                    d = b - a;
                    e = d;
                }
                if (Math.Abs(fc) < Math.Abs(fb)) {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }
                double xm = 0.5 * (c - b);
                if (Math.Abs(xm) <= tol || fb == 0.0) return new RootResult(b, 2 + i);
                if (Math.Abs(e) >= tol && Math.Abs(fa) > Math.Abs(fb)) {
                    double p, q;
                    double s = fb / fa;
                    if (a == c) {
                        p = 2.0 * xm * s;
                        q = 1.0 - s;
                    } else {
                        q = fa / fc;
                        double r = fb / fc;
                        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                    }
                    if (p > 0.0) q = -q;
                    p = Math.Abs(p);
                    double min1 = 3.0 * Math.Abs(xm * q) - Math.Abs(tol * q);
                    double min2 = Math.Abs(e * q);
                    if (2.0 * p < Math.Min(min1, min2)) {
                        e = d;
                        d = p / q;
                    } else {
                        d = xm;
                        e = d;
                    }
                } else {
                    d = xm;
                    e = d;
                }
                a = b;
                fa = fb;
                if (Math.Abs(d) > tol) {
                    b += d;
                } else {
                    b += Math.Sign(xm) * tol;
                }
                fb = f(b);

            }

            throw new NonconvergenceException();
        }

        public static RootResult Ridder (Func<double, double> f, double a, double b) {

            double fa = f(a);
            if (fa == 0.0) return new RootResult(a, 1);

            double fb = f(b);
            if (fb == 0.0) return new RootResult(b, 2);

            if (Math.Sign(fa) == Math.Sign(fb)) throw new InvalidOperationException();

            for (int i = 0; i < 200; i += 2) {

                double m = 0.5 * (a + b);
                if (Math.Abs(a - b) < tol) return new RootResult(m, 2 + i);

                double fm = f(m);
                if (fm == 0.0) return new RootResult(m, 3 + i);

                double s = Math.Sqrt(fm * fm - fa * fb);
                Debug.Assert(s > 0.0); // Since fm * fm > 0 and fa * fb < 0, this must be non-zero.
                double x = m + (m - a) * ((fa > fb) ? fm : -fm ) / s; // Check that a <-> b gives same result. (Yes, same, but few ulp jiggle.)
                Debug.Assert((a <= x && x <= b) || (b <= x && x <= a));

                double fx = f(x);
                if (fx == 0.0) return new RootResult(x, 4 + i);

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

            throw new NonconvergenceException();
        }

        public static RootResult Secant(Func<double, double> f, double a, double b) {


            double fa = f(a);
            if (fa == 0.0) return new RootResult(a, 1);
            double u = a;
            double fu = fa;

            double fb = f(b);
            if (fb == 0.0) return new RootResult(b, 2);
            double v = b;
            double fv = fb;

            if (Math.Sign(fa) == Math.Sign(fb)) throw new InvalidOperationException();

            for (int i = 0; i < 200; i++) {

                double w = (v * fu - u * fv) / (fu - fv);

                // The naive seceant method doesn't maintain a bracket at all. To maintain
                // a bracket, we fall back to the midpoint if the secant would take us
                // outside the bracket. (Should we fall back to false position instead?)
                if (!((a < w && w < b) || (b < w && w < a))) {
                    w = 0.5 * (a + b);
                }

                if (Math.Abs(a - b) <= tol) return new RootResult(w, 2 + i);

                double fw = f(w);
                if (fw == 0.0) return new RootResult(w, 3 + i);

                // Update extrapolation point set
                u = v;
                fu = fv;
                v = w;
                fv = fw;

                // Update bracket
                if (Math.Sign(fw) == Math.Sign(fa)) {
                    a = w;
                    fa = fw;
                } else {
                    Debug.Assert(Math.Sign(fw) == Math.Sign(fb));
                    b = w;
                    fb = fw;
                }

            }

            throw new NonconvergenceException();
        }

        public static RootResult FalsePosition(Func<double, double> f, double a, double b) {

            double fa = f(a);
            if (fa == 0.0) return new RootResult(a, 1);

            double fb = f(b);
            if (fb == 0.0) return new RootResult(b, 2);

            if (Math.Sign(fa) == Math.Sign(fb)) throw new InvalidOperationException();

            int lastSign = 0;
            for (int i = 0; i < 200; i++) {

                double x = (b * fa - a * fb) / (fa - fb);
                //if (x == a || x == b) x = 0.5 * (a + b);
                Debug.Assert((a <= x && x <= b) || (b <= x && x <= a));
                if (Math.Abs(a - b) <= tol) return new RootResult(x, 2 + i);

                double fx = f(x);
                if (fx == 0.0) return new RootResult(x, 3 + i);

                // The "Illinois" improvement to the false position algorithm attempts to
                // address the situation that we keep updating only one side of the bracket.
                // This very much does occur in real life. To avoid it, if we update one side
                // twice, we half the recorded function value of the other side, to move the
                // next test point closer to that side.
                int sign = Math.Sign(fx);
                if (sign == Math.Sign(fa)) {
                    a = x;
                    fa = fx;
                    if (sign == lastSign) fb /= 2.0;
                } else {
                    Debug.Assert(sign == Math.Sign(fb));
                    b = x;
                    fb = fx;
                    if (sign == lastSign) fa /= 2.0;
                }
                lastSign = sign;

            }

            throw new NonconvergenceException();

        }

        public static RootResult Bisection(Func<double, double> f, double a, double b) {

            double fa = f(a);
            if (fa == 0.0) return new RootResult(a, 1);

            double fb = f(b);
            if (fb == 0.0) return new RootResult(b, 2);

            if (Math.Sign(fa) == Math.Sign(fb)) throw new InvalidOperationException();

            for (int i = 0; i < 200; i++) {

                Debug.Assert(Math.Sign(fa) != Math.Sign(fb));

                double m = 0.5 * (a + b);
                Debug.Assert((a < m && m < b) || (b < m && m < a));
                if (Math.Abs(a - b) <= tol) return new RootResult(m, 2 + i);

                double fm = f(m);
                if (fm == 0.0) return new RootResult(m, 3 + i);

                if (Math.Sign(fm) == Math.Sign(fa)) {
                    a = m;
                    fa = fm;
                } else {
                    Debug.Assert(Math.Sign(fm) == Math.Sign(fb));
                    b = m;
                    fb = fm;
                }

            }

            throw new NonconvergenceException();

        }

        private IReadOnlyList<Func<Func<double, double>, double, double, RootResult>> solvers = new List<Func<Func<double, double>, double, double, RootResult>> {
            Bisection, FalsePosition, Ridder, Secant, BrentNR, BrentCPP, Brent
        };

        private IReadOnlyList<(string, Func<double, double>, double, double)> testCases = new List<(string, Func<double, double>, double, double)> {
            ("Cosine", (double x) => Math.Cos(x), 0.0, 2.0),
            ("LogGamma", (double x) => AdvancedMath.LogGamma(x), 0.1, 1.6),
            ("Psi", (double x) => AdvancedMath.Psi(x), 0.1, 3.0),
            ("Ei", (double x) => AdvancedMath.IntegralEi(x), 0.01, 1.0),
            ("J0", (double x) => AdvancedMath.BesselJ(0, x), 0.4, 4.0),
            ("W-", (double w) => w * Math.Exp(w) + 0.25, -1.0 / Math.E, 2.0), // This has very large slope close to endpoint.
            ("W+", (double w) => w * Math.Exp(w) - 1.0, -1.0 / Math.E, 2.0),
            ("P1", (double x) => x * (3.0 + x * (-4.0 + x * 2.0)), -1.0, +1.0), // This polynomial is known to cause naive false position to misbehave.
            ("Ci", (double x) => AdvancedMath.IntegralCi(x), 0.03, 3.0),
            ("Flat", (double x) => (x == 0.0) ? 0.0 : x * Math.Exp(-1.0 / Math.Abs(x)), -0.6, 1.5), // Ridiculously flat function similiar to example by Brent
            ("x9", (double x) => MoreMath.Pow(x, 9), -1.0, 4.0), // Test function used by Brent
            ("Step", (double x) => (x < 0.0) ? -1.0 : +1.0, -3.0, +2.0),
            ("Well", (double u) => Math.Sqrt(1.0 - u * u) - u * Math.Tan(u), 0.0, 1.0), // QM Square Well Bound State
            ("Near", (double x) => (x - 0.1) * (x * x + 0.01), -1.5, +1.0) // Cubic with almost double root near real root
        };

        [TestMethod]
        public void RootFinderTest () {

            foreach (Func<Func<double, double>, double, double, RootResult> solver in solvers) {

                Console.WriteLine(solver.Method.Name);

                int count = 0;
                foreach ((string name, Func<double, double> f, double a, double b) in testCases) {
                    try {
                        RootResult r = solver(f, a, b);
                        Console.WriteLine($"  {name} {r.Root} {r.EvaluationCount}");
                        count += r.EvaluationCount;
                    } catch (NonconvergenceException) {
                        Console.WriteLine($"  {name} Noncovergence");
                    }
                }
                Console.WriteLine(count);

            }
        }


    }
}
