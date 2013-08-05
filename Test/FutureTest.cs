using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.ObjectModel;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.SignalProcessing;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace FutureTest {

    [TestClass]
    public class FutureTest {

        [TestMethod]
        public void ET () {

            SquareMatrix A = new SquareMatrix(4);
            A[0, 0] = 1.0; A[0, 1] = 1.0; A[0, 2] = 1.0; A[0, 3] = 1.0;
            A[1, 0] = 1.0; A[1, 1] = 1.0; A[1, 2] = 1.0; A[1, 3] = 1.0;
            A[2, 0] = 0.0; A[2, 1] = 0.0; A[2, 2] = 1.0; A[2, 3] = 1.0;
            A[3, 0] = 0.0; A[3, 1] = 0.0; A[3, 2] = 1.0; A[3, 3] = 1.0;
            A.Eigensystem();
        }

        // We want to find the rotation that brings a 2 X 2 matrix into triangular form.
        //   (  c  s ) ( a11  a12 ) ( c  -s )
        //   ( -s  c ) ( a21  a22 ) ( s   c )
        // Multiplying out gives
        //   ( c^2 a11 + cs a21 + cs a12 + s^2 a22  c^2 a12 - cs a11 + cs a22 - s^2 a21 )
        //   ( c^2 a21 - cs a11 + cs a22 - s^2 a12  c^2 a22 - cs a21 - cs a12 + s^2 a11 )

        // We want a21' = 0. Since the rotation must vanish when a21 = 0 and grow as a21 grows, make the Ansatz that s ~ a21, i.e.
        //   s = \frac{a_{21}}{\sqrt{a_{21}^2 + b^2}}  c = \frac{b}{\sqrt{a_{21}^2 + b^2}}
        // It is then straightforward to derive
        //   b = \frac{a_{11} - a_{22} \pm q}{2}
        // where q^2 = ( a_{11} - a_{22} )^2 + 4 a12 a21 is the same descriminant as appears in the eigenvalue problem, so this rotation
        // exists iff the eigenvalues are real.

        private void TwoByTwoSchur (ref double a11, ref double a12, ref double a21, ref double a22, out double s, out double c) {

            // compute some quantities we will use 
            double u = a11 + a22;
            double v = a11 - a22;
            double w = a11 * a22 - a12 * a21;
            double q2 = v * v + 4.0 * a12 * a21;
            // note u is the trace and w is the determinant

            if (q2 >= 0.0) {

                // the descriminant is positive so the eigenvalues are real
                double q = Math.Sqrt(q2);
 
                // find the rotation sets a21' = 0
                // in the equation for b, choose the sign so as to minimize cancelation
                double b = (v >= 0.0) ? (v + q) / 2.0 : (v - q) / 2.0;
                double rho = MoreMath.Hypot(a21, b);
                s = a21 / rho;
                c = b / rho;

                // Note that a12' - a21' = a12 - a21, and since a21' = 0, a12' = a12 - a21
                a12 = a12 - a21;
                a21 = 0.0;

                // the eigenvalues are (u \pm q) / 2
                // we avoid cancelation by computing the non-canceling one first and
                // computing the other using the fact that their product equals the determinant
                if (u >= 0.0) {
                    a11 = (u + q) / 2.0;
                    a22 = w / a11;
                } else {
                    a22 = (u - q) / 2.0;
                    a11 = w / a22;
                }
                // we have placed the u + q eigenvalue in the a11 slot. if v >= 0, that is where the rotation puts it
                // if v < 0, the u - q eigenvalue belongs these, so we need to swap them
                if (v < 0) { double t = a11; a11 = a22; a22 = t; }

            } else {

                // In the q2 < 0 case, we can't zero a21. But we can rotate so a11' = a22'

                double r = - (a12 + a21) / v;
                double t = Math.Sign(r) / (Math.Abs(r) + MoreMath.Hypot(1.0, r));
                c = 1.0 / MoreMath.Hypot(1.0, t);
                s = t * c;

                // Since rotations preserve the trace, the equal diagonal elements must equal the average of the previous elements 
                a11 = u / 2.0;
                a22 = a11;

                //double q = Math.Sqrt(-q2);
                //e1 = new Complex(u / 2.0, q / 2.0);
                //e2 = e1.Conjugate;

                //s = 0.0;
                //c = 1.0;

            }
        }

        [TestMethod]
        public void EigenTest () {
            
            Random rng = new Random(1);
            for (int i = 0; i < 100; i++) {

                double a00 = 1.0 - 2.0 * rng.NextDouble();
                double a01 = 1.0 - 2.0 * rng.NextDouble();
                double a10 = 1.0 - 2.0 * rng.NextDouble();
                double a11 = 1.0 - 2.0 * rng.NextDouble();

                SquareMatrix A = new SquareMatrix(2);
                A[0, 0] = a00;
                A[0, 1] = a01;
                A[1, 0] = a10;
                A[1, 1] = a11;

                double s, c;
                TwoByTwoSchur(ref a00, ref a01, ref a10, ref a11, out s, out c);

                if (s > 1.0) continue;

                SquareMatrix T = new SquareMatrix(2);
                T[0, 0] = c;
                T[0, 1] = s;
                T[1, 0] = -s;
                T[1, 1] = c;

                SquareMatrix S = T * A * T.Transpose();

                Console.WriteLine("{0} {1}", a00, S[0,0]);

            }
            
            //A[0, 0] = 1.0; A[0, 1] = 4.0;
            //A[1, 0] = 2.0; A[1, 1] = 3.0;
            //A.Eigenvalues();
        }

        /* ORDER STATISTICS */

        // 24-point Gauss-Hermite integration
        // We chose this because the smallest weight is 10^{-16}

        private static double[] xs = new double[] {
            0.22441454747251558515,
            0.67417110703721223600,
            1.1267608176112450721,
            1.5842500109616941485,
            2.0490035736616989118,
            2.5238810170114269742,
            3.0125461375655648257,
            3.5200068130345247113,
            4.0536644024481495039,
            4.6256627564237872650,
            5.2593829276680443674,
            6.0159255614257397173
        };

        private static double[] ws = new double[] {
            0.42693116386869924965,
            0.28617953534644301790,
            0.12773962178455916065,
            3.744547050323074601E-2,
            7.04835581007267097E-3,
            8.23692482688417458E-4,
            5.68869163640437977E-5,
            2.15824570490233363E-6,
            4.01897117494142968E-8,
            3.04625426998756390E-10,
            6.58462024307817006E-13,
            1.66436849648910887E-16
        };

        private static readonly double SqrtTwo = Math.Sqrt(2.0);
        private static readonly double SqrtPI = Math.Sqrt(Math.PI);

        private static double GaussHermiteIntegrate (Func<double, double> f) {
            double y = 0.0;
            for (int i = 0; i < xs.Length; i++) {
                double x = SqrtTwo * xs[i];
                y += (ws[i] / SqrtPI) * (f(x) + f(-x));
            }
            return (y);
        }

        private static double NormalMeanOrderStatisticExpansion (int i, int n) {

            double p = i / (n + 1.0);
            double q = 1.0 - p;

            double FI = Math.Sqrt(2.0) * AdvancedMath.InverseErf(2.0 * p - 1.0);
            double FI2 = 2.0 * Math.PI * Math.Exp(FI * FI) * FI;

            //double FI3 = Math.Pow(2.0 * Math.PI * Math.Exp(FI * FI), 3.0 / 2.0) * (1.0 + 2.0 * FI * FI);
            //double FI4 = Math.Pow(2.0 * Math.PI * Math.Exp(FI * FI), 2.0) * FI * (7.0 + 6.0 * FI * FI);

            return (FI + p * q / 2.0 * FI2 / (n + 2));

            //return (FI + p * q / 2.0 * FI2 / (n + 2) + p * q / (n + 2) / (n + 2) * ((q - p) / 3.0 * FI3 + p * q / 8.0 * FI4));

        }

        private static double NormalMeanOrderStatisticExpansion2 (int i, int n) {

            double p = i / (n + 1.0);
            double q = 1.0 - p;

            double FI = Math.Sqrt(2.0) * AdvancedMath.InverseErf(2.0 * p - 1.0);
            double FI2 = 2.0 * Math.PI * Math.Exp(FI * FI) * FI;

            double FI3 = Math.Pow(2.0 * Math.PI * Math.Exp(FI * FI), 3.0 / 2.0) * (1.0 + 2.0 * FI * FI);
            double FI4 = Math.Pow(2.0 * Math.PI * Math.Exp(FI * FI), 2.0) * FI * (7.0 + 6.0 * FI * FI);

            return (FI + p * q / 2.0 * FI2 / (n + 2) + p * q / (n + 2) / (n + 2) * ((q - p) / 3.0 * FI3 + p * q / 8.0 * FI4));


        }

        [TestMethod]
        public void TestNormalOrderStatistic () {

            int n = 100;
            //int r = 3 * n / 4;
            int r = 52;
            Distribution d = new NormalDistribution();

            double C = Math.Exp(AdvancedIntegerMath.LogFactorial(n) - AdvancedIntegerMath.LogFactorial(r - 1) - AdvancedIntegerMath.LogFactorial(n - r));

            double m = GaussHermiteIntegrate(x => C * MoreMath.Pow(d.LeftProbability(x), r - 1) * MoreMath.Pow(d.RightProbability(x), n - r) * x);
            //double m = GaussHermiteIntegrate(x => 1.0);

            double m2 = FunctionMath.Integrate(
                //x => 1.0 * Math.Exp(-x * x / 2.0) / Math.Sqrt(2.0 * Math.PI),
                x => C * MoreMath.Pow(d.LeftProbability(x), r - 1) * MoreMath.Pow(d.RightProbability(x), n - r) * x * Math.Exp(-x * x / 2.0) / Math.Sqrt(2.0 * Math.PI),
                Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity)
            );

            Console.WriteLine(m);
            Console.WriteLine(m2);
            Console.WriteLine(NormalMeanOrderStatisticExpansion(r, n));
            Console.WriteLine(NormalMeanOrderStatisticExpansion2(r, n));
            //Console.WriteLine(1.5 / Math.Sqrt(Math.PI));
            

        }

        /* INVERSE KS CDF */

        public static double InverseQ (double Q) {

            if (Q == 0.0) return (Double.PositiveInfinity);

            // This is an inversion technique suggested in Numerical Recipies, 3rd Edition, Section 6.14.12
            // Write the Q expansion as
            //   Q = 2 (y - y^4 + y^9 - y^{16} + y^{25} + \cdots)
            // where y = e^{-2 x^2}.  Rewrite this as
            //   y = \frac{Q}{2} + y^4 - y^9 + y^{16} - y^{25} + \cdots
            // and use this iteratively to generate y starting from y = \frac{Q}{2}. When you have y,
            //   x = \sqrt{-\log(x) / 2}
            // The next term, y^36, is below floating point accuracy for all y <~ 1/3, corresponding to Q <~ 2/3.
            // For Q <~ 1/2, 12 or fewer iterations are required for y to converge.

            double halfQ = Q / 2.0;
            double y = halfQ;
            for (int i = 0; i < 256; i++) {
                double y_old = y;
                double y4, y9, y16, y25;
                ComputePowerSet(y, out y4, out y9, out y16, out y25);
                y = halfQ + y4 - y9 + y16 - y25;
                if (y == y_old) {
                    return (Math.Sqrt(-Math.Log(y) / 2.0));
                }
            }

            throw new NonconvergenceException();

        }

        // This gives all the powers of y required in the Q series with the minimum number of operations.

        private static void ComputePowerSet (double y, out double y4, out double y9, out double y16, out double y25) {
            double y2, y8;
            y2 = y * y;
            y4 = y2 * y2;
            y8 = y4 * y4;
            y9 = y8 * y;
            y16 = y8 * y8;
            y25 = y16 * y9;
        }

        [TestMethod]
        public void TestInverseKS () {

            KolmogorovDistribution d = new KolmogorovDistribution();
            
            double x;
            
            Stopwatch s = Stopwatch.StartNew();
            for (double Q = 0.00001; Q <= 0.5; Q += 0.00001) {
                x = InverseQ(Q);
            }
            s.Stop();
            Console.WriteLine(s.ElapsedMilliseconds);

            s.Restart();
            for (double Q = 0.00001; Q <= 0.5; Q += 0.00001) {
                x = d.InverseRightProbability(Q);
            }
            s.Stop();
            Console.WriteLine(s.ElapsedMilliseconds);
            

            Console.WriteLine(InverseQ(0.5));
            Console.WriteLine(d.InverseRightProbability(0.5));

        }

        /* MINIMIZATION */

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

        private static void CubicHermiteMinimum (double x0, double f0, double m0, double x1, double f1, double m1, out double x, out double fpp) {
            double dx = x1 - x0;
            double t;  CubicHermiteMinimum(f0, m0 * dx, f1, m1 * dx, out t, out fpp);
            x = x0 + dx * t;
            fpp = fpp / (dx * dx);
        }

        private static void CubicHermiteMinimum (double f0, double m0, double f1, double m1, out double t, out double fpp) {

            // Given f,f' = f0,m0 at t = 0 and f,f' = f1,m1 at t = 1, the cubic Hermite interpolating polynomial
            // (http://en.wikipedia.org/wiki/Cubic_Hermite_spline) is
            //   f = (2 t^3 - 3 t^2 + 1) f0 + (t^3 - 2 t^2 + t) m0 + (-2 t^3 + 3 t^2) f1 + (t^3 - t^2) m1
            //     = (1 + 2 t) (1 - t)^2 f0 + t (1 - t)^2 m0 + t^2 (3 - 2 t) f1 + t^2 (t - 1) m1
            //     = [2 (f0 - f1) + (m0 + m1)] t^3 - [3 (f0 - f1) + 2 m0 + m1] t^2 + m0 t + p0

            // The derivative of this function is
            //   f' = 3 [2 (f0 - f1) + (m0 + m1)] t^2 - 2 [3 (f0 - f1) + 2 m0 + m1] t + m0 = a t^2 - b t + c
            // Set this equal to zero to get a quadratic equation a t^2 - b t + c = 0 for minimum (and maximum).
            // The descriminant q^2 = b^2 - 4 a c = ?
            // If a > 0 the cubic is upward sloping, so minimum is rightmost solution. If a < 0, minimum is leftmost solution.

            // The second derivative is f'' = 2 a t - b.

            double df = f0 - f1;
            double sm = m0 + m1;
            double q2 = 4.0 * MoreMath.Pow(3.0 * df + sm, 2) - 4.0 * m0 * m1;

            // If q^2 < 0, there is no minimum

            if (q2 < 0.0) {
                t = Double.NaN;
                fpp = Double.NaN;
                return;
            }

            double q = Math.Sqrt(q2);

            double a = 3.0 * (2.0 * df + sm);
            double b = 2.0 * (3.0 * df + 2.0 * m0 + m1);
            double c = m0;

            Console.WriteLine("a={0} b={1} c={2}", a, b, c);

            // If a is very small or zero, our cubic becomes a parabola.
            // This happens, for example, given two points with equal function values and opposite slopes.
            if (Math.Abs(a) <= 1.0E-7 * Math.Abs(b) && Math.Abs(a * c) <= 1.0E-14 * b * b) {
                Console.WriteLine("parabola");
                if (b < 0.0) {
                    t = c / b;
                    fpp = -b;
                } else {
                    t = Double.NaN;
                    fpp = Double.NaN;
                }
                return;
            }

            double t1, t2;
            if (b > 0) {
                t1 = (b + q) / (2.0 * a);
                t2 = (2.0 * c) / (b + q);
            } else {
                t1 = (b - q) / (2.0 * a);
                t2 = (2.0 * c) / (b - q);
            }
            if (a > 0) {
                t = Math.Max(t1, t2);
            } else {
                t = Math.Min(t1, t2);
            }

            fpp = 2.0 * a * t - b;


            /*
            if (a > 0) {
                // rightmost solution is minimum
                // use expression with no cancelation
                if (b > 0) {
                    return ((b + q) / (2.0 * a));
                } else {
                    return ((2.0 * c) / (b - q));
                }
            } else {
                // leftmost solution is minimum
                // again, use expression with no cancelation
                if (b > 0) {
                    return ((2.0 * c) / (b + q));
                } else {
                    return ((b - q) / (2.0 * a));
                }
            }
            */
        }

        [TestMethod]
        public void TestCubicHermite () {

            double x, fpp; CubicHermiteMinimum(3.1428858, -0.99999916, 0.0012931460, 3.2067118, -0.997880, 0.0650731, out x, out fpp);

            //double x = CubicHermiteMinimum(1.0, -1.0, 1.0, +1.0);
            Console.WriteLine(x);

        }

        [TestMethod]
        public void TestParabolicFit () {

            double x0, fpp;
            ParabolicFit(2, 2, 1, 5, 0, 10, out x0, out fpp);

            Console.WriteLine(x0);
            Console.WriteLine(fpp);

        }

        private static double FindMinimum (
            Func<double, double> f,
            double a, double b
        ) {

            // evaluate three points within the bracket
            double u = (3.0 * a + b) / 4.0;
            double v = (a + b) / 2.0;
            double w = (a + 3.0 * b) / 4.0;

            double fu = f(u); double fv = f(v); double fw = f(w);

            Console.WriteLine("f({0})={1}  f({2})={3}  f({4})={5}", u, fu, v, fv, w, fw);

            // move in the bracket boundaries, if possible
            if (fv < fu) { a = u; if (fw < fv) a = v; }
            if (fv < fw) { b = w; if (fu < fv) b = v; }

            Console.WriteLine("a={0} b={1}", a, b);

            // sort u, v, w by fu, fv, fw values
            // these three comparisons are the most efficient three-item sort
            if (fv < fu) { double t = v; v = u; u = t; t = fv; fv = fu; fu = t; }
            if (fw < fu) { double t = w; w = u; u = t; t = fw; fw = fu; fu = t; }
            if (fw < fv) { double t = w; w = v; v = t; t = fw; fw = fv; fv = t; }

            // An evaluation budget of 32 is sufficient for all our test cases except for |x|, which requires 82 (!) evaluations to converge. Parabolic fitting just does a very poor job
            // for this function (at all scales, since it is scale invariant). We should look into cubic fitting.

            EvaluationSettings settings = new EvaluationSettings() { EvaluationBudget = 128, AbsolutePrecision = 0.0, RelativePrecision = 0.0 };
            return (FindMinimum(f, a, b, u, fu, v, fv, w, fw, settings, 3));

        }

        private static double FindMinimum (
            Func<double,double> f,
            double a, double b,
            double u, double fu, double v, double fv, double w, double fw,
            EvaluationSettings settings, int count
        ) {

            double tol = 0.0;

            while (count < settings.EvaluationBudget) {

                Console.WriteLine("n={0} tol={1}", count, tol);
                Console.WriteLine("[{0}  f({1})={2}  f({3})={4}  f({5})={6}  {7}]", a, u, fu, v, fv, w, fw, b);

                Debug.Assert(a < b);
                Debug.Assert((a <= u) && (u <= b));
                Debug.Assert((fu <= fv) && (fv <= fw));

                // Expected final situation is a<tol><tol>u<tol><tol>b, leaving no point left to evaluate that is not within tol of an existing point.

                if ((b - a) <= 4.0 * tol) return (u);

                // While a < u < b is guaranteed, a < v, w < b is not guaranteed, since the bracket can sometimes be made tight enough to exclude v or w.
                // For example, if u < v < w, then we can set b = v, placing w outside the bracket.

                double x, fpp;
                ParabolicFit(u, fu, v, fv, w, fw, out x, out fpp);
                Console.WriteLine("parabolic x={0} f''={1}", x, fpp);
 
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

                    Console.WriteLine("golden section x={0}", x);

                }

                // ensure we don't evaluate within tolerance of an existing point
                if (Math.Abs(x - u) < tol) { Console.WriteLine("shift from u (x={0})", x); x = (x > u) ? u + tol : u - tol; }
                if ((x - a) < tol) { Console.WriteLine("shift from a (x={0})", x); x = a + tol; }
                if ((b - x) < tol) { Console.WriteLine("shift from b (x={0})", x); x = b - tol; }

                count++;
                double fx = f(x);
                Console.WriteLine("f({0}) = {1}", x, fx);

                Console.WriteLine("delta={0}", fu - fx);

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

                        Console.WriteLine("bad point");
                        //throw new NotImplementedException();
                    }

                }

                // if the user has specified a tollerance, use it
                if ((settings.RelativePrecision > 0.0 || settings.AbsolutePrecision > 0.0)) {
                    tol = Math.Max(Math.Abs(u) * settings.RelativePrecision, settings.AbsolutePrecision);
                } else {
                    // otherwise, try to get the tollerance from the curvature
                    if (fpp > 0.0) {
                        tol = Math.Sqrt(2.0 * (Math.Abs(fu) * 1.0E-14 + 1.0E-28) / fpp);
                    } else {
                        // but if we don't have a useable curvature either, wing it
                    }
                }

            }

            throw new NonconvergenceException();

        }

        public delegate void FuncWithDerivative (double x, out double f, out double fp);

        private double FindMinimumWithDerivative (
            FuncWithDerivative f,
            double a, double b
        ) {

            // pick two points in the interval
            double u = 2.0 / 3.0 * a + 1.0 / 3.0 * b;
            double v = 1.0 / 3.0 * a + 2.0 / 3.0 * b;

            // evalue the function there
            double fu, fpu, fv, fpv;
            f(u, out fu, out fpu);
            f(v, out fv, out fpv);

            // move in the bound at the higher side
            if (fu > fv) {
                a = u;
            } else {
                b = v;
            }

            // if f(v) < f(u), swap the points to ensure that u and v are ordered as required
            if (fv < fu) {
                double t;
                t = u; u = v; v = t;
                t = fu; fu = fv; fv = t;
                t = fpu; fpu = fpv; fpv = t;
            }

            // An evaluation budget of 32 is sufficient for all our test cases except for |x|, which requires 82 (!) evaluations to converge. Parabolic fitting just does a very poor job
            // for this function (at all scales, since it is scale invariant). We should look into cubic fitting.

            return (FindMinimumWithDerivative(f, a, b, u, fu, fpu, v, fv, fpv, new EvaluationSettings() { EvaluationBudget = 64, AbsolutePrecision = 0.0, RelativePrecision = 0.0 }));
        }

        private double FindMinimumWithDerivative (
            FuncWithDerivative f,
            double a, double b,
            double u, double fu, double fpu,
            double v, double fv, double fpv,
            EvaluationSettings settings
        ) {

            double tol = 0.0;

            int count = 0;
            while (count < settings.EvaluationBudget) {

                Console.WriteLine("n = {0}, tol = {1}", count, tol);
                Console.WriteLine("[{0} f({1})={2}({3}) f({4})={5}({6}) {7}]", a, u, fu, fpu, v, fv, fpv, b);

                // a, b bracket minimum a < u, v < b and f(u), f(v) <= f(a), f(b)
                Debug.Assert(a < b);
                Debug.Assert((a <= u) && (u <= b));
                //Debug.Assert((a <= v) && (v <= b));
                Debug.Assert(fu <= fv);

                if ((b - a) <= 4.0 * tol) return (u);

                // compute the minimum of the interpolating Hermite cubic
                double x, fpp; CubicHermiteMinimum(u, fu, fpu, v, fv, fpv, out x, out fpp);

                Console.WriteLine("cubic x = {0}, fpp = {1}", x, fpp);

                // if the cubic had no minimum, or the minimum lies outside our bounds, fall back to bisection
                if (Double.IsNaN(x) || (x <= a) || (x >= b)) {

                    // the derivative tells us which side to choose
                    if (fpu > 0.0) {
                        x = (a + u) / 2.0;
                    } else {
                        x = (u + b) / 2.0;
                    }

                    Console.WriteLine("bisection x = {0}", x);

                }

                // ensure we don't evaluate within tolerance of an existing point
                if (Math.Abs(x - u) < tol) { Console.WriteLine("shift from u (x={0})", x); x = (x > u) ? u + tol : u - tol; }
                if ((x - a) < tol) { Console.WriteLine("shift from a (x={0})", x); x = a + tol; }
                if ((b - x) < tol) { Console.WriteLine("shift from b (x={0})", x); x = b - tol; }

                // evaluate the function plus derivative at the predicted minimum
                double fx, fpx;
                f(x, out fx, out fpx);
                count++;

                Console.WriteLine("f({0}) = {1}({2})", x, fx, fpx);

                // check if we have converged
                double df = fu - fx;
                Console.WriteLine("df={0}", df);
                if ((Math.Abs(df) < settings.AbsolutePrecision) || (2.0 * Math.Abs(df) < settings.RelativePrecision * (Math.Abs(fu) + Math.Abs(fx)))) {
                    Console.WriteLine("count = {0}", count);
                    return (x);
                }

                if (fx < fu) {

                    // x is the new lowest point: f(x) < f(u) < f(v)
                    // this is the expected outcome

                    // move the bracket
                    if (x < u) {
                        b = u;
                    } else {
                        a = u;
                    }

                    // x -> u -> v
                    v = u; fv = fu; fpv = fpu;
                    u = x; fu = fx; fpu = fpx;

                } else {

                    // move the bracket
                    if (x < u) {
                        a = x;
                    } else {
                        b = x;
                    }

                    if (fx < fv) {

                        // x lies between other two known points: f(u) < f(x) < f(v)

                        // x -> v
                        v = x; fv = fx; fpv = fpx;

                    } else {

                        // x is higher than both other points: f(u) < f(v) < f(x)
                        // this is a really poor outcome; we expected to get a point lower than our other two and we got a point higher than both
                        // next time we should bisect
                        Console.WriteLine("bad point");
                        //throw new NotImplementedException();

                        //v = x; fv = fx; fpv = fpx;

                    }

                }

                // if the user has specified a tollerance, use it
                if ((settings.RelativePrecision > 0.0 || settings.AbsolutePrecision > 0.0)) {
                    tol = Math.Max(Math.Abs(u) * settings.RelativePrecision, settings.AbsolutePrecision);
                } else {
                    // otherwise, try to get the tollerance from the curvature
                    if (fpp > 0.0) {
                        tol = Math.Sqrt(2.0 * 1.0E-14 * (Math.Abs(fu) + 1.0E-14) / fpp);
                    } else {
                        // but if we don't have a useable curvature either, wing it
                        if (tol == 0.0) tol = 1.0E-7;
                    }
                }


            }

            throw new NonconvergenceException();

        }

        public double FindMinimum (Func<double, double> f, double x, double d, EvaluationSettings settings) {

            // evaluate at x and x + d
            double fx = f(x);
            double y = x + d;
            double fy = f(y);
            int count = 2;

            // if we stepped uphill, reverse direction of steps and exchange x & y
            if (fy > fx) {
                double t = x; x = y; y = t;
                t = fx; fx = fy; fy = t;
                d = -d;
            }

            // we now know f(x) >= f(y) and we are stepping downhill
            // continue stepping until we step uphill
            double z, fz;
            while (true) {

                if (count >= settings.EvaluationBudget) throw new NonconvergenceException();

                z = y + d;
                fz = f(z);
                count++;

                Console.WriteLine("f({0})={1} f({2})={3} f({4})={5} d={6}", x, fx, y, fy, z, fz, d);

                if (fz > fy) break;

                // increase the step size each time
                d = AdvancedMath.GoldenRatio * d;

                // x <- y <- z
                x = y; fx = fy; y = z; fy = fz;


            }

            // we x and z now bracket a local minimum, with y the lowest point evaluated so far
            double a = Math.Min(x, z); double b = Math.Max(x, z);
            if (fz < fx) { double t = x; x = z; z = t; t = fx; fx = fz; fz = t; }

            return (FindMinimum(f, a, b, y, fy, x, fx, z, fz, settings, count));

        }

        [TestMethod]
        public void TestMinimizationWithoutDerivative () {

            //double t = FindMinimum(x => x * Math.Log(x), 0.1, 10.0); // 13
            //double t = FindMinimum(x => Math.Cos(x), 0.0, 5.0); // 7
            //double t = FindMinimum(x => AdvancedMath.Gamma(x), 0.0, 2.0); // 8
            //double t = FindMinimum(x => Math.Exp(-x), 0.0, 1.0); // 29 (not a minimum)
            //double t = FindMinimum(x => Math.Abs(x), -2.0, 3.0); // 82!
            //double t = FindMinimum(x => Math.Cosh(x), -4.0, 3.0); // 7
            //double t = FindMinimum(x => -1.0 / Math.Cosh(x), -1.5, 2.5); // 3
            //double t = FindMinimum(x => Math.Pow(x, 4.0), -2.0, 1.0); // 19

            //double t = FindMinimum(x => x * Math.Log(x), 1.0, -0.03125, new EvaluationSettings() { EvaluationBudget = 128, AbsolutePrecision = 0.0, RelativePrecision = 0.0 });
            double t = FindMinimum(AdvancedMath.Gamma, 1.0, 0.03125, new EvaluationSettings() { EvaluationBudget = 128, AbsolutePrecision = 0.0, RelativePrecision = 0.0 });
            Console.WriteLine(t);
        }

        [TestMethod]
        public void TestOldMinimization () {
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (x * Math.Log(x)); }, Interval.FromEndpoints(0.1, 10.0));
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (Math.Cos(x)); }, Interval.FromEndpoints(0.0, 5.0));
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (AdvancedMath.Gamma(x)); }, Interval.FromEndpoints(0.0, 2.0));
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (Math.Exp(-x)); }, Interval.FromEndpoints(0.0, 1.0));
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (Math.Abs(x)); }, Interval.FromEndpoints(-2.0, 3.0));
            FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (Math.Pow(x,4)); }, Interval.FromEndpoints(-2.0, 1.0));

        }

        [TestMethod]
        public void TestMinimizationWithDerivative () {

            //double x = FindMinimumWithDerivative(XLogXWithDerivative, 0.1, 10.0); // 9
            //double x = FindMinimumWithDerivative(CosineWithDerivative, 0.0, 5.0); // 5
            //double x = FindMinimumWithDerivative(GammaWithDerivative, 0.0, 2.0); // 9
            //double x = FindMinimumWithDerivative(ExponentialDecayWithDerivative, 0.0, 1.0); // 29 (not a minimum)
            //double x = FindMinimumWithDerivative(AbsoluteValueWithDerivative,-2.0, 3.0); // 19
            double x = FindMinimumWithDerivative(CoshWithDerivative, -4.0, 3.0); // 7
            //double x = FindMinimumWithDerivative(SechWithDerivative, -1.5, 2.5); // 5
            //double x = FindMinimumWithDerivative(FourthPowerWithDerivative, -2.0, 1.0); // 17
            Console.WriteLine(x);

        }

        private void CosineWithDerivative (double x, out double f, out double fp) {
            f = Math.Cos(x);
            fp = -Math.Sin(x);
        }

        private void GammaWithDerivative (double x, out double f, out double fp) {
            f = AdvancedMath.Gamma(x);
            fp = f * AdvancedMath.Psi(x);
        }

        private void XLogXWithDerivative (double x, out double f, out double fp) {
            f = x * Math.Log(x);
            fp = 1.0 + Math.Log(x);
        }

        private void ExponentialDecayWithDerivative (double x, out double f, out double fp) {
            f = Math.Exp(-x);
            fp = -f;
        }

        private void AbsoluteValueWithDerivative (double x, out double f, out double fp) {
            f = Math.Abs(x);
            fp = Math.Sign(x);
        }

        private void CoshWithDerivative (double x, out double f, out double fp) {
            f = Math.Cosh(x);
            fp = Math.Sinh(x);
        }

        private void SechWithDerivative (double x, out double f, out double fp) {
            f = -1.0 / Math.Cosh(x);
            fp = Math.Tanh(x) / Math.Cosh(x);
        }

        private void FourthPowerWithDerivative (double x, out double f, out double fp) {
            f = Math.Pow(x, 4);
            fp = 4.0 * Math.Pow(x, 3);
        }

        /* PRIME FACTORIZATION */

        private void FermatFactor (int n) {

            int c = 0;
            int s = (int)Math.Floor(Math.Sqrt(n));
            int x = 2 * s + 1; int y = 1; int r = s * s - n;
            while (r != 0) {
                r += x; x += 2;
                c++;
                do {
                    r -= y; y += 2;
                    c++;
                } while (r > 0);

            }
            int p = (x - y) / 2;
            int q = (x + y - 2) / 2;
            if (p == 1) {
                //Console.WriteLine(q);
            } else {
                FermatFactor(p);
                FermatFactor(q);
            }
        }

        private void PollardRhoFactor (int n) {

            int x = 5; int y = 2; int k = 1; int l = 1;

            for (int c = 0; c < 10000; c++) {
            //while (true) {
                int g = (int) AdvancedIntegerMath.GCF(Math.Abs(y - x), n);
                if (g == n) {
                    Console.WriteLine("n={0}", n);
                    return;
                } else if (g == 1) {
                    k--;
                    if (k == 0) {
                        y = x;
                        l = 2 * l;
                        k = l;
                    }
                    //Console.WriteLine("{0}^2 mod {1}", x, n);
                    x = AdvancedIntegerMath.PowMod(x, 2, n) + 1;
                    if (x == n) x = 0;
                } else {
                    //Console.WriteLine("g={0}", g);
                    n = n / g;
                    x = x % n;
                    y = y % n;
                }
            }



        }

#if FUTURE
        [TestMethod]
        public void TestFactor () {
            /*
            Stopwatch s1 = Stopwatch.StartNew();
            FermatFactor(1157625);
            s1.Stop(); Console.WriteLine(s1.ElapsedMilliseconds);
            
            Stopwatch s2 = Stopwatch.StartNew();
            PollardRhoFactor(37);
            s2.Stop(); Console.WriteLine(s2.ElapsedMilliseconds);
            */
            // for (3*5*7)^3 = 1157625, Pollard's Rho fails to factor 125 = 5^3
            
            int c = 0;
            Stopwatch s1 = Stopwatch.StartNew();

            //for (int i = 1; i < 1000000; i+=2) {
            int i = 220211;
                List<Factor> factors = AdvancedIntegerMath.Factor(i);
                //FermatFactor(i);
                
                int m = 1;
                foreach (Factor factor in factors) {
                    Console.WriteLine(factor);
                    if (!AdvancedIntegerMath.IsPrime(factor.Value)) {
                        c++;
                        Console.WriteLine("for {0}, factor {1} is not prime", i, factor.Value);
                    }
                    //Console.WriteLine("{0} {1}", factor.Value, factor.Multiplicity);
                    m *= (int)MoreMath.Pow(factor.Value, factor.Multiplicity);
                }
                if (m != i) {
                    Console.WriteLine("for {0}, factors do not multiply to number", i);
                }
                //Assert.IsTrue(m == i);
                //Console.WriteLine(m);
                
            // }
            s1.Stop(); Console.WriteLine(s1.ElapsedMilliseconds);

            Console.WriteLine(c);
            
        }
#endif
#if FUTURE
        [TestMethod]
        public void BigAdd () {

            BigFloat a = new BigFloat(new byte[] { 2, 2, 5 }, 0);
            Console.WriteLine(a);
            BigFloat b = new BigFloat(new byte[] { 3, 3, 3 }, 0);
            Console.WriteLine(b);
            BigFloat c = a + b;
            Console.WriteLine(c);

        }
#endif
        [TestMethod]
        public void EinTest () {

            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 3.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 2.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 2.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 1.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 1.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 0.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 0.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -0.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -1.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -1.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -2.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -2.5)));

        }

        [TestMethod]
        public void MultiRootTest () {

            // here is a system with a zeros at (3,4) and (1,0)
            Func<double[], double[]> f = delegate (double[] u) {
                double x = u[0]; double y = u[1];
                double a = 2.0 * y - x;
                double b = (x * x + x * (y * y - 2.0) - 4.0 * y) / (x + 4.0);
                double c = Math.Sqrt(x * x + y * y);
                return (new double[] { b - a, c - a });
            };

            FunctionMath.FindZero(f, new double[] { 1.0, 1.0 });

        }

        /* INVERSE BETA */

        public static double ApproximateInverseBetaSeries (double a, double b, double P) {

            double bigB = AdvancedMath.Beta(a, b);

            double z = Math.Pow(a * P * bigB, 1.0 / a);

            double z2 = (b - 1.0) / (a + 1.0) * z;

            double z3 = z2 * (a * a + 3.0 * b * a - a + 5.0 * b - 4.0) / (a + 1.0) / (a + 2.0) / 2.0 * z;

            double z4 = z2 * (a * a * a * a + (6.0 * b - 1.0) * a * a * a + (b + 2.0) * (8.0 * b - 5.0) * a * a +
                (33.0 * b * b - 30.0 * b + 4.0) * a + b * (31.0 * b - 47.0) + 18.0) / MoreMath.Pow(a + 1.0, 2) / (a + 2.0) / (a + 3.0) / 3.0 * z * z;

            Console.WriteLine("z={0} z2={1} z3={2} z4={3}", z, z2, z3, z4);

            return (z * (1.0 + z2 + z3 + z4));

        }

        public static double ApproximateInverseBeta (double a, double b, double P) {

            if (P > 0.5) {
                return (1.0 - ApproximateInverseBeta(b, a, 1.0 - P));
            } else {
                double bigB = AdvancedMath.Beta(a, b);
                double z = Math.Pow(a * P * bigB, 1.0 / a);
                double z2 = (b - 1.0) / (a + 1.0) * z;
                if (z2 < 0.25) {
                    double z3 = z2 * (a * a + 3.0 * b * a - a + 5.0 * b - 4.0) / (a + 1.0) / (a + 2.0) / 2.0 * z;
                    return (z * (1.0 + z2 + z3));
                } else {
                    throw new NotImplementedException();
                }
            }

        }

        public static double RefineInverseBeta (double a, double b, double P, double x) {

            for (int i = 0; i < 8; i++) {
                double x_old = x;
                double y = AdvancedMath.LeftRegularizedBeta(a, b, x) - P;
                double yp = Math.Pow(x, a - 1.0) * Math.Pow(1.0 - x, b - 1.0) / AdvancedMath.Beta(a, b);
                double dx = -y / yp;
                x += dx;
                if (x == x_old) return (x);
            }
            return (x);
        }

        [TestMethod]
        public void TestBeta () {

            double a = 200.0; double b = 200.0; double P = 1.0E-5;
            double x1 = ApproximateInverseBetaSeries(a, b, P);
            if ((0.0 < x1) && (x1 < 1.0)) {
                Console.WriteLine("x1 {0} {1}", x1, AdvancedMath.LeftRegularizedBeta(a, b, x1));
            }

            double x2 = 1.0 - ApproximateInverseBetaSeries(b, a, 1.0 - P);
            if ((0.0 < x2) && (x2 < 1.0)) {
                Console.WriteLine("x2 {0} {1}", x2, AdvancedMath.LeftRegularizedBeta(a, b, x2));
            }

            //x1 = RefineInverseBeta(a, b, P, x1);
            //Console.WriteLine("{0} {1}", x1, AdvancedMath.LeftRegularizedBeta(a, b, x1));

            NormalDistribution N = new NormalDistribution();
            double m = a / (a + b); double s = Math.Sqrt(a * b / (a + b + 1.0)) / (a + b);
            double x3 = m + s * N.InverseLeftProbability(P);
            if ((0.0 < x3) && (x3 < 1.0)) {
                Console.WriteLine("x3 {0} {1}", x3, AdvancedMath.LeftRegularizedBeta(a, b, x3));
            }

            //Console.WriteLine(AdvancedMath.Beta(a, b, 0.35) / AdvancedMath.Beta(a, b));
            //Console.WriteLine(AdvancedMath.Beta(a, b, 0.40) / AdvancedMath.Beta(a, b));
            //Console.WriteLine(AdvancedMath.Beta(a, b, 0.45) / AdvancedMath.Beta(a, b));

        }

    }

}
