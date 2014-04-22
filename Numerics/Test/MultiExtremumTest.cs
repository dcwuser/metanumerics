using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Test {

    // There are a number of sources of multidimensional test optimization problems, e.g.
    //   http://en.wikipedia.org/wiki/Test_functions_for_optimization
    //   http://www.sfu.ca/~ssurjano/optimization.html
    //   http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf  
    // Many test problems have local minima. Demanding it find the global minimum is not
    // a fair test of a local minimizer, unless it is started definitively closest to
    // the global minimum.

    [TestClass]
    public class MultiExtremumTest {

        // Test functions for local minizer

        [TestMethod]
        public void Vardim () {

            // This function is described by Powell in his NEWOUA paper as one that caused problems for the simplest formulation of that minimizer.

            Func<IList<double>, double> f = (IList<double> x) => {
                double s1 = 0.0;
                double s2 = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    s1 += i * (x[i] - 1.0);
                    s2 += MoreMath.Sqr(x[i] - 1.0);
                }
                return (s2 + s1 * s1 + s1 * s1 * s1 * s1);
            };

            for (int n = 2; n < 8; n++) {

                ColumnVector start = new ColumnVector(n);
                ColumnVector solution = new ColumnVector(n);
                solution.Fill((i, j) => 1.0);

                MultiExtremum min = MultiFunctionMath.FindMinimum(f, start);

                Console.WriteLine(min.Value);
                Console.WriteLine(min.EvaluationCount);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Value, 0.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Location, solution, Math.Sqrt(TestUtilities.TargetPrecision)));

            }

        }

        [TestMethod]
        public void GoldsteinPrice () {

            // Goldstein-Price has a valley with a complicated shape and a minimum value of 3 at (0,-1).
            // Start at 2,2, coverges to some weird value. Investigate.

            Func<IList<double>, double> fGoldsteinPrice = (IList<double> v) => {
                double x = v[0];
                double y = v[1];
                return (
                    (1 + MoreMath.Pow(x + y + 1, 2) * (19 - 14 * x + 3 * x * x - 14 * y + 6 * x * y + 6 * y * y)) *
                    (30 + MoreMath.Pow(2 * x - 3 * y, 2) * (18 - 32 * x + 12 * x * x + 48 * y - 36 * x * y + 27 * y * y))
                );
            };

            ColumnVector start = new ColumnVector(1.0, 1.0);
            //ColumnVector start = new ColumnVector(2.0, 2.0);

            MultiExtremum min = MultiFunctionMath.FindMinimum(fGoldsteinPrice, start);

            Console.WriteLine(min.EvaluationCount);
            Console.WriteLine(min.Value);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Value, 3.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Location, new ColumnVector(0.0, -1.0), Math.Sqrt(TestUtilities.TargetPrecision)));

        }

        [TestMethod]
        public void Beale () {

            // Beale is a very interesting function.
            // The only local minimum is at (3,1/2) where the function value is 0.
            // But along the lines (0,+\infty) and (-\infty,1) the function value decreases toward 0 at infinity.
            // From the point of view of a local minimizer those are perfectly valid downhill directions and
            // if the minimzer gets caught in them it will move toward infinity until its evaluation budget is
            // exhausted. Starting from y > 0 will probably keep us safe, or we can do bounded optimization.

            Func<IList<double>, double> fBeale = (IList<double> x) =>
                MoreMath.Sqr(1.5 - x[0] + x[0] * x[1]) +
                MoreMath.Sqr(2.25 - x[0] + x[0] * x[1] * x[1]) +
                MoreMath.Sqr(2.625 - x[0] + x[0] * x[1] * x[1] * x[1]);
            ColumnVector start = new ColumnVector(1.0, 2.0);

            MultiExtremum min = MultiFunctionMath.FindMinimum(fBeale, start);

            Console.WriteLine(min.EvaluationCount);
            Console.WriteLine(min.Value);
            ColumnVector solution = new ColumnVector(3.0, 0.5);
            Assert.IsTrue(min.Dimension == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Value, 0.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Location, solution, Math.Sqrt(TestUtilities.TargetPrecision)));

        }

        [TestMethod]
        public void StylblinskiTang () {

            Func<IList<double>, double> fStyblinskiTang = (IList<double> x) => {
                double fst = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    double x1 = x[i];
                    double x2 = MoreMath.Sqr(x1);
                    fst += x2 * (x2 - 16.0) + 5.0 * x1;
                }
                return (fst / 2.0);
            };

            // solution coordinate is root of 5 - 32 x + 4 x^3 = 0 with negative second derivative.
            // There are two such roots. 
            double root1 = -2.9035340277711770951;
            double root2 = 2.7468027709908369925;

            for (int n = 2; n < 8; n++) {

                ColumnVector start = new ColumnVector(n);
                //ColumnVector start = new ColumnVector(-1.0, -2.0, -3.0, -4.0, -5.0, -6.0);

                MultiExtremum minimum = MultiFunctionMath.FindMinimum(fStyblinskiTang, start);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine(minimum.Value);
                for (int i = 0; i < minimum.Dimension; i++) {
                    Console.WriteLine(minimum.Location[i]);
                    Assert.IsTrue(
                        TestUtilities.IsNearlyEqual(minimum.Location[i], root1, Math.Sqrt(TestUtilities.TargetPrecision)) ||
                        TestUtilities.IsNearlyEqual(minimum.Location[i], root2, Math.Sqrt(TestUtilities.TargetPrecision))
                    );
                }
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, fStyblinskiTang(minimum.Location)));

            }

        }

        [TestMethod]
        public void ThreeHumpCamel () {

            // This function has three local minima, so not at all starting points should be expected to bring us to the global minimum at the origin.

            Func<IList<double>, double> function = (IList<double> x) => 2.0 * MoreMath.Pow(x[0], 2) - 1.05 * MoreMath.Pow(x[0], 4) + MoreMath.Pow(x[0], 6) / 6.0 + x[0] * x[1] + MoreMath.Pow(x[1], 2);
            ColumnVector start = new ColumnVector(1.0, 1.0);

            MultiExtremum minimum = MultiFunctionMath.FindMinimum(function, start);

            Console.WriteLine(minimum.Value);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.0));
            Console.WriteLine("{0} {1}", minimum.Location[0], minimum.Location[1]);
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, new ColumnVector(2), Math.Sqrt(TestUtilities.TargetPrecision)));

        }

        [TestMethod]
        public void McCormick () {

            Func<IList<double>, double> fMcCormick = (IList<double> x) => Math.Sin(x[0] + x[1]) + MoreMath.Sqr(x[0] - x[1]) - 1.5 * x[0] + 2.5 * x[1] + 1.0;

        }

        [TestMethod]
        public void Rosenbrock () {

            // This is a multi-dimensional generalization of the famous Rosenbrock, aka banana function,
            // which has a narrow parabolic valley whose floor slopes only gently to the minimum.

            Func<IList<double>, double> fRosenbrock = delegate(IList<double> x) {
                double s = 0.0;
                for (int i = 0; i < (x.Count - 1); i++) {
                    s += 100.0 * MoreMath.Pow(x[i + 1] - x[i] * x[i], 2) + MoreMath.Pow(1.0 - x[i], 2);
                }
                return (s);
            };

            for (int n = 2; n < 8; n++) {

                double[] start = new double[n];
                MultiExtremum minimum = MultiFunctionMath.FindMinimum(fRosenbrock, start);
                
                Console.WriteLine(minimum.EvaluationCount);
                ColumnVector solution = new ColumnVector(n);
                for (int i = 0; i < solution.Dimension; i++) solution[i] = 1.0;

                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.0));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, solution, Math.Sqrt(TestUtilities.TargetPrecision)));

            }

        }

        [TestMethod]
        public void SmoothedEasom () {

            // This function is mostly flat except very near (\pi, \pi).
            // For (1,1) or (0,0) "converges" to minimum of 0 at (1.30, 1.30). This is probably a local minimum of the cosine product.

            //Func<IList<double>, double> function = (IList<double> x) => -Math.Exp(-(MoreMath.Sqr(x[0] - Math.PI) + MoreMath.Sqr(x[1] - Math.PI)));
            Func<IList<double>, double> function = (IList<double> x) => -Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Exp(-(MoreMath.Sqr(x[0] - Math.PI) + MoreMath.Sqr(x[1] - Math.PI)));
            ColumnVector start = new ColumnVector(2.0, 2.0);

            MultiExtremum minimum = MultiFunctionMath.FindMinimum(function, start);
            Console.WriteLine(minimum.EvaluationCount);
            Console.WriteLine(minimum.Value);
            Console.WriteLine("{0} {1}", minimum.Location[0], minimum.Location[1]);

        }

        [TestMethod]
        public void Perm () {

            Func<IList<double>, double> fPerm = (IList<double> x) => {
                double s = 0.0;
                for (int i = 1; i <= x.Count; i++) {
                    double t = 0.0;
                    for (int j = 0; j < x.Count; j++) {
                        t += (j + 1) * (MoreMath.Pow(x[j], i) - 1.0 / MoreMath.Pow(j + 1, i));
                    }
                    s += MoreMath.Sqr(t);
                }
                return (s);
            };

            int n = 4;

            ColumnVector start = new ColumnVector(n);

            MultiExtremum minimum = MultiFunctionMath.FindMinimum(fPerm, start);

            Console.WriteLine(minimum.EvaluationCount);
            Console.WriteLine(minimum.Value);
            for (int i = 0; i < minimum.Dimension; i++) Console.WriteLine(minimum.Location[i]);

        }

        [TestMethod]
        public void SumOfPowers () {

            Func<IList<double>, double> function = (IList<double> x) => {
                double s = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    s += MoreMath.Pow(Math.Abs(x[i]), i + 2);
                }
                return (s);
            };

            for (int n = 2; n < 8; n++) {
                ColumnVector start = new ColumnVector(n);
                for (int i = 0; i < n; i++) start[i] = 1.0;

                MultiExtremum minimum = MultiFunctionMath.FindMinimum(function, start);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine(minimum.Value);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.0, 1.0E-10));
                //Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, new ColumnVector(n), Math.Sqrt(TestUtilities.TargetPrecision)));

            }

        }

        // Global minimization tests

        [TestMethod]
        public void Ackley () {

            // Ackley's function has many local minima, and a global minimum at (0, 0) -> 0.

            Func<IList<double>, double> function = (IList<double> x) => {
                double s = 0.0;
                double c = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    s += x[i] * x[i];
                    c += Math.Cos(2.0 * Math.PI * x[i]);
                }
                return (-20.0 * Math.Exp(-0.2 * Math.Sqrt(s / x.Count)) - Math.Exp(c / x.Count) + 20.0 + Math.E);
            };

            EvaluationSettings settings = new EvaluationSettings() { AbsolutePrecision = 1.0E-8, EvaluationBudget = 10000000 };

            for (int n = 2; n < 32; n = (int) Math.Round(AdvancedMath.GoldenRatio * n)) {
                Console.WriteLine("n={0}", n);

                Interval[] box = new Interval[n];
                for (int i = 0; i < box.Length; i++) box[i] = Interval.FromEndpoints(-32.0, 32.0);

                MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box, settings);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine(minimum.Value);
                foreach (double coordinate in minimum.Location) Console.WriteLine(coordinate);
            }
        }

        [TestMethod]
        public void Schwefel () {

            // Test over [-500,500], minimum at (420.969,...) -> -418.983*d, many local minima

            Func<IList<double>, double> function = (IList<double> x) => {
                double s = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    s += x[i] * Math.Sin(Math.Sqrt(Math.Abs(x[i])));
                }
                return (-s);
            };

            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-6, EvaluationBudget = 10000000 };

            for (int n = 2; n < 32; n = (int) Math.Round(AdvancedMath.GoldenRatio * n)) {
                Console.WriteLine("n={0}", n);

                Interval[] box = new Interval[n];
                for (int i = 0; i < box.Length; i++) box[i] = Interval.FromEndpoints(-500.0, 500.0);

                MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box, settings);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine(minimum.Value);
                foreach (double coordinate in minimum.Location) Console.WriteLine(coordinate);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, -418.983 * n, settings.RelativePrecision));
            }

        }

        [TestMethod]
        public void Bukin () {

            // Burkin has a narrow valley, not aligned with any axis, punctuated with many tiny "wells" along its bottom.
            // The deepest well is at (-10,1)-> 0.
            Func<IList<double>, double> function = (IList<double> x) => 100.0 * Math.Sqrt(Math.Abs(x[1] - 0.01 * x[0] * x[0])) + 0.01 * Math.Abs(x[0] + 10.0);

            IList<Interval> box = new Interval[] { Interval.FromEndpoints(-15.0, -5.0), Interval.FromEndpoints(-3.0, 3.0) };

            EvaluationSettings settings = new EvaluationSettings() { AbsolutePrecision = 1.0E-8, EvaluationBudget = 1000000 };
            settings.Update += (object result) => {
                MultiExtremum e = (MultiExtremum) result;
                Console.WriteLine("After {0} evaluations, best value {1}", e.EvaluationCount, e.Value);
            };
            MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box, settings);

            Console.WriteLine(minimum.EvaluationCount);
            Console.WriteLine(minimum.Value);
            Console.WriteLine("{0} {1}", minimum.Location[0], minimum.Location[1]);

        }

        [TestMethod]
        public void Easom () {

            Func<IList<double>, double> function = (IList<double> x) => -Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Exp(-(MoreMath.Sqr(x[0] - Math.PI) + MoreMath.Sqr(x[1] - Math.PI)));
            IList<Interval> box = new Interval[] { Interval.FromEndpoints(-10.0, 10.0), Interval.FromEndpoints(-10.0, 10.0) };

            MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box);

            Console.WriteLine(minimum.EvaluationCount);
            Console.WriteLine(minimum.Value);


        }

        [TestMethod]
        public void PackCirclesInSquare () {

            // See http://en.wikipedia.org/wiki/Circle_packing_in_a_square
            // This is a great test problem because it is simple to understand, hard to solve,
            // known solutions are proven, and it covers many dimensions.

            double[] solutions = new double[] { 0.0, 0.0,
                Math.Sqrt(2.0),
                Math.Sqrt(6.0) - Math.Sqrt(2.0),
                1.0,
                Math.Sqrt(2.0) / 2.0,
                Math.Sqrt(13.0) / 6.0,
                4.0 - 2.0 * Math.Sqrt(3.0),
                (Math.Sqrt(6.0) - Math.Sqrt(2.0)) / 2.0,
                1.0 / 2.0,
                Double.NaN
            };

            for (int n = 2; n < 8; n++) {

                if (n == 6 || n == 7) continue;

                Console.WriteLine("n={0}", n);

                Func<IList<double>, double> function = (IList<double> x) => {
                    // Intrepret coordinates as (x_1, y_1, x_2, y_2, \cdots, x_n, y_n)
                    // Iterate over all pairs of points, finding smallest distance betwen any pair of points.
                    double sMin = Double.MaxValue;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < i; j++) {
                            double s = MoreMath.Hypot(x[2 * i] - x[2 * j], x[2 * i + 1] - x[2 * j + 1]);
                            if (s < sMin) sMin = s;
                        }
                    }
                    return (-sMin);
                };

                IList<Interval> box = new Interval[2 * n];
                for (int i = 0; i < box.Count; i++) box[i] = Interval.FromEndpoints(0.0, 1.0);

                EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-6, AbsolutePrecision = 1.0E-8, EvaluationBudget = 10000000 };

                MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box, settings);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine("{0} {1}", solutions[n], -minimum.Value);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(-minimum.Value, solutions[n], settings.RelativePrecision)); 

            }

        }

        [TestMethod]
        public void PackSpheresInCube () {

            // See http://www.combinatorics.org/ojs/index.php/eljc/article/download/v11i1r33/pdf

            double[] solutions = new double[] { 0.0, 0.0,
                Math.Sqrt(3.0),
                Math.Sqrt(2.0),
                Math.Sqrt(2.0),
                Math.Sqrt(5.0) / 2.0,
                3.0 * Math.Sqrt(2.0) / 4.0,
                1.001089824548965412749095030329,
                1.0,
                Math.Sqrt(3.0) / 2.0,
                3.0 / 4.0
            };

            for (int n = 2; n < 7; n++) {

                Func<IList<double>, double> function = (IList<double> x) => {
                    // Intrepret coordinates as (x_1, y_1, z_1, x_2, y_2, z_2, \cdots, x_n, y_n, z_n)
                    // Iterate over all pairs of points, finding smallest distance betwen any pair of points.
                    double sMin = Double.MaxValue;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < i; j++) {
                            double s = Math.Sqrt(
                                MoreMath.Sqr(x[3 * i] - x[3 * j]) +
                                MoreMath.Sqr(x[3 * i + 1] - x[3 * j + 1]) +
                                MoreMath.Sqr(x[3 * i + 2] - x[3 * j + 2])
                            );
                            if (s < sMin) sMin = s;
                        }
                    }
                    return (-sMin);
                };

                IList<Interval> box = new Interval[3 * n];
                for (int i = 0; i < box.Count; i++) box[i] = Interval.FromEndpoints(0.0, 1.0);

                EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-6, AbsolutePrecision = 1.0E-8, EvaluationBudget = 10000000 };

                MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box, settings);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine("{0} {1}", solutions[n], -minimum.Value);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(-minimum.Value, solutions[n], settings.RelativePrecision));

            }
        }
    }
}
