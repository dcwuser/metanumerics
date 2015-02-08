using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;

namespace Test {

    // There are a number of sources of multidimensional test optimization problems, e.g.
    //   http://en.wikipedia.org/wiki/Test_functions_for_optimization
    //   http://www.sfu.ca/~ssurjano/optimization.html
    //   http://www.zsd.ict.pwr.wroc.pl/files/docs/functions.pdf  
    // Many test problems have multiple local minima. Demanding it find the global minimum is not
    // a fair test of a local minimizer, unless it is started definitively closest to
    // the global minimum.

    [TestClass]
    public class MultiExtremumTest {

        // Test functions for local minizer

        [TestMethod]
        public void MinimizeQuadratic1D () {

            Func<IList<double>, double> f = (IList<double> x) => 1.0 + MoreMath.Sqr(x[0] - 3.0);
            MultiExtremum m = MultiFunctionMath.FindLocalMinimum(f, new double[] { 1.0 });

            Assert.IsTrue(TestUtilities.IsNearlyEqual(1.0, m.Value));

        }

        [TestMethod]
        public void MinimizeQuadratic2D () {

            Func<IList<double>, double> f = (IList<double> x) => 1.0 + 2.0 * MoreMath.Sqr(x[0] - 3.0) + 4.0 * (x[0] - 3.0) * (x[1] - 5.0) + 6.0 * MoreMath.Sqr(x[1] - 5.0);
            MultiExtremum m = MultiFunctionMath.FindLocalMinimum(f, new double[] { 1.0, 1.0 });

            Console.WriteLine(m.Value);

            SymmetricMatrix H = new SymmetricMatrix(2);
            H[0, 0] = 4.0;
            H[0, 1] = 4.0;
            H[1, 1] = 12.0;

            Assert.IsTrue(m.Dimension == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(m.Value, 1.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(m.Location, new ColumnVector(3.0, 5.0), new EvaluationSettings() { RelativePrecision = Math.Sqrt(TestUtilities.TargetPrecision) }));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(m.HessianMatrix, H, new EvaluationSettings() { RelativePrecision = 8.0 * Math.Sqrt(Math.Sqrt(TestUtilities.TargetPrecision)) }));

        }

        [TestMethod]
        public void MinimizePerturbedQuadratic2D () {

            EvaluationSettings s = new EvaluationSettings() { EvaluationBudget = 100, RelativePrecision = 1.0E-10 };

            Func<IList<double>, double> f = (IList<double> x) =>
                1.0 + 2.0 * MoreMath.Sqr(x[0] - 3.0) + 4.0 * (x[0] - 3.0) * (x[1] - 5.0) + 6.0 * MoreMath.Sqr(x[1] - 5.0) +
                7.0 * MoreMath.Pow(x[0] - 3.0, 4) + 8.0 * MoreMath.Pow(x[1] - 5.0, 4);

            MultiExtremum m = MultiFunctionMath.FindLocalMinimum(f, new double[] { 1.0, 1.0 }, s);

            Assert.IsTrue(m.EvaluationCount <= s.EvaluationBudget);
            Assert.IsTrue(m.Dimension == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(m.Value, 1.0, s));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(m.Location, new ColumnVector(3.0, 5.0), new EvaluationSettings() { RelativePrecision = Math.Sqrt(s.RelativePrecision) }));

        }

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

                Console.WriteLine("n = {0}", n);

                ColumnVector start = new ColumnVector(n);
                ColumnVector solution = new ColumnVector(n);
                solution.Fill((i, j) => 1.0);

                MultiExtremum min = MultiFunctionMath.FindLocalMinimum(f, start);

                Console.WriteLine(min.EvaluationCount);
                Console.WriteLine("{0} ({1}) ?= 0.0", min.Value, min.Precision);

                Assert.IsTrue(min.Dimension == n);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Value, 0.0, new EvaluationSettings() { AbsolutePrecision = 4.0 * min.Precision }));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Location, solution, new EvaluationSettings() { RelativePrecision = Math.Sqrt(4.0 * min.Precision) }));

            }

        }

        [TestMethod]
        public void GoldsteinPrice () {

            // Goldstein-Price has a valley with a complicated shape and a global minimum value of 3 at (0,-1).
            // It also has local minima, so we have to start close to this minimum if we expect to end at it.

            Func<IList<double>, double> fGoldsteinPrice = (IList<double> v) => {
                double x = v[0];
                double y = v[1];
                return (
                    (1 + MoreMath.Pow(x + y + 1, 2) * (19 - 14 * x + 3 * x * x - 14 * y + 6 * x * y + 6 * y * y)) *
                    (30 + MoreMath.Pow(2 * x - 3 * y, 2) * (18 - 32 * x + 12 * x * x + 48 * y - 36 * x * y + 27 * y * y))
                );
            };

            ColumnVector start = new ColumnVector(0.5, -0.5);
            MultiExtremum min = MultiFunctionMath.FindLocalMinimum(fGoldsteinPrice, start);

            Console.WriteLine(min.EvaluationCount);
            Console.WriteLine(min.Value);

            Assert.IsTrue(min.Dimension == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Value, 3.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Location, new ColumnVector(0.0, -1.0), Math.Sqrt(TestUtilities.TargetPrecision)));

            MultiExtremum min2 = MultiFunctionMath.FindGlobalMinimum(fGoldsteinPrice, new Interval[] { Interval.FromEndpoints(-2.0, 2.0), Interval.FromEndpoints(-2.0, 2.0) });

            Assert.IsTrue(min2.Dimension == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(min2.Value, 3.0, new EvaluationSettings() { AbsolutePrecision = min2.Precision }));
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(min2.Location, new ColumnVector(0.0, -1.0), Math.Sqrt(TestUtilities.TargetPrecision)));

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

            ColumnVector start = new ColumnVector(2.0, 1.0);

            MultiExtremum min = MultiFunctionMath.FindLocalMinimum(fBeale, start);

            Console.WriteLine(min.EvaluationCount);
            Console.WriteLine(min.Value);

            ColumnVector solution = new ColumnVector(3.0, 0.5);
            Assert.IsTrue(min.Dimension == 2);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Value, 0.0, new EvaluationSettings() { AbsolutePrecision = min.Precision }));
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

            // tested up to n=16, works with slightly decreasing accuracy of Location
            for (int n = 2; n < 8; n++) {

                Console.WriteLine(n);

                ColumnVector start = new ColumnVector(n);
                //ColumnVector start = new ColumnVector(-1.0, -2.0, -3.0, -4.0, -5.0, -6.0);

                MultiExtremum minimum = MultiFunctionMath.FindLocalMinimum(fStyblinskiTang, start);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine(minimum.Value);
                for (int i = 0; i < minimum.Dimension; i++) {
                    Console.WriteLine(minimum.Location[i]);
                    Assert.IsTrue(
                        TestUtilities.IsNearlyEqual(minimum.Location[i], root1, Math.Sqrt(Math.Sqrt(TestUtilities.TargetPrecision))) ||
                        TestUtilities.IsNearlyEqual(minimum.Location[i], root2, Math.Sqrt(Math.Sqrt(TestUtilities.TargetPrecision)))
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

            MultiExtremum minimum = MultiFunctionMath.FindLocalMinimum(function, start);

            Console.WriteLine("{0} ({1}) ?= 0.0", minimum.Value, minimum.Precision);
            Console.WriteLine("{0} {1}", minimum.Location[0], minimum.Location[1]);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.0, new EvaluationSettings() { AbsolutePrecision = 2.0 * minimum.Precision }));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, new ColumnVector(2), new EvaluationSettings { AbsolutePrecision = 2.0 * Math.Sqrt(minimum.Precision) }));

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


            for (int n = 2; n < 8; n++)
            {
                Console.WriteLine("n={0}", n);

                double[] start = new double[n];
                MultiExtremum minimum = MultiFunctionMath.FindLocalMinimum(fRosenbrock, start);
                
                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine("{0} {1} ?= 0.0", minimum.Value, minimum.Precision);

                ColumnVector solution = new ColumnVector(n);
                for (int i = 0; i < solution.Dimension; i++) solution[i] = 1.0;

                Assert.IsTrue(minimum.Dimension == n);
                Assert.IsTrue(minimum.Precision > 0.0);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.0, new EvaluationSettings() { AbsolutePrecision = 2.0 * minimum.Precision }));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, solution, new EvaluationSettings() { RelativePrecision = 2.0 * Math.Sqrt(minimum.Precision) }));

            }

        }

        [TestMethod]
        public void SmoothedEasom () {

            // This function is mostly flat except very near (\pi, \pi).
            // For (1,1) or (0,0) "converges" to minimum of 0 at (1.30, 1.30). This is probably a local minimum of the cosine product.

            //Func<IList<double>, double> function = (IList<double> x) => -Math.Exp(-(MoreMath.Sqr(x[0] - Math.PI) + MoreMath.Sqr(x[1] - Math.PI)));
            Func<IList<double>, double> function = (IList<double> x) => -Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Exp(-(MoreMath.Sqr(x[0] - Math.PI) + MoreMath.Sqr(x[1] - Math.PI)));
            ColumnVector start = new ColumnVector(2.0, 2.0);

            MultiExtremum minimum = MultiFunctionMath.FindLocalMinimum(function, start);
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

            MultiExtremum minimum = MultiFunctionMath.FindLocalMinimum(fPerm, start);

            Console.WriteLine(minimum.EvaluationCount);
            Console.WriteLine(minimum.Value);
            for (int i = 0; i < minimum.Dimension; i++) Console.WriteLine(minimum.Location[i]);

        }

        [TestMethod]
        public void SumOfPowers () {
            
            // This test is difficult because the minimum is emphatically not quadratic.
            // We do get close to the minimum but we massivly underestimate our error.

            Func<IList<double>, double> function = (IList<double> x) => {
                double s = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    s += MoreMath.Pow(Math.Abs(x[i]), i + 2);
                }
                return (s);
            };

            for (int n = 2; n < 8; n++) {

                Console.WriteLine(n);

                ColumnVector start = new ColumnVector(n);
                for (int i = 0; i < n; i++) start[i] = 1.0;

                EvaluationSettings settings = new EvaluationSettings() { AbsolutePrecision = 1.0E-8, EvaluationBudget = 32 * n * n * n };

                MultiExtremum minimum = MultiFunctionMath.FindLocalMinimum(function, start, settings);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine("{0} {1}", minimum.Value, minimum.Precision);
                Console.WriteLine("|| {0} {1} ... || = {2}", minimum.Location[0], minimum.Location[1], minimum.Location.FrobeniusNorm());

                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.0, new EvaluationSettings() { AbsolutePrecision = 1.0E-4 }));
                //Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, new ColumnVector(n), new EvaluationSettings() { AbsolutePrecision = 1.0E-2 }));

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

            for (int n = 2; n < 16; n = (int) Math.Round(AdvancedMath.GoldenRatio * n)) {
                Console.WriteLine("n={0}", n);

                Interval[] box = new Interval[n];
                for (int i = 0; i < box.Length; i++) box[i] = Interval.FromEndpoints(-32.0, 32.0);

                MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box, settings);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine(minimum.Value);
                foreach (double coordinate in minimum.Location) Console.WriteLine(coordinate);

                Assert.IsTrue(minimum.Dimension == n);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, 0.0, new EvaluationSettings() { AbsolutePrecision = 2.0 * minimum.Precision }));

                ColumnVector solution = new ColumnVector(n);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, solution, new EvaluationSettings() { AbsolutePrecision = 2.0 * Math.Sqrt(minimum.Precision) }));

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

            for (int n = 2; n < 32; n = (int) Math.Round(AdvancedMath.GoldenRatio * n)) {
                Console.WriteLine("n={0}", n);

                Interval[] box = new Interval[n];
                for (int i = 0; i < box.Length; i++) box[i] = Interval.FromEndpoints(-500.0, 500.0);

                MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box);

                Assert.IsTrue(minimum.Dimension == n);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine("{0} ({1}) ?= {2}", minimum.Value, minimum.Precision, -418.983 * n);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, -418.983 * n, new EvaluationSettings() { AbsolutePrecision = minimum.Precision }));

                foreach (double coordinate in minimum.Location) Console.WriteLine(coordinate);
                ColumnVector solution = new ColumnVector(n);
                solution.Fill((int i, int j) => 420.969);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Location, solution, new EvaluationSettings() { RelativePrecision = 2.0 * Math.Sqrt(Math.Abs(minimum.Precision / minimum.Value)) }));
            }

        }

        [TestMethod]
        public void Griewank () {

            // See http://mathworld.wolfram.com/GriewankFunction.html

            for (int n = 2; n < 8; n++) {

                Console.WriteLine(n);

                Func<IList<double>, double> function = (IList<double> x) => {
                    double s = 0.0;
                    double p = 1.0;
                    for (int i = 0; i < x.Count; i++) {
                        s += MoreMath.Sqr(x[i]);
                        p *= Math.Cos(x[i] / Math.Sqrt(i + 1.0));
                    }
                    return (1.0 + s / 4000.0 - p);
                };

                Interval[] box = new Interval[n];
                for (int i = 0; i < n; i++) box[i] = Interval.FromEndpoints(-100.0, 100.0);

                EvaluationSettings settings = new EvaluationSettings() { AbsolutePrecision = 1.0E-6, EvaluationBudget = 1000000 };

                MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box, settings);

                Console.WriteLine(minimum.Dimension);
                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine("{0} ({1}) ?= 0.0", minimum.Value, minimum.Precision);
                Console.WriteLine("{0} {1} ...", minimum.Location[0], minimum.Location[1]);

                // We usually find the minimum at 0, but sometimes land a valley or two over.

            }

        }

        [TestMethod]
        public void Bukin () {

            // Burkin has a narrow valley, not aligned with any axis, punctuated with many tiny "wells" along its bottom.
            // The deepest well is at (-10,1)-> 0.
            Func<IList<double>, double> function = (IList<double> x) => 100.0 * Math.Sqrt(Math.Abs(x[1] - 0.01 * x[0] * x[0])) + 0.01 * Math.Abs(x[0] + 10.0);

            IList<Interval> box = new Interval[] { Interval.FromEndpoints(-15.0, 0.0), Interval.FromEndpoints(-3.0, 3.0) };

            EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 0.0, AbsolutePrecision = 1.0E-4, EvaluationBudget = 1000000 };
            /*
            settings.Update += (object result) => {
                MultiExtremum e = (MultiExtremum) result;
                Console.WriteLine("After {0} evaluations, best value {1}", e.EvaluationCount, e.Value);
            };
            */
            MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(function, box, settings);

            Console.WriteLine(minimum.EvaluationCount);
            Console.WriteLine("{0} ({1})", minimum.Value, minimum.Precision);
            Console.WriteLine("{0} {1}", minimum.Location[0], minimum.Location[1]);

            // We do not end up finding the global minimum.

        }

        [TestMethod]
        public void EasomLocal () {

            Func<IList<double>, double> function = (IList<double> x) => Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Exp(-(MoreMath.Sqr(x[0] - Math.PI) + MoreMath.Sqr(x[1] - Math.PI)));

            // We can't start too far from minimum, since cosines introduce many local minima.
            ColumnVector start = new ColumnVector(1.5, 2.0);

            MultiExtremum maximum = MultiFunctionMath.FindLocalMaximum(function, start);

            Console.WriteLine(maximum.Value);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Value, 1.0, new EvaluationSettings() { AbsolutePrecision = 2.0 * maximum.Precision }));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Location, new ColumnVector(Math.PI, Math.PI), new EvaluationSettings { AbsolutePrecision = 2.0 * Math.Sqrt(maximum.Precision) }));

        }

        [TestMethod]
        public void EasomGlobal () {

            Func<IList<double>, double> function = (IList<double> x) => Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Exp(-(MoreMath.Sqr(x[0] - Math.PI) + MoreMath.Sqr(x[1] - Math.PI)));
            IList<Interval> box = new Interval[] { Interval.FromEndpoints(-10.0, 10.0), Interval.FromEndpoints(-10.0, 10.0) };

            MultiExtremum maximum = MultiFunctionMath.FindGlobalMaximum(function, box);

            Console.WriteLine(maximum.EvaluationCount);
            Console.WriteLine(maximum.Value);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Value, 1.0, new EvaluationSettings() { AbsolutePrecision = 2.0 * maximum.Precision }));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Location, new ColumnVector(Math.PI, Math.PI), new EvaluationSettings { AbsolutePrecision = 2.0 * Math.Sqrt(maximum.Precision) }));


        }

        [TestMethod]
        public void PackCirclesInSquare () {

            // Put n points into the unit square. Place them so as to maximize the minimum distance between them.
            // See http://en.wikipedia.org/wiki/Circle_packing_in_a_square, http://hydra.nat.uni-magdeburg.de/packing/csq/csq.html
            // This is a great test problem because it is simple to understand, hard to solve,
            // solutions are known/proven for many n, and it covers many dimensions.

            double[] solutions = new double[] {
                Double.NaN,
                Double.NaN,
                Math.Sqrt(2.0), /* at opposite corners */
                Math.Sqrt(6.0) - Math.Sqrt(2.0),
                1.0, /* at each corner */
                Math.Sqrt(2.0) / 2.0 /* at each corner and in center */,
                Math.Sqrt(13.0) / 6.0,
                4.0 - 2.0 * Math.Sqrt(3.0),
                (Math.Sqrt(6.0) - Math.Sqrt(2.0)) / 2.0,
                1.0 / 2.0, /* 3 X 3 grid */
                0.421279543983903432768821760651,
                0.398207310236844165221512929748,
                Math.Sqrt(34.0) / 15.0,
                0.366096007696425085295389370603,
                2.0 / (4.0 + Math.Sqrt(3.0)),
                2.0 / (Math.Sqrt(6.0) + 2.0 + Math.Sqrt(2.0)),
                1.0 / 3.0
            };

            //int n = 11;
            for (int n = 2; n < 8; n++)
            {

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
                    return (sMin);
                };

                IList<Interval> box = new Interval[2 * n];
                for (int i = 0; i < box.Count; i++) box[i] = Interval.FromEndpoints(0.0, 1.0);

                EvaluationSettings settings = new EvaluationSettings() { RelativePrecision = 1.0E-4, AbsolutePrecision = 1.0E-6, EvaluationBudget = 10000000 };

                MultiExtremum maximum = MultiFunctionMath.FindGlobalMaximum(function, box, settings);

                Console.WriteLine(maximum.EvaluationCount);
                Console.WriteLine("{0} {1}", solutions[n], maximum.Value);

                Assert.IsTrue(maximum.Dimension == 2 * n);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Value, solutions[n], new EvaluationSettings() { AbsolutePrecision = 2.0 * maximum.Precision })); 

            }

        }

        [TestMethod]
        public void PackCirclesInCircle () {

            // Put n points into the unit circle. Place them so as to maximize the minimum distance between them.
            // http://en.wikipedia.org/wiki/Circle_packing_in_a_circle, http://hydra.nat.uni-magdeburg.de/packing/csq/csq.html
            // This can also be ecoded via a box constraint if we use polar coordinates.
            // In that case, 0 < r < 1 and -\pi < \theta < +\pi are the coordinate bounds.

            double[] solutions = {
                Double.NaN,
                Double.NaN,
                2.0, /* uniformly along edge, i.e. opposite sides */
                Math.Sqrt(3.0), /* uniformly along edge, i.e. triangle */
                Math.Sqrt(2.0), /* uniformly along edge, i.e. square */
                Math.Sqrt(2.0 / (1.0 + 1.0 / Math.Sqrt(5.0))), /* uniformly along edge, i.e. pentagon */
                1.0, /* uniformly along edge, i.e. hexagon */
                1.0, /* six uniformly along edge plus one in center */
                2.0 * Math.Sin(Math.PI / 7.0),
                Math.Sqrt(2.0 / (2.0 + Math.Sqrt(2.0))),
                0.710978235561246550830725976904,
                2.0 * Math.Sin(Math.PI / 9.0)
            };

            for (int n = 2; n < 8; n++)
            {
                // We have a problem with n = 5. Chooses to put one in middle and others at 90 degree angles around it instead of distributing uniformly along edge.
                if (n == 5) continue;
                
                Console.WriteLine(n);

                // We assume that the point r=1, t=0 exists and encode only the remaining n-1 points in a 2(n-1)-length vector.

                Func<IList<double>, double> f = (IList<double> u) => {
                    double sMin = Double.MaxValue;
                    for (int i = 0; i < (n - 1); i++) {
                        double ri = u[2 * i];
                        double ti = u[2 * i + 1];
                        double xi = ri * Math.Cos(ti);
                        double yi = ri * Math.Sin(ti);
                        for (int j = 0; j < i; j++) {
                            double rj = u[2 * j];
                            double tj = u[2 * j + 1];
                            double xj = rj * Math.Cos(tj);
                            double yj = rj * Math.Sin(tj);
                            double s = MoreMath.Hypot(xi - xj, yi - yj);
                            if (s < sMin) sMin = s;
                        }
                        // also compare distance to fixed point (1, 0)
                        double s1 = MoreMath.Hypot(xi - 1.0, yi);
                        if (s1 < sMin) sMin = s1;
                    }
                    return (sMin);
                };

                Interval[] box = new Interval[2 * (n - 1)];
                for (int i = 0; i < (n - 1); i++) {
                    box[2 * i] = Interval.FromEndpoints(0.0, 1.0);
                    box[2 * i + 1] = Interval.FromEndpoints(-Math.PI, Math.PI);
                }

                MultiExtremum maxmimum = MultiFunctionMath.FindGlobalMaximum(f, box);

                Console.WriteLine(maxmimum.EvaluationCount);
                Console.WriteLine("{0} {1}", maxmimum.Value, maxmimum.Precision);
                for (int i = 0; i < (n - 1); i++) Console.WriteLine("{0} {1}", maxmimum.Location[2 * i], maxmimum.Location[2 * i + 1]);
                Assert.IsTrue(maxmimum.Dimension == 2 * (n - 1));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(maxmimum.Value, solutions[n], new EvaluationSettings() { AbsolutePrecision = 2.0 * maxmimum.Precision }));

            }

        }

        [TestMethod]
        public void PackSpheresInCube () {

            // See http://www.combinatorics.org/ojs/index.php/eljc/article/download/v11i1r33/pdf

            double[] solutions = new double[] {
                Double.NaN,
                Double.NaN,
                Math.Sqrt(3.0), /* opposite corners */
                Math.Sqrt(2.0), /* three non-adjacent corners*/
                Math.Sqrt(2.0), /* on opposite faces in opposite corners */
                Math.Sqrt(5.0) / 2.0,
                3.0 * Math.Sqrt(2.0) / 4.0,
                Math.Sqrt(2.0 / 3.0 * (8.0 + Math.Sqrt(3.0) - 2.0 * Math.Sqrt(10.0 + 4.0 * Math.Sqrt(3.0)))),
                1.0, /* in each corner */
                Math.Sqrt(3.0) / 2.0, /* in each corner and at center */
                3.0 / 4.0,
                0.710116382462,
                0.707106806467,
                1.0 / Math.Sqrt(2.0),
                1.0 / Math.Sqrt(2.0),
                5.0 / 8.0
            };

            //int n = 11;
            for (int n = 2; n < 8; n++)
            {

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
                    return (sMin);
                };

                IList<Interval> box = new Interval[3 * n];
                for (int i = 0; i < box.Count; i++) box[i] = Interval.FromEndpoints(0.0, 1.0);

                MultiExtremum maximum = MultiFunctionMath.FindGlobalMaximum(function, box);

                Console.WriteLine(maximum.EvaluationCount);
                Console.WriteLine("{0} {1} ({2})", solutions[n], maximum.Value, maximum.Precision);

                Assert.IsTrue(maximum.Dimension == 3 * n);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Value, solutions[n], new EvaluationSettings() { AbsolutePrecision = 2.0 * maximum.Precision }));

            }
        }

        // The Thompson problem asks for the minimum-energy configuration of n unit charges confined to the surface of a sphere
        // See http://en.wikipedia.org/wiki/Thomson_problem

        private double[] thompsonSolutions = new double[] {
                Double.NaN,
                Double.NaN,
                1.0 / 2.0, /* n = 2: antipodal */
                1.732050808, /* n = 3: equilateral triangle in great circle */
                6.0 * Math.Sqrt(3.0 / 8.0), /* n = 4: tetrahedron */
                1.0 / 2.0 + Math.Sqrt(3.0) + 3.0 * Math.Sqrt(2.0), /* n = 5: triangular bipyramid */
                12.0 / Math.Sqrt(2.0) + 3.0 / 2.0, /* n = 6: octahedron */
                14.452977414,
                19.675287861,
                25.759986531,
                32.716949460,
                40.596450510,
                49.165253058 /* n = 12: isosahedron */
        };

        private Func<IList<double>, double> GetThompsonFunction (int n) {

            // Define a function that assumes one point at spherical coordinates (0, 0) and others at given 2 * (n - 1) spherical coordinates

            Func<IList<double>, double> f = (IList<double> u) => {

                // add a point at 0,0
                double[] v = new double[2 * n];
                u.CopyTo(v, 0);

                double e = 0.0;

                // iterate over all distinct pairs of points
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < i; j++) {
                        // compute the chord length between points i and j
                        double dx = Math.Cos(v[2 * j]) * Math.Cos(v[2 * j + 1]) - Math.Cos(v[2 * i]) * Math.Cos(v[2 * i + 1]);
                        double dy = Math.Cos(v[2 * j]) * Math.Sin(v[2 * j + 1]) - Math.Cos(v[2 * i]) * Math.Sin(v[2 * i + 1]);
                        double dz = Math.Sin(v[2 * j]) - Math.Sin(v[2 * i]);
                        double d = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                        e += 1.0 / d;
                    }
                }

                return (e);

            };

            return (f);

        }


        [TestMethod]
        public void ThompsonGlobal () {

            for (int n = 2; n < 8; n++)
            {
                Console.WriteLine("n={0}", n);

                Func<IList<double>, double> f = GetThompsonFunction(n);

                Interval[] box = new Interval[2 * (n - 1)];
                for (int i = 0; i < (n - 1); i++) {
                    box[2 * i] = Interval.FromEndpoints(-Math.PI, Math.PI);
                    box[2 * i + 1] = Interval.FromEndpoints(-Math.PI / 2.0, Math.PI / 2.0);
                }

                MultiExtremum minimum = MultiFunctionMath.FindGlobalMinimum(f, box);

                Console.WriteLine(minimum.EvaluationCount);
                Console.WriteLine("{0} ({1})", minimum.Value, minimum.Precision);
                for (int i = 0; i < (n - 1); i++) {
                    Console.WriteLine("{0} {1}", minimum.Location[2 * i], minimum.Location[2 * i + 1]);
                }
                Console.WriteLine("0.0 0.0");

                Assert.IsTrue(minimum.Dimension == 2 * (n - 1));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(minimum.Value, thompsonSolutions[n], new EvaluationSettings() { AbsolutePrecision = 2.0 * minimum.Precision }));

            }

        }

        /*
        [TestMethod]
        public void ThomsonLocalOld () {


            for (int n = 2; n < 9; n++) {
                Console.WriteLine(n);

                // define the thompson metric
                int count = 0;
                Func<double[], double> f2 = (double[] u) => {

                    // add a point at 0,0
                    double[] v = new double[2 * n];
                    u.CopyTo(v, 0);

                    double e = 0.0;

                    // iterate over all distinct pairs of points
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < i; j++) {
                            // compute the chord length between points i and j
                            double dx = Math.Cos(v[2 * j]) * Math.Cos(v[2 * j + 1]) - Math.Cos(v[2 * i]) * Math.Cos(v[2 * i + 1]);
                            double dy = Math.Cos(v[2 * j]) * Math.Sin(v[2 * j + 1]) - Math.Cos(v[2 * i]) * Math.Sin(v[2 * i + 1]);
                            double dz = Math.Sin(v[2 * j]) - Math.Sin(v[2 * i]);
                            double d = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                            e += 1.0 / d;
                        }
                    }

                    count++;
                    return (e);

                };

                // random distribution to start
                // using antipodal pairs gives us a better starting configuration
                Random r = new Random(1001110000);
                double[] start = new double[2 * (n - 1)];
                for (int i = 0; i < (n - 1) / 2; i++) {
                    int j = 4 * i;
                    start[j] = -Math.PI + 2.0 * r.NextDouble() * Math.PI;
                    start[j + 1] = Math.Asin(2.0 * r.NextDouble() - 1.0);
                    start[j + 2] = -(Math.PI - start[j]);
                    start[j + 3] = -start[j + 1];
                }
                // add one more point if necessary
                if (n % 2 == 0) {
                    start[2 * n - 4] = -Math.PI + 2.0 * r.NextDouble() * Math.PI;
                    start[2 * n - 3] = Math.Asin(2.0 * r.NextDouble() - 1.0);
                }

                //EvaluationSettings s = new EvaluationSettings() { RelativePrecision = 1.0E-8, AbsolutePrecision = 1.0E-12, EvaluationBudget = 10000 };

                EvaluationSettings set = new EvaluationSettings() { RelativePrecision = 1.0E-9, EvaluationBudget = 50000 };
                SpaceExtremum min = FunctionMath.FindMinimum(f2, start, set);
                double tol = set.RelativePrecision * min.Value;

                Console.WriteLine(min.Dimension);
                Console.WriteLine(count);
                Console.WriteLine("{0} ({1}) ?= {2}", min.Value, tol, thompsonSolutions[n]);

                Assert.IsTrue(min.Dimension == 2 * (n - 1));
                Console.WriteLine(TestUtilities.IsNearlyEqual(min.Value, thompsonSolutions[n], new EvaluationSettings() { AbsolutePrecision = 4.0 * tol }));

            }
                
        }
        */

        [TestMethod]
        public void ThomsonLocal () {

            for (int n = 2; n < 8; n++)
            {
                Console.WriteLine(n);

                // define the thompson metric
                Func<IList<double>, double> f = GetThompsonFunction(n);

                // random distribution to start
                // using antipodal pairs gives us a better starting configuration
                Random r = new Random(1001110000);
                double[] start = new double[2 * (n - 1)];
                for (int i = 0; i < (n - 1) / 2; i++) {
                    int j = 4 * i;
                    start[j] = -Math.PI + 2.0 * r.NextDouble() * Math.PI;
                    start[j + 1] = Math.Asin(2.0 * r.NextDouble() - 1.0);
                    start[j + 2] = -(Math.PI - start[j]);
                    start[j + 3] = -start[j + 1];
                }
                // add one more point if necessary
                if (n % 2 == 0) {
                    start[2 * n - 4] = -Math.PI + 2.0 * r.NextDouble() * Math.PI;
                    start[2 * n - 3] = Math.Asin(2.0 * r.NextDouble() - 1.0);
                }

                EvaluationSettings set = new EvaluationSettings() { RelativePrecision = 1.0E-9 };
                MultiExtremum min = MultiFunctionMath.FindLocalMinimum(f, start, set);

                Console.WriteLine(min.Dimension);
                Console.WriteLine(min.EvaluationCount);
                Console.WriteLine("{0} ({1}) ?= {2}", min.Value, min.Precision, thompsonSolutions[n]);

                Assert.IsTrue(min.Dimension == 2 * (n - 1));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(min.Value, thompsonSolutions[n], new EvaluationSettings() { AbsolutePrecision = 4.0 * min.Precision }));

            }

        }

        // See also
        // http://en.wikipedia.org/wiki/Tammes_problem, http://mathworld.wolfram.com/SphericalCode.html

        // Known solutions are d_2 = 2, d_3 = \sqrt{3}, d_4 = \sqrt{8/3}, d_5 = d_6 = \sqrt{2}
        // http://neilsloane.com/packings/ lists best-known solutions. They give \alpha, the angle between closest points, which can be converted
        // to distance via d = 2 \sin(\alpha / 2).

        // Local minimization does this very poorly, probably because there are discontinutities where the minimum-distance pair changes

        [TestMethod]
        public void TammesGlobal () {

            double[] tammesSolutions = new double[] { Double.NaN, Double.NaN, 2.0, Math.Sqrt(3.0), Math.Sqrt(8.0 / 3.0), Math.Sqrt(2.0), Math.Sqrt(2.0) };

            for (int n = 2; n < 8; n++) {
                Console.WriteLine("n = {0}", n);

                // Define a function that returns the minimum distance between points
                Func<IList<double>, double> f = (IList<double> u) => {
                    // add a point at 0,0
                    double[] v = new double[2 * n];
                    u.CopyTo(v, 0);

                    double dmin = Double.PositiveInfinity;

                    // iterate over all distinct pairs of points
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < i; j++) {
                            // compute the chord length between points i and j
                            double dx = Math.Cos(v[2 * j]) * Math.Cos(v[2 * j + 1]) - Math.Cos(v[2 * i]) * Math.Cos(v[2 * i + 1]);
                            double dy = Math.Cos(v[2 * j]) * Math.Sin(v[2 * j + 1]) - Math.Cos(v[2 * i]) * Math.Sin(v[2 * i + 1]);
                            double dz = Math.Sin(v[2 * j]) - Math.Sin(v[2 * i]);
                            double d = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                            if (d < dmin) dmin = d;
                        }
                    }

                    return (dmin);
                };

                /*
                // Start with a random arrangement
                Random r = new Random(1001110000);
                double[] start = new double[2 * n];
                for (int i = 0; i < n; i++) {
                    start[2 * i] = -Math.PI + 2.0 * r.NextDouble() * Math.PI;
                    start[2 * i + 1] = Math.Asin(2.0 * r.NextDouble() - 1.0);
                }

                // Maximize the minimum distance
                MultiExtremum maximum = MultiFunctionMath.FindLocalMinimum(f, start);
                */
                
                Interval[] limits = new Interval[2 * (n - 1)];
                for (int i = 0; i < (n - 1); i++) {
                    limits[2 * i] = Interval.FromEndpoints(-Math.PI, +Math.PI);
                    limits[2 * i + 1] = Interval.FromEndpoints(-Math.PI / 2.0, +Math.PI / 2.0);
                }
                MultiExtremum maximum = MultiFunctionMath.FindGlobalMaximum(f, limits);

                Console.WriteLine(maximum.EvaluationCount);
                Console.WriteLine("{0} ({1})", maximum.Value, maximum.Precision);

                Assert.IsTrue(maximum.Dimension == 2 * (n - 1));

                if (n < tammesSolutions.Length) Assert.IsTrue(TestUtilities.IsNearlyEqual(maximum.Value, tammesSolutions[n], new EvaluationSettings() { AbsolutePrecision = 2.0 * maximum.Precision }));

            }

        }

    }
}
