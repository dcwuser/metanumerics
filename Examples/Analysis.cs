using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;

namespace Examples {
    
    public static class Analysis {

        [ExampleMethod]
        public static void Integration () {

            // This is the area of the upper half-circle, so the value should be pi/2.
            IntegrationResult i = FunctionMath.Integrate(x => Math.Sqrt(1.0 - x * x), -1.0, +1.0);
            Console.WriteLine($"i = {i.Estimate} after {i.EvaluationCount} evaluations.");

            Interval zeroToPi = Interval.FromEndpoints(0.0, Math.PI);
            IntegrationResult watson = MultiFunctionMath.Integrate(
                x => 1.0 / (1.0 - Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Cos(x[2])),
                new Interval[] { zeroToPi, zeroToPi, zeroToPi }
            );

            IntegrationResult i2 = FunctionMath.Integrate(
                x => Math.Sqrt(1.0 - x * x), -1.0, +1.0,
                new IntegrationSettings() { RelativePrecision = 1.0E-4 }
            );

            // This integrand has a log singularity at x = 0.
            // The value of this integral is the Catalan constant.
            IntegrationResult soft = FunctionMath.Integrate(
                x => -Math.Log(x) / (1.0 + x * x), 0.0, 1.0
            );

            // This integral has a power law singularity at x = 0.
            // The value of this integral is 4.
            IntegrationResult hard = FunctionMath.Integrate(
                x => Math.Pow(x, -3.0 / 4.0), 0.0, 1.0,
                new IntegrationSettings() { RelativePrecision = 1.0E-6 }
            );

            // The value of this infinite integral is sqrt(pi)
            IntegrationResult infinite = FunctionMath.Integrate(
                x => Math.Exp(-x* x), Double.NegativeInfinity, Double.PositiveInfinity
            );

            Func<IReadOnlyList<double>,double> distance = z => {
                ColumnVector x = new ColumnVector(z[0], z[1], z[2]);
                ColumnVector y = new ColumnVector(z[3], z[4], z[5]);
                ColumnVector d = x - y;
                return(d.Norm());
            };

            Interval oneBox = Interval.FromEndpoints(0.0, 1.0);
            Interval[] sixBox = new Interval[] { oneBox, oneBox, oneBox, oneBox, oneBox, oneBox };

            IntegrationSettings settings= new IntegrationSettings() {
                RelativePrecision = 1.0E-4,
                AbsolutePrecision = 0.0,
                Listener = r => {
                    Console.WriteLine($"Estimate {r.Estimate} after {r.EvaluationCount} evaluations.");
                }
            };

            IntegrationResult numeric = MultiFunctionMath.Integrate(distance, sixBox, settings);
            Console.WriteLine($"The numeric result is {numeric.Estimate}.");

            double analytic = 4.0 / 105.0 + 17.0 / 105.0 * Math.Sqrt(2.0) - 2.0 / 35.0 * Math.Sqrt(3.0)
                + Math.Log(1.0 + Math.Sqrt(2.0)) / 5.0 + 2.0 / 5.0 * Math.Log(2.0 + Math.Sqrt(3.0))
                - Math.PI / 15.0;
            Console.WriteLine($"The analytic result is {analytic}.");

        }

        [ExampleMethod]
        public static void Optimization () {

            Interval range = Interval.FromEndpoints(0.0, 6.28);
            Extremum min = FunctionMath.FindMinimum(x => Math.Sin(x), range);
            Console.WriteLine($"Found minimum of {min.Value} at {min.Location}.");
            Console.WriteLine($"Required {min.EvaluationCount} evaluations.");

            MultiExtremum rosenbrock = MultiFunctionMath.FindLocalMinimum(
                x => MoreMath.Sqr(2.0 - x[0]) + 100.0 * MoreMath.Sqr(x[1] - x[0] * x[0]),
                new ColumnVector(0.0, 0.0)
            );
            ColumnVector xm = rosenbrock.Location;
            Console.WriteLine($"Found minimum of {rosenbrock.Value} at ({xm[0]},{xm[1]}).");
            Console.WriteLine($"Required {rosenbrock.EvaluationCount} evaluations.");
         
            Func<IReadOnlyList<double>, double> leviFunction = z => {
                double x = z[0];
                double y = z[1];
                return(
                    MoreMath.Sqr(MoreMath.Sin(3.0 * Math.PI * x)) +
                    MoreMath.Sqr(x - 1.0) * (1.0 + MoreMath.Sqr(MoreMath.Sin(3.0 * Math.PI * y))) +
                    MoreMath.Sqr(y - 1.0) * (1.0 + MoreMath.Sqr(MoreMath.Sin(2.0 * Math.PI * y)))
                );
            };

            Interval[] leviRegion = new Interval[] {
                Interval.FromEndpoints(-10.0, +10.0),
                Interval.FromEndpoints(-10.0, +10.0)
            };

            MultiExtremum levi = MultiFunctionMath.FindGlobalMinimum(leviFunction, leviRegion);
            
            ColumnVector zm = levi.Location;
            Console.WriteLine($"Found minimum of {levi.Value} at ({zm[0]},{zm[1]}).");
            Console.WriteLine($"Required {levi.EvaluationCount} evaluations.");

            // Select a dimension
            int n = 6;

            // Define a function that takes 2n polar coordinates in the form
            // phi_1, theta_1, phi_2, theta_2, ..., phi_n, theta_n, computes
            // the sum of the potential energy 1/d for all pairs.
            Func<IReadOnlyList<double>, double> thompson = (IReadOnlyList<double> v) => {
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

            // Limit all polar coordinates to their standard ranges.
            Interval[] box = new Interval[2 * n];
            for (int i = 0; i < n; i++) {
                box[2 * i] = Interval.FromEndpoints(-Math.PI, Math.PI);
                box[2 * i + 1] = Interval.FromEndpoints(-Math.PI / 2.0, Math.PI / 2.0);
            }

            // Use settings to monitor proress toward a rough solution.
            MultiExtremumSettings roughSettings = new MultiExtremumSettings() {
                RelativePrecision = 0.05, AbsolutePrecision = 0.0,
                Listener = r => {
                    Console.WriteLine($"Minimum {r.Value} after {r.EvaluationCount} evaluations.");
                }
            };
            MultiExtremum roughThompson = MultiFunctionMath.FindGlobalMinimum(thompson, box, roughSettings);

            // Use settings to monitor proress toward a refined solution.
            MultiExtremumSettings refinedSettings = new MultiExtremumSettings() {
                RelativePrecision = 1.0E-5, AbsolutePrecision = 0.0,
                Listener = r => {
                    Console.WriteLine($"Minimum {r.Value} after {r.EvaluationCount} evaluations.");
                }
            };
            MultiExtremum refinedThompson = MultiFunctionMath.FindLocalMinimum(thompson, roughThompson.Location, refinedSettings);

            Console.WriteLine($"Minimum potential energy {refinedThompson.Value}.");
            Console.WriteLine($"Required {roughThompson.EvaluationCount} + {refinedThompson.EvaluationCount} evaluations.");

            /*
            // Define a function that takes 2n coordinates x1, y1, x2, y2, ... xn, yn
            // and finds the smallest distance between two coordinate pairs.
            Func<IReadOnlyList<double>, double> function = (IReadOnlyList<double> x) => {
                double sMin = Double.MaxValue;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < i; j++) {
                        double s = MoreMath.Hypot(x[2 * i] - x[2 * j], x[2 * i + 1] - x[2 * j + 1]);
                        if (s < sMin) sMin = s;
                    }
                }
                return (sMin);
            };

            // Limit all coordinates to the unit box.
            Interval[] box = new Interval[2 * n];
            for (int i = 0; i < box.Length; i++) box[i] = Interval.FromEndpoints(0.0, 1.0);

            // Use settings to monitor proress toward a rough solution.
            MultiExtremumSettings roughSettings = new MultiExtremumSettings() {
                RelativePrecision = 1.0E-2, AbsolutePrecision = 0.0,
                Listener = r => {
                    Console.WriteLine($"Minimum {r.Value} after {r.EvaluationCount} evaluations.");
                }
            };
            MultiExtremum roughMaximum = MultiFunctionMath.FindGlobalMaximum(function, box, roughSettings);

            // Use settings to monitor proress toward a rough solution.
            MultiExtremumSettings refinedSettings = new MultiExtremumSettings() {
                RelativePrecision = 1.0E-8, AbsolutePrecision = 0.0,
                Listener = r => {
                    Console.WriteLine($"Minimum {r.Value} after {r.EvaluationCount} evaluations.");
                }
            };
            MultiExtremum refinedMaximum = MultiFunctionMath.FindLocalMaximum(function, roughMaximum.Location, refinedSettings);
            */
        }

        [ExampleMethod]
        public static void IntegrateOde () {

            Func<double, double, double> rhs = (x, y) => - x * y;
            OdeResult sln = FunctionMath.IntegrateOde(rhs, 0.0, 1.0, 2.0);
            Console.WriteLine($"Numeric solution y({sln.X}) = {sln.Y}.");
            Console.WriteLine($"Required {sln.EvaluationCount} evaluations.");
            Console.WriteLine($"Analytic solution y({sln.X}) = {Math.Exp(-MoreMath.Sqr(sln.X) / 2.0)}");

            // Lotka-Volterra equations
            double A = 0.1;
            double B = 0.02;
            double C = 0.4;
            double D = 0.02;
            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> lkRhs = (t, y) => {
                return(new double[] {
                    A * y[0] - B * y[0] * y[1], D * y[0] * y[1] - C * y[1]
                });
            };
            MultiOdeSettings lkSettings = new MultiOdeSettings() {
                Listener = r =>  {Console.WriteLine($"t={r.X} rabbits={r.Y[0]}, foxes={r.Y[1]}"); }
            };
            MultiFunctionMath.IntegrateOde(lkRhs, 0.0, new double[] {20.0, 10.0}, 50.0, lkSettings);

            Func<double, IReadOnlyList<double>, IReadOnlyList<double>> rhs1 = (x, u) => {
                return new double[] { u[1], -u[0] };
            };
            MultiOdeSettings settings1 = new MultiOdeSettings() { EvaluationBudget = 100000 };
            MultiOdeResult result1 = MultiFunctionMath.IntegrateOde(
                rhs1, 0.0, new double[] { 0.0, 1.0}, 500.0, settings1
            );
            double s1 = MoreMath.Sqr(result1.Y[0]) + MoreMath.Sqr(result1.Y[1]);
            Console.WriteLine($"y({result1.X}) = {result1.Y[0]}, (y)^2 + (y')^2 = {s1}");
            Console.WriteLine($"Required {result1.EvaluationCount} evaluations.");

            Func<double, double, double> rhs2 = (x, y) => -y;
            OdeSettings settings2 = new OdeSettings() { EvaluationBudget = 100000 };
            OdeResult result2 = FunctionMath.IntegrateConservativeOde(
                rhs2, 0.0, 0.0, 1.0, 500.0, settings2
            );
            double s2 = MoreMath.Sqr(result2.Y) + MoreMath.Sqr(result2.YPrime);
            Console.WriteLine($"y({result2.X}) = {result2.Y}, (y)^2 + (y')^2 = {s2}");
            Console.WriteLine($"Required {result2.EvaluationCount} evaluations");

            Console.WriteLine(MoreMath.Sin(500.0));
        }

    }

}