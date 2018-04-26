using System;
using System.Diagnostics;

using Meta.Numerics;

namespace Examples {
    
    public static class Core {

        [ExampleMethod]
        public static void MoreMathFunctions () {

            double x = 1.0E-12;

            Console.WriteLine($"Exp({x}) - 1.0 = {Math.Exp(x) - 1.0}");
            Console.WriteLine($"ExpMinusOne({x}) = {MoreMath.ExpMinusOne(x)}");

            Console.WriteLine($"Log(1.0 + {x}) = {Math.Log(1.0 + x)}");
            Console.WriteLine($"LogOnePlus({x}) = {MoreMath.LogOnePlus(x)}");

            double z = 1.0E199;
            Console.WriteLine($"Math.Sin({z}) = {Math.Sin(z)}");
            Console.WriteLine($"MoreMath.Sin({z}) = {MoreMath.Sin(z)}");

            Console.WriteLine($"sqrt({z}^2 + 1.0) = {Math.Sqrt(z * z + 1.0)}");
            Console.WriteLine($"Hypot({z}, 1) = {MoreMath.Hypot(z, 1.0)}");

            int n = 128;
            double y = 1.0 - 1.0 / n;
            double s1 = 0.0;
            Stopwatch t1 = Stopwatch.StartNew();
            for (int i = 0; i < 10000000; i++) {
                s1 += Math.Pow(y, n);
            }
            t1.Stop();
            Console.WriteLine($"{s1} in {t1.ElapsedMilliseconds}ms");

            double s2 = 0.0;
            Stopwatch t2 = Stopwatch.StartNew();
            for (int i = 0; i < 10000000; i++) {
                s2 += MoreMath.Pow(y, n);
            }
            t2.Stop();
            Console.WriteLine($"{s2} in {t2.ElapsedMilliseconds}ms");

        }

    }


}