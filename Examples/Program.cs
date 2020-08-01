using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace Examples
{
    // Add this attribute to static methods to make
    // them available as examples.
    internal class ExampleMethodAttribute : Attribute { }

    class Program
    {
        private static MethodInfo[] GetExampleMethods () {
            Assembly assembly = Assembly.GetExecutingAssembly();
            MethodInfo[] methods = assembly.GetTypes()
                .SelectMany(t => t.GetMethods())
                .Where(m => m.GetCustomAttributes(typeof(ExampleMethodAttribute), false).Length > 0)
                .ToArray();
            return(methods);
        }

        private static double EllipticPi(double n, double k) {

            if (n > 1.0) throw new ArgumentOutOfRangeException(nameof(n));
            if (k < 1.0 || k > 1.0) throw new ArgumentOutOfRangeException(nameof(k));

            if (n == 1.0 || k == 1.0) return Double.PositiveInfinity;

            // DLMF 19.8 describes how to compute \Pi(n, k) via AGM plus some auxiluary
            // calculations. Here a and g, along with p,  converge to AGM(1,k') in the
            // usual way; the sum of auxiluary variables q computed along the way gives \Pi(n, k).
            // This method appears to have been derived by Carlson in "Three Improvements in
            // Reduction and Computation of Elliptic Integrals", Journal of Research of the
            // National Institute of Standards and Technology, 2002 Sep-Oct 107(5): 413-418
            // (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4861378/)
            // as a specialization of his method to calcuate R_J.

            double a = 1.0;
            double g = Math.Sqrt((1.0 - k) * (1.0 + k));
            double n1 = 1.0 - n;
            double pSquared = n1;
            double p = Math.Sqrt(pSquared);
            double q = 1.0;
            double s = q;

            for (int i = 1; i < 100; i++) {

                double s_old = s;
                double ag = a * g;
                double e = (pSquared - ag) / (pSquared + ag);
                q = 0.5 * q * e;
                s += q;

                double p_old = p;
                p = (pSquared + ag ) / (2.0 * p);

                if (p == p_old && s == s_old) {
                    return Math.PI / 4.0 / p * (2.0 + n / n1 * s);
                }

                pSquared = p * p;
                a = 0.5 * (a + g);
                g = Math.Sqrt(ag);

            }

            throw new Exception();

        }


        static void Main(string[] args)
        {
            double e3 = EllipticPi(-0.25,0.99);

            MethodInfo[] methods = GetExampleMethods();
            Dictionary<string, MethodInfo> index = new Dictionary<string, MethodInfo>();
            foreach (MethodInfo method in methods) {
                index.Add(method.Name, method);
            }

            if (args.Length == 0) {
                foreach(string key in index.Keys) {
                    Console.WriteLine(key);
                }
            } else {
                foreach (string arg in args) {
                    MethodInfo method = index[arg];
                    method.Invoke(null, null);
                }
            }

        }
    }
}
