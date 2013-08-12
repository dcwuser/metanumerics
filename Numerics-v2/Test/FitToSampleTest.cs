using System;
using System.Diagnostics;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {




    [TestClass]
    public class FitToSampleTest {

        [TestMethod]
        public void TestSpearmanViaTruncated () {

            int n = 12;
            double x = Math.Sqrt((n - 1) / 2.0);
            double lambda = 1.0 / Math.Sqrt((n - 1) * (1.0 - 2.0 / Math.Sqrt(Math.PI) * x * Math.Exp(-x * x) / AdvancedMath.Erf(x)));
            double D = 1.0 / lambda;

            Console.WriteLine("n={0} lambda={1} D={2}", n, lambda, D);

            double D0 = Math.Sqrt(n - 1);
            Func<double, double> fd = delegate(double d) {
                return (1.0 - Math.Sqrt(2.0 / Math.PI) * d * Math.Exp(-d * d / 2.0) / AdvancedMath.Erf(d / Math.Sqrt(2.0)) - d * d / (n - 1));
            };
            double D1 = FunctionMath.FindZero(fd, D0);
            lambda = 1.0 / D1;

            Console.WriteLine("n={0} lambda={1} D={2}", n, lambda, D1);

            TruncatedNormalDistribution t = new TruncatedNormalDistribution(D1);
            Console.WriteLine("M1 = {0}", t.Mean);
            Console.WriteLine("V={0}", t.Variance * MoreMath.Pow(lambda, 2));
            Console.WriteLine("C2={0}", t.MomentAboutMean(2) * MoreMath.Pow(lambda, 2));
            Console.WriteLine("C4={0}", t.MomentAboutMean(4) * MoreMath.Pow(lambda, 4));
            Console.WriteLine("K4={0}", (t.MomentAboutMean(4) - 3.0 * MoreMath.Pow(t.MomentAboutMean(2), 2)) * MoreMath.Pow(lambda, 4));
            Console.WriteLine("support {0}", t.Support.RightEndpoint * lambda);

        }


        public static long ISqrt (long n) {
            long op, res, one;

            op = n;
            res = 0;

            /* "one" starts at the highest power of four <= than the argument. */
            one = 1 << 30;  /* second-to-top bit set */
            while (one > op) one >>= 2;

            while (one != 0) {
                if (op >= res + one) {
                    op -= res + one;
                    res += one << 1;  // <-- faster than 2 * one
                    //res += 2 * one;
                }
                res >>= 1;
                one >>= 2;
            }
            return res;

        }

        public static long sqrt32(long n) {
            long c = 0x8000;
            long g = 0x8000;

            for(;;) {
                if(g*g > n) g ^= c;
                c >>= 1;
                if(c == 0) return g;
                g |= c;
            }
        }


        [TestMethod]
        public void TestIsqrt () {

            int n = 10000000;
            long xi = 0; double xf = 0.0;

            Stopwatch s1 = Stopwatch.StartNew();
            for (int i = 1; i <= n; i++) {
                xi = ISqrt(i);
            }
            s1.Stop();
            Console.WriteLine(s1.ElapsedMilliseconds);

            Stopwatch s2 = Stopwatch.StartNew();
            for (int i = 1; i <= n; i++) {
                xf = Math.Sqrt(i);
            }
            s2.Stop();
            Console.WriteLine(s2.ElapsedMilliseconds);

        }

        // Normal deviate review article http://www.cse.cuhk.edu.hk/~phwl/mt/public/archives/papers/grng_acmcs07.pdf

        private static readonly double cl = Math.Sqrt(2.0 / Math.E);

        public static double LevaNormalRng (Random rng) {

            const double s = 0.449871;
            const double t = -0.386595;
            const double a = 0.19600;
            const double b = 0.25472;
            const double r1 = 0.27597;
            const double r2 = 0.27846;

            while (true) {
                double u = rng.NextDouble();
                double v = rng.NextDouble();
                double x = u - s;
                double y = cl * v - t;
                double Q = x * x + y * (a * y - b * x);
                if (Q < r1) {
                    return (v / u);
                } else if ((Q < r2) && (v * v < -4.0 * u * u * Math.Log(u))) {
                    return (v / u);
                }
            }

        }

        private static double savedDeviate;

        public static double BoxMullerAccept (Random rng) {
            if (savedDeviate != 0.0) {
                double t = savedDeviate;
                savedDeviate = 0.0;
                return (t);
            }
            while (true) {
                // generate random numbers in a box
                // and reject unless they fall in the unit disc
                double x = 2.0 * rng.NextDouble() - 1.0;
                double y = 2.0 * rng.NextDouble() - 1.0;
                double r2 = x * x + y * y;
                if ((r2 < 1.0) && (r2 != 0.0)) {
                    // in the unit disc, r is random length,
                    // x/r and y/r are random sine and cosine,
                    // so we can use Box-Mueller to generate
                    // random normal deviate
                    double f = Math.Sqrt(-2.0 * Math.Log(r2) / r2);
                    savedDeviate = x * f;
                    return(y * f);
                }
            } 
        }

        [TestMethod]
        public void NormalRngTime () {

            Random rng = new Random(1);
            double x = 0.0;
            int n = 10000000;

            Stopwatch s1 = Stopwatch.StartNew();
            for (int i = 0; i < n; i++) {
                x = LevaNormalRng(rng);
            }
            s1.Stop();
            Console.WriteLine(s1.ElapsedMilliseconds);

            NormalDistribution d = new NormalDistribution();
            Stopwatch s2 = Stopwatch.StartNew();
            for (int i = 0; i < n; i++) {
                x = BoxMullerAccept(rng);
            }
            s2.Stop();
            Console.WriteLine(s2.ElapsedMilliseconds);


        }

        private double[] c = {
            1.6449340668482264365, -4.3030772285491516199,
            11.718339177218881145, -26.531419164640120037,
            57.676112859609746010, -120.62540774769279835,
            247.65840041981135684, -502.64981061634897817,
            1013.6553090068470808, -2036.6542815160842067,
            4083.6549298829502902, -8178.6547668453805814,
            16369.654755395653115, -32752.654672804177704,
            65519.654613867397313, -131054.65456242378508
        };

        [TestMethod]
        public void TestMethod() {

            double x = 1.0 / 16.0;

            Console.WriteLine(AdvancedMath.Gamma(1.0 + 2.0 * x) - MoreMath.Pow(AdvancedMath.Gamma(1.0 + x), 2));

            double f = 0.0;
            double xx = x;

            double[][] w = new double[c.Length][];
            w[0] = new double[c.Length];

            for (int i = 0; i < c.Length; i++) {
                xx = xx * x;
                f += c[i] * xx;
                Console.WriteLine(f);
                w[0][i] = f;
                if (i > 3) {
                    double s = w[0][i - 2] - MoreMath.Pow(w[0][i - 1] - w[0][i - 2], 2) / (w[0][i] - 2.0 * w[0][i - 1] + w[0][i - 2]);
                    Console.WriteLine(" -> {0}", s);
                }
            }

            for (int j = 1; j < c.Length; j++) {
                w[j] = new double[c.Length - j];
                for (int k = 0; k < w[j - 1].Length - 1; k++) {
                    w[j][k] = (w[j - 1][k] + w[j - 1][k + 1]) / 2.0;
                    Console.Write("{0}  ", w[j][k]);
                }
                Console.WriteLine();
            }

            Console.WriteLine(w[c.Length - 1][0]);

        }

    }
}
