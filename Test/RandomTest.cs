using System;
using System.Text;
using System.Diagnostics;
using System.IO;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    public interface IDeviateGenerator {

        double GetNext (Random rng);

    }
    

    public class BoxMullerNormalGenerator : IDeviateGenerator {

        private bool haveNextDeviate = false;
        private double nextDeviate;

        public double GetNext (Random rng) {

            if (haveNextDeviate) {

                haveNextDeviate = false;
                return (nextDeviate);

            } else {

                // pick a point in the unit disc
                double u = rng.NextDouble();
                double t = 2.0 * Math.PI * rng.NextDouble();

                double a = Math.Sqrt(-2.0 * Math.Log(u));

                // store one deviate
                nextDeviate = a * Math.Sin(t);
                haveNextDeviate = true;

                // return the other
                return (a * Math.Cos(t));

            }

        }
    }

    public class PolarRejectionNormalDeviateGenerator : IDeviateGenerator {

        private bool haveNextDeviate = false;
        private double nextDeviate;

        public double GetNext (Random rng) {

            if (haveNextDeviate) {
                haveNextDeviate = false;
                return (nextDeviate);
            } else {

                // just as in Box-Mueller we need a point in the unit disc
                // but in order to avoid trig functions, we find it by
                // picking a point in the circumscribing square, then
                // rejecting it if it lies outside the unit circle
                double x, y, r2;
                do {
                    // pick a point in the square (-1,1)^2
                    x = 2.0 * rng.NextDouble() - 1.0;
                    y = 2.0 * rng.NextDouble() - 1.0;
                    // determine if it lies in the unit disc
                    r2 = x * x + y * y;
                } while ((r2 > 1.0) || (r2 == 0.0));

                double f = Math.Sqrt(-2.0 * Math.Log(r2)/r2);

                // store one deviate
                nextDeviate = f * x;
                haveNextDeviate = true;

                // return the other
                return (f * y);

            }

        }

    }


    public class RatioOfUniformsNormalGenerator : IDeviateGenerator {

        private static readonly double c0 = Math.Sqrt(2.0 / Math.E);
        private static readonly double c1 = 4.0 * Math.Pow(Math.E, 0.25);
        private static readonly double c2 = 4.0 * Math.Pow(Math.E, -1.35);

        public double GetNext (Random rng) {

            while (true) {

                // pic a point in (0,1) X (-\sqrt{2/e},\sqrt{2/2})
                double u = rng.NextDouble();
                double v = c0 * (2.0 * rng.NextDouble() - 1.0);

                // form ratio; this is our candidate deviate
                double x = v / u;

                // determine if the candidate is in acceptance or rejection region
                double x2 = x * x;
                if (x2 < 5.0 - c1 * u) {
                   // inside squeeze, accept immediately
                    return(x);
                } else if (x2 > c2 / u + 1.4) {
                    // outside squeeze, reject immediately
                    continue;
                }

                // between squeezes, do costly border evaluation
                if (v * v < - 4.0 * u * u * Math.Log(u)) return(x);
            }
        }
    }

    public class LevaNormalGenerator : IDeviateGenerator {

        private static readonly double c = Math.Sqrt(2.0 / Math.E);

        private const double s = 0.449871, t = -0.386595;
        private const double a = 0.19600, b = 0.25472;
        private const double r1 = 0.27597, r2 = 0.27846;

        public double GetNext (Random rng) {

            while (true) {

                double u = rng.NextDouble();
                double v = c * (2.0 * rng.NextDouble() - 1.0);

                double x = u - s;
                double y = Math.Abs(v) - t;
                double q = x * x + y * (a * y - b * x);

                if (q < r1) {
                    // inside squeeze, accept
                    return(v / u);
                } else if (q > r2) {
                    // outside squeeze, reject
                    continue;
                }

                // evaluate exact border
                if (v * v < -4.0 * u * u * Math.Log(u)) return (v / u);

            }

        }

    }

    public class CauchyGenerator : IDeviateGenerator {

        public double GetNext (Random rng) {

            // pick a random point within the unit semicircle using rejection
            double x, y;
            do {
                x = 2.0 * rng.NextDouble() - 1.0;
                y = rng.NextDouble();
            } while ((x * x + y * y > 1.0) || (y == 0.0));

            // its tangent is the tangent of a random angle and is thus Cauchy distributed
            return (x / y);

        }

    }

    public class AhrensDieterLowAlphaGammaGenerator : IDeviateGenerator {

        public AhrensDieterLowAlphaGammaGenerator (double alpha) {
            if (alpha > 1.0) throw new ArgumentOutOfRangeException("alpha");
            a = alpha;
            b = (Math.E + a) / Math.E;
        }

        double a, b;

        public double GetNext (Random rng) {

            while (true) {

                double P = b * rng.NextDouble();

                if (P < 1.0) {
                    double x = Math.Pow(P, 1.0 / a);
                    if (rng.NextDouble() > Math.Exp(-x)) continue;
                    return (x);
                } else {
                    double x = -Math.Log((b - P) / a);
                    if (rng.NextDouble() > Math.Pow(x, a - 1.0)) continue;
                    return (x);
                }

            }

        }

    }

    public class MarsagliaTsangGammaGenerator : IDeviateGenerator {

        public MarsagliaTsangGammaGenerator (IDeviateGenerator normalGenerator, double alpha) {
            if (alpha < 1.0) throw new ArgumentOutOfRangeException("alpha");
            this.normalGenerator = normalGenerator;
            a = alpha;
            a1 = a - 1.0 / 3.0;
            a2 = 1.0 / Math.Sqrt(9.0 * a1);
        }

        IDeviateGenerator normalGenerator;

        double a, a1, a2;

        public double GetNext (Random rng) {

            double u, v, z;
            while (true) {

                // generate a candidate point
                do {
                    z = normalGenerator.GetNext(rng);
                    v = 1.0 + a2 * z;
                } while (v <= 0.0);
                v = v * v * v;
                u = rng.NextDouble();

                // check whether the point is outside the acceptance region
                // first via a simple interrior squeeze
                double z2 = z * z; double z4 = z2 * z2;
                if (u > 1.0 - 0.331 * z4) {
                    // then, if necessary, via the exact rejection boundary
                    if (Math.Log(u) > z2 / 2.0 + a1 * (1.0 - v + Math.Log(v))) continue;
                }

                return (a1 * v);

            }

        }
    }


    [TestClass]
    public class RandomTest {

        [TestMethod]
        public void TimeNormalGenerators () {

            Random rng = new Random(1);
            //IDeviateGenerator nRng = new BoxMullerNormalGenerator();
            //IDeviateGenerator nRng = new PolarRejectionNormalDeviateGenerator();
            //IDeviateGenerator nRng = new RatioOfUniformsNormalGenerator();
            IDeviateGenerator nRng = new LevaNormalGenerator();

            //Sample sample = new Sample();
            ContinuousDistribution nrm = new NormalDistribution();

            Stopwatch timer = Stopwatch.StartNew();
            double sum = 0.0;
            for (int i = 0; i < 10000000; i++) {
                sum += nrm.InverseLeftProbability(rng.NextDouble());
                //sum += nRng.GetNext(rng);
                //sample.Add(nRng.GetNext(rng));
            }
            timer.Stop();

            //Console.WriteLine(sample.KolmogorovSmirnovTest(new NormalDistribution()).RightProbability);
            Console.WriteLine(sum);
            Console.WriteLine(timer.ElapsedMilliseconds);
        }

        [TestMethod]
        public void TimeCauchyGenerators () {

            Random rng = new Random(1);
            IDeviateGenerator nRng = new CauchyGenerator();
            ContinuousDistribution d = new CauchyDistribution();

            Sample sample = new Sample();

            Stopwatch timer = Stopwatch.StartNew();
            double sum = 0.0;
            for (int i = 0; i < 10000000; i++) {
                sum += nRng.GetNext(rng);
                //sum += d.InverseLeftProbability(rng.NextDouble());
                //sample.Add(nRng.GetNext(rng));
            }
            timer.Stop();

            //Console.WriteLine(sample.KolmogorovSmirnovTest(d).RightProbability);
            Console.WriteLine(sum);
            Console.WriteLine(timer.ElapsedMilliseconds);

        }

        [TestMethod]
        public void TimeGammaGenerators () {

            double alpha = 1.0;

            Random rng = new Random(1);
            //IDeviateGenerator nRng = new AhrensDieterGammaGenerator(alpha);
            IDeviateGenerator nRng = new MarsagliaTsangGammaGenerator(new PolarRejectionNormalDeviateGenerator(), alpha);
            ContinuousDistribution d = new GammaDistribution(alpha);

            //double sum = 0.0;
            Sample sample = new Sample();

            Stopwatch timer = Stopwatch.StartNew();
            for (int i = 0; i < 1000000; i++) {
                //double x = nRng.GetNext(rng);
                double x = d.InverseLeftProbability(rng.NextDouble());
                //sum += x;
                sample.Add(x);
            }
            timer.Stop();

            Console.WriteLine(sample.KolmogorovSmirnovTest(d).RightProbability);
            //Console.WriteLine(sum);
            Console.WriteLine(timer.ElapsedMilliseconds);

        }

        //[TestMethod]
        public void TestMethod1 () {

            FileStream s = File.OpenWrite(@"C:\Users\dawright\downloads\diehard\random.bin");

            Random rng = new Random(1);
            byte[] buffer = new byte[1024];
            for (int i = 0; i < 16384; i++) {
                rng.NextBytes(buffer);
                s.Write(buffer, 0, buffer.Length);
            }

            s.Close();
        }

    }
}
