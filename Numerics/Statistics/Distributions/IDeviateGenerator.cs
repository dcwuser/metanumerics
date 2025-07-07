using System;
using System.Diagnostics;

namespace Meta.Numerics.Statistics.Distributions {

    internal interface IDeviateGenerator<T> {

        T GetNext(Random rng);

    }

    // Revisit generators. See https://arxiv.org/html/2411.01415v1#S5.

    internal static class DeviateGeneratorFactory {

        public static IDeviateGenerator<double> GetNormalGenerator() {
            return (new BoxMullerRejectionNormalDeviateGenerator());
        }

        public static IDeviateGenerator<double> GetCauchyGenerator() {
            return (new CauchyGenerator());
        }

        public static IDeviateGenerator<double> GetGammaGenerator(double alpha) {
            // choose a gamma generator depending on whether the shape parameter is less than one or greater than one
            if (alpha < 1.0) {
                return (new BestGammaGenerator(alpha));
            } else {
                return (new MarsagliaTsangGammaGenerator(alpha));
            }
        }

        public static IDeviateGenerator<double> GetBetaGenerator(double alpha, double beta) {
            return (new BetaFromGammaGenerator(GetGammaGenerator(alpha), GetGammaGenerator(beta)));
        }

    }

    // Many normal generators are described and compared in Thomas, et al., "Gaussian Random Number Generators", ACM Computing Surveys 39 (2007)
    // Most of our normal generators are written based on the description there.

    // Timing for generating 10^7 deviates in ms:
    //   Inverse CDF            12600 ms
    //   Box-Muller             1060 ms
    //   Box-Muller Rejection   920 ms
    //   Ratio of Uniforms      1210 ms
    //   Leva's R of U          1130 ms
    // Our measurements indicate that simple Box-Muller with rejection is slightly faster than Leva's ratio of uniforms algorithm, and
    // the review article indicates that they are equally good, so we stick with Box-Muller for now. We have not yet attempted to implement
    // the Zigurraut algorithm.

    // One problem with Box-Muller is that the saving of the next value makes it stateful and therefore not thread-safe (at least not without adding
    // locking logic that would slow it down considerably). On the other hand the Random class itself isn't thread-safe, so this seems pretty minor.

    internal class BoxMullerRejectionNormalDeviateGenerator : IDeviateGenerator<double> {

        private bool haveNextDeviate;
        private double nextDeviate;

        public double GetNext(Random rng) {

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

                double f = Math.Sqrt(-2.0 * Math.Log(r2) / r2);

                // store one deviate
                nextDeviate = f * x;
                haveNextDeviate = true;

                // return the other
                return (f * y);

            }

        }

    }

#if PAST

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
                    return (x);
                } else if (x2 > c2 / u + 1.4) {
                    // outside squeeze, reject immediately
                    continue;
                }

                // between squeezes, do costly border evaluation
                if (v * v < -4.0 * u * u * Math.Log(u)) return (x);
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
                    return (v / u);
                } else if (q > r2) {
                    // outside squeeze, reject
                    continue;
                }

                // evaluate exact border
                if (v * v < -4.0 * u * u * Math.Log(u)) return (v / u);

            }

        }

    }

#endif

    // This Cauchy generator uses the same basic idea as Box-Muller with rejection. It is suggested in Numerical Recipies.
    // Our timing indicates that it is about 25% faster than simply evaluating the inverse tangent function of the inverse CDF.

    internal class CauchyGenerator : IDeviateGenerator<double> {

        public double GetNext(Random rng) {

            // Pick a random point within the unit semicircle using rejection.
            double x, y;
            do {
                x = 2.0 * rng.NextDouble() - 1.0;
                y = rng.NextDouble();
            } while ((x * x + y * y > 1.0) || (y == 0.0));

            // Its tangent is the tangent of a random angle and is thus Cauchy distributed.
            return (x / y);

        }

    }

    // Using these rejection-based gamma generators is about six times faster than evaluating the inverse CDF.

    // The Ahrens-Dieter generator for \alpha < 1 is very well described by Best
    // Best also describes two improvements

    // Best's generator is an improvement on Ahrens-Dieter with a tighter, multi-stage squeeze.
    // Best DJ, "A Note on Gamma Variate Generators with Shape Parameter less than Unity", Computing 30 (1983) 185-88
    // In my tests, Best is 10-20% faster than Ahrens-Dieter for a < 1.

    internal class BestGammaGenerator : IDeviateGenerator<double> {

        private readonly double a, b, z, ai, am1;

        public BestGammaGenerator(double alpha) {
            Debug.Assert(0.0 < alpha && alpha < 1.0);
            this.a = alpha;
            this.z = 0.07 + 0.75 * Math.Sqrt(1.0 - a);
            this.b = 1.0 + Math.Exp(-z) * a / z;
            this.ai = 1.0 / a;
            this.am1 = a - 1.0;
        }

        public double GetNext(Random rng) {
            while (true) {
                double u = rng.NextDouble();
                double P = b * u;
                if (P <= 1.0) {
                    double x = z * Math.Pow(P, ai);
                    double v = rng.NextDouble();
                    if (v <= (2.0 - x) / (2.0 + x)) return x;
                    if (v <= Math.Exp(-x)) return x;
                } else {
                    double x = -Math.Log(z * (b - P) / a);
                    double y = x / z;
                    double v = rng.NextDouble();
                    if (v * (a + y - a * y) < 1.0) return x;
                    if (v <= Math.Pow(y, am1)) return x;
                }
            }
        }
    }

    // The Ahrens-Dieter generator is described by Best in his article above.
    // It is an acceptance-rejection generator for a < 1.

#if PAST
    internal class AhrensDieterLowAlphaGammaGenerator : IDeviateGenerator<double> {

        public AhrensDieterLowAlphaGammaGenerator(double alpha) {
            Debug.Assert(0.0 < alpha && alpha <= 1.0);
            a = alpha;
            b = (Math.E + a) / Math.E;
            ai = 1.0 / a;
            am1 = a - 1.0;
        }

        private readonly double a, b, ai, am1;

        public double GetNext(Random rng) {

            while (true) {

                double P = b * rng.NextDouble();

                if (P < 1.0) {
                    double x = Math.Pow(P, ai);
                    if (rng.NextDouble() > Math.Exp(-x)) continue;
                    return (x);
                } else {
                    double x = -Math.Log((b - P) / a);
                    if (rng.NextDouble() > Math.Pow(x, am1)) continue;
                    return (x);
                }

            }

        }

    }

#endif

    // Marsaglia, G and Tsang W W, "A Simple Method for GEnerating Gamma Variables",
    // ACM Transactions on Mathematical Software 26 (2000) 363-372

    internal class MarsagliaTsangGammaGenerator : IDeviateGenerator<double> {

        public MarsagliaTsangGammaGenerator(double alpha) {
            Debug.Assert(alpha >= 1.0);
            d = alpha - 1.0 / 3.0;
            c = 1.0 / Math.Sqrt(9.0 * d);
            zGenerator = DeviateGeneratorFactory.GetNormalGenerator();
        }

        double c, d;
        IDeviateGenerator<double> zGenerator;

        public double GetNext(Random rng) {

            while (true) {
                double x, v;
                do {
                    x = zGenerator.GetNext(rng);
                    v = 1.0 + c * x;
                } while (v <= 0.0);
                v = v * v * v;
                double u = rng.NextDouble();
                double x2 = x * x;
                double x4 = x2 * x2;
                if (u < 1.0 - 0.0331 * x4) return d * v;
                if (Math.Log(u) < 0.5 * x2 + d * (1.0 - v + Math.Log(v))) return d * v;
            }
        }

    }

#if PAST

    // Marsaglia and Tsang also describe a "boost" prodcedue to generate deviates for
    // a < 1 from a generator that works for a > 1, but we don't use this because
    // the dedicated a < 1 generators are faster.

    public class BoostedMarsagliaTsangGamma : IDeviateGenerator<double> {

        public BoostedMarsagliaTsangGamma(double alpha) {
            Debug.Assert(0.0 < alpha && alpha < 1.0);
            ai = 1.0 / alpha;
            boostedGammaGenerator = new MarsagliaTsangGammaGenerator(alpha + 1.0);
        }

        private double ai;
        private IDeviateGenerator<double> boostedGammaGenerator;

        public double GetNext(Random rng) {

            double g = boostedGammaGenerator.GetNext(rng);
            double u = rng.NextDouble();
            return g * Math.Pow(u, ai);

        }

    }

    internal class MarsagliaTsangGammaGeneratorOld : IDeviateGenerator<double> {

        public MarsagliaTsangGammaGenerator (IDeviateGenerator<double> normalGenerator, double alpha) {
            Debug.Assert(alpha >= 1.0);
            this.normalGenerator = normalGenerator;
            a = alpha;
            a1 = a - 1.0 / 3.0;
            a2 = 1.0 / Math.Sqrt(9.0 * a1);
        }

        private readonly IDeviateGenerator<double> normalGenerator;

        private readonly double a, a1, a2;

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
                double z2 = z * z;
                double z4 = z2 * z2;
                if (u > 1.0 - 0.331 * z4) {
                    // then, if necessary, via the exact rejection boundary
                    if (Math.Log(u) > z2 / 2.0 + a1 * (1.0 - v + Math.Log(v))) continue;
                }

                return (a1 * v);

            }

        }
    }

#endif

    internal class BetaFromGammaGenerator : IDeviateGenerator<double> {

        public BetaFromGammaGenerator (IDeviateGenerator<double> alphaGenerator, IDeviateGenerator<double> betaGenerator) {
            Debug.Assert(alphaGenerator != null);
            Debug.Assert(betaGenerator != null);
            this.alphaGenerator = alphaGenerator;
            this.betaGenerator = betaGenerator;
        }

        private readonly IDeviateGenerator<double> alphaGenerator, betaGenerator;

        public double GetNext (Random rng) {
            double x = alphaGenerator.GetNext(rng);
            double y = betaGenerator.GetNext(rng);
            return (x / (x + y));
        }
    }

}
