using System;
using System.Collections.Generic;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Extended;

namespace Meta.Numerics.Functions
{
    internal static class BesselMath
    {

        public static double AiryAiZero (int k) {

            if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

            // The very 1st root requires 3 Halley cycles, so special-case it.
            if (k == 1) return -2.33810741045976704;


            // DLMF 9.9 gives an asymptotic expansion for Airy roots.
            double t = 3.0 / 8.0 * Math.PI * (4 * k - 1);
            double x = -T(t);

            // We have f=Ai(x), f'= Ai'(x), f'' = Ai''(x) = x Ai(x).
            // Halley's method triples digits of accuracy with each iteration.
            // Since we know it works, hard-code a small, sufficient number of
            // iterations to avoid messsing around with termination criteria.

            // Passing into a root-finder, even one that uses derivatives, requires 5-30
            // function evaluations, so this is faster, as I verified via direct comparison.

            for (int i = 0; i < 2; i++) {
                SolutionPair a = AdvancedMath.Airy(x);
                double r = a.FirstSolutionValue / a.FirstSolutionDerivative;
                x -= r / (1.0 - 0.5 * r * x);
            }

            return x;

        }

        public static double AiryBiZero (int k) {

            if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

            // The 1st root requires 5 Hally iterations, so special-case it.
            if (k == 1) return -1.17371322270912792;

            double t = 3.0 / 8.0 * Math.PI * (4 * k - 3);
            double x = -T(t);

            for (int i = 0; i < 2; i++) {
                SolutionPair b = AdvancedMath.Airy(x);
                double r = b.SecondSolutionValue / b.SecondSolutionDerivative;
                x -= r / (1.0 - 0.5 * r * x);
            }

            return x;

        }

        private static double T (double t) {
            // Even for t_1 ~ 3.53, terms decrease up to 1/t^6, yielding a relative error of ~2.2E-4
            // for higher roots, even more accurate.
            double u = 1.0 / (t * t);
            return Math.Pow(t, 2.0 / 3.0) * (
                1.0 + u * (5.0 / 48.0 + u * (-5.0 / 36.0 + u * (77125.0 / 82944.0)))
            );
        }

    }
}
