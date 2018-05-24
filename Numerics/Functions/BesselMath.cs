using System;
using System.Collections.Generic;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Analysis;

namespace Meta.Numerics.Functions
{
    internal static class BesselMath
    {

        public static double AiryAiZero (int k) {
            if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

            double tm = 3.0 / 8.0 * Math.PI * (4 * k - 2);
            //double t0 = 3.0 / 8.0 * Math.PI * (4 * k - 1);
            double tp = 3.0 / 8.0 * Math.PI * (4 * k - 0);

            double am = -T(tm);
            //double a0 = -T(t0);
            double ap = -T(tp);

            double a1 = FunctionMath.FindZero(AdvancedMath.AiryAi, Interval.FromEndpoints(ap, am));

            return (a1);
        }

        public static double AiryBiZero (int k) {
            if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

            double tm = 3.0 / 8.0 * Math.PI * (4 * k - 4);
            double tp = 3.0 / 8.0 * Math.PI * (4 * k - 2);

            double bm = (tm == 0.0) ? 0.0 : -T(tm);
            double bp = -T(tp);

            double b1 = FunctionMath.FindZero(AdvancedMath.AiryBi, Interval.FromEndpoints(bp, bm));

            return (b1);
        }

        private static double T (double t) {
            double u = 1.0 / (t * t);
            return ( Math.Pow(t, 2.0 / 3.0) * (
                1.0 + u * (5.0 / 48.0 + u * (-5.0/36.0 + 77125.0 / 82944.0 * u))
            ));
        }

    }
}
