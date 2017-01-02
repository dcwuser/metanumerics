using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

#if FUTURE

    public sealed class NoncentralChiSquaredDistribution : Distribution {

        public NoncentralChiSquaredDistribution (double nu, double lambda) {
            if (nu <= 0.0) throw new ArgumentOutOfRangeException("nu");
            if (lambda < 0.0) throw new ArgumentOutOfRangeException("lambda");
            this.nu = nu;
            this.lambda = lambda;
        }

        private readonly double lambda;
        private readonly double nu;

        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        public override double Mean {
            get {
                return (nu + lambda);
            }
        }

        public override double Variance {
            get {
                return (2.0 * (nu + 2.0 * lambda));
            }
        }

        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (0.0);
            } else {
                return (MoreMath.Pow(2.0, r - 1) * AdvancedIntegerMath.Factorial(r - 1) * (nu + r * lambda));
            }
        }

        public override double ProbabilityDensity (double x) {
            if (x < 0.0) {
                return (0.0);
            } else if (x == 0.0) {
                if (nu < 1.0) {
                    return (Double.PositiveInfinity);
                } else if (nu == 1.0) {
                    return (1.0);
                } else {
                    return (0.0);
                }
            } else {
                return (ProbabilityDensity_Series(x));
            }
        }

        private double ProbabilityDensity_Series (double x) {

            double z = lambda * x / 4.0;
            double halfNu = nu / 2.0;

            double t = 1.0 / AdvancedMath.Gamma(halfNu);
            double s = t;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                t *= z / k / (k - 1 + halfNu);
                s += t;
                if (s == s_old) {
                    return (s);
                }
            }
            throw new NonconvergenceException();
        }
    }

#endif

}
