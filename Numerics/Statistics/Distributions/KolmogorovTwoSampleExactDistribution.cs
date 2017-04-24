using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    internal class KolmogorovTwoSampleExactDistribution : DiscreteDistribution {

        public KolmogorovTwoSampleExactDistribution (int n, int m) {
            if (n < 1) throw new ArgumentOutOfRangeException("n");
            if (m < 1) throw new ArgumentOutOfRangeException("m");

            // Ensure n >= m
            if (m > n) Global.Swap(ref n, ref m);

            this.n = n;
            this.m = m;
        }

        private readonly int n, m;

        // We want to count all lattice paths that lead from 0,0 to i,j while satisfying some condition. For the 2-sample KS test, that condition
        // is | i/n - j/m | < c/mn. This count satisfies
        //   N(i,j) = N(i, j-1) + N(i-1, j)

        public long LatticePathSum (int c) {

            Func<int, int, int> indicator = (int i, int j) => {
                if (Math.Abs(m * i - n * j) <= c) {
                    return (1);
                } else {
                    return (0);
                }
            };

            long[] column = new long[m + 1];

            for (int j = 0; j <= m; j++) column[j] = indicator(0, j);

            for (int i = 1; i <= n; i++) {
                column[0] = indicator(i, 0);
                for (int j = 1; j <= m; j++) {
                    if (indicator(i, j) == 0) {
                        column[j] = 0;
                    } else {
                        column[j] = column[j] + column[j - 1];
                    }
                }
            }

            return (column[m]);

        }

        public override DiscreteInterval Support {
            get {
                // It may be possible to improve these bounds with a little more work.
                return (new DiscreteInterval(0, (int) AdvancedIntegerMath.LCM(n, m)));
            }
        }

        public override double LeftExclusiveProbability (int k) {
            return (LeftInclusiveProbability(k - 1));
        }

        public override double LeftInclusiveProbability (int k) {
            if (k < Support.LeftEndpoint) {
                return (0.0);
            } else if (k >= Support.RightEndpoint) {
                return (1.0);
            } else {
                return (LatticePathSum(k * ((int) AdvancedIntegerMath.GCF(n, m))) / AdvancedIntegerMath.BinomialCoefficient(n + m, n));
            }
        }

        public override double ProbabilityMass (int k) {
            if ((k < Support.LeftEndpoint) || (k > Support.RightEndpoint)) {
                return (0.0);
            } else {
                int c = k * ((int) AdvancedIntegerMath.GCF(n, m));
                int c1 = (k - 1) * ((int) AdvancedIntegerMath.GCF(n, m));
                return ((LatticePathSum(c) - LatticePathSum(c1)) / AdvancedIntegerMath.BinomialCoefficient(n + m, n));
            }
        }

        // had to do this because base call causes loop: Mean -> ExpectationValue -> Mean

        public override double Mean {
            get {
                double M1 = 0.0;
                for (int k = Support.LeftEndpoint; k <= Support.RightEndpoint; k++) {
                    M1 += ProbabilityMass(k) * k;
                }
                return (M1);
            }
        }

        public override double Variance {
            get {
                double M1 = Mean;
                double C2 = 0.0;
                for (int k = Support.LeftEndpoint; k <= Support.RightEndpoint; k++) {
                    C2 += ProbabilityMass(k) * MoreMath.Sqr(k - M1);
                }
                return (C2);
            }
        }

        /*
        public override double ExpectationValue (Func<int, double> f) {
            double s = 0.0;
            for (int k = 0; k <= Maximum; k++) {
                s += ProbabilityMass(k) * f(k);
            }
            return (s);
        }
        */
    }

}
