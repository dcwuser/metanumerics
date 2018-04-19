using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Statistics.Distributions {

    // Kendall's tau is a non-parametric measure of correlation in a bivariate sample.

    // Let K be the number of discordant pairs, i.e. pairs i & j where x_i > x_j but y_i < y_j.
    // The minimum number of discordant pairs is none, i.e. 0; the maximum number is all pairs, i.e. n (n - 1) / 2.
    // Since the distribution is symmetric, the mean is n (n - 1) / 4.

    // One way to find the distribution of K is to enumerate all n! possible relative orderings,
    // e.g. for n = 3, the possible pair are 1 & 2, 2 & 3, 1 & 3. Suppose the x-values are ordered abc.
    // Counting the number of discordant (!) pairs for each possible ordering of y-values:
    //   abc => 1 & 2 = ab & ab, 2 & 3 = bc & bc, 1 & 3 = ac & ac => K = 0
    //   bac => 1 & 2 = ab & ba (!), 2 & 3 = bc & ac, 1 & 3 = ac & bc => K = 1
    //   acb => 1 & 2 = ab & ac, 2 & 3 = bc & cb (!), 1 & 3 = ac & ab => K = 1
    //   cba => 1 & 2 = ab & cb (!), 2 & 3 = bc & ba (!), 1 & 3 = ac & ca (!) => K = 3
    //   cab => 1 & 2 = ab & ca (!), 2 & 3 = bc & ab, 1 & 3 = ac & cb (!) => K = 2
    //   bca => 1 & 2 = ab & bc, 2 & 3 = bc & ca (!), 1 & 3 = ac & ba (!) => K = 2
    // So the distribution is:
    //   P(K=0)=1/6, P(K=1)=2/6, P(K=2)=2/6, P(K=3)=1/6

    // These distributions can also be produced via the generating function
    //   \frac{1}{n!} \prod_{j=1}^{n} \frac{x^j - 1}{x - 1} = \sum_{k=0}^{n(n-1)/2} P_k x^k
    // Note that (x-1) is always a factor of (x^j - 1), with quotient 1 + x + \cdots + x^{j-1},
    // i.e. the polynomial of degree j-1 with all unit coefficients.
    // So for n = 3, the relevant polynomial is:
    //   1 X (1 + x) X (1 + x + x^2) = 1 + 2 x + 2 x^2 + x^3
    // whose coefficients, as claimed, give the same distribution of K derived above.

    // In David, Kendall, & Stuart, "Some Questions of Distribution in the Theory of Rank Correlation",
    // Biometrika 38 (1951) 131, there are formulas for the cumulants of the distribution of K:
    //    K1 = M1 = n (n - 1) / 4
    //    K2 = C2 = n (n - 1) (2 n + 5) / 72
    //    K4 = - n * (6 n^4 + 15 n^3 + 10 n^2 - 31) / 3600
    //    K6 = n (6 n^6 + 21 n^5 + 21 n^4 - 7 n^2 - 41) / 10584
    // and of course odd cumulants above the first vanish due to the symmetry of the distribution.
    // The counts we produce agree with these formulas.

    // Linear mapping of K = 0 to \tau = 1 and K = n(n-1)/2 to \tau = -1 implies \tau = 1 - \frac{4 K}{n (n - 1)}

    internal sealed class KendallExactDistribution : DiscreteDistribution {

        public KendallExactDistribution (int n) {
            if (n < 2) throw new ArgumentOutOfRangeException(nameof(n));

            // since n! > Int64.MaxValue for n > 20, we can only go up to n = 20
            if (n > 20) throw new ArgumentOutOfRangeException(nameof(n));

            this.n = n;
            ComputeCounts();

        }

        private int n;
        private long[] counts;
        private long total;

        private void ComputeCounts () {

            // Use the generating function method to generate counts.
            // Multiply factors (1 + x + \cdots + x^{k-1}) for all k up to n to get the coefficients.
            // Multiply factors k for all k up to n to get n!, the total number of counts.

            counts = new long[] { 1 };
            total = 1;
            for (int k = 2; k <= n; k++) {
                counts = MultiplyByOnes(counts, k - 1);
                total *= k;
            }

        }

        // A simple utility method to multiply ( a_0 + a_1 x + ... + a_n x^ad ) X ( 1 + 1 x + ... + 1 x^bd )
        // We use the standard O(N^2) polynomial multiplication algorithm, but knowing that the coefficients of b
        // are all 1 allows us to avoid multiplications and do pure additions instead.

        private static long[] MultiplyByOnes (long[] a, int bDegree) {
            long[] ab = new long[a.Length + bDegree];
            for (int i = 0; i < a.Length; i++) {
                for (int j = 0; j <= bDegree; j++) {
                    ab[i + j] += a[i];
                }
            }
            return (ab);
        }

        public override double ProbabilityMass (int k) {
            if ((k < 0) || (k >= counts.Length)) {
                return (0.0);
            } else {
                return (((double) counts[k]) / total);
            }
        }

        public override double LeftExclusiveProbability (int k) {
            if (k < 1) {
                return (0.0);
            } else if (k >= counts.Length) {
                return (1.0);
            } else {
                return (((double) CountSum(0, k - 1)) / total);
            }
        }

        public override double RightExclusiveProbability (int k) {
            if (k < 0) {
                return (1.0);
            } else if (k >= counts.Length - 1) {
                return (0.0);
            } else {
                return (((double) CountSum(k + 1, counts.Length - 1)) / total);
            }
        }

        private long CountSum (int kmin, int kmax) {
            Debug.Assert(kmin >= 0);
            Debug.Assert(kmax < counts.Length);
            Debug.Assert(kmin <= kmax);
            long sum = 0;
            for (int k = kmin; k <= kmax; k++) sum += counts[k];
            return (sum);
        }

        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return (new DiscreteInterval(0, n * (n - 1) / 2));
            }
        }

        public override double Mean {
            get {
                return (n * (n - 1) / 4.0);
            }
        }

        public override double Variance {
            get {
                return (n * (n - 1) * (2 * n + 5) / 72.0);
            }
        }

        public override double Skewness {
            get {
                return (0.0);
            }
        }

        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentNullException(nameof(r));
            } else {
                double M = 0.0;
                for (int i = 0; i < counts.Length; i++) {
                    M += counts[i] * Math.Pow(i, r);
                }
                return (M / total);
            }
        }

        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                double M1 = Mean;
                double C = 0.0;
                for (int i = 0; i < counts.Length; i++) {
                    C += counts[i] * Math.Pow(i - M1, r);
                }
                return (C / total);
            }
        }

    }

}
