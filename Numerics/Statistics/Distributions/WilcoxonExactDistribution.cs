using System;
using System.Diagnostics;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of the (transformed) Wilcoxon rank sum statistic.
    /// </summary>
    internal class WilcoxonDistribution : DiscreteDistribution {

        public WilcoxonDistribution (int n) {
            Debug.Assert(n > 1);
            Debug.Assert(n < 32);
            this.n = n;
            this.max = n * (n + 1) / 2;
            this.counts = ComputeCounts(n);
        }

        private readonly int n;
        private readonly int max;
        private readonly int[] counts;

        public override double ProbabilityMass (int k) {
            double norm = 1.0 / (1 << n);
            if (k <= max / 2) {
                if (k < 0) {
                    return (0.0);
                } else {
                    return (norm * counts[k]);
                }
            } else {
                if (k <= max) {
                    return (norm * counts[max - k]);
                } else {
                    return (0.0);
                }
            }
        }

        public override DiscreteInterval Support {
            get {
                return (new DiscreteInterval(0, max));
            }
        }

        public override double Mean {
            get {
                return (max / 2.0);
            }
        }

        public override double Variance {
            get {
                return (n * (n + 1.0) * (2.0 * n + 1.0) / 24.0);
            }
        }

        // At some point, I found a formula for the rth cumulant which had it
        // proportional to S_{r}(n) = \sum_{i=1}^{n} i^r. This appears to be
        // basically right, but with an n-independent factor I don't know. Also
        // I can no longer find the paper.

        private int LeftInclusiveSum (int k) {
            int sum = 0;
            for (int j = 0; j <= k; j++) {
                sum += counts[j];
            }
            return (sum);
        }

        public override double LeftInclusiveProbability (int k) {
            double norm = 1.0 / (1 << n);
            if (k <= max / 2) {
                if (k < 0) {
                    return (0.0);
                } else {
                    return (norm * LeftInclusiveSum(k));
                }
            } else {
                if (k <= max) {
                    return (1.0 - norm * LeftInclusiveSum(max - k));
                } else {
                    return (1.0);
                }
            }
        }

        private static int[] ComputeCounts (int n) {

            // The standard Wilcoxon W is the sum of the signed ranks 1 to n. The maximum is thus
            // n(n+1)/2 and the minimum is the negative of this number. This is a discrete
            // distribution, but it is problematic because every other integer has probability zero.
            // This screws up our smoothing techniques. So instead of summing 1 to n multiplied by
            // +1 and -1, we will sum 1 to n multiplied by +1 and 0. This produces exactly the same
            // counts, but ranging from 0 to n(n+1)/2 and with no gaps.

            // Under the null hypothesis, both signs are equally possible for all ranks. That is
            // 2^n possible sign combinations. For example, for n = 3,
            //    3 2 1
            //    0 0 0    W = 0 + 0 + 0 = 0
            //    0 0 1    W = 0 + 0 + 1 = 1
            //    0 1 0    W = 0 + 2 + 0 = 2
            //    0 1 1    W = 0 + 2 + 1 = 3
            //    1 0 0    W = 3 + 0 + 0 = 3
            //    1 0 1    W = 3 + 0 + 1 = 4
            //    1 1 0    W = 3 + 2 + 0 = 5
            //    1 1 1    W = 3 + 2 + 1 = 6
            // W has the promised range. The distribution is
            //    W  Count
            //    0  1
            //    1  1
            //    2  1
            //    3  2
            //    4  1
            //    5  1
            //    6  1
            // Note the distribution is symmetric around its midpoint. This is true generally.
            // Note the counts sum to 2^n. This is true generally.
            // Note in our example W increased monotonically with the sign pattern intrepreted
            // as a binary number. This is not true generally.

            // The example points to an obvious algorithm: Go through the numbers 0 to 2^n - 1.
            // Intepret the binary representation of each number as a sign pattern. For each sign
            // pattern, compute W. Accumulate the counts. This algorithm is correct and workable.
            // You can tweak it a bit by messing with Gray codes and such. But it is much
            // slower than the better algorithm we present now.
            
            // Note that the counts F obey the following recurrence:
            //   F_{n}(k) = F_{n-1}(k) - F_{n-1}(k - n)
            // where F for arguments outside the allowed range is taken to be 0. 
         
            Debug.Assert(n >= 2);

            int max = 3;
            int mid = max / 2;
            int[] counts = new int[] { 1, 1 };

            for (int m = 3; m <= n; m++) {
                int previousMax = max;
                int previousMid = mid;
                int[] previousCounts = counts;
                max = m * (m + 1) / 2;
                mid = max / 2;
                counts = new int[mid + 1];

                // F_{m}(k) = F_{m-1}(k) + F_{m-1}(k - m)

                for (int k = 0; k <= mid; k++) {
                    int k1 = k;
                    if (k1 > previousMid) k1 = previousMax - k1;
                    int count = previousCounts[k1];

                    int k2 = k - m;
                    if (k2 >= 0) {
                        if (k2 > previousMid) k2 = previousMax - k2;
                        count += previousCounts[k2];
                    }

                    counts[k] = count;
                }

            }

            return (counts);

        }

    }
}
