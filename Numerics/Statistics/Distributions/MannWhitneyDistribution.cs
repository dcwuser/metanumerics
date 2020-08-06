using System;
using System.Diagnostics;

using Meta.Numerics.Extended;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of the Mann-Whitney test statistic.
    /// </summary>
    internal sealed class MannWhitneyExactDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new Mann-Whitney distribution.
        /// </summary>
        /// <param name="m">The number of elements in the first sample.</param>
        /// <param name="n">The number of elements in the second sample.</param>
        public MannWhitneyExactDistribution (int m, int n) {

            if (m < 1) throw new ArgumentOutOfRangeException(nameof(m));
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n));

            this.m = m;
            this.n = n;
            this.counts = GaussianBinomialCoefficients(m + n, n);

            total = 0.0;
            for (int i = 0; i < n * m; i++) {
                total += (double) counts[i];
            }

        }

        private readonly int m, n;
        private readonly Int128[] counts;
        private readonly double total;

        // This routine is based on the recursion
        //   [ m n ] = ( 1 - q^m ) / ( 1 - q^(m-n) ) [ m-1 n ]
        // and the starting point [ n n ] = 1.

        // The coefficients are integers and get large quickly as m and n increase.
        // First we used decimal because it was the has the largest precise integer
        // values of any built-in, but now that we have Int128, we use that instead.

        private static Int128[] GaussianBinomialCoefficients (int m, int n) {

            Debug.Assert(m >= 0);
            Debug.Assert(n >= 0);
            Debug.Assert(m >= n);

            // create  an array to hold our coefficients
            Int128[] c = new Int128[(m - n) * n + 1];

            // start with [n n] = 1 * q^0
            c[0] = 1;

            // keep track of current degree of our polynomial
            int d = 0;

            // create a scratch array for intermediate use
            // it needs to be larger than the previous array by (m-n) to hold intermediate polynomials
            Int128[] b = new Int128[c.Length + (m - n)];

            // iterate from [n n] up to [m n]
            for (int k = n + 1; k <= m; k++) {

                // multiply by (1-q^k)
                for (int i = 0; i <= d; i++) {
                    b[i] = c[i];
                }
                d = d + k;
                for (int i = k; i <= d; i++) {
                    b[i] = b[i] - c[i - k];
                }

                // divide by (1-q^(k-n))
                for (int i = d - (k - n); i >= 0; i--) {
                    c[i] = -b[k - n + i];
                    b[k - n + i] = b[k - n + i] + c[i];
                    b[i] = b[i] - c[i];
                }
                d = d - (k - n);

            }

            // we're done
            return (c);

        }

        /// <inheritdoc />
        public override DiscreteInterval Support {
	        get {
                return(new DiscreteInterval(0, m * n));
            }
        }

        /// <inheritdoc />
        public override double ProbabilityMass(int k) {
            if ((k < 0) || (k > m * n)) {
                return(0.0);
            } else {
                return(((double) counts[k]) / total);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (m * n / 2.0);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (m * n * (m + n + 1) / 12.0);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                return (12.0 / 5.0 * (m + m * m + m * n + n * n + n) / (m * n * (m + n + 1)));
            }
        }

        /// <inheritdoc />
        public override double LeftExclusiveProbability (int k) {
            if (k < 0) {
                return (0.0);
            } else if (k >= m * n) {
                return (1.0);
            } else {
                double P = 0;
                for (int i = 0; i < k; i++) {
                    P += (double) counts[i];
                }
                return (P / total);
            }

        }

       /// <inheritdoc />
       public override double RightExclusiveProbability (int k) {
           if (k < 0) {
               return (1.0);
           } else if (k >= counts.Length) {
               return (0.0);
           } else {
               double P = 0;
               for (int i = k+1; i < counts.Length; i++) {
                   P += (double) counts[i];
               }
               return (P / total);
           }

       }

       /// <inheritdoc />
       public override double CentralMoment (int r) {
           if (r < 0) {
               throw new ArgumentOutOfRangeException(nameof(r));
           } else if (r == 0) {
               return (1.0);
           } else if (r % 2 != 0) {
               return (0.0);
           } else {
               return (base.CentralMoment(r));
           }
       }

    }


    // Froda & van Eeden, Canadian Journal of Statistics 28 (2000) quote a closed formula for cumulants
    // from a dutch manuscript van Dantzig, "Kader cursus Mathematische Statistiek" (1947-1950)
    //   K_{2\nu} = \frac{B_{2\nu}{2\nu} \sum_{j=1}^{m} \left[ (n + j)^{2v} - j^{2v} \right]
    // Using values for the Bernoulli numbers and formulas for sums of powers, this implies
    //   K_2 = \frac{1}{12} m n ( m + n + 1 )
    //   K_4 = -\frac{1}{60} m n ( m + n + 1 ) ( m + m^2 + m n + n^2 + n )
    // Testing against our distribution, these appear to be correct.
}
