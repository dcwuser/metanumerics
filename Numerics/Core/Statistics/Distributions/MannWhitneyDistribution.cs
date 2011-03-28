using System;
using System.Diagnostics;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of the Mann-Whitney test statistic.
    /// </summary>
    internal sealed class MannWhitneyExactDistribution : DiscreteDistribution {

        /// <summary>
        /// Instantiates a new Mann-Whitney distribution.
        /// </summary>
        /// <param name="m">The number of elements in the first sample.</param>
        /// <param name="n">The number of elements in the second sample.</param>
        public MannWhitneyExactDistribution (int m, int n) {

            if (m < 1) throw new ArgumentOutOfRangeException("m");
            if (n < 1) throw new ArgumentOutOfRangeException("n");

            this.m = m;
            this.n = n;
            this.counts = GaussianBinomialCoefficients(m + n, n);

            total = 0.0;
            for (int i = 0; i < n * m; i++) {
                total += (double) counts[i];
            }

        }

        private int m, n;
        private double total;
        private decimal[] counts;

        // this routine is based on the recurrsion
        // [ m n ] = ( 1 - q^m ) / ( 1 - q^(m-n) ) [ m-1 n ]
        // and the starting point [ n n ] = 1

        // the coefficients are integers and get large quickly as m and n increase
        // we use decimal because it handles larger integers than long
        // we can't use double because the calculation requires delicate cancelations
        // among large intermediate values, thus necessicating exact integer arithmetic
        // look into using an arbitrary-sized integer structure in the future

        private static decimal[] GaussianBinomialCoefficients (int m, int n) {

            if (m < 0) throw new ArgumentOutOfRangeException("m");
            if (n < 0) throw new ArgumentOutOfRangeException("n");

            Debug.Assert(m >= n);

            // create  an array to hold our coefficients
            decimal[] c = new decimal[(m - n) * n + 1];

            // start with [n n] = 1 * q^0
            c[0] = 1;

            // keep track of current degree of our polynomial
            int d = 0;

            // create a scratch array for intermediate use
            // it needs to be larger than the previous array by (m-n) to hold intermediate polynomials
            decimal[] b = new decimal[c.Length + (m - n)];

            // interate from [n n] up to [m n]
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

#if PAST
        /// <inheritdoc />
        public override DiscreteInterval  Support {
	        get {
                return(DiscreteInterval.FromEndpoints(0, m * n));
            }
        }
#endif

        /// <inheritdoc />
        public override int Minimum {
            get {
                return (0);
            }
        }

        /// <inheritdoc />
        public override int Maximum {
            get {
                return (m * n);
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
       public override double LeftProbability (int k) {
            if (k < 0) {
                return (0.0);
            } else if (k >= m * n) {
                return (1.0);
            } else {
                double P = 0;
                for (int i = 0; i <= k; i++) {
                    P += (double) counts[i];
                }
                return (P / total);
            }

        }

       /// <inheritdoc />
       public override double RightProbability (int k) {
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
       public override double MomentAboutMean (int n) {
           if (n < 0) {
               throw new ArgumentOutOfRangeException("n");
           } else if (n == 0) {
               return (1.0);
           } else if (n % 2 != 0) {
               return (0.0);
           } else {
               return (base.MomentAboutMean(n));
           }
       }

    }

}
