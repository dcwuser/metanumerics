using System;
using System.Diagnostics;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a geometric distribution.
    /// </summary>
    /// <remarks>
    /// <para>The probability of obtaining each integer value in a geometric distribution is lower than the probability of obtaining
    /// the previous value by a constant factor. The probability thus decreases geometricly with the value, giving the distribution its name.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Geometric_distribution"/>
    /// <seealso href="https://mathworld.wolfram.com/GeometricDistribution.html"/>
    public sealed class GeometricDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new geometric distribution.
        /// </summary>
        /// <param name="p">The probabability of the value zero.</param>
        public GeometricDistribution (double p) {
            if ((p <= 0.0) || (p >= 1.0)) throw new ArgumentOutOfRangeException(nameof(p));
            this.p = p;
            this.q = 1.0 - p;
            this.lnq = MoreMath.LogOnePlus(-p);
        }

        private readonly double p, q;

        private readonly double lnq;


        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return DiscreteInterval.Semiinfinite;
            }
        }

        /// <inheritdoc />
        public override double ProbabilityMass (int k) {
            if (k < 0) {
                return 0.0;
            } else {
                return p * MoreMath.Pow(q, k);
            }
        }

        /// <inheritdoc />
        public override double LeftExclusiveProbability (int k) {
            if (k <= 0) {
                return 0.0;
            } else {
                return 1.0 - MoreMath.Pow(q, k);
            }
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
            if (k < 0) {
                return 1.0;
            } else {
                return MoreMath.Pow(q, k + 1);
            }
        }

        /// <inheritdoc />
        public override int InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return ((int) Math.Min(Math.Floor(MoreMath.LogOnePlus(-P) / lnq), Int32.MaxValue));
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return q / p;
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return q / (p * p);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return Math.Sqrt(q) / p;
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (2.0 - p) / Math.Sqrt(q);
            }
        }

        // Raw moments M_r = \frac{q E_r(q)}{p^r} where E_n(q) is the nth Eulerian polynomial, i.e.
        // the polynomial with coefficient < n m > of x^m. (https://en.wikipedia.org/wiki/Eulerian_number)

        // Note that the Eulerian polynomials and their coefficients the Eulerian numbers are not
        // the same as the Euler polynomials and Euler numbers.

        // The (falling) factorial moments of the geometric distribution are also simply expressible
        // in closed form as F_r \mu^r r!, so we could also use the trick of using the Stirling
        // numbers to convert these to raw moments. But since we are doing the Eulerian numbers
        // anyway for cumulants, we will use them for raw moments, too.

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else {
                return q * EulerianPolynomial(r, q) / MoreMath.Pow(p, r);
            }
        }

        // Cumulants K_r = \frac{q E_{r-1}(q)}{p^r} (also Eulerian polynomials but not quite same
        // as raw moments). See Shenton & Bowman, "The Geometric Distribution's Central Moments
        // and Eulerian Numbers of the Second Kind", Far East J. Theo. Stat. 7 (1) 2002 1-17
        // (http://www.csm.ornl.gov/~bowman/fjts7.pdf)

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 0.0;
            } else if (r == 1) {
                return q / p;
            } else {
                return q * EulerianPolynomial(r - 1, q) / MoreMath.Pow(p, r);
            }
        }

        // Shenton & Bowman also relate central moments to Eulerian numbers of the second kind.
        // We should implement this, too.

        // Raw and central moments can also be related to the negative-order polylog functions
        // and the Lerch transcendent, but we don't compute those yet.

        private static double EulerianPolynomial (int n, double x) {

            Debug.Assert(n > 0);

            // Determine the Eulerian numbers <n m> using
            //   < n m > = (n - m) < (n-1) (m-1) > + (m + 1) < (n-1) m >
            double[] c = new double[n];
            c[0] = 1.0;           
            for (int j = 2; j <= n; j++) {
                c[j - 1] = 1.0;
                for (int k = j - 2; k > 0; k--) {
                    c[k] = (k + 1) * c[k] + (j - k) * c[k - 1];
                }
            }

            // Evaluate the polynomial \sum_{m} <n m> x^m.
            double f = 1.0;
            double xj = 1.0; // tracks x^j
            for (int j = 1; j < n; j++) {
                xj *= x;
                f += c[j] * xj;
            }
            return (f);

        }

        /// <inheritdoc />
        public override int GetRandomValue (Random rng) {
            if (rng is null) throw new ArgumentNullException(nameof(rng));
            double u = rng.NextDouble();
            return ((int) Math.Min(Math.Floor(Math.Log(u) / lnq), Int32.MaxValue));
        }

    }

}
