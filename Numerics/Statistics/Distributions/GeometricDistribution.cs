using System;
using System.Collections.Generic;
using System.Text;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a geometric distribution.
    /// </summary>
    /// <remarks>
    /// <para>The probability of obtaining each integer value in a geometric distribution is lower than the probability of obtaining
    /// the previous value by a constant factor. The probability thus decreases geometricly with the value, giving the distribution its name.</para>
    /// </remarks>
    public sealed class GeometricDistribution : DiscreteDistribution {

        /// <summary>
        /// Initializes a new geometric distribution.
        /// </summary>
        /// <param name="p">The probabability of the value zero.</param>
        public GeometricDistribution (double p) {
            if ((p <= 0.0) || (p >= 1.0)) throw new ArgumentOutOfRangeException("p");
            this.p = p;
            this.q = 1.0 - p;
        }

        private readonly double p, q;

#if PAST
        /// <inheritdoc />
        public override DiscreteInterval Support {
            get {
                return (DiscreteInterval.FromEndpoints(0, Int32.MaxValue));
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
                return (Int32.MaxValue);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityMass (int k) {
            if (k < 0) {
                return (0.0);
            } else {
                return (p * MoreMath.Pow(q, k));
            }
        }

        /// <inheritdoc />
        public override double LeftExclusiveProbability (int k) {
            if (k <= 0) {
                return (0.0);
            } else {
                return (1.0 - MoreMath.Pow(q, k));
            }
        }

        /// <inheritdoc />
        public override double RightExclusiveProbability (int k) {
            if (k < 0) {
                return (1.0);
            } else {
                return (MoreMath.Pow(q, k + 1));
            }
        }

        /// <inheritdoc />
        public override int InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            if (P < p) {
                return (0);
            } else {
                double Q = 1.0 - P;
                if (Q == 0.0) {
                    return (Int32.MaxValue);
                } else {
                    double kd = Math.Log(Q) / Math.Log(q) - 1.0;
                    if (kd > Int32.MaxValue) {
                        return (Int32.MaxValue);
                    } else {
                        return ((int)Math.Ceiling(kd));
                    }
                }
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (q / p);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (q / (p * p));
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.Sqrt(q) / p);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return ((2.0 - p) / Math.Sqrt(q));
            }
        }

        // raw moments involve the negative polylog function, central moments the Lerch transcendent

        // Raw moments M_r = E_n(q) / p^r where E_n(q) is the nth Eulerian polynomial, i.e. the polynomial with coefficient < n m > of x^m.
        // Note that the Eulerian polynomials and their coefficients the Eulerian numbers are not the same as the Euler polynomials and Euler numbers.

        /// <inheritdoc />
        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException("r");
            } else if (r == 0) {
                return (1.0);
            } else {
                return (q * EulerianPolynomial(r, q) / MoreMath.Pow(p, r));
            }
        }

        public static double EulerianPolynomial (int n, double x) {

            // Determine coefficients using < n m > = (n - m) < (n-1) (m-1) > + (m + 1) < (n-1) m >
            double[] c0 = new double[n];
            double[] c1 = new double[n];
            c0[0] = 1.0;
            for (int i = 2; i <= n; i++) {
                c1[0] = 1.0;
                for (int j = 1; j < i; j++) {
                    c1[j] = (i - j) * c0[j - 1] + (j + 1) * c0[j];
                }
                double[] t = c0; c0 = c1; c1 = t;
            }

            // Evaluate the polynomial \sum_{m = 1} <n m> x^m.
            double f = 1.0;
            double xj = 1.0; // tracks x^j
            for (int j = 1; j < n; j++) {
                xj *= x;
                f += c0[j] * xj;
            }
            return (f);

        }


    }

}
