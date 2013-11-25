using System;
using System.Collections.Generic;
using System.Text;

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

        private double p, q;

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
                return (q / (p*p));
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

    }

}
