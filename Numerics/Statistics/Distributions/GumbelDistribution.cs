using System;
using System.Collections.Generic;

using Meta.Numerics.Matrices;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Gumbel distribution.
    /// </summary>
    /// <seealso href="http://en.wikipedia.org/wiki/Gumbel_distribution"/>
    public sealed class GumbelDistribution : ContinuousDistribution {

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity));
            }
        }

        /// <summary>
        /// Initializes a new standard Gumbel distribution.
        /// </summary>
        public GumbelDistribution () : this(0.0, 1.0) {
        }

        /// <summary>
        /// Initializes a new Gumbel distribution with the given parameters.
        /// </summary>
        /// <param name="location">The location parameter.</param>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="scale"/> is negative or zero.</exception>
        public GumbelDistribution (double location, double scale) {
            if (scale <= 0.0) throw new ArgumentOutOfRangeException(nameof(scale));
            this.m = location;
            this.s = scale;
        }

        private readonly double m, s;


        /// <summary>
        /// Gets the location parameter of the distribution.
        /// </summary>
        public double Location {
            get {
                return (m);
            }
        }

        /// <summary>
        /// Gets the scale parameter of the distribution.
        /// </summary>
        public double Scale {
            get {
                return (s);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            // check for infinite e; when e = PositiveInfinity, Exp(-e) = 0, but the computer doesn't recognize that the
            // latter is even smaller than the former is big, and gives NaN rather than zero when they are multiplied
            // for some reason e ~ Infinity rather than PositiveInfinity, so we check for that; i actually thought that
            // infinities were always either PositiveInfinity or NegativeInfinity, but that appears not to be the case
            if (Double.IsInfinity(e)) {
                return (0.0);
            } else {
                return (e * Math.Exp(-e) / s);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            return (Math.Exp(-e));
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            return (-MoreMath.ExpMinusOne(-e));
        }

        /// <inheritdoc />
        public override double Hazard (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            return (e / MoreMath.ExpMinusOne(e) / s);
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (m - s * Math.Log(-Math.Log(P)));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (m - s * Math.Log(-MoreMath.LogOnePlus(-Q)));
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (m - s * Math.Log(Global.LogTwo));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (m + s * AdvancedMath.EulerGamma);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (MoreMath.Sqr(Math.PI * s) / 6.0);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.PI / Math.Sqrt(6.0) * s);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (12.0 * Math.Sqrt(6.0) * AdvancedMath.Apery / MoreMath.Pow(Math.PI, 3));
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                return (12.0 / 5.0);
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                return (base.CentralMoment(r));
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                double[] M = RawMoments(r);
                return (M[r]);
            }
        }

        // Substituting u = e^{-z}, the raw moment integral for the standard Gumbel distribution can be re-written as 
        //   M_r = \int_{0}^{\infty} \! e^{-u} ( - \log u )^r \, du
        // Koelbig, "On the Integral \int_{0}^{\infty} \!e^{-\mu t} t^{\nu - 1} \log^m t \, dt", Mathematics of Computation 41 (1983) 171
        // derives formulas for integrals of this form. In particular, his equation (23) implies
        //   M_j = \sum_{k=0}{j-1} \frac{(j-1)!}{k!} \hat{\zeta}(j-k) M_k
        // where \hat{\zeta}(n) = \zeta(n) s^n for n > 1 and m + \gamma s for n = 1. This expresses a given raw moment in terms
        // of all the lower raw moments, making computation of the rth moment O(r^2).

        internal override double[] RawMoments (int rMax) {

            // Create an array to hold the moments.
            double[] M = new double[rMax + 1];
            M[0] = 1.0;
            if (rMax == 0) return (M);

            // Pre-compute the zeta-hat values, since we will use them repeatedly.
            double[] Z = new double[rMax + 1];
            double sj = s; // tracks s^j
            Z[1] = m + s * AdvancedMath.EulerGamma;
            for (int j = 2; j < Z.Length; j++) {
                sj *= s;
                Z[j] = AdvancedMath.RiemannZeta(j) * sj;
            }

            // Compute higher moments in turn.
            for (int j = 1; j <= rMax; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += Z[j - k] * M[k] / AdvancedIntegerMath.Factorial(k);
                }
                M[j] = AdvancedIntegerMath.Factorial(j - 1) * sum;
            }

            return (M);
            
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (0.0);
            } else if (r == 1) {
                return (Mean);
            } else {
                double C = AdvancedIntegerMath.Factorial(r - 1) * AdvancedMath.RiemannZeta(r) * MoreMath.Pow(s, r);
                return (C);
            }
        }

        /// <summary>
        /// Find the parameters of a Gumbel distribution that best fit the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The fit result.</returns>
        public static GumbelFitResult FitToSample (Sample sample) {
            return (Univariate.FitToGumbel(sample.data));
        }

    }

}
