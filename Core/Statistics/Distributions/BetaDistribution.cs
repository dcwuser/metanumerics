using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a &#x3B2; distribution.
    /// </summary>
    /// <remarks>
    /// <para>The &#x3B2; distribution is defined on the interval [0,1]. Depending on its two shape parameters, it can take on a wide
    /// variety of forms on this interval.</para>
    /// <para>The &#x3B2; distribution appears in a few statistical applications. Because of its versatility, it is also sometimes
    /// used as an ad hoc model for a distribution on a finite interval.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Beta_distribution"/>
    public class BetaDistribution : Distribution {

        /// <summary>
        /// Instantiates a new &#x3B2; distribution.
        /// </summary>
        /// <param name="alpha">The left shape parameter, which controls the form of the distribution near x=0.</param>
        /// <param name="beta">The right shape parameter, which controls the form of the distribution near x=1.</param>
        /// <remarks>
        /// <para>The <paramref name="alpha"/> shape parameter controls the form of the distribution near x=0. The
        /// <paramref name="beta"/> shape parameter controls the form of the distribution near z=1. If a shape parameter
        /// is less than one, the PDF diverges on the side of the distribution it controls. If a shape parameter
        /// is greater than one, the PDF goes to zero on the side of the distribution it controls. If the left and right
        /// shapre parameters are equal, the distribution is symmetric about x=1/2.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Beta_distribution" />
        public BetaDistribution (double alpha, double beta) {
            if (alpha <= 0.0) throw new ArgumentOutOfRangeException("alpha");
            if (beta <= 0.0) throw new ArgumentOutOfRangeException("beta");
            this.alpha = alpha;
            this.beta = beta;
        }

        double alpha, beta;

        /// <summary>
        /// Gets the left shape parameter.
        /// </summary>
        public double Alpha {
            get {
                return (alpha);
            }
        }

        /// <summary>
        /// Gets the right shape parameter.
        /// </summary>
        public double Beta {
            get {
                return (beta);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, 1.0));
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if ((x < 0.0) || (x > 1.0)) {
                return (0.0);
            } else {
                return (Math.Pow(x, alpha - 1.0) * Math.Pow(1.0 - x, beta - 1.0) / AdvancedMath.Beta(alpha, beta));
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x >= 1.0) {
                return (1.0);
            } else {
                return (AdvancedMath.Beta(alpha, beta, x) / AdvancedMath.Beta(alpha, beta));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (alpha / (alpha + beta));
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                double ab = alpha + beta;
                return (alpha * beta / (ab + 1.0) / (ab * ab));
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                double ab = alpha + beta;
                return (2.0 * (beta - alpha) / (ab + 2.0) * Math.Sqrt((ab + 1.0) / (alpha * beta)));
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else {
                return (AdvancedMath.Beta(alpha + n, beta) / AdvancedMath.Beta(alpha, beta));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else {
                return (CentralMomentFromRawMoment(n));
            }
        }

        // IParameterizedDistribution implementation
        // can't make this IParameterized yet, because fitting routine doesn't respect limits
        /*

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { alpha, beta });
        }

        void IParameterizedDistribution.SetParameters(IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 2) throw new InvalidOperationException();
            //if ((parameters[0] <= 0.0) || (parameters[1] <= 0.0)) throw new InvalidOperationException();
            alpha = parameters[0];
            beta = parameters[1];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }
        */

    }

}