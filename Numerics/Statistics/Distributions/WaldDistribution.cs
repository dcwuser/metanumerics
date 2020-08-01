using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Wald (Inverse Gaussian) distribution.
    /// </summary>
    /// <remakrs>
    /// <para>The Wald distribution, also called the inverse Gaussian distribution, is the distribution of first
    /// passage times for a random walk.</para>
    /// <para>The Wald distribution is often called the inverse Gaussian distribution, but it is not the
    /// distribution of 1/x when x is Gaussian-distributed.</para>
    /// <para>If a system exhibits one-dimensional Brownian motion with mean displacement mu t and variance
    /// sigma t-squared.</para>
    /// <para>This can be phrased in terms of the Gambler's ruin problem: given an initial endowment x, a gambler
    /// repeatedly plays a game in which he wins 1 dollar with probability p and looses one dollar with probability
    /// q = 1 - p. If q > p, he will eventually loose all his endowment. What is the probability distribution that
    /// he will do so after exactly t games?</para>
    /// </remakrs>
    /// <seealso href="http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution"/>
    public sealed class WaldDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new Wald distribution.
        /// </summary>
        /// <param name="mean">The mean value, which must be positive.</param>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        public WaldDistribution (double mean, double shape) {
            if (mean <= 0.0) throw new ArgumentOutOfRangeException(nameof(mean));
            if (shape <= 0.0) throw new ArgumentOutOfRangeException(nameof(shape));
            this.mu = mean;
            this.lambda = shape;
        }

        private readonly double mu;
        private readonly double lambda;
        private readonly IDeviateGenerator<double> rngGenerator = new BoxMullerRejectionNormalDeviateGenerator();

        /// <summary>
        /// Gets the shape parameter of the distribution.
        /// </summary>
        public double Shape {
            get {
                return (lambda);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return(mu);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (mu * mu * mu / lambda);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (3.0 * Math.Sqrt(mu / lambda));
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                return (15.0 * mu / lambda);
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return(1.0);
            } else {
                double M1 = mu;
                if (r == 1) return(M1);
                double mu2 = mu * mu;
                double M2 = mu2 * (lambda + mu) / lambda;
                // use recursion M_{k+1} = (2k-1) M_{k} \mu^2 / \lambda + \mu^2 M_{k-1} 
                for (int k = 2; k < r; k++) {
                    double M3 = ((2 * k - 1) * M2 / lambda + M1) * mu2;
                    M1 = M2;
                    M2 = M3;
                }
                return (M2);
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
                double[] K = Cumulants(r);
                return (MomentMath.CumulantToCentral(K, r));
            }
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
                // Mathworld (http://mathworld.wolfram.com/InverseGaussianDistribution.html) says:
                //   K_{r+1} = \frac{(2 r)!}{2^r r! \mu^{2r+1} \lambda^r}
                //           = \frac{(2r - 1)!! \mu^{2r+1}}{\lambda^r}
                return (AdvancedIntegerMath.DoubleFactorial(2 * r - 3) * MoreMath.Pow(mu * mu / lambda, r - 1) * mu);
            }
        }

        internal double[] Cumulants (int rMax) {
            // This is just a recursive version of the cumulant formula.
            double[] K = new double[rMax + 1];
            K[0] = 0.0;
            K[1] = mu;
            double t = mu * mu / lambda;
            for (int r = 2; r <= rMax; r++) {
                K[r] = (2 * r - 3) * t * K[r - 1];
            }
            return (K);
        }


        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double z = (x - mu) / mu;
                return (Math.Sqrt(lambda / Global.TwoPI / (x * x * x)) * Math.Exp(-lambda * z * z / (2.0 * x)));
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double s = Math.Sqrt(lambda / x);
                double z1 = s * (x - mu) / mu;
                double z2 = s * (x + mu) / mu;
                return (NormalDistribution.Phi(z1) + Math.Exp(2.0 * lambda / mu) * NormalDistribution.Phi(-z2));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else {
                // Start from the left probability formula and use 1-\Phi(z) = \Phi(-z)
                // This formula accurately reproduces very small probabilities in the right tail.
                double s = Math.Sqrt(lambda / x);
                double z1 = s * (x - mu) / mu;
                double z2 = s * (x + mu) / mu;
                return (NormalDistribution.Phi(-z1) - Math.Exp(2.0 * lambda / mu) * NormalDistribution.Phi(-z2));
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <summary>
        /// Determines the parameters of the Wald distribution that best fits a sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The fit.</returns>
        /// <remarks>
        /// <para>The returned fit parameters are the <see cref="Mean"/> and <see cref="Shape"/>, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="WaldDistribution(double,double)"/> constructor to
        /// specify a new Wald distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static WaldFitResult FitToSample (Sample sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            return (Univariate.FitToWald(sample.data));
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {

            if (rng == null) throw new ArgumentNullException(nameof(rng));

            // This is a rather weird transformation generator described in Michael et al, "Generating Random Variates
            // Using Transformations with Multiple Roots", The American Statistician 30 (1976) 88-90.

            // u ~ U(0,1), v ~ ChiSquare(1), i.e. square of standard normal deviate

            double v = MoreMath.Sqr(rngGenerator.GetNext(rng));
            double w = mu * v;
            double x = mu * (1.0 + (w - Math.Sqrt((4.0 * lambda + w ) * w)) / (2.0 * lambda));
            double u = rng.NextDouble();
            if (u <= mu / (mu + x)) {
                return (x);
            } else {
                return (mu * mu / x);
            }
        }

    }

}