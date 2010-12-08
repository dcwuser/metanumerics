using System;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Wald distribution.
    /// </summary>
    /// <remakrs>
    /// <para>The Wald distribution, also called the inverse Gaussian distribution, is the distribution of first
    /// passage times for a random walk.</para>
    /// <para>This can be phrased in terms of the Gambler's ruin problem: given an initial endowment x, a gambler
    /// repeatedly plays a game in which he wins 1 dollar with probability p and looses one dollar with probability
    /// q = 1 - p. If q > p, he will eventually loose all his endowment. What is the probability distribution that
    /// he will do so after exactly t games?</para>
    /// </remakrs>
    /// <seealso href="http://en.wikipedia.org/wiki/Inverse_Gaussian_distribution"/>
    public class WaldDistribution : Distribution {

        /// <summary>
        /// Instantiates a new Wald distribution.
        /// </summary>
        /// <param name="mean">The mean value, which must be positive.</param>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        public WaldDistribution (double mean, double shape) {
            if (mean <= 0.0) throw new ArgumentOutOfRangeException("mean");
            if (shape <= 0.0) throw new ArgumentOutOfRangeException("shape");
            this.mu = mean;
            this.lambda = shape;
        }

        private double mu;
        private double lambda;

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
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return(1.0);
            } else {
                double M1 = mu;
                if (n == 1) return(M1);
                double mu2 = mu * mu;
                double M2 = mu2 * (lambda + mu) / lambda;
                // use recursion M_{k+1} = (2k-1) M_{k} \mu^2 / \lambda + \mu^2 M_{k-1} 
                for (int k = 2; k < n; k++) {
                    double M3 = ((2 * k - 1) * M2 / lambda + M1) * mu2;
                    M1 = M2;
                    M2 = M3;
                }
                return (M2);
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
                return(CentralMomentFromRawMoment(n));
            }
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
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

    }


}