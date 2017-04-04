using System;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Wald distribution.
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
        private readonly IDeviateGenerator rngGenerator = new BoxMullerRejectionNormalDeviateGenerator();

        /// <summary>
        /// Gets the shape parameter of the distribution.
        /// </summary>
        public double ShapeParameter {
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
        /// <para>The returned fit parameters are the <see cref="Mean"/> and <see cref="ShapeParameter"/>, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="WaldDistribution(double,double)"/> constructor to
        /// specify a new Wald distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        public static FitResult FitToSample (Sample sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            // For maximum likelihood estimation, take logs of pdfs and sum:
            //    \ln p = \frac{1}{2} \ln \lambda - \frac{1}{2} \ln (2\pi) - \frac{3}{2} \ln x
            //            - \frac{\lambda x}{2\mu^2} + \frac{\lambda}{\mu} - \frac{\lambda}{2x}
            //    \ln L = \sum_i p_i
            
            // Take derivative wrt \mu
            //    \frac{\partial \ln L}{\partial \mu} = \sum_i \left[ \frac{\lambda x_i}{\mu^3} - \frac{\lambda}{\mu^2} \right]
            // and set equal to zero to obtain
            //    \mu = \frac{1}{n} \sum_i x_i = <x>
            // which agrees with method of moments.

            // Take derivative wrt \lambda
            //    \frac{\partial \ln L}{\partial \lambda} = \sum_i \left[ \frac{1}{2 \lambda} -\frac{x_i}{2\mu^2} + \frac{1}{\mu} - \frac{1}{2 x_i} \right]
            // Set equal to zero, plug in our expression for \mu, and solve for \lambda to get
            //    \frac{n}{\lambda} = \sum_i \left( \frac{1}{x_i} - \frac{1}{\mu} \right)
            //  i.e. \lambda^{-1} = <(x^{-1} - \mu^{-1})>

            double mu = sample.Mean;
            double mui = 1.0 / mu;
            double lambda = 0.0;
            foreach (double value in sample) {
                if (value <= 0.0) throw new InvalidOperationException();
                lambda += (1.0 / value - mui);
            }
            lambda = (sample.Count - 3) / lambda;

            // If x ~ IG(\mu, \lambda), then \sum_i x_i ~ IG(n \mu, n^2 \lambda), so \hat{\mu} ~ IG (\mu, n \lambda). This gives us
            // not just the exact mean and variance of \hat{\mu}, but its entire distribution. Since its mean is \mu, \hat{\mu} is an
            // unbiased estimator. And by the variance formula for IG, the variance of \hat{\mu} is \frac{\mu^3}{n \lambda}.

            // Tweedie, "Statistical Properties of Inverse Gaussian Distributions" (http://projecteuclid.org/download/pdf_1/euclid.aoms/1177706964)
            // showed that \frac{n \lambda}{\hat{\lambda}} ~ \chi^2_{n-1}. Since the mean of \chi^2_{k} is k, the MLE estimator of
            // \frac{1}{\lambda} can be made unbiased by replacing n by (n-1). However, we are estimating \lambda, not \frac{1}{\lambda}.
            // By the relation between chi squared and inverse chi squared distributions, \frac{\hat{\lambda}}{n \lambda} ~ I\chi^2_{n-1}.
            // The mean of I\chi^2_{k} is \frac{1}{n-2}, so to get an unbiased estimator of \lambda, we need to replace n by (n-3). This is
            // what we have done above. Furthermore, the variance of I\chi^2_{k} is \frac{2}{(k-2)^2 (k-4)}, so the variance of \hat{\lambda}
            // is \frac{2 \lambda^2}{(n-5)}.

            // We can also get covariances from the MLE approach. To get a curvature matrix, take additional derivatives
            //   \frac{\partial^2 \ln L}{\partial \mu^2} = \sum_i \left[ -\frac{3 \lambda x_i}{\mu^4} + \frac{2 \lambda}{\mu^3} \right]
            //   \frac{\partial^2 \ln L}{\partial \mu \partial \lambda} = \sum_i \left[ \frac{x_i}{\mu^3} - \frac{1}{\mu^2} \right]
            //   \frac{\partial^2 \ln L}{\partial \lambda^2} =\sum_i \left[ - \frac{1}{2 \lambda^2} \right]
            // and substitutue in best-fit values of \mu and \lambda
            //   \frac{\partial^2 \ln L}{\partial \mu^2} = - \frac{n \lambda}{\mu^3}
            //   \frac{\partial^2 \ln L}{\partial \mu \partial \lambda} = 0
            //   \frac{\partial^2 \ln L}{\partial \lambda^2} = - \frac{n}{2 \lambda^2}
            // Mixed derivative vanishes, so matrix is trivially invertible to obtain covariances. These results agree with the
            // results from exact distributions in the asymptotic regime.

            double v_mu_mu = mu * mu * mu / lambda / sample.Count;
            double v_lambda_lambda = 2.0 * lambda * lambda / (sample.Count - 5);
            double v_mu_lambda = 0.0;

            ContinuousDistribution dist = new WaldDistribution(mu, lambda);
            TestResult test = sample.KolmogorovSmirnovTest(dist);
            return (new FitResult(mu, Math.Sqrt(v_mu_mu), lambda, Math.Sqrt(v_lambda_lambda), v_mu_lambda, test));
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {

            // This is a rather weird transformation generator described in Michael et al, "Generating Random Variates
            // Using Transformations with Multiple Roots", The American Statistician 30 (1976) 88-90.

            double u = rngGenerator.GetNext(rng);
            double y = MoreMath.Sqr(rngGenerator.GetNext(rng));
            double muy = mu * y;
            double x = mu * (1.0 + (muy - Math.Sqrt((4.0 * lambda + muy ) * muy)) / (2.0 * lambda));
            double z = rng.NextDouble();
            if (z <= mu / (mu + x)) {
                return (x);
            } else {
                return (mu * mu / x);
            }
        }

    }

    /*
    public sealed class WaldFitResult : FitResult {

        public UncertainValue Mean {
            get {
                return (this.Parameter(0));
            }
        }

        public UncertainValue Shape {
            get {
                return (this.Parameter(1));
            }
        }

        public WaldDistribution Distribution {
            get {
                return (new WaldDistribution(this.Parameter(0).Value, this.Parameter(1).Value));
            }
        }

    }
    */
}