using System;

using Meta.Numerics.Functions;


namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a non-central chi squared distribution.
    /// </summary>
    /// <seealso href="https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution"/>
    public sealed class NoncentralChiSquaredDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new instance of the non-central chi squared distribution.
        /// </summary>
        /// <param name="nu">The number of degrees of freedom, which must be positive.</param>
        /// <param name="lambda">The non-centrality parameter, which must be greater than zero.</param>
        public NoncentralChiSquaredDistribution (int nu, double lambda) {
            if (nu <= 0) throw new ArgumentOutOfRangeException(nameof(nu));
            if (lambda <= 0.0) throw new ArgumentOutOfRangeException(nameof(lambda));
            this.nu = nu;
            this.lambda = lambda;
        }

        private readonly int nu;
        private readonly double lambda;
        private readonly IDeviateGenerator<double> zRng = new BoxMullerRejectionNormalDeviateGenerator();

        /// <summary>
        /// Gets the number of degrees of freedom of the distribution.
        /// </summary>
        public int DegreesOfFreedom {
            get {
                return nu;
            }
        }

        /// <summary>
        /// Gets the non-centrality parameter of the distribution.
        /// </summary>
        public double Noncentrality {
            get {
                return lambda;
            }
        }

        /// <inheritdoc/>

        public override Interval Support {
            get {
                return Interval.Semiinfinite;
            }
        }

        /// <inheritdoc/>

        public override double ProbabilityDensity (double x) {

            if (x < 0.0) {
                return 0.0;
            } else {
                double sqrt_x = Math.Sqrt(x);
                double sqrt_lambda = Math.Sqrt(lambda);
                return 0.5 * Math.Pow(x / lambda, 0.25 * nu - 0.5) * Math.Exp(-0.5 * MoreMath.Sqr(sqrt_x - sqrt_lambda)) * AdvancedMath.ScaledModifiedBesselI(0.5 * nu - 1.0, sqrt_x * sqrt_lambda);
                /*
                double I = AdvancedMath.ModifiedBesselI(nu / 2.0 - 1.0, Math.Sqrt(lambda * x));
                double E = Math.Exp(-(lambda + x) / 2.0);
                double P = Math.Pow(x / lambda, nu / 4.0 - 0.5);
                return (I * E * P / 2.0);
                */
                // Consider implementing asymptotic series and/or power series explicitly to reduce overhead
                // and noise from sqrt, exponent, and power multiplies.
            }
        }

        // 0F1(k/2, lambda x/4)

        /*
        private double Series (double x) {
            double t = lambda * x / 2.0;

            double ds = Math.Exp(-(lambda + x) / 2.0) * Math.Pow(x / 2.0, nu / 2.0 - 1.0) / AdvancedMath.Gamma(nu / 2.0) / 2.0;
            double s = ds;
            for (int j = 1; j < 100; j++) {
                double s_old = s;
                ds *= t / j / (nu + 2 * (j - 1));
                s += ds;
                if (s == s_old) {
                    return (s);
                }
            }
            throw new NonconvergenceException();
        }
        */

        /// <inheritdoc/>
        public override double LeftProbability (double x) {
            if (x < 0.0) {
                return 0.0;
            } else {
                AdvancedMath.Marcum(0.5 * nu, Math.Sqrt(lambda), Math.Sqrt(x), out double P, out _);
                return P;
                //CumulativeDistributionFunction(x, out double P, out double Q);
                //return (P);
            }
        }

        /// <inheritdoc/>
        public override double RightProbability (double x) {
            if (x < 0.0) {
                return 1.0;
            } else {
                AdvancedMath.Marcum(0.5 * nu, Math.Sqrt(lambda), Math.Sqrt(x), out _, out double Q);
                return Q;
                //CumulativeDistributionFunction(x, out double P, out double Q);
                //return (Q);
            }
        }

        private void CumulativeDistributionFunction(double x, out double P, out double Q) {

            double nu2 = nu / 2.0;
            double x2 = x / 2.0;
            double lambda2 = lambda / 2.0;

            if (x < Mean) {

                double t = 1.0;
                double s = AdvancedMath.LeftRegularizedGamma(nu2, x2);
                for (int k = 1; k < Global.SeriesMax; k++) {
                    double s_old = s;
                    t *= lambda2 / k;
                    s += t * AdvancedMath.LeftRegularizedGamma(nu2 + k, x2);
                    if (s == s_old) {
                        double e = Math.Exp(-lambda2);
                        P = e * s;
                        Q = 1.0 - P;
                        return;
                    }
                }
                throw new NonconvergenceException();

            } else {
                
                // The Q-case is exactly analogous, but the recurrence for the right regularized
                // Gamma function is upward stable, so we can use it to cut down on incomplete
                // Gamma function evaluations.

                double t = 1.0;
                double q = AdvancedMath.RightRegularizedGamma(nu2, x2);
                double s = q;
                double dq = AdvancedMath.PowerOverFactorial(x2, nu2) * Math.Exp(-x2);
                for (int k = 1; k < Global.SeriesMax; k++) {
                    double s_old = s;
                    t *= lambda2 / k;
                    q += dq;
                    s += t * q;
                    if (s == s_old) {
                        double e = Math.Exp(-lambda2);
                        Q = e * s;
                        P = 1.0 - Q;
                        return;
                    }
                    dq *= x2 / (nu2 + k);
                }
                throw new NonconvergenceException();

            }

        }


        /// <inheritdoc/>
        public override double Mean {
            get {
                return nu + lambda;
            }
        }

        /// <inheritdoc/>
        public override double Variance {
            get {
                return 2.0 * (nu + 2.0 * lambda);
            }
        }

        /// <inheritdoc/>
        public override double Skewness {
            get {
                return Math.Pow(2.0 / (nu + 2.0 * lambda), 3.0 / 2.0) * (nu + 3.0 * lambda);
            }
        }


        /// <inheritdoc/>
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else {

                // Integrating the pdf gives: 
                //   M_r = < x^r > = 2^r Hypergeometric1F1(-r, \nu / 2, -\lambda / 2) ( \nu / 2 )_r
                // The confluent hypergeometric recurrence then implies:
                //   M_{r+1| = \left[ ( \nu + \lambda ) 4 r \right] M_{r} - 2 n \left[ \nu + 2(n - 1) \right] M_{r-1}
                // This appears to be stable for the \nu, \lambda regions I have tried.

                double mu = nu + lambda;

                double M0 = 1.0;
                double M1 = mu;
                for (int k = 1; k < r; k++) {
                    double M2 = (mu + 4 * k) * M1 - (2 * k) * (nu + 2 * (k - 1)) * M0;
                    M0 = M1;
                    M1 = M2;
                }

                return M1;
            }
        }

        /// <inheritdoc/>
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 0.0;
            } else {
                // This expression for the cumulant appears in the Wikipedia article.
                return MoreMath.Pow(2.0, r - 1) * (nu + r * lambda) * AdvancedIntegerMath.Factorial(r - 1);
            }
        }

        private double[] Cumulants (int rMax) {
            // This is just a recursive reformulation of the cumulant expression.
            double[] K = new double[rMax + 1];
            K[0] = 0.0;
            K[1] = nu + lambda;
            double t = 1.0;
            for (int r = 2; r <= rMax; r++) {
                t *= 2.0 * (r - 1);
                K[r] = t * (nu + r * lambda);
            }
            return K;
        }

        /// <inheritdoc/>
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else {
                double[] K = Cumulants(r);
                return MomentMath.CumulantToCentral(K, r);
            }
        }

        /// <inheritdoc/>
        public override double GetRandomValue (Random rng) {
            if (rng is null) throw new ArgumentNullException(nameof(rng));
            double mu = Math.Sqrt(lambda / nu);
            double x = 0.0;
            for (int k = 0; k < nu; k++) {
                x += MoreMath.Sqr(mu + zRng.GetNext(rng));
            }
            return x;
        }

    }
}

