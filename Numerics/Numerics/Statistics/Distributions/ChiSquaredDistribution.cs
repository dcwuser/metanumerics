using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {


    /// <summary>
    /// Represents a &#x3C7;<sup>2</sup> distribution.
    /// </summary>
    /// <remarks>
    /// <para>A chi squared distribution is an asymmetrical distribution ranging from zero to infinity with a peak near its
    /// number of degrees of freedom &#x3BD;. It is a one-parameter distribution determined entirely by the parameter nu.</para>
    /// <img src="../images/ChiSquaredPlot.png" />
    /// <para>The figure above shows the &#x3C7;<sup>2</sup> distribution for &#x3BD; = 6, as well as the normal distribution
    /// with equal mean and variance for reference.</para>
    /// <para>The sum of the squares of &#x3BD; independent standard-normal distributed variables is distributed as &#x3C7;<sup>2</sup>
    /// with &#x3BD; degrees of freedom.</para>
    /// <img src="../images/ChiSquaredFromNormal.png" />
    /// <para>The &#x3C7;<sup>2</sup> distribution appears in least-squares fitting as the distribution of the sum-of-squared-deviations
    /// under the null hypothesis that the model explains the data. For example, the goodness-of-fit statistic returned by the
    /// model our model fitting methods (<see cref="UncertainMeasurementSample{T}.FitToFunction"/>, <see cref="UncertainMeasurementSample{T}.FitToLinearFunction"/>,
    /// <see cref="UncertainMeasurementSample.FitToLine"/>, and others) follows a &#x3C7;<sup>2</sup> distribution.</para>
    /// </remarks>
    /// <seealso cref="ContingencyTable.PearsonChiSquaredTest"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Chi-square_distribution" />
    public sealed class ChiSquaredDistribution : Distribution {

        // internally, we use our Gamma distribution machinery to do our heavy lifting

        private int nu;
        private GammaDistribution gamma;

        /// <summary>
        /// Initializes a new &#x3C7;<sup>2</sup> distribution.
        /// </summary>
        /// <param name="nu">The number of degrees of freedom, which must be positive.</param>
        public ChiSquaredDistribution (int nu) {
            if (nu < 1) throw new ArgumentOutOfRangeException("nu");
            this.nu = nu;
            this.gamma = new GammaDistribution(nu / 2.0);
        }

        /// <summary>
        /// Gets the number of degrees of freedom &#x3BD; of the distribution.
        /// </summary>
        public int DegreesOfFreedom {
            get {
                return (nu);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            return (gamma.ProbabilityDensity(x / 2.0) / 2.0);
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            return (gamma.LeftProbability(x / 2.0));
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            return (gamma.RightProbability(x / 2.0));
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            return (2.0 * gamma.InverseLeftProbability(P));
        }

        // improve this
        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else if (n == 2) {
                return (nu * (nu + 2.0));
            } else {
                // nu ( nu + 2 ) ( nu + 4 ) ... (nu + 2n - 2 )
                double nu2 = nu / 2.0;
                return (Math.Exp(n * Global.LogTwo + AdvancedMath.LogGamma(nu2 + n) - AdvancedMath.LogGamma(nu2)));
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
                // use C_{n} = 2^n U(-n, 1-n-\mu/2, -\mu/2) where U is irregular confluent hypergeometric
                // use recursion U(a-1,b-1,z) = (1-b+z) U(a,b,z) + z a U(a+1,b+1,z) to derive
                // C_{n+1} = 2n (C_{n} + \nu C_{n-1})
                double C1 = 0.0;
                double C2 = 2.0 * nu;
                for (int k = 2; k < n; k++) {
                    double C3 = (2*k) * (C2 + nu * C1);
                    C1 = C2;
                    C2 = C3;
                }
                return(C2);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (nu);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (2.0 * nu);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                // start with an approximation
                double xm = nu - (2.0 / 3.0) + (4.0 / 27.0) / nu - (8.0 / 729.0) / (nu * nu);
                // polish it using Newton's method
                while (true) {
                    double dx = (0.5 - LeftProbability(xm)) / ProbabilityDensity(xm);
                    xm += dx;
                    if (Math.Abs(dx) <= Global.Accuracy * xm) break;
                }
                return (xm);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (Math.Sqrt(8.0 / nu));
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