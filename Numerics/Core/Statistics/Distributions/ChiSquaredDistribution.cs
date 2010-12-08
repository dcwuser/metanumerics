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
    /// model our model fitting methods (<see cref="DataSet{T}.FitToFunction"/>, <see cref="DataSet{T}.FitToLinearFunction"/>,
    /// <see cref="DataSet.FitToLine"/>, and others) follows a &#x3C7;<sup>2</sup> distribution.</para>
    /// </remarks>
    /// <seealso cref="ContingencyTable.PearsonChiSquaredTest"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Chi-square_distribution" />
    public class ChiSquaredDistribution : Distribution {

        private int nu;

        /// <summary>
        /// Initializes a new &#x3C7;<sup>2</sup> distribution.
        /// </summary>
        /// <param name="nu">The number of degrees of freedom, which must be positive.</param>
        public ChiSquaredDistribution (int nu) {
            if (nu < 1) throw new ArgumentOutOfRangeException("nu");
            this.nu = nu;
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
            if (x < 0.0) {
                return (0.0);
            } else {
                double x2 = x / 2.0;
                double nu2 = nu / 2.0;
                return (Math.Exp((nu2 - 1.0) * Math.Log(x2) - x2 - AdvancedMath.LogGamma(nu2)) / 2.0);
                //return (0.5 * Math.Pow(x2, nu2 - 1.0) * Math.Exp(-x2) / AdvancedMath.Gamma(nu2));
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x < 0.0) {
                return (0.0);
            } else {
                return (AdvancedMath.LeftRegularizedGamma(0.5 * nu, 0.5 * x));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < 0.0) {
                return (1.0);
            } else {
                return (AdvancedMath.RightRegularizedGamma(0.5 * nu, 0.5 * x));
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {

            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");

            // limit cases
            if (P == 0.0) {
                return (0.0);
            } else if (P == 1.0) {
                return (Double.PositiveInfinity);
            }

            if (nu == 1) {
                // analytic formula for nu=1
                double e = AdvancedMath.InverseErf(P);
                return (2.0 * e * e);
            } else if (nu == 2) {
                // analytic formula for nu=2
                return (-2.0 * Math.Log(1.0 - P));
            } else {
                // try a normal approximation
                double z = Global.SqrtTwo * AdvancedMath.InverseErf(2.0 * P - 1.0);
                double x = Mean + StandardDeviation * z;
                // if we have run off end due to small P, use inverted small-x expansion instead
                if (x <= 0.0) {
                    double nu2 = nu / 2.0;
                    x = 2.0 * Math.Exp((Math.Log(P) + AdvancedMath.LogGamma(nu2 + 1.0)) / nu2);
                }
                // start a search from our initial guess
                Function<double, double> f = delegate(double y) {
                    return (LeftProbability(y) - P);
                };
                x = FunctionMath.FindZero(f, x);
                return (x);
                // use newton's method; this failed to converge -- look into this
                /*
                for (int k = 0; k < 10; k++ ) {
                    Console.WriteLine("x={0}", x);
                    double x_old = x;
                    Console.WriteLine("P={0}", LeftProbability(x));
                    Console.WriteLine("p={0}", ProbabilityDensity(x));
                    double dx = (P - LeftProbability(x)) / ProbabilityDensity(x);
                    Console.WriteLine("dx={0}", dx);
                    x += dx;
                    if (x == x_old) return (x);
                }
                throw new NonconvergenceException();
                */
                // we need to be careful; what happens if we get x <= 0 in this loop?
            }
            //return (base.InverseLeftProbability(P));
            // change to deal with Inverse so we can't run off the end
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