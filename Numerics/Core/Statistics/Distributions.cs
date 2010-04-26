using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a probability distribution.
    /// </summary>
	public abstract class Distribution {

        /// <summary>
        /// Returns the probability density at the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The probability density p(x1).</returns>
		public abstract double ProbabilityDensity (double x);

        /// <summary>
        /// Returns the cumulative probability to the left of (below) the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The integrated probability P(x1) to obtain a result below the reference point.</returns>
		public abstract double LeftProbability (double x);

        /// <summary>
        /// Return the cumulative probability to the right of (above) the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The integrated probability 1-P(x1) to obtain a result above the reference point.</returns>
		public virtual double RightProbability (double x) {
			return( 1.0 - LeftProbability(x) );
		}

        /// <summary>
        /// Returns the point at which the cumulative distribution function attains a given value. 
        /// </summary>
        /// <param name="P">The left cumulative probability P, which must lie between 0 and 1.</param>
        /// <returns>The point x1 at which the left cumulative probability attains the value P.</returns>
		public virtual double InverseLeftProbability (double P) {
            // a quick little implementation using Newton's method 
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            Function<double,double> f = delegate(double x) {
                return(LeftProbability(x) - P);
            };
            double y = FunctionMath.FindZero(f, Mean);
            return(y);
            // since we have the PDF = CDF', change to a method using Newton's method
		}

        /// <summary>
        /// Returns the given moment of the distribution.
        /// </summary>
        /// <param name="n">The order of the moment to determine.</param>
        /// <returns>The moment M<sub>n</sub> about the origin.</returns>
        /// <seealso cref="MomentAboutMean"/>
		public virtual double Moment (int n) {
			// inheritors may implement
			throw new NotImplementedException();
		}

        /// <summary>
        /// Returns the given moment of the distribution, about the mean. 
        /// </summary>
        /// <param name="n">The order of the moment to determine.</param>
        /// <returns>The moment of order n about the mean.</returns>
        /// <seealso cref="Moment" />
		public virtual double MomentAboutMean (int n) {
			// inheritors may implement
			throw new NotImplementedException();
		}

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
		public virtual double Mean {
			get {
				return(Moment(1));
			}
		}

        /// <summary>
        /// Gets the standard deviation of the distribution.
        /// </summary>
		public virtual double StandardDeviation {
			get {
				return( Math.Sqrt(Variance) );
			}
		}

        /// <summary>
        /// Gets the variance of the distribution.
        /// </summary>
        public virtual double Variance {
            get {
                return (MomentAboutMean(2));
            }
        }

        /// <summary>
        /// Ges the skewness of the distribution.
        /// </summary>
        /// <remarks>
        /// <para>The skewness of a distribution is a measurement of its asymmetry about its mean.
        /// It is the third moment about the mean, measured in units of the cubed standard deviation.</para>
        /// </remarks>
        public virtual double Skewness {
            get {
                return (MomentAboutMean(3) / Math.Pow(MomentAboutMean(2), 3.0 / 2.0));
            }
        }

        /// <summary>
        /// Gets the median of the distribution.
        /// </summary>
        /// <remarks>The median is the point with equal integrated probability above and below, i.e. with P(x1) = 0.5.</remarks>
        public virtual double Median {
            get {
                return (InverseLeftProbability(0.5));
            }
        }

        /// <summary>
        /// Gets the interval over which the distribution is nonvanishing.
        /// </summary>
        public virtual Interval Support {
            get {
                return (Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity));
            }
        }

        // compute central moments from raw moments
        // this is subject to loss of precision from cancelation, so be careful

        internal virtual double CentralMomentFromRawMoment (int n) {

            double m = Mean;

            double mm = 1.0;
            double C = Moment(n);
            for (int k = 1; k <= n; k++) {
                mm = mm * (-m);
                C += AdvancedIntegerMath.BinomialCoefficient(n,k) * mm * Moment(n - k);
            }

            return (C);

        }

        // compute raw moments from central moments
        // this doesn't suffer from the cancelation problem of the reverse calculation

        internal virtual double RawMomentFromCentralMoments (int n) {

            double m = Mean;

            double mm = 1.0;
            double M = MomentAboutMean(n);
            for (int k = 1; k <= n; k++) {
                mm = mm * m;
                M += AdvancedIntegerMath.BinomialCoefficient(n, k) * mm * MomentAboutMean(n - k);
            }

            return (M);

        }

    }


    /// <summary>
    /// Represents a uniform distribution over an interval.
    /// </summary>
	public class UniformDistribution : Distribution {

		private Interval range;

        /// <summary>
        /// Gets the range of the uniform distribution.
        /// </summary>
		public Interval Range {
			get {
				return(range);
			}
		}

        /// <summary>
        /// Initializes a new uniform distribution on the given interval.
        /// </summary>
        /// <param name="range">The range of the distribution.</param>
		public UniformDistribution (Interval range) {
			this.range = range;
		}

        /// <inheritdoc />
		public override double ProbabilityDensity (double x) {
			if (range.ClosedContains(x)) {
				return( 1.0 / (range.Width) );
			} else {
				return(0.0);
			}
		}

        /// <inheritdoc />
        public override double LeftProbability (double x) {
			if (x < range.LeftEndpoint) {
				return(0.0);
			} else if (x > range.RightEndpoint) {
				return(1.0);
			} else {
				return( (x - range.LeftEndpoint) / range.Width );
			}
		}

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
			if ((P<0.0) || (P>1.0)) throw new ArgumentOutOfRangeException("P");
			return( range.LeftEndpoint * (1.0-P) + range.RightEndpoint * P );
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

                double m = range.Midpoint;
                double w = range.Width;

                if (Math.Abs(m) > w) {

                    // the midpoint is greater than the width
                    // start from the approximate value m^n and compute a correction factor in terms of (w/m)

                    double f = Math.Pow(m, n);

                    double r = w / m;
                    double rr = 4.0 * r * r;

                    double dg = 1.0;
                    double g = dg;
                    for (int k = 2; k <= n; k += 2) {
                        dg = dg * rr;
                        g += g * AdvancedIntegerMath.BinomialCoefficient(n, k) / (k + 1.0);
                    }

                    return (f * g);

                } else {
                    // the width is large compared to the midpoint
                    // it should be safe to do a simple subtraction of endpoint powers

                    return ((Math.Pow(range.RightEndpoint, n + 1) - Math.Pow(range.LeftEndpoint, n + 1)) / range.Width / (n + 1));
                }
            }
		}

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if ((n % 2) != 0) {
                return (0.0);
            } else {
                return (Math.Pow(range.Width / 2.0, n) / (n + 1));
            }
		}

        /// <inheritdoc />
        public override double Mean {
			get {
				return(range.Midpoint);
			}
		}

        /// <inheritdoc />
        public override double StandardDeviation {
			get {
				return( range.Width / Math.Sqrt(12.0) );
			}
		}

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (Mean);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (range);
            }
        }

	}

    /// <summary>
    /// Represents a normal (Gaussian) distribution.
    /// </summary>
    /// <remarks>A normal distribution is a bell-shaped curve centered at its mean and falling off symmetrically on each side. It is
    /// a two-parameter distribution determined by giving its mean and standard deviation, i.e. its center and width. Its range is the
    /// entire real number line, but the tails, i.e. points more than a few standard deviations from the means, fall off extremely rapidly.
    /// <para>A normal distribution with mean zero and standard deviation one is called a standard normal distribution. Any normal distribution
    /// can be converted to a standard normal distribution by reparameterzing the data in terms of "standard deviations from the mean",
    /// i.e. z = (x - &#x3BC;) / &#x3C3;.</para>
    /// <para>Normal distribution appear in many contexts. In practical work, the normal distribution is often used as a crude
    /// model for the distribution of any continuous parameter that tends to cluster near its average, for example human height
    /// and weight. In more refined theoretical work, the normal distribution often emerges as a limiting distribution. For example,
    /// it can be shown that, if a large number of errors affect a measurement, then for nearly any underlying distribution
    /// of error terms, the distribution of total error tends to a normal distribution.</para>
    /// <para>The normal distribution is sometimes called a Gaussian distribtuion, after the mathematician Friedrich Gauss.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Normal_distribution"/>
	public class NormalDistribution : Distribution, IParameterizedDistribution {

		private double mu, sigma;

		private static readonly double c = 1.0/Math.Sqrt(2.0 * Math.PI);

        /// <summary>
        /// Initializes a new normal distribution with the given mean and standard deviation.
        /// </summary>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The standard deviation, which must be positive.</param>
		public NormalDistribution (double mu, double sigma) {
            if (sigma <= 0.0) throw new ArgumentOutOfRangeException("sigma");
            this.mu = mu;
			this.sigma = sigma;
		}

        /// <summary>
        /// Initializes a new standard normal distribution.
        /// </summary>
        /// <remarks>The standard normal distribution has mean zero and standard deviation one.</remarks>
		public NormalDistribution () : this(0.0,1.0) {}

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
			double z = (x - mu)/sigma;
			return( (c/sigma) * Math.Exp(-z*z/2.0) );
		}

		private static readonly double Sqrt2 = Math.Sqrt(2.0);

        /// <inheritdoc />
        public override double LeftProbability (double x) {
			double z = (x - mu)/sigma;
			if (z<0.0) {
				return( 0.5 * AdvancedMath.Erfc(-z/Sqrt2) );
			} else {
				return( 0.5 * (1.0 + AdvancedMath.Erf(z/Sqrt2)) );
			}
		}

        /// <inheritdoc />
        public override double RightProbability (double x) {
			double z = (x - mu)/sigma;
			if (z<0.0) {
				return( 0.5 * (1.0 + AdvancedMath.Erf(-z/Sqrt2)) );
			} else {
				return( 0.5 * AdvancedMath.Erfc(z/Sqrt2) );
			}
		}

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else {
                // use an expansion in powers of the mean and moments about the mean
                // i think this still has problems
                //double mu = Mean;
                double M = 0.0;
                for (int k = 0; k <= n; k=k+2) {
                    M += AdvancedIntegerMath.BinomialCoefficient(n, k) * Math.Pow(mu, n-k) * MomentAboutMean(k);
                }
                return (M);
            }
		}

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
			if (n < 0) throw new ArgumentOutOfRangeException("n");
			if ((n % 2) == 0) {
                // compute (n-1)!! sigma^n
                double sigma2 = sigma * sigma;
                double m = 1.0;
                for (int j = 1; j < n; j = j + 2) {
                    m *= sigma2 * j;
                }
				return(m);
			} else {
				return(0.0);
			}
		}

        /// <inheritdoc />
        public override double Mean {
			get {
				return(mu);
			}
		}

        /// <inheritdoc />
        public override double StandardDeviation {
			get {
				return(sigma);
			}
		}

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (mu);
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");

            double z = Math.Sqrt(2.0) * AdvancedMath.InverseErf(2.0 * P - 1.0);

            return (mu + sigma * z);
        }

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { mu, sigma });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 2) throw new DimensionMismatchException();
            if (parameters[1] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            mu = parameters[0];
            sigma = parameters[1];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }

	}

    /// <summary>
    /// Represents an exponential distribution.
    /// </summary>
    /// <remarks>
    /// <para>An exponential distribution falls off exponentially in the range from zero to infinity. It is a one-parameter
    /// distribution, determined entirely by its rate of fall-off.</para>
    /// <para>The exponential distribution describes the distribution of decay times of radioactive particles.</para>
    /// </remarks>
    /// <seealso href="WeibullDistribution"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Exponential_distribution"/>
	public class ExponentialDistribution : Distribution, IParameterizedDistribution {

		private double mu;

        /// <summary>
        /// Initializes a new exponential distribution with the given mean.
        /// </summary>
        /// <param name="mu">The mean, which must be positive.</param>
		public ExponentialDistribution (double mu) {
			if (mu <= 0.0) throw new ArgumentOutOfRangeException("mu");
			this.mu = mu;
		}

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x < 0.0) {
                return (0.0);
            } else {
                return (Math.Exp(-x / mu) / mu);
            }
		}

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x < 0.0) {
                return (0.0);
            } else {
                if (x < 0.25 * mu) {
                    // if x1 << mu, use the series expansion
                    double xx = x / mu;
                    double df = xx;
                    double f = df;
                    for (int k = 2; k < 250; k++) {
                        double f_old = f;
                        df = -df * xx / k;
                        f += df;
                        if (f == f_old) return (f);
                    }
                    throw new NonconvergenceException();
                } else {
                    // otherwise, just use the analytic expression
                    return (1.0 - Math.Exp(-x / mu));
                }
            }
		}

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < 0.0) {
                return (1.0);
            } else {
                return (Math.Exp(-x / mu));
            }
		}

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
			if ((P<0.0) || (P>1.0)) throw new ArgumentOutOfRangeException("P");
            // if P is very small, use the inverted series
			return(-mu * Math.Log(1.0-P));
		}

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else {
                return (Math.Pow(mu, n) * AdvancedIntegerMath.Factorial(n));
            }
		}

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else {
                double f = 0.0;
                double df = -1.0;
                for (int k = 2; k <= n; k++) {
                    df = -df / k;
                    f += df;
                }
                return (Math.Pow(mu, n) * AdvancedIntegerMath.Factorial(n) * f);
            }
		}

        /// <inheritdoc />
        public override double Mean {
			get {
				return(mu);
			}
		}

        /// <inheritdoc />
        public override double StandardDeviation {
			get {
				return(mu);
			}
		}

        /// <inheritdoc />
        public override double Median {
            get {
                return (mu * Math.Log(2.0));
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (2.0);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { mu });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 1) throw new DimensionMismatchException();
            if (parameters[0] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            mu = parameters[0];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }

	}

    /// <summary>
    /// Represents a &#x3C7;<sup>2</sup> distribution.
    /// </summary>
    /// <remarks>
    /// <para>A chi squared distribution is an asymmetrical distribution ranging from zero to infinity with a peak near its
    /// number of degrees of freedom &#x3BD; (nu). It is a one-parameter distribution determined entirely by the parameter nu.</para>
    /// <para>Technically, the &#x3C7;<sup>2</sup> distribution is the distribution of a sum of the squares of &#x3BD; independent
    /// standard-normally distributed variables. In practice it is most often encountered in data fitting or contingency testing.</para>
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
				return(nu);
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
                return (AdvancedMath.LeftGamma(0.5 * nu, 0.5 * x));
            }
		}

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < 0.0) {
                return (1.0);
            } else {
                return (AdvancedMath.RightGamma(0.5 * nu, 0.5 * x));
            }
		}

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {

			if ((P<0.0) || (P>1.0)) throw new ArgumentOutOfRangeException("P");

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
                double z = Math.Sqrt(2.0) * AdvancedMath.InverseErf(2.0 * P - 1.0);
                double x = Mean + StandardDeviation * z;
                // if we have run off end due to small P, use inverted small-x expansion instead
                if (x <= 0.0) {
                    double nu2 = nu / 2.0;
                    x = 2.0 * Math.Exp((Math.Log(P) + AdvancedMath.LogGamma(nu2 + 1.0)) / nu2);
                }
                // start a search from our initial guess
                Function<double,double> f = delegate(double y) {
                    return(LeftProbability(y) - P);
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
                //Console.WriteLine("n = {0}, nu2+n={1}, L1={2}, L2={3}, n ln 2 = {4}", n, nu2 + n, AdvancedMath.LogGamma(nu2+n), AdvancedMath.LogGamma(nu2), n * Math.Log(2.0));
                //Console.WriteLine("arg = {0}", n * Math.Log(2.0) + AdvancedMath.LogGamma(nu2 + n) - AdvancedMath.LogGamma(nu2));
                //Console.WriteLine("exp = {0}", Math.Exp(n * Math.Log(2.0) + AdvancedMath.LogGamma(nu2 + n) - AdvancedMath.LogGamma(nu2)));
                return (Math.Exp(n * Math.Log(2.0) + AdvancedMath.LogGamma(nu2 + n) - AdvancedMath.LogGamma(nu2)));
            }
		}

        // implement this
        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else if (n == 2) {
                return (2.0 * nu);
            } else if (n == 3) {
                return (8.0 * nu);
            } else {

                // uses confluent hypergeometric function
                // 2^n U(-n, 1-n-mu/2,-mu/2)
                // this is a temporary fix until we get a better one

                double C;
                switch (n) {
                    case 4:
                        C = 12.0 * nu * (4.0 + nu);
                        break;
                    case 5:
                        C = 32.0 * nu * (12.0 + 5.0 * nu);
                        break;
                    case 6:
                        C = 40.0 * nu * (96.0 + 52.0 * nu + 3.0 * nu * nu);
                        break;
                    case 7:
                        C = 96.0 * nu * (480.0 + 308.0 * nu + 35.0 * nu * nu);
                        break;
                    default:
                        throw new NotImplementedException();
                }
                return(C);
            }
		}

        /// <inheritdoc />
        public override double Mean {
			get {
				return(nu);
			}
		}

        /// <inheritdoc />
        public override double StandardDeviation {
			get {
				return(Math.Sqrt(2.0*nu));
			}
		}

        /// <inheritdoc />
        public override double Median {
            get {
                // start with an approximation
                double xm = nu - (2.0 / 3.0) + (4.0 / 27.0) / nu - (8.0 / 729.0) / (nu * nu);
                // polish it using Newton's method
                while (true) {
                    double dx = (0.5-LeftProbability(xm))/ProbabilityDensity(xm);
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

    /// <summary>
    /// Represents the distribution of Student't t statistic.
    /// </summary>
    /// <remarks><para>The mean of n independent standard-normal distributed variables, divided by their root mean square,
    /// is distributed according a Student distribution with n degrees of freedom. Since this is the form of the expression
    /// for the mean of a sample divided by its standard deviation, the Student distribution expresses the distribution of
    /// sample means arround the population mean, for a normally distributed population.</para></remarks>
    /// <seealso cref="Sample.StudentTTest(Double)"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Student_t_distribution" />
	public class StudentDistribution : Distribution {

		private double nu;

        /// <summary>
        /// Initializes a new Student distribution.
        /// </summary>
        /// <param name="nu">The number of degrees of freedom, which must be positive.</param>
		public StudentDistribution (double nu) {
			if (nu <= 0.0) throw new ArgumentOutOfRangeException("nu", nu, "Degrees of freedom must be positive.");
			this.nu = nu;
		}

        /// <summary>
        /// Gets the number of degrees of freedom.
        /// </summary>
		public double DegreesOfFreedom {
			get {
				return(nu);
			}
		}

        /// <inheritdoc />
        public override double ProbabilityDensity (double t) {
			return( Math.Pow(1.0+t*t/nu, -0.5*(nu+1.0)) / AdvancedMath.Beta(0.5,0.5*nu) / Math.Sqrt(nu) );
		}

        /// <inheritdoc />
        public override double LeftProbability (double t) {
			double t2 = t*t;
			if (t<0.0) {
				return( 0.5 * AdvancedMath.Beta(0.5*nu, 0.5, nu/(nu+t2))/AdvancedMath.Beta(0.5, 0.5*nu) );
			} else {
				return( 0.5 *(1.0 + AdvancedMath.Beta(0.5, 0.5*nu,  t2/(nu+t2))/AdvancedMath.Beta(0.5, 0.5*nu)) );
			}
		}

        /// <inheritdoc />
        public override double RightProbability (double t) {
			double t2 = t*t;
			if (t<0.0) {
				return( 0.5*(1.0 + AdvancedMath.Beta(0.5, 0.5*nu, t2/(nu+t2))/AdvancedMath.Beta(0.5, 0.5*nu)) );
			} else {
				return( 0.5 * AdvancedMath.Beta(0.5*nu, 0.5, nu/(nu+t2))/AdvancedMath.Beta(0.5, 0.5*nu) );
			}
		}

        /// <inheritdoc />
        public override double Moment (int n) {
			if (n < 0) throw new ArgumentOutOfRangeException("n");
 
                if ((n % 2) == 0) {
                    if (n < nu) {
                        double m = 1.0;
                        for (int j = 1; j <= n/2; j++) {
                            m *= (2 * j - 1) * nu / (nu - 2 * j);
                        }
                        return (m);
                    } else {
                        return (Double.PositiveInfinity);
                    }
                } else {
                    return (0.0);
                }
 
		}

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
			return( Moment(n) );
		}

        /// <inheritdoc />
        public override double Mean {
			get {
				return(0.0);
			}
		}

        /// <inheritdoc />
        public override double StandardDeviation {
			get {
                if (nu <= 2.0) {
                    return (Double.PositiveInfinity);
                } else {
                    return (Math.Sqrt(nu / (nu - 2.0)));
                }
			}
		}

        /// <inheritdoc />
        public override double Skewness {
            get {
                return(0.0);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (0.0);
            }
        }

	}

    /// <summary>
    /// Represents a Fischer distribution.
    /// </summary>
    /// <remarks><para>The Fisher distribution is the distribution of the F statistic (under the null hypothesis) in a F-test.</para></remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/F_distribution"/>
	public class FisherDistribution : Distribution {

		private double nu1;
		private double nu2;

        /// <summary>
        /// Gets the number of degrees of freedom in the numerator.
        /// </summary>
		public double NumeratorDegreesOfFreedom {
			get {
				return(nu1);
			}
		}

        /// <summary>
        /// Gets the number of degrees of freedom in the denominator.
        /// </summary>
		public double DenominatorDegreesOfFreedom {
			get {
				return(nu2);
			}
		}

        /// <summary>
        /// Instantiates a new Fisher distribution.
        /// </summary>
        /// <param name="nu1">The number of degrees of freedom in the numerator, which must be positive.</param>
        /// <param name="nu2">The number of degrees of freedom in the denominator, which must be positive.</param>
        /// <seealso cref="Sample.FisherFTest"/>
		public FisherDistribution (double nu1, double nu2) {
			if (nu1 <= 0.0) throw new ArgumentOutOfRangeException("nu1");
			if (nu2 <= 0.0) throw new ArgumentOutOfRangeException("nu2");
			this.nu1 = nu1;
			this.nu2 = nu2;
		}

        /// <inheritdoc />
        public override double ProbabilityDensity (double F) {
            if (F <= 0.0) {
                return (0.0);
            } else {
                double N = Math.Pow(nu1, 0.5 * nu1) * Math.Pow(nu2, 0.5 * nu2) /
                    AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2);
                return (N * Math.Pow(F, 0.5 * nu1 - 1.0) * Math.Pow(nu2 + nu1 * F, -0.5 * (nu1 + nu2)));
            }
		}

        /// <inheritdoc />
        public override double LeftProbability (double F) {
            if (F <= 0.0) {
                return (0.0);
            } else {
                return (AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2, nu1 * F / (nu2 + nu1 * F)) / AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2));
            }
		}

        /// <inheritdoc />
        public override double RightProbability (double F) {
            if (F <= 0.0) {
                return (1.0);
            } else {
                return (AdvancedMath.Beta(0.5 * nu2, 0.5 * nu1, nu2 / (nu2 + nu1 * F)) / AdvancedMath.Beta(0.5 * nu2, 0.5 * nu1));
            }
		}

        /// <inheritdoc />
        public override double Mean {
			get {
                if (nu2 < 2.0) {
                    return (System.Double.PositiveInfinity);
                } else {
                    return (nu2 / (nu2 - 2.0));
                }
			}
		}

        /// <inheritdoc />
        public override double StandardDeviation {
			get {
				if (nu2 <= 4.0) {
                    return(System.Double.PositiveInfinity);
                } else {
				    double m = Mean;
				    return( Math.Sqrt(2.0 * m * m * (nu1 + nu2 - 2.0) / nu1 / (nu2 - 4.0)) );
			    }
            }
		}

        /// <inheritdoc />
        public override double Variance {
            get {
                if (nu2 <= 4.0) {
                    return (System.Double.PositiveInfinity);
                } else {
                    double m = Mean;
                    return (2.0 * m * m * (nu1 + nu2 - 2.0) / nu1 / (nu2 - 4.0));
                }
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
                if (nu2 <= 2.0 * n) {
                    return (System.Double.PositiveInfinity);
                } else {

                    // (nu2)^n  Gamma(nu1/2 + n)  Gamma(nu2/2 - n)
                    // (---)    ----------------  ----------------
                    // (nu1)      Gamma(nu1/2)      Gamma(nu2/2)

                    // this can be computed using the recursion relation for the Gamma function

                    double r = nu2 / nu1;
                    double M = 1.0;
                    for (int k = 0; k < n; k++) {
                        M = M * r * (nu1 + 2*k) / (nu2 - 2*(k+1));
                    }

                    return (M);

                }
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
            } else if (n == 2) {
                return (Variance);
            } else if (n == 3) {
                if (nu2 <= 6.0) {
                    return (Double.PositiveInfinity);
                } else {
                    double m = Mean;
                    return (8 * m * m * m * (2 * nu1 + nu2 - 2.0) * (nu1 + nu2 - 2.0) / (nu1 * nu1) / (nu2 - 4.0) / (nu2 - 6.0));
                }
            } else {
                if (nu2 < 2.0 * n) {
                    return (Double.PositiveInfinity);
                } else {
                    return (CentralMomentFromRawMoment(n));
                }
            }
		}
		
	}


    /// <summary>
    /// Represents the distribution of the Kolmogorov-Smirnov D statistic.
    /// </summary>
    /// <remarks><para>The D statistic in a Kolmogorov-Smirnov test is distributed (under the null hypothesis) according to a Kolmogorov disribution.</para></remarks>
    /// <seealse cref="Sample.KolmogorovSmirnovTest(Meta.Numerics.Statistics.Distribution)" />
    public class KolmogorovDistribution : Distribution {

        /// <summary>
        /// Instantiates a new asymptotic Kolmogorov distribution.
        /// </summary>
        public KolmogorovDistribution () {}

        // the sample size; when N=0 we will report the asymptotic distribution

        internal KolmogorovDistribution (double scale) {
            this.scale = scale;
        }

        private double scale = 1.0;

        /// <inheritdoc />
        public override double ProbabilityDensity (double d) {

            if (d < 1.0) {
                return (AsymptoticPPrime(d));
            } else {
                return (AsymptoticQPrime(d));
            }

        }

        // the asymptotic PDF for x <~ 1

        private static double AsymptoticPPrime (double x) {

            if (x <= 0.0) return (0.0);
            
            double p = 0.0;
            for (int k = 1; k < Global.SeriesMax; k += 2) {
                double p_old = p;
                double z = k * Math.PI / x / 2.0;
                double dp = Math.Exp(-z * z / 2.0) * (z * z - 1.0);
                p += dp;
                if (p == p_old) return (Math.Sqrt(2.0 * Math.PI) / (x * x) * p);
            }

            throw new NonconvergenceException();
        }

        // the asymptotic PDF for x >~ 1

        private static double AsymptoticQPrime (double x) {

            double p = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double p_old = p;
                double kx = k * x;
                double dp = k * k * Math.Exp(-2.0 * kx * kx);
                if (k % 2 == 0) dp = -dp;
                p += dp;
                if (p == p_old) return (8.0 * p * x);
            }

            throw new NonconvergenceException();
        }


        /// <inheritdoc />
        public override double LeftProbability (double d) {

            if (d < scale) {
                return (AsymptoticP(d/scale));
            } else {
                return (1.0 - AsymptoticQ(d/scale));
            }

        }

        /// <inheritdoc />
        public override double RightProbability (double d) {

            if (d < scale) {
                return (1.0 - AsymptoticP(d/scale));
            } else {
                return (AsymptoticQ(d/scale));
            }

        }

        // implements \frac{\sqrt{2\pi}}{x1} \sum{k=0}{\infty} e^{ \frac{(2k+1)^2 \pi^2}{8 x1^2} }
        // convergence is rapid; 4 terms at x~1 and still just 10 terms at x~3

        private static double AsymptoticP (double x) {

            if (x <= 0.0) {
                return (0.0);
            } else {

                double p = 0.0;
                for (int k = 1; k < Global.SeriesMax; k += 2) {
                    double p_old = p;
                    double z = k * Math.PI / x / 2.0;
                    double dp = Math.Exp(-z * z / 2.0);
                    p += dp;
                    //Console.WriteLine("{0} {1} {2} {3} {4}", k, z, Math.Exp(-z * z / 2.0), dp, p);
                    if (p == p_old) return (Math.Sqrt(2.0 * Math.PI) / x * p);
                }

                throw new NonconvergenceException();
            }
        }

        // implements \sum_{k=-\infty}^{\infty} (-1)^k e^{-2 k^2 x1^2}
        // convergence is very rapid; 5 terms at x~1 and just 2 terms at x~3

        private static double AsymptoticQ (double x) {
            double xx = x * x;
            double f = 0.0;
            int sign = -1;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                sign = - sign;
                double df = sign * Math.Exp(-(2 * k * k) * xx);
                f = f_old + df;
                if (f == f_old) {
                    return (2.0 * f);
                }
            }
            throw new NonconvergenceException();
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (Math.Sqrt(Math.PI / 2.0) * Math.Log(2.0) * scale);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                double ln2 = Math.Log(2.0);
                return (Math.PI / 2.0 * (Math.PI / 6.0 - ln2 * ln2) * scale * scale);

            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.Sqrt(Variance));
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (0.82757355518991 * scale);
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
            } else if (n == 2) {
                return (Math.PI * Math.PI / 12.0 * scale * scale);
            } else {
                return (AdvancedMath.Gamma(n / 2.0 + 1.0) * AdvancedMath.DirichletEta(n) / Math.Pow(2.0, n / 2.0 - 1.0) * Math.Pow(scale, n));
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

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /*

        // an exact formula for 1/2 <= t <= 1

        private static double Smallest_P (int n, double t) {

            Debug.Assert((0.5 <= t) && (t <= 1.0));

            return (Math.Exp(n * Math.Log((2.0 * t - 1.0) / n) + AdvancedIntegerMath.LogFactorial(n)));

        }

        // an exact formula for n-1 <= t <= n

        private static double Smallest_Q (int n, double t) {

            Debug.Assert((n - 1) <= t);

            return (2.0 * Math.Pow(1.0 - t / n, n));

        }
        */

    }

#if FUTURE

    public class FiniteKolmogorovDistribution : KolmogorovDistribution {

        public FiniteKolmogorovDistribution (int size) {
            if (size < 1) throw new ArgumentOutOfRangeException("size");
            N = size;
            sqrtN = Math.Sqrt(N);
        }

        int N;
        double sqrtN;

        int maxN = 256;

        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.5 / N, 1.0));
            }
        }

        public override double ProbabilityDensity (double d) {

            if (N < maxN) {

                double t = d * N;
                if (2.0 * t < N) {
                    return (DurbinPPrime(N, t) * N);
                } else {
                    return (DurbinQPrime(N, t) * N);
                }

            } else {

                return (base.ProbabilityDensity(d * sqrtN) * sqrtN);

            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double d) {

            if (d <= 1.0 / (2.0 * N)) {
                return (0.0);
            } else if (d >= 1.0) {
                return (1.0);
            } else {

                if (N < maxN) {

                    // use Durbin formulas for small N

                    double t = d * N;
                    if (2.0 * t < N) {
                        // the Durbin series formula is faster than the durbin matrix formula,
                        // but it has alternating sign terms and and can suffer a catastrophic
                        // loss of accuracy
                        return (MatrixP(N, t));
                        //return (DurbinP(N, t));
                    } else {
                        return (1.0 - DurbinQ(N, t));
                    }

                } else {
                    return (base.LeftProbability(d * sqrtN));
                }

            }

        }

        /// <inheritdoc />
        public override double RightProbability (double d) {

            if (d <= 1.0 / (2.0 * N)) {
                return (1.0);
            } else if (d >= 1.0) {
                return (0.0);
            } else {

                if (N < maxN) {

                    // use Durbin formulas for small N

                    double t = d * N;
                    if (2.0 * t < N) {
                        return (1.0 - MatrixP(N, t));
                        //return (1.0 - DurbinP(N, t));
                    } else {
                        return (DurbinQ(N, t));
                    }

                } else {

                    return (base.RightProbability(d * sqrtN));

                }
            }

        }

            
        // Durbin's formula for exact P_n(t) that holds for n > Truncate(2 t), i.e. small t
        // we holding an array of P_m(t) for m < n; this keeps us having to reevaluate repeatedly

        // the formula unfortunately involves canceling terms of increasing size;
        // it breaks down when n gets large

        private static double DurbinP (int n, double t) {

            if (t <= 0.5) return (0.0);

            int t2 = (int) Math.Truncate(2.0 * t);

            //double s = 0.0;
            //for (int j = 1; j <= t2; j++) {
            //    double ds = Math.Exp(
            //        j * Math.Log(2.0 * t - j) - AdvancedIntegerMath.LogFactorial(j) +
            //        (n - j) * Math.Log(n - j) - AdvancedIntegerMath.LogFactorial(n - j)
            //    ) * DurbinP1(n - j, t);
            //    if (j % 2 == 0) ds = - ds;
            //    s += ds;
            //}
            //s = Math.Exp(AdvancedIntegerMath.LogFactorial(n) - n * Math.Log(n)) * s;
            //return (s);


            // an array for P_m(t)
            double[] P = new double[n+1];
            P[0] = 1.0;

            // populate up to m = t2 using the Durbin Q formula
            for (int m = 1; m <= t2; m++) {
                P[m] = 1.0 - DurbinQ(m, t);
            }

            // populate higher m using the recurrsion relation
            for (int m = t2 + 1; m <= n; m++) {
                double B = 1.0; // binomial coefficient (n j)
                double s = 0.0;
                for (int j = 1; j <= t2; j++) {
                    B = B * (m - (j - 1)) / j;
                    double C = Math.Pow((2.0 * t - j) / m, j) * Math.Pow(1.0 * (m - j) / m, m - j);
                    double ds = B * C * P[m - j];
                    if (j % 2 == 0) ds = -ds;
                    s += ds;
                }
                P[m] = s;
            }

            /*
            for (int i = 0; i <= n; i++) {
                Console.WriteLine("p[{0}]={1}", i, P[i]);
            }
            */

            return (P[n]);

        }

        // Durbin's formula for exact Q_n(t) that holds for 2t > n

        private static double DurbinQ (int n, double t) {

            if (t >= n) return (0.0);

            double s = Math.Pow(1.0 - t / n, n); // j = 0 term
            int jmax = (int) Math.Truncate(n - t);
            double B = 1.0; // binomial coefficient ( n j )
            for (int j = 1; j <= jmax; j++) {
                B = B * (n - (j-1)) / j;
                double C = Math.Pow((t + j) / n, j - 1) * Math.Pow((n - j - t) / n, n - j) * t / n;
                s += B * C;
            }
            return (2.0 * s);

        }

        private static double DurbinQPrime (int n, double t) {

            if (t >= n) return (0.0);

            double s = 2.0 * Math.Pow(1.0 - t / n, n - 1);
            int jmax = (int) Math.Truncate(n - t);
            if (jmax == (n - t)) jmax--;
            for (int j = 1; j <= jmax; j++) {
                double B = Math.Exp(AdvancedIntegerMath.LogFactorial(n-1) - AdvancedIntegerMath.LogFactorial(j) - AdvancedIntegerMath.LogFactorial(n - j));
                double C = Math.Pow((t + j) / n, j - 1) * Math.Pow((n - j - t) / n, n - j);
                double D = 1.0 + (j-1) * t / (t + j) - (n-j) * t / (n - j - t);
                s += -2.0 * B * C * D;
            }
            return (s);
        }

        private static double DurbinPPrime (int n, double t) {

            int t2 = (int) Math.Truncate(2.0 * t);

            // arrays for PDF and CDF for m <= n
            double[] p = new double[n+1];
            double[] P = new double[n+1];
            p[0] = 0.0;
            P[0] = 1.0;

            // populate up to m = t2 using the Durbin Q formula
            for (int m = 1; m <= t2; m++) {
                P[m] = 1.0 - DurbinQ(m, t);
                p[m] = DurbinQPrime(m, t);
            }

            // compute higher p and P using the recursion formula
            for (int m = t2 + 1; m <= n; m++) {
                double B = 1.0; // binomial coefficient (m j)
                double s = 0.0;
                double sp = 0.0;
                double jmax = t2;
                if (t2 == 2.0 * t) jmax--;
                for (int j = 1; j <= jmax; j++) {
                    B = B * (m - (j - 1)) / j;
                    double C = Math.Pow((2.0 * t - j) / m, j) * Math.Pow(1.0 * (m - j) / m, m - j);
                    double ds = B * C * P[m - j];
                    double dsp = B * C * (p[m - j] + 2.0 * j / (2.0 * t - j) * P[m - j]);
                    if (j % 2 == 0) {
                        ds = - ds;
                        dsp = - dsp;
                    }
                    s += ds;
                    sp += dsp;
                }
                P[m] = s;
                p[m] = sp;
            }

            /*
            for (int i = 0; i <= n; i++) {
                Console.WriteLine("P[{0}]={1} p[{0}]={2}", i, P[i], p[i]);
            }
            */

            // return the desired PDF value
            return (p[n]);

        }

        // Durbin's matrix form, also programed by Marsaglia
        // all number are positive, so this does not suffer from the cancelation probmlems of Durbin's recursion

        private static double MatrixP (int n, double t) {

            // compute stuff used in matrix entries
            int tp = (int) Math.Truncate(t) + 1;
            double h = tp - t;
            int p = 2 * tp - 1;


            // construct the matrix
            SquareMatrix H = new SquareMatrix(p);

            // superdiagonal
            for (int j = 1; j < p; j++) {
                H[j-1,j] = 1.0;
            }

            // diagonal and subdiagonals
            double F = 1.0; // factorial
            double hh = h; // power of h
            for (int i = 1; i < p; i++) {
                H[i - 1, 0] = (1.0 - hh) / F;
                H[p-1, p-i] = H[i-1,0];
                for (int j = i+1; j < p; j++) {
                    H[j - 1, j - i] = 1.0 / F;
                }
                hh = hh * h;
                F = F * (i+1);
            }

            // lower-left element
            double g = 1.0 - 2.0 * hh;
            if (h > 0.5) g = g + Math.Pow(2.0 * h - 1.0, p);
            g = g / F;
            H[p-1,0] = g;

            // raise the matrix to the nth power
            SquareMatrix HN = MatrixPower(H, n);

            // return the appropriate element
            double hf = Math.Exp(AdvancedIntegerMath.LogFactorial(n) - n * Math.Log(n));
            return (hf * HN[tp-1, tp-1]);


        }

        private static SquareMatrix MatrixPower (SquareMatrix A, int n) {

            SquareMatrix B = null;

            SquareMatrix D = A.Clone();

            while (true) {
                if (n % 2 != 0) {
                    if (B == null) {
                        B = D.Clone();
                    } else {
                        B = B * D;
                    }
                }
                n = n / 2;
                if (n == 0) break;
                D = D * D;
            }

            return (B);


        }

        ///<inheritdoc />
        public override double Mean {
            get {
                if (N < maxN) {

                    if (N == 1) {
                        return ((3.0 / 4.0) / N);
                    } else if (N == 2) {
                        return ((13.0 / 12.0) / N);
                    } else if (N == 3) {
                        return ((293.0 / 216.0) / N);
                    }

                    throw new NotImplementedException();
                } else {
                    return (base.Mean / sqrtN);
                }
            }
        }

        ///<inheritdoc />
        public override double Variance {
            get {
                if (N < maxN) {
                    if (N == 1) {
                        return (1.0 / 48.0 / N);
                    } else if (N == 2) {
                        return ((7.0 / 72.0) / N);
                    }
                    throw new NotImplementedException();
                } else {
                    return (base.Variance / N);
                }
            }
        }

        ///<inheritdoc />
        public override double Moment (int n) {
            throw new NotImplementedException();
        }

        ///<inheritdoc />
        public override double MomentAboutMean (int n) {
            throw new NotImplementedException();
        }

    }

#endif

    /// <summary>
    /// Represents the asymptotic distribution of Kuiper's V statistic.
    /// </summary>
    public class KuiperDistribution : Distribution {

        /// <summary>
        /// Instantiates a new Kuiper distribution.
        /// </summary>
        public KuiperDistribution () {
        }

        internal KuiperDistribution (double scale) {
            this.scale = scale;
        }

        private double scale = 1.0;

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x < 1.0) {
                return (AsymptoticP(x/scale));
            } else {
                return (1.0 - AsymptoticQ(x/scale));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < 1.0) {
                return (1.0 - AsymptoticP(x/scale));
            } else {
                return (AsymptoticQ(x/scale));
            }
        }

        // series \sqrt{2\pi}{x^3} \sum_{k=1}^{\infty} k^2 \pi^2 e^{-k^2 \pi^2 / 2 x^2}
        // useful for small x

        private static double AsymptoticP (double x) {

            if (x <= 0.0) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;
                double z = k * Math.PI / x;
                double z2 = z * z;
                double ds = z2 * Math.Exp(-z2 / 2.0);
                s += ds;

                if (s == s_old) return (Math.Sqrt(2.0 * Math.PI) / x * s);

            }

            throw new NonconvergenceException();

        }

        // series \sum_{k=1}^{\infty} (4 k^2 x^2 - 1) e^{-2 k^2 x^2}
        // useful for large x

        private static double AsymptoticQ (double x) {

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (4.0 * z2 - 1.0) * Math.Exp(-2.0 * z2);
                s += ds;

                if (s == s_old) return (2.0 * s);
            }

            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x < 1.0) {
                return (AsymptoticPPrime(x/scale)/scale);
            } else {
                return (AsymptoticQPrime(x/scale)/scale);
            }
        }

        private static double AsymptoticQPrime (double x) {

            double s = 0.0;
            for (int k = 1; k < 20; k++) {

                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (k * k) * (4.0 * z2 - 3.0) * Math.Exp(-2.0 * z2);
                s += ds;

                if (s == s_old) return (8.0 * x * s);
            }

            throw new NonconvergenceException();

        }

        private static double AsymptoticPPrime (double x) {

            if (x <= 0.0) return (0.0);

            double s = 0.0;
            for (int k = 1; k < 20; k++) {

                double s_old = s;
                double z = Math.PI * k / x;
                double z2 = z * z;
                double ds = z2 * (z2 - 3.0) * Math.Exp(-z2 / 2.0);
                s += ds;

                if (s == s_old) return (Math.Sqrt(2.0 * Math.PI) / (x * x) * s);

            }

            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                // from numerical integral; would be nice to get an analytic result
                return (1.2533141373155 * scale);
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
            } else if (n == 2) {
                return (Math.PI * Math.PI / 6.0 * scale * scale);
            } else {
                return (AdvancedMath.RiemannZeta(n) * AdvancedMath.Gamma(1 + n / 2.0) * (n - 1) / Math.Pow(2.0, n / 2.0 - 1) * Math.Pow(scale, n));
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

    }

    /// <summary>
    /// Represents a log-normal distribution.
    /// </summary>
    /// <remarks>
    /// <para>A logrithm of a log-normal distributed variable is distributed normally.</para>
    /// <para>The log-normal distribution is commonly used in financial engineering as a model of stock prices.
    /// If the rate of return on an asset is distributed normally, then the price will be distribution log-normally.</para>
    /// </remarks>
    /// <seealso cref="NormalDistribution"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Log-normal_distribution" />
    public class LognormalDistribution : Distribution, IParameterizedDistribution {

        /// <summary>
        /// Initializes a log normal distribution.
        /// </summary>
        /// <param name="mu">The mean of the underlying normal distribution.</param>
        /// <param name="sigma">The standard deviation of the underlying normal distribution.</param>
        /// <remarks>
        /// <para>Note that the <paramref name="mu"/> and <paramref name="sigma"/> parameters are
        /// <em>not</em> the mean and standard deviation of the log-normal deviates; they are the
        /// mean and standard deviation of their logrithms z = ln x.</para>
        /// </remarks>
        public LognormalDistribution (double mu, double sigma) {
            if (sigma <= 0.0) throw new ArgumentOutOfRangeException("sigma");
            this.mu = mu;
            this.sigma = sigma;
        }

        private double mu = 0.0;
        private double sigma = 1.0;

        /// <inheritdoc />
        public override double ProbabilityDensity (double y) {
            if (y <= 0.0) return(0.0);
            double z = (Math.Log(y) - mu) / sigma;
            return( Math.Exp(-z*z / 2.0) / y / (Math.Sqrt(2.0 * Math.PI) * sigma));
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (Math.Exp(mu + sigma * sigma / 2.0));
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (Math.Exp(mu));
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.Sqrt(Es2m1) * Mean);
            }
        }

        // compute e^(sigma^2) - 1, accounting for cancelation with sigma is small
        // this expression appears in several contexts for this distribution
        private double Es2m1 {
            get {
                double sigma2 = sigma * sigma;
                if (sigma2 < 1.0E-4) {
                    return(sigma2 * (1.0 + sigma2 / 2.0 + sigma2 * sigma2 / 6.0 + sigma2 * sigma2 * sigma2 / 24.0));
                } else {
                    return(Math.Exp(sigma2) - 1.0);
                }
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double y) {
            if (y <= 0.0) return (0.0);
            double z = (Math.Log(y) - mu) / sigma;
            if (z < 0.0) {
                return (0.5 * AdvancedMath.Erfc(-z / Math.Sqrt(2.0)));
            } else {
                return (0.5 * (1.0 + AdvancedMath.Erf(z / Math.Sqrt(2.0))));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double y) {
            if (y <= 0.0) return (1.0);
            double z = (Math.Log(y) - mu) / sigma;
            if (z < 0.0) {
                return (0.5 * (1.0 + AdvancedMath.Erf(-z / Math.Sqrt(2.0))));
            } else {
                return (0.5 * AdvancedMath.Erfc(z / Math.Sqrt(2.0)));
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            double z = Math.Sqrt(2.0) * AdvancedMath.InverseErf(2.0 * P - 1.0);
            return (Math.Exp(mu + sigma * z));
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            return (Math.Exp(n * mu + n * n * sigma * sigma / 2.0));
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else if (n == 2) {
                return (Es2m1 * Mean * Mean);
            } else {
                // this isn't great, but it does the job
                // expand in terms of moments about the origin
                // there is likely to be some cancelation, but the distribution is wide enough that it may not matter
                double m = -Mean;
                double C = 0.0;
                for (int k = 0; k <= n; k++) {
                    C += AdvancedIntegerMath.BinomialCoefficient(n, k) * Moment(k) * Math.Pow(m, n - k);
                }
                return (C);
            }

        }

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { mu, sigma });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 2) throw new DimensionMismatchException();
            if (parameters[1] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            mu = parameters[0];
            sigma = parameters[1];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }


    }

    /// <summary>
    /// Represents a Weibull distribution.
    /// </summary>
    /// <remarks>
    /// <para>The Weibull distribution is a generalized form of the exponential distriubtion
    /// for which the decay probability is not constant, but instead increases with time
    /// (for shape parameters greater than one) or, less commonly, decreases with time (for
    /// shape parameters less than one). When the shape parameter is one, the Weibull
    /// distribution reduced to the exponential distribution.</para>
    /// <para>The Weibull distribution is commonly used in engineering applications to
    /// model the time-to-failure of industrial componets.</para>
    /// </remarks>
    public class WeibullDistribution : Distribution, IParameterizedDistribution {

        /// <summary>
        /// Initializes a new Weibull distribution.
        /// </summary>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        /// <param name="shape">The shape parameter, which must be positive.</param>
        public WeibullDistribution (double scale, double shape) {
            if (scale <= 0.0) throw new ArgumentOutOfRangeException("scale");
            if (shape <= 0.0) throw new ArgumentOutOfRangeException("shape");
            this.scale = scale;
            this.shape = shape;
        }

        private double scale;

        private double shape;

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            double z = x / scale;
            double zz = Math.Pow(z, shape);
            return (Math.Exp(-zz) * zz * shape / x);
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double z = x / scale;
                double zz = Math.Pow(z, shape);
                if (zz < 0.5) {
                    // for small zz, use a series to avoid cancelation (and avoid doing exponentiation)
                    double dP = zz;
                    double P = dP;
                    for (int k = 2; k < Global.SeriesMax; k++) {
                        double P_old = P;
                        dP = - dP * zz / k;
                        P += dP;
                        if (P == P_old) return (P);
                    }
                    throw new NonconvergenceException();
                } else {
                    return (1.0 - Math.Exp(-zz));
                }
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else {
                double z = x / scale;
                double zz = Math.Pow(z, shape);
                return (Math.Exp(-zz));
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (scale * Math.Pow(-Math.Log(1.0 - P), 1.0 / shape));
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (scale * AdvancedMath.Gamma(1.0 + 1.0 / shape));
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (scale * Math.Pow(Math.Log(2.0), 1.0 / shape));
            }
        }

        // standard deviation inherits from base; nothing special to say about it

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0.0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                return (Math.Pow(scale, n) * AdvancedMath.Gamma(1.0 + n / shape));
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

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { scale, shape });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 2) throw new DimensionMismatchException();
            if (parameters[0] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            if (parameters[1] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            scale = parameters[0];
            shape = parameters[1];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }

    }

    /// <summary>
    /// Represents the distribution of the Mann-Whitney statistic.
    /// </summary>
    /// <remarks>
    /// <para>The expected use of this class is as the TestResult.Distribution of a
    /// Mann-Whitney U-test. You probably don't want to use this distribution directly
    /// unless you are studying the Mann-Whitney U-test.</para>
    /// </remarks>
    public class MannWhitneyDistribution : Distribution {

        /// <summary>
        /// Instantiates a new Mann-Whitney distribution.
        /// </summary>
        /// <param name="m">The number of elements in the first sample.</param>
        /// <param name="n">The number of elements in the second sample.</param>
        public MannWhitneyDistribution (int m, int n) {

            if (m < 1) throw new ArgumentOutOfRangeException("m");
            if (n < 1) throw new ArgumentOutOfRangeException("n");

            mc = m;
            nc = n;

            total = Math.Exp(AdvancedIntegerMath.LogFactorial(n + m) - AdvancedIntegerMath.LogFactorial(n) - AdvancedIntegerMath.LogFactorial(m));

            if (0.99 * total < ((double) Decimal.MaxValue)) {

                // for a total small enough to fit into a decimal, compute the exact count of partitions
                // contributing to each u; this operation is m^2 n, but it's integer math and we are
                // only doing it up to m~n~50, so it's still pretty fast

                counts = GaussianBinomialCoefficients(m + n, m);

                decimal sum = 0M;
                for (int u = 0; u <= m * n; u++) {
                    sum += counts[u];
                }
                total = (double) sum;

            } else {

                // for larger values, we will use a normal approximation

                normal = new NormalDistribution(Mean, StandardDeviation);

            }

            
            // i'd like to use an edgeworth approximation, and I have formulas for the required higher
            // cumulants, but edgeworth breaks down for extreme values (giving negative probabilities
            // and probabilities larger than one); i could transition back to a normal approximation
            // at those extremes, but i don't know how to do this and keep all my quantities self-consistent,
            // i.e. p is derivative of P, moment functions return actual moment integrals, etc.

            /*
            double M1 = m * n / 2.0;
            double C2 = 0.0;
            double C4 = 0.0;
            double C6 = 0.0;
            for (int i = 0; i <= m * n; i++) {
                double z = (i - M1);
                double z2 = z * z;
                double z4 = z2 * z2;
                double z6 = z4 * z2;
                C2 += ((double) counts[i]) * z2;
                C4 += ((double) counts[i]) * z4;
                C6 += ((double) counts[i]) * z6;
            }
            C2 = C2 / ((double) total);
            C4 = C4 / ((double) total);
            C6 = C6 / ((double) total);

            double m2 = m * m;
            double n2 = n * n;
            double m3 = m2 * m;
            double n3 = n2 * n;
            double m4 = m2 * m2;
            double n4 = n2 * n2;

            double S1 = m * n / 2.0;
            double S2 = m * n * (m + 1.0 + n) / 12.0;
            double S3 = 0.0;
            double S4 = - m * n * (m + n + 1.0) * (m2 + m + m * n + n + n2) / 120.0;            
            double S6 = m * n * (m + n + 1.0) * (
                2.0 * m4 + 4.0 * m3 * n + 4.0 * m3 + 7.0 * m2 * n + m2 - m +
                2.0 * n4 + 4.0 * m * n3 + 4.0 * n3 + 7.0 * m * n2 + n2 - n +
                6.0 * m2 * n2 + 2.0 * m * n ) / 504.0;

            Console.WriteLine("M1 {0} {1}", M1, S1);
            Console.WriteLine("C2 {0} {1}", C2, S2);
            Console.WriteLine("C4 {0} {1}", C4, S4 + 3.0 * S2 * S2);
            Console.WriteLine("C6 {0} {1}", C6, S6 + 15.0 * S4 * S2 + 10.0 * S3 * S3 + 15.0 * S2 * S2 * S2);
            */

        }

        private int mc, nc;
        private double total;
        private decimal[] counts;
        private NormalDistribution normal;

        public override double  ProbabilityDensity(double x) {
 	        throw new NotImplementedException();
        }

        private double Probability (int u) {
            if ((u < 0) || (u > mc*nc)) return (0.0);

            if (counts != null) {
                return ((double) counts[u] / total);
            } else {
                return (normal.ProbabilityDensity(u));
            }

        }

        /// <inheritdoc />
        public override double  LeftProbability(double x) {
            return( LeftInclusiveProbability((int) Math.Truncate(x)) );
        }   

        private double LeftInclusiveProbability (int u) {

            if (u < 0) return (0.0);
            if (u > mc*nc) return (1.0);

            if (counts != null) {
                double P = 0;
                for (int i = 0; i <= u; i++) {
                    P += (double) counts[i];
                }
                return (P / total);
            } else {
                return( normal.LeftProbability(u) );
            }


        }

        /// <inheritdoc />
        public override double RightProbability (double u) {
            return( RightExclusiveProbability((int) Math.Truncate(u)) );
        }

        private double RightExclusiveProbability (int u) {

            if (u < 0) return (1.0);
            if (u > mc*nc) return (0.0);

            if (counts != null) {
                double Q = 0;
                for (int i = u + 1; i <= mc * nc; i++) {
                    Q += (double) counts[i];
                }
                return (Q / total);
            } else {
                return( normal.RightProbability(u) );
            }
        }

        
        // this routine is based on the recurrsion
        // [ m n ] = ( 1 - q^m ) / ( 1 - q^(m-n) ) [ m-1 n ]
        // and the starting point [ n n ] = 1

        // the coefficients are integers and get large quickly as m and n increase
        // we use decimal because it handles larger integers than long
        // we can't use double because the calculation requires delicate cancelations
        // among large intermediate values, thus necessicating exact integer arithmetic
        // look into using an arbitrary-sized integer structure in the future

        private decimal[] GaussianBinomialCoefficients (int m, int n) {

            if (m < 0) throw new ArgumentOutOfRangeException("m");
            if (n < 0) throw new ArgumentOutOfRangeException("n");

            Debug.Assert(m >= n);

            // create  an array to hold our coefficients
            decimal[] c = new decimal[(m - n) * n + 1];

            // start with [n n] = 1 * q^0
            c[0] = 1;

            // keep track of current degree of our polynomial
            int d = 0;

            // create a scratch array for intermediate use
            // it needs to be larger than the previous array by (m-n) to hold intermediate polynomials
            decimal[] b = new decimal[c.Length + (m-n)];

            // interate from [n n] up to [m n]
            for (int k = n + 1; k <= m; k++) {

                // multiply by (1-q^k)
                for (int i = 0; i <= d; i++) {
                    b[i] = c[i];
                }
                d = d + k;
                for (int i = k; i <= d; i++) {
                    b[i] = b[i] - c[i - k];
                }

                // divide by (1-q^(k-n))
                for (int i = d - (k - n); i >= 0; i--) {
                    c[i] = -b[k - n + i];
                    b[k - n + i] = b[k - n + i] + c[i];
                    b[i] = b[i] - c[i];
                }
                d = d - (k - n);

            }

            // we're done
            return (c);


        }

        // these are some other ways we tried to generate counts; they are not as efficient as
        // the gaussian binomial coefficient technique, but I keep them around for the record

        /*
        public static int[] ComputeProbabilities (int n, int m) {

            int[] counts = new int[n * m + 1];

            int nm = n + m;

            // initialze to the lexographically first sequence
            // we could save memory by doing this witha BitArray, but we loose
            // about a factor of 2 in performance
            bool[] sequence = new bool[nm];
            for (int i = 0; i < n; i++) {
                sequence[i] = false;
            }
            for (int i = n; i < nm; i++) {
                sequence[i] = true;
            }
            int W = 0;

            while (true) {

                counts[W] = counts[W] + 1;

                // find the descending-order string at the end (the first 0<-1 transition)
                int j = nm - 1;
                int n0 = 0;
                int n1 = 0;
                while (j > 0) {

                    if (!sequence[j]) {
                        n0++;
                    } else {
                        n1++;
                        if (!sequence[j - 1]) break;
                    }

                    j--;

                }

                // check whether we're done
                if (j == 0) break;

                // flip the 0 to a 1
                sequence[j - 1] = true;

                // put the remainder of the string in ascending-order (0s, then 1s)
                for (int i = 0; i < (n0 + 1); i++) {
                    sequence[j] = false;
                    j++;
                }
                for (int i = 0; i < (n1 - 1); i++) {
                    sequence[j] = true;
                    j++;
                }

                // change W
                int dW = (n0 + 1) * (n1 - 1) - n1;
                W -= dW;

            }

            int sum = 0;
            for (int i = 0; i < counts.Length; i++) {
                sum += counts[i];
                Console.Write(counts[i] + " ");
            }
            Console.WriteLine();
            Console.WriteLine(sum);


            return(counts);

        }

        public static int CountOrderings (int m, int n, int u) {
            
            if ((m < 0) || (n < 0) || (u < 0) || (u > m*n)) return(0);

            if (u == 0) return (1);

            return (CountOrderings(m - 1, n, u - n) + CountOrderings(m, n - 1, u));

        }
        */

        /// <inheritdoc />
        public override double Mean {
            get {
                return (mc * nc / 2.0);
            }
        }

        /// <inheritdoc />
        public override double Median {
	        get { 
		        return (Mean);
	        }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (mc * nc * (mc + nc + 1) / 12.0);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        // add moments

    }

    /// <summary>
    /// Represents a logistic distribution.
    /// </summary>
    /// <remarks>
    /// <para>Like the normal distribution, the logistic distribution is a symmetric, unimodal distribution
    /// distribution with exponentially supressed tails.</para>
    /// </remarks>
    public class LogisticDistribution : Distribution, IParameterizedDistribution {

        /// <summary>
        /// Initializes a new standard logistic distribution.
        /// </summary>
        /// <remarks>
        /// <para>The standard logistic distribution has zero mean and a unit width parameter.</para>
        /// </remarks>
        public LogisticDistribution () {
            m = 0.0;
            s = 1.0;
        }

        /// <summary>
        /// Initializes a new logistic distribution with the given mean and width parameters.
        /// </summary>
        /// <param name="m">The mean.</param>
        /// <param name="s">The width parameter.</param>
        public LogisticDistribution (double m, double s) {
            if (s <= 0.0) throw new ArgumentOutOfRangeException("s");
            this.m = m;
            this.s = s;
        }

        // mean and width parameter

        private double m, s;

        /// <inheritdoc />
        public override double Mean {
            get {
                return (m);
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (m);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.PI * s / Math.Sqrt(3.0));
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (Math.PI * Math.PI * s * s / 3.0);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            double z = (x - m) / s;
            double cosh = Math.Cosh(z / 2.0);
            return (1.0 / (4.0 * s * cosh * cosh));
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double z = (x - m) / s;
            return (1.0 / (1.0 + Math.Exp(-z)));
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double z = (x - m) / s;
            return (1.0 / (1.0 + Math.Exp(z)));
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || ( P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (m + s * Math.Log(P / (1.0 - P)));
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return(1.0);
            } else if (n == 1) {
                return (Mean);
            } else {
                return (RawMomentFromCentralMoments(n));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return(1.0);
            } else {
                if (n % 2 == 0) {
                    return (2.0 * AdvancedIntegerMath.Factorial(n) * AdvancedMath.DirichletEta(n) * Math.Pow(s, n));
                } else {
                    return (0.0);
                }
            }
        }

        double[] IParameterizedDistribution.GetParameters () {
            return (new double[] { m, s });
        }

        void IParameterizedDistribution.SetParameters (IList<double> parameters) {
            if (parameters == null) throw new ArgumentNullException("parameters");
            if (parameters.Count != 2) throw new DimensionMismatchException();
            if (parameters[1] <= 0.0) throw new ArgumentOutOfRangeException("parameters");
            m = parameters[0];
            s = parameters[1];
        }

        double IParameterizedDistribution.Likelihood (double x) {
            return (ProbabilityDensity(x));
        }

    }

    /// <summary>
    /// Represents an parameterized likelihood distribution.
    /// </summary>
    public interface IParameterizedDistribution {

        /// <summary>
        /// Gets the parameter values of the distribution.
        /// </summary>
        /// <returns></returns>
        double[] GetParameters ();

        /// <summary>
        /// Sets the parameter values of the distribution.
        /// </summary>
        /// <param name="parameters">A list of parameter values.</param>
        void SetParameters (IList<double> parameters);

        /// <summary>
        /// Gets the likelihood of a value, given the current parameters.
        /// </summary>
        /// <param name="x">The value.</param>
        /// <returns>The likelihood of the value.</returns>
        double Likelihood (double x);

    }

    /// <summary>
    /// Represents a triangular distribution.
    /// </summary>
    /// <remarks>
    /// <para>Like a uniform distribution, a triangular distribution is confined to a finite interval. Unlike a
    /// uniform distribution, a triangular distribution is not uniform across the interval.</para>
    /// <para>Triangular distributions are often used in project planning, where a maximum, minimum, and most
    /// likely value for some quantity is known or supposed.</para>
    /// </remarks>
    public class TriangularDistribution : Distribution {

        /*
        public TriangularDistribution (Interval range) {
            r = range;
            c = range.Midpoint;
        }

        private TriangularDistribution (Interval range, double peak) {
            if (!range.ClosedContains(peak)) throw new ArgumentOutOfRangeException("peak");
            r = range;
            c = peak;
        }
        */

        /// <summary>
        /// Initializes a new triangular distribution.
        /// </summary>
        /// <param name="a">One point of the distribution.</param>
        /// <param name="b">One point of the distribution.</param>
        /// <param name="c">One point of the distribution.</param>
        public TriangularDistribution (double a, double b, double c) {

            // check for validity
            if ((a == b) && (b == c)) throw new ArgumentException();

            // order points

            // record values
            this.a = a;
            this.b = b;
            this.c = c;

            this.ab = b - a;
            this.bc = c - b;
            this.ac = c - a;

            this.h = 2.0 / ac;

        }

        private double a, b, c;
        private double ab, bc, ac;
        double h;

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if ((x <= a) || (x >= c)) {
                return (0.0);
            } else {
                if (x < b) {
                    return (h * (x - a) / ab);
                } else {
                    return (h * (c - x) / bc);
                }
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= a) {
                return (0.0);
            } else if (x >= c) {
                return (1.0);
            } else {
                if (x < b) {
                    return(LeftTriangleArea(x));
                } else {
                    return(1.0 - RightTrigangleArea(x));
                }
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= a) {
                return (1.0);
            } else if (x >= c) {
                return (0.0);
            } else {
                if (x < b) {
                    return (1.0 - LeftTriangleArea(x));
                } else {
                    return (RightTrigangleArea(x));
                }
            }
        }

        private double LeftTriangleArea (double x) {
            Debug.Assert(x >= a); Debug.Assert(x <= b);
            double ax = x - a;
            return (ax * ax / ab / ac);
        }

        private double RightTrigangleArea (double x) {
            Debug.Assert(x >= b); Debug.Assert(x <= c);
            double xc = c - x;
            return (xc * xc / bc / ac);
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            double Pb = ab / ac;
            if (P < Pb) {
                return (a + Math.Sqrt(ab * ac * P));
            } else {
                return(c - Math.Sqrt(bc * ac * (1.0 - P)));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return ((a + b + c) / 3.0);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return ((ab * ab + bc * bc + ac * ac) / 36.0);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(a, c));
                //return (r);
            }
        }

        private double MomentAboutMode (int n) {
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return ((bc - ab) / 3.0);
            } else {
                // ac = ab + bc is always a factor of the numerator; find a way to divide it out analytically
                if (n % 2 == 0) {
                    return (2.0 * (Math.Pow(bc, n + 1) + Math.Pow(ab, n + 1)) / (n + 1) / (n + 2) / ac);
                } else {
                    return (2.0 * (Math.Pow(bc, n + 1) - Math.Pow(ab, n + 1)) / (n + 1) / (n + 2) / ac);
                }
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
                double M = MomentAboutMode(n);
                double t = 1.0;
                for (int k = n - 1; k >= 0; k--) {
                    t *= b;
                    M += AdvancedIntegerMath.BinomialCoefficient(n, k) * MomentAboutMode(k) * t;
                }
                return (M);
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
                double M = MomentAboutMode(n);
                double s = -MomentAboutMode(1);
                double t = 1.0;
                for (int k = n - 1; k >= 0; k--) {
                    t *= s;
                    M += AdvancedIntegerMath.BinomialCoefficient(n, k) * MomentAboutMode(k) * t;
                }
                //Console.WriteLine("n={0}, M={1}", n, M);
                return (M);
            }
        }
    }

    public class BetaDistribution : Distribution {

        /// <summary>
        /// Instantiates a new &#x3B2; distribution.
        /// </summary>
        /// <param name="alpha">The left shape parameter, which controls the form of the distribution near x=0.</param>
        /// <param name="beta">The right shape parameter, which controls the form of the distribution near x=1.</param>
        public BetaDistribution (double alpha, double beta) {
            if (alpha <= 0.0) throw new ArgumentOutOfRangeException("alpha");
            if (beta <= 0.0) throw new ArgumentOutOfRangeException("beta");
            this.alpha = alpha;
            this.beta = beta;
        }

        double alpha, beta;

        public double Alpha {
            get {
                return (alpha);
            }
        }

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
                return (AdvancedMath.Beta(alpha, beta, x));
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
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                return (AdvancedMath.Beta(alpha + n, beta) / AdvancedMath.Beta(alpha, beta));
            }
        }

    }

    /*
    public class WaldDistribution : Distribution {


        public WaldDistribution (double mu, double lambda) {
            if (mu <= 0.0) throw new ArgumentOutOfRangeException("mu");
            if (lambda <= 0.0) throw new ArgumentOutOfRangeException("lambda");
            this.mu = mu;
            this.lambda = lambda;
        }

        double mu;
        double lambda;

        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        public override double ProbabilityDensity (double x) {
            return (Math.Sqrt(lambda / (2.0 * Math.PI * x * x * x)) * Math.Exp(-lambda * (x - mu) / (2.0 * x * mu * mu)));
        }

        public override double Mean {
            get {
                return (mu);
            }
        }


        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                // use recursion relation M_{n+1} = ((2n-1) (M_n / lambda) + M_{n-1}) mu^2
                double MM = mu;
                if (n == 1) return (MM);
                double mu2 = mu * mu;
                double M0 = m2 * (lambda + mu) / lambda;
                if (n == 2) return (M0);
                for (int k = 2; k < n; k++) {
                    double MP = ((2 * k - 1) * (M0 / lambda) + MM) * mu2;
                    MM = M0;
                    M0 = MP;
                }
                return (M0);
            }
        }

    }
    */

    // Deviates
    // Maximum likelyhood estimation
    // Skewness
    // Cumulants

	// Gamma Distribution
	// Beta Distribution
    // Wald (inverse guassian) Distribution

}
