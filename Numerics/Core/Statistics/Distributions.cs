
using System;
using Meta.Numerics;
using Meta.Numerics.Functions;

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
				return( Math.Sqrt( MomentAboutMean(2) ) );
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
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (n == 0) {
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
	public class NormalDistribution : Distribution {

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
                double mu = Mean;
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
	public class ExponentialDistribution : Distribution {

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
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
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
                return (0.5 * Math.Pow(x2, nu2 - 1.0) * Math.Exp(-x2) / AdvancedMath.Gamma(nu2));
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
            return (base.InverseLeftProbability(P));
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
                // this is a temporary fix until be get a better one

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
        public override double Median {
            get {
                return (0.0);
            }
        }

	}

    /*
    /// <summary>
    /// Represents a Fischer distribution.
    /// </summary>
    /// <remarks><para>The Fisher distribution is the distribution of the F statistic (under the null hypothesis) in a F-test.</para></remarks>
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
            if (F < 0.0) {
                return (0.0);
            } else {
                double N = Math.Pow(nu1, 0.5 * nu1) * Math.Pow(nu2, 0.5 * nu2) /
                    AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2);
                return (N * Math.Pow(F, 0.5 * nu1 - 1.0) * Math.Pow(nu2 + nu1 * F, -0.5 * (nu1 + nu2)));
            }
		}

        /// <inheritdoc />
        public override double LeftProbability (double F) {
            if (F < 0.0) {
                return (0.0);
            } else {
                return (AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2, nu1 * F / (nu2 + nu1 * F)) / AdvancedMath.Beta(0.5 * nu1, 0.5 * nu2));
            }
		}

        /// <inheritdoc />
        public override double RightProbability (double F) {
            if (F < 0.0) {
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
				if (nu2 < 4.0) {
                    return(System.Double.PositiveInfinity);
                } else {
				    double m = Mean;
				    return( Math.Sqrt(2.0 * m * m * (nu1 + nu2 - 2.0) / nu1 / (nu2 - 4.0)) );
			    }
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
                if (nu2 <= 2.0 * n) {
                    return (System.Double.PositiveInfinity);
                } else {
                    return (Math.Exp(n * Math.Log(nu2 / nu1) +
                        AdvancedMath.LogGamma(0.5 * nu1 + n) +
                        AdvancedMath.LogGamma(0.5 * nu1 - n) -
                        AdvancedMath.LogGamma(0.5 * nu1) -
                        AdvancedMath.LogGamma(0.5 * nu2)));
                }
            }
		}

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
			throw new NotImplementedException();
		}
		
	}
    */


    /*
    public class WeibullDistribution : Distribution {

        private int k;
        private double lambda;

        public WeibullDistribution (int k, double lambda) {
            this.k = k;
            this.lambda = lambda;
        }

        public override double PDF (double x1) {
            return( (k/lambda) * Math.Pow(x1/lambda, k-1) * Math.Exp(-Math.Pow(x1/lambda, k)) );
        }

        public override double CDF (double x1) {
            throw new Exception("The method or operation is not implemented.");
        }

    }
    */

    /*
    public class LogNormalDistribution : Distribution {

        public LogNormalDistribution (double mu, double sigma) {
            this.mu = mu;
            this.sigma = sigma;
        }

        private double mu;
        private double sigma;

        public override double PDF (double x1) {
            if (x1 <= 0.0) {
                return (0.0);
            } else {
                double x3 = (Math.Log(x1) - mu)/sigma;
                return (Math.Exp(-x3 * x3 / 2.0) / x1 / sigma / (2.0 * Math.PI)); 
            }
        }

        public override double CDF (double x1) {
            if (x1 <= 0.0) {
                return (0.0);
            } else {

            }
            throw new Exception("The method or operation is not implemented.");
        }

        public override double Median {
            get {
                return (Math.Exp(mu));
            }
        }

        public override double Mean {
            get {
                return (Math.Exp(mu + sigma * sigma / 2.0));
            }
        }

        public override double StandardDeviation {
            get {
                return base.StandardDeviation;
            }
        }

    }
    */

    /// <summary>
    /// Represents a Kolmogorov distribution.
    /// </summary>
    /// <remarks><para>The D statistic in a Kolmogorov-Smirnov test is distributed (under the null hypothesis) according to a Kolmogorov disribution.</para></remarks>
    /// <seealse cref="Sample.KolmogorovSmirnovTest(Meta.Numerics.Statistics.Distribution)" />
    public class KolmogorovDistribution : Distribution {

        /// <summary>
        /// Instantiates a new Kolmogorov distribution.
        /// </summary>
        public KolmogorovDistribution () {
        }

        internal KolmogorovDistribution (double scale) {
            this.scale = scale;
        }

        private double scale = 1.0;

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x < scale) {
                return (p_SmallX(x / scale) / scale);
            } else {
                return (p_LargeX(x / scale) / scale);
            }
        }

        private static double p_SmallX (double x) {

            /*
            double x1 = Math.PI / x / 2.0;
            double x2 = x1 * x1;
            double r = Math.Exp(- x2 / 2.0);
            double rr = r * r;
            double f = 0.0;
            for (int k = 1; k < Global.SeriesMax; k += 2) {
                double f_old = f;
                double df = r * (k * x2 - 1.0);
                f += df;
                if (f == f_old) {
                    return (Math.Sqrt(2.0 * Math.PI) / x / x * f);
                }
                r = r * rr;
            }
            */

            
            double p = 0.0;
            for (int k = 1; k < Global.SeriesMax; k += 2) {
                double p_old = p;
                double z = k * Math.PI / x / 2.0;
                double dp = Math.Exp(-z * z / 2.0) * (z * z - 1.0);
                p += dp;
                //Console.WriteLine("{0} {1} {2} {3} {4}", k, z, Math.Exp(-z * z / 2.0), dp, p);
                if (p == p_old) return (Math.Sqrt(2.0 * Math.PI) / (x * x) * p);
            }
            
            
            /*
            for (int k = 0; k < Global.SeriesMax; k++) {
                double p_old = p;
                double kp = Math.PI * (2 * k + 1);
                double t = kp * kp / 4.0 / x;
                double dp = (t - 1.0) * Math.Exp(-t / 2.0 / x);
                p = p_old + dp;
                if (p == p_old) {
                    return (Math.Sqrt(2.0 * Math.PI) / x * p);
                }
            }
            */

            throw new NonconvergenceException();
        }

        private static double p_LargeX (double x) {

            double p = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double p_old = p;
                double kx = k * x;
                double dp = k * k * Math.Exp(-2.0 * kx * kx);
                if (k % 2 == 0) dp = -dp;
                p += dp;
                if (p == p_old) return (8.0 * p * x);
            }

            /*
            double xx = 2.0 * x * x;
            double f = 0.0;
            int sign = 1;
            for (int k = 1; k < 100; k++) {
                double f_old = f;
                int kk = k * k;
                double df = sign * kk * x * Math.Exp(-kk * xx);
                f = f_old + df;
                if (f == f_old) return (8.0 * f);
                sign = -sign;
            }
            */

            throw new NonconvergenceException();
        }


        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x < scale) {
                return (P_SmallX(x/scale));
            } else {
                return (1.0 - Q_LargeX(x/scale));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else if (x < scale) {
                return (1.0 - P_SmallX(x/scale));
            } else {
                return (Q_LargeX(x/scale));
            }
        }


        // implements \frac{\sqrt{2\pi}}{x1} \sum{k=0}{\infty} e^{ \frac{(2k+1)^2 \pi^2}{8 x1^2} }
        // convergence is rapid; 4 terms at x~1 and still just 10 terms at x~3

        private static double P_SmallX (double x) {
            if (x == 0) {
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

                /*
                double xx = Math.PI / x;
                double df = Math.Exp(-xx * xx / 8.0);
                double ff = df * df;
                double f = df;
                for (int k = 1; k < Global.SeriesMax; k++) {
                    double f_old = f;
                    df = df * ff;
                    f += df;
                    if (f == f_old) {
                        return (Math.Sqrt(2.0 * Math.PI) / x * f);
                    }
                }
                */

                /*
                double xi = Math.PI / x;
                double xx = xi * xi / 8.0;
                double f = 0.0;
                for (int k = 1; k < Global.SeriesMax; k++) {
                    double f_old = f;
                    int n = 2*k - 1;
                    double df = Math.Exp(- n * n * xx);
                    f = f_old + df;
                    if (f == f_old) {
                        return (Math.Sqrt(2.0 * Math.PI) / x * f);
                    }
                }
                */

                throw new NonconvergenceException();
            }
        }

        // implements \sum_{k=-\infty}^{\infty} (-1)^k e^{-2 k^2 x1^2}
        // convergence is very rapid; 5 terms at x~1 and just 2 terms at x~3

        private static double Q_LargeX (double x) {
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
        public override double StandardDeviation {
            get {
                double ln2 = Math.Log(2.0);
                return (Math.Sqrt(Math.PI / 2.0 * (Math.PI / 6.0 - ln2*ln2)) * scale);
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
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else if (n == 2) {
                return (Math.PI * Math.PI / 12.0 * scale * scale);
            } else {
                return (AdvancedMath.Gamma(n / 2.0 + 1.0) * AdvancedMath.DirichletEta(n) / Math.Pow(2.0, n / 2.0 - 1.0) * Math.Pow(scale,n));
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
                // this isn't great, but it does the job
                // expand in terms of moments about the origin
                // there is likely to be some cancelation, but the distribution is wide enough that it may not matter
                double mm = -Mean;
                double C = 0.0;
                for (int k = 0; k <= n; k++) {
                    C += AdvancedIntegerMath.BinomialCoefficient(n, k) * Moment(k) * Math.Pow(mm, n - k);
                }
                return (C);
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
    /// Represents a log-normal distribution.
    /// </summary>
    /// <remarks>
    /// <para>A logrithm of a log-normal distributed variable is distributed normally.</para>
    /// <para>The log-normal distribution is commonly used in financial engineering as a model of stock prices.
    /// If the rate of return on an asset is distributed normally, then the price will be distribution log-normally.</para>
    /// </remarks>
    /// <seealso cref="NormalDistribution"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Log-normal_distribution" />
    public class LognormalDistribution : Distribution {

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
    public class WeibullDistribution : Distribution {

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
                // this isn't great, but it does the job
                // expand in terms of moments about the origin
                // there is likely to be some cancelation, but the distribution is wide enough that it may not matter
                double mm = -Mean;
                double C = 0.0;
                for (int k = 0; k <= n; k++) {
                    C += AdvancedIntegerMath.BinomialCoefficient(n, k) * Moment(k) * Math.Pow(mm, n - k);
                }
                return (C);
            }
        }

    }

    // Deviates
    // Maximum likelyhood estimation

	// Weibull Distribution
	// Gamma Distribution
	// Beta Distribution
	// Kuiper Distribution
    // Logistic Distribution
    // Wald (inverse guassian) distribution

}

