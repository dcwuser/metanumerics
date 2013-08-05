using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a continuous probability distribution.
    /// </summary>
	public abstract class Distribution {

        /// <summary>
        /// Returns the probability density at the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The probability density p(x).</returns>
        /// <remarks>
        /// <para>The probability density function (PDF) gives the relative probability of obtaining different values.</para>
        /// </remarks>
		public abstract double ProbabilityDensity (double x);

        /// <summary>
        /// Returns the cumulative probability to the left of (below) the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The integrated probability to obtain a result below the reference point.</returns>
        /// <remarks>
        /// <para>The left probability function is commonly called the cumulative distribution function (CDF).</para>
        /// </remarks>
        public virtual double LeftProbability (double x) {
            if (x <= Support.LeftEndpoint) {
                return (0.0);
            } else if (x >= Support.RightEndpoint) {
                return (1.0);
            } else {
                return (FunctionMath.Integrate(ProbabilityDensity, Interval.FromEndpoints(Support.LeftEndpoint, x)));
            }
        }

        /// <summary>
        /// Return the cumulative probability to the right of (above) the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The integrated probability 1-P(x1) to obtain a result above the reference point.</returns>
        /// <remarks>
        /// <para>In survival analysis, the right probability function is commonly called the survival function, because it gives the
        /// fraction of the population remaining after the given time.</para>
        /// </remarks>
		public virtual double RightProbability (double x) {
			return( 1.0 - LeftProbability(x) );
		}

        /// <summary>
        /// Returns the point at which the cumulative distribution function attains a given value. 
        /// </summary>
        /// <param name="P">The left cumulative probability P, which must lie between 0 and 1.</param>
        /// <returns>The value x at which <see cref="LeftProbability"/> equals P.</returns>
        /// <remarks>
        /// <para>The inverse left probability is commonly called the quantile function. Given a quantile,
        /// it tells which variable value is the lower border of that quantile.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="P"/> lies outside [0,1].</exception>
		public virtual double InverseLeftProbability (double P) {
            // find x where LeftProbability(x) = P 
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            Func<double,double> f = delegate(double x) {
                return(LeftProbability(x) - P);
            };
            double y = FunctionMath.FindZero(f, Mean);
            return(y);
            // since we have the PDF = CDF', change to a method using Newton's method
		}

        /// <summary>
        /// Returns the point at which the right probability function attains the given value.
        /// </summary>
        /// <param name="Q">The right cumulative probability, which must lie between 0 and 1.</param>
        /// <returns>The value x for which <see cref="RightProbability"/> equals Q.</returns>
        public virtual double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException("Q");
            Func<double, double> f = delegate(double x) {
                return (RightProbability(x) - Q);
            };
            double y = FunctionMath.FindZero(f, Mean);
            return (y);
        }

        /// <summary>
        /// Returns the given moment of the distribution.
        /// </summary>
        /// <param name="n">The order of the moment to determine.</param>
        /// <returns>The moment M<sub>n</sub> about the origin.</returns>
        /// <seealso cref="MomentAboutMean"/>
		public virtual double Moment (int n) {
            return (ExpectationValue(delegate(double x) { return (ProbabilityDensity(x) * MoreMath.Pow(x, n)); }));
		}

        /// <summary>
        /// Returns the given moment of the distribution, about the mean. 
        /// </summary>
        /// <param name="n">The order of the moment to determine.</param>
        /// <returns>The moment of order n about the mean.</returns>
        /// <seealso cref="Moment" />
		public virtual double MomentAboutMean (int n) {
			// inheritors may implement
            double mu = Mean;
            return (ExpectationValue(delegate(double x) {
                return (ProbabilityDensity(x) * MoreMath.Pow(x - mu, n));
            }));
		}

        internal virtual double Cumulant (int n) {
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

        /// <summary>
        /// Computes the expectation value of the given function.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <returns>The expectation value of the function.</returns>
        public virtual double ExpectationValue (Func<double, double> f) {
            return (FunctionMath.Integrate(f, Support));
        }

        /// <summary>
        /// Returns a random value.
        /// </summary>
        /// <param name="rng">A random number generator.</param>
        /// <returns>A number distributed according to the distribution.</returns>
        /// <remarks>
        /// <para>Note that the random number generator <paramref name="rng"/> will be advanced by this method. The next call to its
        /// generator methods will not give the same value as it would had it not been passed to this method.</para>
        /// </remarks>
        public virtual double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException("rng");
            return (InverseLeftProbability(rng.NextDouble()));
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
            int B = 1; // binomial coefficient; use recurrence B(n,k) = (n-k+1/k) B(n,k-1)
            for (int k = 1; k <= n; k++) {
                B = B * (n - k + 1) / k;
                mm = mm * m;
                M += B * mm * MomentAboutMean(n - k);
            }

            return (M);

        }

    }


    // Moment formulas used

    // Distribution             Cumulant            Central             Raw
    // ChiSquared               ?                   Safe Recurrence     Closed
    // Exponential              Closed              Closed              Sum
    // Normal                   Closed              Closed              Sum
#if FUTURE
    public static class MomentMath {

        // C_n = \sum_{k=0}^{n} \binom{n}{k} M_{n-k} (-\mu)^k

        public static double CentralMomentFromRawMoments (double[] rawMoments) {

            Debug.Assert(rawMoments[0] == 1.0);

            int n = rawMoments.Length - 1;

            if (n == 0) return (1.0);
            if (n == 1) return (0.0);

            // the first term
            double MM1 = - rawMoments[1];
            int k = n;
            IEnumerator<double> B = AdvancedIntegerMath.BinomialCoefficients(n).GetEnumerator(); B.MoveNext();
            double MP = 1.0;
            double C = rawMoments[k];

            // the remaining terms
            while (B.MoveNext()) {
                k--;
                MP *= MM1;
                C += B.Current * MP * rawMoments[k];
            }
            return (C);
        }

        // M_n = \sum_{k=0}^{n} \binom{n}{k} C_{n-k} \mu^k

        public static double RawMomentFromCentralMoments (double[] centralMoments, double mean) {
            int n = centralMoments.Length - 1;
            double MP = 1.0;
            double M = 0.0;
            IEnumerator<double> B = AdvancedIntegerMath.BinomialCoefficients(n).GetEnumerator();
            for (int i = n; i >= 0; i--) {
                B.MoveNext();
                M += MP * B.Current * centralMoments[i];
                MP *= mean;
            }
            return (M);
        }

        // Faa di Bruno's formula related raw moments to cumulants

        // M_n = \sum_{partitions of n} \frac{n!}{m_1 \cdots m_k} C_1 \cdots C_k

        public static double RawMomentFromCumulants (double[] cumulants) {
            int n = cumulants.Length - 1;
            double M = 0.0;
            // each partition of n corresponds to a term, with parition integer -> cumulant
            foreach (int[] partition in AdvancedIntegerMath.Partitions(n)) {
                double dM = AdvancedIntegerMath.Factorial(n);
                int m = 1; int v1 = 0;
                foreach (int v in partition) {
                    if (v == v1) {
                        m++;
                    } else {
                        m = 1;
                    }
                    dM *= cumulants[v] / AdvancedIntegerMath.Factorial(v) / m;
                    v1 = v;
                }
                M += dM;
            }
            return (M);
        }

        public static double CentralMomentFromCumulants (double[] cumulants) {
            int n = cumulants.Length - 1;
            double M = 0.0;
            // each partition of n corresponds to a term, with parition integer -> cumulant
            foreach (int[] partition in AdvancedIntegerMath.Partitions(n)) {
                double dM = AdvancedIntegerMath.Factorial(n);
                int m = 1; int v1 = 0;
                foreach (int v in partition) {
                    if (v == 0) {
                        dM = 0.0;
                        break;
                    } else if (v == v1) {
                        m++;
                    } else {
                        m = 1;
                    }
                    dM *= cumulants[v] / AdvancedIntegerMath.Factorial(v) / m;
                    v1 = v;
                }
                M += dM;
            }
            return (M);
        }

    }
#endif

    /// <summary>
    /// Represents an parameterized likelihood distribution.
    /// </summary>
    public interface IParameterizedDistribution {

        /// <summary>
        /// Gets the parameter values of the distribution.
        /// </summary>
        /// <returns>The parameter values characterizing the distribution.</returns>
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

#if FUTURE
    internal static class DistributionMath {

        public static double ComputeCentralMomentFromRawMoments (double[] rawMoments, int n) {

            throw new NotImplementedException();
        }

        public static double ComputeRawMomentFromCentralMoments (double mean, double[] centralMoments, int n) {

            if (centralMoments == null) throw new ArgumentNullException("C");
            if (centralMoments.Length <= n) throw new InvalidOperationException();


            throw new NotImplementedException();
        }

        public static double ComputeRawMomentFromCulumants (double[] cumulants, int n) {

            if (cumulants == null) throw new ArgumentNullException("cumulants");
            if (n < 0) throw new ArgumentNullException("n");
            if (cumulants.Length <= n) throw new InvalidOperationException();

            double K = 0.0;
            IntegerPartitionEnumerator e = new IntegerPartitionEnumerator(n);
            while (e.MoveNext()) {
                double K1 = 1.0;
                int[] fs = e.Current;
                foreach (int f in fs) {
                    K1 *= cumulants[f];
                }
                K += K1;
            }

            return (K);

        }

        public static double ComputeCentralMomentFromCumulants (double[] cumulants, int n) {
            throw new NotImplementedException();
        }

        public static double ComputeCumulantFromRawMoments (double[] rawMoments, int n) {

            if (rawMoments == null) throw new ArgumentNullException("rawMoments");
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (rawMoments.Length <= n) throw new InvalidOperationException();

            double K = rawMoments[n];
            for (int k = n - 1; k >= 0; k++) {
            }

            return (K);
        }

        public static double CumputeCumulantFromCentralMoments (int mean, double[] centralMoments, int n) {
            throw new NotImplementedException();
        }

    }
#endif

    // Deviates
    // Maximum likelyhood estimation
    // Cumulants

}

