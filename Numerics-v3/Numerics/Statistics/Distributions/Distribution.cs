using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents all continuous, univariate probability distribution.
    /// </summary>
    public abstract class Distribution : UnivariateDistribution {

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
            return (1.0 - LeftProbability(x));
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
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            Func<double, double> f = delegate(double x) {
                return (LeftProbability(x) - P);
            };
            double y = FunctionMath.FindZero(f, Mean);
            return (y);
            // since we have the PDF = CDF', change to a method using Newton's method
        }

        /// <summary>
        /// Returns the point at which the right probability function attains the given value.
        /// </summary>
        /// <param name="Q">The right cumulative probability, which must lie between 0 and 1.</param>
        /// <returns>The value x for which <see cref="RightProbability"/> equals Q.</returns>
        public virtual double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            Func<double, double> f = delegate(double x) {
                return (RightProbability(x) - Q);
            };
            double y = FunctionMath.FindZero(f, Mean);
            return (y);
        }

        /// <inheritdoc />
        public override double Moment (int r) {
            if (r == 0) {
                return (1.0);
            } else {
                return (ExpectationValue(x => MoreMath.Pow(x, r)));
            }
        }

        /// <inheritdoc />
        public override double MomentAboutMean (int r) {
            if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                double mu = Mean;
                return (ExpectationValue(x => MoreMath.Pow(x - mu, r)));
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
            if (f == null) throw new ArgumentNullException(nameof(f));
            return (FunctionMath.Integrate(x => f(x) * ProbabilityDensity(x), Support));
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
        /// <exception cref="ArgumentNullException"><paramref name="rng"/> is null.</exception>
        public virtual double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException(nameof(rng));
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
                C += AdvancedIntegerMath.BinomialCoefficient(n, k) * mm * Moment(n - k);
            }

            return (C);

        }


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

