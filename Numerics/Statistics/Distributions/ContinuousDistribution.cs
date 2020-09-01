using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents all continuous, univariate probability distributions.
    /// </summary>
    public abstract class ContinuousDistribution : UnivariateDistribution {

        /// <summary>
        /// Returns the probability density at the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The probability density p(x).</returns>
        /// <remarks>
        /// <para>The probability density function (PDF) gives the relative probabilities of obtaining different values.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Probability_density_function"/>
        public abstract double ProbabilityDensity (double x);

        /// <summary>
        /// Returns the cumulative probability to the left of (below) the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The integrated probability P(x) to obtain a result below the reference point.</returns>
        /// <remarks>
        /// <para>The left probability function is commonly called the cumulative distribution function (CDF).</para>
        /// <para>If you want a right-tailed value, you should use <see cref="RightProbability(double)"/>
        /// instead.</para>
        /// </remarks>
        /// <seealso href="https://en.wikipedia.org/wiki/Cumulative_distribution_function"/>
        public virtual double LeftProbability (double x) {
            if (x <= Support.LeftEndpoint) {
                return 0.0;
            } else if (x >= Support.RightEndpoint) {
                return 1.0;
            } else {
                return FunctionMath.Integrate(ProbabilityDensity, Support.LeftEndpoint, x).Estimate.Value;
            }
        }

        /// <summary>
        /// Returns the cumulative probability to the right of (above) the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The integrated probability Q(x) = 1 - P(x) to obtain a result above the reference point.</returns>
        /// <remarks>
        /// <para>In survival analysis, the right probability function is commonly called the survival function, because it gives the
        /// fraction of the population remaining after the given time.</para>
        /// <para>If you want a right-tailed probability, it is better to call this
        /// method directly instead of computing 1.0 - LeftProbability(x), because the
        /// latter will loose accuracy as P(x) gets close to 1. (In fact, in the far right tail
        /// where Q(x) &lt; 1.0E-16, it will give 0.0, whereas this method is likely to give
        /// an accurate tiny value.)</para>
        /// </remarks>
        /// <seealso href="http://mathworld.wolfram.com/SurvivalFunction.html"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Survival_function"/>
        public virtual double RightProbability (double x) {
            if (x <= Support.LeftEndpoint) {
                return 1.0;
            } else if (x >= Support.RightEndpoint) {
                return 0.0;
            } else {
                return FunctionMath.Integrate(ProbabilityDensity, x, Support.RightEndpoint).Estimate.Value;
            }
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
        /// <seealso href="https://en.wikipedia.org/wiki/Quantile_function"/>
        public virtual double InverseLeftProbability (double P) {
            // find x where LeftProbability(x) = P 
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            Func<double, double> f = (double x) => (LeftProbability(x) - P);
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
            Func<double, double> f = (double x) => (RightProbability(x) - Q);
            double y = FunctionMath.FindZero(f, Mean);
            return (y);
        }

        /// <summary>
        /// Computes the hazard function.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The hazard function p(x)/Q(x).</returns>
        /// <remarks>
        /// <para>Also known as the failure rate or force of mortality.</para>
        /// </remarks>
        /// <seealso href="http://mathworld.wolfram.com/HazardFunction.html"/>
        public virtual double Hazard (double x) {
            return (ProbabilityDensity(x) / RightProbability(x));
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r == 0) {
                return (1.0);
            } else {
                // If there will be no cancelation in the integral, use a pure relative accuracy target.
                // Uh, that's not true, because x can be negative.
                IntegrationSettings settings = new IntegrationSettings();
                if ((r % 2 == 0) || !this.Support.OpenContains(0.0)) settings.AbsolutePrecision = 0.0;
                IntegrationResult result = FunctionMath.Integrate(x => this.ProbabilityDensity(x) * MoreMath.Pow(x, r), this.Support, settings);
                return (result.Estimate.Value);
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r == 0) {
                return 1.0;
            } else if (r == 1) {
                return 0.0;
            } else {
                double mu = Mean;
                return ExpectationValue(x => MoreMath.Pow(x - mu, r));
            }
        }

        /// <summary>
        /// Gets the median of the distribution.
        /// </summary>
        /// <remarks>The median is the point with equal integrated probability above and below, i.e. with P(x1) = 0.5.</remarks>
        public virtual double Median {
            get {
                return InverseLeftProbability(0.5);
            }
        }

        /// <summary>
        /// Gets the interval over which the distribution is non-vanishing.
        /// </summary>
        public virtual Interval Support {
            get {
                return Interval.Infinite;
            }
        }

        /// <summary>
        /// Computes the expectation value of the given function.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <returns>The expectation value of the function.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="f"/> is <see langword="null"/>.</exception>
        public virtual double ExpectationValue (Func<double, double> f) {
            if (f is null) throw new ArgumentNullException(nameof(f));
            return FunctionMath.Integrate(x => f(x) * ProbabilityDensity(x), Support).Estimate.Value;
        }

        /// <summary>
        /// Generates a random variate.
        /// </summary>
        /// <param name="rng">A random number generator.</param>
        /// <returns>A value distributed according to the distribution.</returns>
        /// <remarks>
        /// <para>Note that the random number generator <paramref name="rng"/> will be advanced by this method. The next call to its
        /// generator methods will not give the same value as it would had it not been passed to this method.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="rng"/> is <see langword="null"/>.</exception>
        public virtual double GetRandomValue (Random rng) {
            if (rng is null) throw new ArgumentNullException(nameof(rng));
            return InverseLeftProbability(rng.NextDouble());
        }

        /// <summary>
        /// Generates the given number of random variates.
        /// </summary>
        /// <param name="rng">A random number generator.</param>
        /// <param name="n">The number of variates to generate.</param>
        /// <returns>An iterator that returns the requested number
        /// of variates distributed according to the distribution.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="rng"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is negative.</exception>
        public virtual IEnumerable<double> GetRandomValues (Random rng, int n)
        {
            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
            for (int i = 0; i < n; i++)
            {
                yield return GetRandomValue(rng);
            }

        }

        // compute central moments from raw moments
        // this is subject to loss of precision from cancelation, so be careful

        internal virtual double CentralMomentFromRawMoment (int n) {

            double m = Mean;

            double mm = 1.0;
            double C = RawMoment(n);
            for (int k = 1; k <= n; k++) {
                mm = mm * (-m);
                C += AdvancedIntegerMath.BinomialCoefficient(n, k) * mm * RawMoment(n - k);
            }

            return (C);

        }


    }

}

