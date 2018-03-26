using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Statistics {

    public static partial class Univariate {

        /// <summary>
        /// Estimates the mean of the underlying population.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>An estimate, with uncertainty, of the mean of the distribution from which the sample was drawn.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than 2 values.</exception>
        /// <remarks>
        /// <para>In contrast to the <see cref="Mean(IReadOnlyCollection{double})"/> method, this method estimates
        /// the mean of the underlying population from which the sample was drawn, and provides a error estimate
        /// for that value.</para>
        /// </remarks>
        public static UncertainValue PopulationMean (this IReadOnlyCollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();
            return (EstimateFirstCumulant(sample));
        }

        /// <summary>
        /// Estimates of the variance of the underlying population.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>An estimate, with uncertainty, of the variance of the distribution from which the sample was drawn.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than 3 values.</exception>
        /// <remarks>
        /// <para>In contrast to the <see cref="Variance(IReadOnlyCollection{double})"/> method, this method estimates
        /// the variance of the underlying population from which the sample was drawn, and provides a error estimate
        /// for that value.</para>
        /// </remarks>
        public static UncertainValue PopulationVariance (this IReadOnlyCollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();
            return (EstimateSecondCumulant(sample));
        }

        /// <summary>
        /// Estimates of the standard deviation of the underlying population.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>An estimate, with uncertainty, of the standard deviation of the distribution from which the sample was drawn.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than 3 values.</exception>
        /// <remarks>
        /// <para>In contrast to the <see cref="StandardDeviation(IReadOnlyCollection{double})"/> method, this method
        /// estimates the standard deviation of the underlying population from which the sample was drawn, and provides
        /// an error estimate for that value.</para>
        /// <para>Note that the value returned by this method is not exactly the square root
        /// of the value returned by <see cref="PopulationVariance"/> (nor is it exactly the same as
        /// the value of <see cref="CorrectedStandardDeviation(IReadOnlyCollection{double})"/>, which is
        /// exactly the square root of the variance best estimate). This is not an error. The variance
        /// estimator has a distribution. The mean of that distribution is equal to the variance of
        /// the underlying population, which is what makes the estimator unbiased. An unbiased estimator
        /// of the standard deviation will have that same property -- its mean will equal the standard
        /// deviation of the underlying population. But the square root of the mean of a distributed
        /// quantity is not the mean of the square root of that quantity, so the square root of the
        /// unbiased estimator of the variance will not an be unbiased estimator of the standard deviation.
        /// This method estimates that bias and corrects for it. The estimation is not exact, so the
        /// value returned by this method will not be perfectly unbiased, but it is likely to be less
        /// biased than <see cref="CorrectedStandardDeviation(IReadOnlyCollection{double})"/>.</para>
        /// </remarks>
        public static UncertainValue PopulationStandardDeviation (this IReadOnlyCollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();
            return (EstimateStandardDeviation(sample));
        }


        /// <summary>
        /// Estimates the given raw moment of the underlying population.
        /// </summary>
        /// <param name="sample">The sample data.</param>
        /// <param name="r">The order of the moment.</param>
        /// <returns>An estimate, with uncertainty, of the <paramref name="r"/>th raw moment of the
        /// distribution from which the sample was drawn.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than 2 values.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> is negative.</exception>
        public static UncertainValue PopulationRawMoment (this IReadOnlyCollection<double> sample, int r) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();
            if (r < 0) {
                // we don't do negative moments
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                // the zeroth moment is exactly one for any distribution
                return (new UncertainValue(1.0, 0.0));
            } else if (r == 1) {
                // the first moment is just the mean
                return (EstimateFirstCumulant(sample));
            } else {
                // moments of order two and higher
                int n = sample.Count;
                double M_r = sample.RawMoment(r);
                double M_2r = sample.RawMoment(2 * r);
                return (new UncertainValue(M_r, Math.Sqrt((M_2r - M_r * M_r) / n)));
            }
        }

        /// <summary>
        /// Estimates the given central moment of the underlying population.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="r">The order of the moment.</param>
        /// <returns>An estimate, with uncertainty, of the <paramref name="r"/>th central moment of the distribution
        /// from which the sample was drawn.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than 3 values.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="r"/> is negative.</exception>
        public static UncertainValue PopulationCentralMoment (this IReadOnlyCollection<double> sample, int r) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();
            if (r < 0) {
                // we don't do negative moments
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                // the zeroth moment is exactly one for any distribution
                return (new UncertainValue(1.0, 0.0));
            } else if (r == 1) {
                // the first moment about the mean is exactly zero by the definition of the mean
                return (new UncertainValue(0.0, 0.0));
            } else if (r == 2) {
                return (EstimateSecondCumulant(sample));
            } else if (r == 3) {
                return (EstimateThirdCumulant(sample));
            } else {
                return (EstimateCentralMoment(sample, r));
            }
        }

        // First three cumulants are M1, C2, and C3. Unbiased estimators of these are called k-statistics.
        // See http://mathworld.wolfram.com/k-Statistic.html for k-statistic formulas in terms of sample central
        // moments and expressions for their variances in terms of population cumulants.

        // These formulas are also found in Dodge & Rousson, "The Complications of the Fourth Central Moment",
        // The American Statistician 53 (1999) 267. Also Stuart & Ord, "Kendall's Advanced Theory of Statistics",
        // Volume I, Chapter 12

        // The first cumulant, the mean:
        //   k_1 = m_1
        //   V(k_1) = \frac{K_2}{n}
        // To get V(k_1), we need K_2. We have an unbiased estimate k_2 = \frac{n}{n-1} c_2, and if we plug that in, we get
        //   V(k_1) = \frac{c_2}{n-1}
        // None of this is surprising. We already knew that the m_1 is an unbiased estimator of M_1, with variance C_2 / n.

        private static UncertainValue EstimateFirstCumulant (IEnumerable<double> sample) {
            Debug.Assert(sample != null);

            int n;
            double mean, sumOfSquares;
            ComputeMomentsUpToSecond(sample, out n, out mean, out sumOfSquares);

            double m1 = mean;
            double c2 = sumOfSquares / n;

            return (new UncertainValue(m1, Math.Sqrt(c2 / (n - 1))));
        }

        // The second cumulant, the variance:
        //   k_2 = \frac{n}{n-1} c_2
        //   V(k_2) = \frac{K_4}{n} + \frac{2 K_2^2}{n-1}
        // Evaluating V(k_2) is not trivial. Just as we used k_2 for K_2 above, we could just use k_2 and k_4 for K_2 and K_4.
        // But that isn't exactly right, because an even though k_2 is an unbiased estimator of K_2, k_2^2 isn't an unbiased
        // estimator of K_2^2. Using polykays, Stuart & Ord derive an unbiased estimator (formula 12.131), which is also
        // quoted by Mathworld:
        //   V(k_2) = \frac{1}{n+1} \left[ \frac{(n-1)}{n} k_4 + 2 k_2^2 \right]
        // Plugging in the formulas for k-statistics in terms of central moments and doing a lot of algebra gets us to
        //   V(k_2) = \frac{n}{(n-2)(n-3)} \left[ c_4 - \frac{(n^2 - 3)}{(n - 1)^2} c_2^2 \right]
        // Numerical simulations show this to indeed be unbiased, but it has a problem. The coefficient of c_2^2 is
        // greater than unity, so the quantity in brackets can be negative. Therefore we fall back on the asymptotic
        // limit
        //   V(k_2) = \frac{c_4 - c_2^2}{n}
        // which cannot be negative because of the inequality (verified in numerical simulations) c_4 > c_2^2. 

        // See if we can do better, either by adjusting denominator and/or making exact for normal distributions.
        // Being a little off here isn't terrible, since this is just an error estimate.

        // For normal, V(k_2) = \frac{2 \sigma^4}{n-1}, so change denominator to n-1?

        private static UncertainValue EstimateSecondCumulant (IEnumerable<double> sample) {
            Debug.Assert(sample != null);

            int n;
            double m1, c2, c3, c4;
            ComputeMomentsUpToFourth(sample, out n, out m1, out c2, out c3, out c4);
            Debug.Assert(c2 >= 0.0);
            Debug.Assert(c4 >= 0.0);
            /*
            double m1p = Mean(sample);
            double c2p = CentralMoment(sample, 2);
            double c3p = CentralMoment(sample, 3);
            double c4p = CentralMoment(sample, 4);
            c4 = c4p;
            c2 = c2p;
            */
            double k2 = c2 * n / (n - 1);
            double v = c4 - c2 * c2;
            Debug.Assert(v >= 0.0);

            return (new UncertainValue(k2, Math.Sqrt(v / (n - 1))));
        }

        // The third cumulant, which characterizes skew:
        //   k_3 = \frac{n^2}{(n-1)(n-2)} c_3
        //   V(k_3) = \frac{K_6}{n} + \frac{9 K_4 K_2}{n-1} + \frac{9 K_3^2}{n-1} + \frac{6 n K_2^3}{(n-1)(n-2)}    
        // Evaluating V(k_3) has the same issues as V(k_2) above. In this case, we won't even attempt an unbiased
        // estimator, but instead just go directly to the asymptotic form. In terms of central moments, this is:
        //   V(k_3) = \frac{c_6 - 6 c_4 c_2 - c_3^2 + 9 c_2^3}{n}
        // Numerical simulations suggest this is guaranteed positive, but it would be nice to find a proof.

        private static UncertainValue EstimateThirdCumulant (IReadOnlyCollection<double> sample) {
            Debug.Assert(sample != null);

            int n;
            double mean, sumOfSquares;
            ComputeMomentsUpToSecond(sample, out n, out mean, out sumOfSquares);
            double c2 = sumOfSquares / n;
            double c3 = CentralMoment(sample, 3, mean);
            double c4 = CentralMoment(sample, 4, mean);
            double c6 = CentralMoment(sample, 6, mean);

            double k3 = c3 * n * n / (n - 1) / (n - 2);
            double v = c6 - c3 * c3 - 3.0 * (2.0 * c4 - 3.0 * c2 * c2) * c2;
            Debug.Assert(v >= 0.0);
            return (new UncertainValue(k3, Math.Sqrt(v / n)));
        }

        // An obvious estimate of the standard deviation is just the square root of the estimated variance, computed
        // with our usual error propagation methods. This is fine as a first approximation, and is all
        // we did in previous versions, but we can do better.

        // The reason it's not exactly right is that the square root of the mean of a distributed quantity
        // is not the mean of the square root of the quantity. This means that, even though we have found
        // in k_2 an unbiased estimator of K_2 = C_2, its square root is not an unbiased estimator of
        // \sqrt{C_2}. To discover its bias, use Taylor expansion
        //   E(f(x)) = \int \! dp \, p(x) \, f(x)
        //           = \int \! dp \, p(x) \left[ f(x_0) + f'(x_0) (x - x_0) + \half f''(x_0) (x - x_0)^2 + \cdots \right]
        //           = f(x_0) \int \! dp \, p(x) + f'(x_0) \int \! dp \, p(x) (x - x_0) + \half f''(x_0) \int \! dp \, p(x) (x - x_0)^2 + \cdots
        //           = f(x_0) + 0 + \half f''(x_0) V(x) + \cdots
        //  In our case f(x) = \sqrt{k_2} so
        //    E(\sqrt{k_2}) = C_2^{1/2} - \frac{1}{8} C_2^{-3/2} V(k_2) + \cdots
        //  So a less biased estimator of \sqrt{C_2} is
        //    \hat{\sigma} = k_2^{1/2} + \frac{1}{8} k_2^{-3/2} V(k_2)
        //                 = \sqrt{k_2} \left[ 1 + \frac{c_4 - c_2^2}{8 n c_2^2} 
        // This is covered in Stuart & Ord in exercise 10.20.

        private static UncertainValue EstimateStandardDeviation (IEnumerable<double> sample) {

            Debug.Assert(sample != null);

            int n;
            double m1, c2, c3, c4;
            ComputeMomentsUpToFourth(sample, out n, out m1, out c2, out c3, out c4);

            double v = (c4 - c2 * c2) / (n - 1);
            double s = Math.Sqrt(c2 * n / (n - 1)) * (1.0 + v / (8.0 * c2 * c2));
            double ds = 0.5 * Math.Sqrt(v) / s;

            return (new UncertainValue(s, ds));
        }

        // Stuart & Ord, "Kendall's Advanced Theory of Statistics", Volume I, Chapter 10 derives
        // formulas for the the expectation and variance of arbitrary c_r in terms of C_r as a series in 1/n.
        //   E(c_r) = C_r - \frac{r C_r - r(r-1) / 2 C_{r-2} C_2}{n} + \cdots
        //   V(c_r) = \frac{C_{2r} - C_{r}^2 + r^2 C_{r-1}^2 C_2 - 2 r C_{r+1} C_{r-1}}{n} + \cdots
        // We can turn these around to get an estimator for C_r and its variance in terms of c_r.
        // We have verified that these agree, to these orders, with the exact expressions for M1, C2, and C3 above.

        private static UncertainValue EstimateCentralMoment (IReadOnlyCollection<double> sample, int r) {
            Debug.Assert(sample != null);
            Debug.Assert(r >= 2);

            // The formula for the best estimate of the r'th population moment involves C_r, C_{r-2}, and C_2.
            // The formula for the uncertainty in this estimate involves C_{2r}, C_{r+1}, C_{r-1}, and C_2.
            // We need all these sample moments.
            int n;
            double m, sumOfSquaredDeviations;
            ComputeMomentsUpToSecond(sample, out n, out m, out sumOfSquaredDeviations);
            Debug.Assert(sumOfSquaredDeviations >= 0);

            double c_2 = sumOfSquaredDeviations / n;
            double c_rm2 = CentralMoment(sample, r - 2, m);
            double c_rm1 = CentralMoment(sample, r - 1, m);
            double c_r = CentralMoment(sample, r, m);
            double c_rp1 = CentralMoment(sample, r + 1, m);
            double c_2r = CentralMoment(sample, 2 * r, m);

            double c = c_r + r * (c_r - 0.5 * (r - 1) * c_rm2 * c_2) / n;
            double v = c_2r - c_r * c_r + r * r * c_rm1 * c_rm1 * c_2 - 2 * r * c_rp1 * c_rm1;
            Debug.Assert(v >= 0.0);

            return (new UncertainValue(c, Math.Sqrt(v / n)));
        }

    }
}
