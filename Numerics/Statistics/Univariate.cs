using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Statistics {


    /// <summary>
    /// Contains methods for analyzing univariate samples.
    /// </summary>
    public static partial class Univariate {

        /// <summary>
        /// Computes the sample mean.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The sample mean.</returns>
        /// <remarks>
        /// <para>Note this method can be called either as a static method of the <see cref="Univariate"/> class, or as an
        /// extension method of any collection that implements <see cref="IReadOnlyCollection{Double}"/>.</para>
        /// </remarks>
        public static double Mean (this IReadOnlyCollection<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) return (Double.NaN);

            int n = 0;
            double mean = 0.0;
            foreach (double value in sample) {
                n++;
                double delta = value - mean;
                mean += delta / n;
            }

            return (mean);
        }

        /// <summary>
        /// Computes the sample variance.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The sample variance.</returns>
        /// <remarks>
        /// <para>Note that a few authors use the term "sample variance" to mean "the population variance as estimated by the sample". We do not
        /// adopt this aberant convention, which is contrary to all other conventional meanings of sample moments. In other words, what is returned
        /// by this methods is the sum of squared deviations divided by n, not divided by (n-1). If you want an estimate of the population
        /// variance, use the <see cref="PopulationVariance"/> method.</para>
        /// <para>Note this method can be called either as a static method of the <see cref="Univariate"/> class, or as an
        /// extension method of any collection that implements <see cref="IReadOnlyCollection{Double}"/>.</para>
        /// </remarks>
        public static double Variance (this IReadOnlyCollection<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) return (Double.NaN);

            int n;
            double mean, sumOfSquares;
            ComputeMomentsUpToSecond(sample, out n, out mean, out sumOfSquares);
            Debug.Assert(sample.Count == n);

            return (sumOfSquares / n);
        }

        /// <summary>
        /// Computes the sample skewness.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The sample skewness.</returns>
        /// <remarks>
        /// <para>Note this method can be called either as a static method of the <see cref="Univariate"/> class, or as an
        /// extension method of any collection that implements <see cref="IReadOnlyCollection{Double}"/>.</para>
        /// </remarks>
        public static double Skewness (this IReadOnlyCollection<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) return (Double.NaN);

            int n;
            double m1, c2, c3;
            ComputeMomentsUpToThird(sample, out n, out m1, out c2, out c3);
            Debug.Assert(sample.Count == n);
            Debug.Assert(c2 >= 0.0);

            return (c3 / Math.Pow(c2, 3.0 / 2.0));

        }

        // We can compute variance (and higher central moments) in the same pass through the data that we use to compute
        // the mean. This approach is not only faster but also more stable than simpler alternatives. Formulae can be
        // derived by forming the difference between moments of the (n-1)-sized sample and the full n-sized sample.

        // For the mean:
        //   n m' = \sum_{k=1}^n x_k  \qquad  (n-1) m = \sum_{k=1}^{n-1} x_k
        //   n m' = (n-1) m + x_n
        //     m' = m + \frac{x_n - m}{n} = m + d
        // For the variance:
        //   n v' = \sum_{k=1}^n ( x_k - m')^2  \qquad (n-1) v = \sum_{k=1}^{n-1} (x_k - m)^2  
        //   n v' = \sum_{k=1}^n ( x_k - m - d)^2 = \sum_{k=1}^n [ ( x_k - m )^2 - 2 ( x_k - m ) d + d^2 ]
        //   n v' = (n-1) v + (n d)^2 - 0 + n d^2 = (n-1) v + n (n + 1) d^2 = (n-1) v + e
        // These formulas, for central moments up to 4, are shown in https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics

        // The second moment update formula can also be written using the difference between x_n and both the pre-update and post-update mean.
        //   (x_n - m) (x_n - m') = (x_n - m - d)(x_n - m') = (d n - d) d n = n (n - 1) d^2 = e

        internal static void ComputeMomentsUpToSecond (IEnumerable<double> sample, out int n, out double mean, out double sumOfSquaredDeviations) {

            // In some cases (e.g. t-tests) it's helpful to access the sum of squared deviations directly. If we were to return
            //   c_2 = sumOfSquaredDeviations / n
            // these callers would need to multiply by N just to get back the number they want. Therefore we return sumOfSquaredDeviations.

            Debug.Assert(sample != null);

            n = 0;
            mean = 0.0;
            sumOfSquaredDeviations = 0.0;
            foreach (double value in sample) {
                n++;
                double delta = value - mean;
                mean += delta / n;
                sumOfSquaredDeviations += (value - mean) * delta;
            }

            Debug.Assert(n > 0);
            Debug.Assert(sumOfSquaredDeviations >= 0.0);
        }

        private static void ComputeMomentsUpToThird (IEnumerable<double> values, out int n, out double m1, out double c2, out double c3) {

            Debug.Assert(values != null);

            n = 0;
            m1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;

            foreach (double x in values) {
                double c2_old = c2;

                n++;
                double d = (x - m1) / n;
                m1 += d;
                double e = d * d * n * (n - 1);
                c2 += e;
                c3 += d * (e * (n - 2) - 3.0 * c2_old);
            }

            c2 /= n;
            c3 /= n;

            Debug.Assert(n > 0);
            Debug.Assert(c2 >= 0.0);
        }

        private static void ComputeMomentsUpToFourth(IEnumerable<double> values, out int n, out double m1, out double c2, out double c3, out double c4) {

            Debug.Assert(values != null);

            n = 0;
            m1 = 0.0;
            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;

            foreach (double x in values) {
                double c2_old = c2;
                double c3_old = c3;

                n++;
                double d = (x - m1) / n;
                m1 += d;
                double e = d * d * n * (n - 1);
                c2 += e;
                c3 += d * (e * (n - 2) - 3.0 * c2_old);
                c4 += d * (d * (e * (n * n - 3 * n + 3) + 6.0 * c2_old) - 4.0 * c3_old);
            }

            c2 /= n;
            c3 /= n;
            c4 /= n;

            Debug.Assert(n > 0);
            Debug.Assert(c2 >= 0.0);
            Debug.Assert(c4 >= 0.0);

        }

        /// <summary>
        /// Computes the sample standard deviation.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The sample standard deviation.</returns>
        /// <remarks>
        /// <para>Note that a few authors use the term "sample standard deviation" to mean "the square root of the population variance as estimated by the sample".
        /// We do not adopt this aberant convention, which is contrary to all other conventional meanings of sample moments. In other words, what is returned
        /// by this methods is the square root of the sum of squared deviations divided by n, not divided by (n-1). If you want an estimate of the population
        /// standard deviation, use the <see cref="PopulationStandardDeviation"/> method.</para>
        /// <para>This method can be called either as a static method of the <see cref="Univariate"/> class, or as an
        /// extension method of any collection that implements <see cref="IReadOnlyCollection{Double}"/>.</para>
        /// </remarks>
        public static double StandardDeviation (this IReadOnlyCollection<double> sample) {
            return (Math.Sqrt(Variance(sample)));
        }

        /// <summary>
        /// Computes the given sample raw moment.
        /// </summary>
        /// <param name="sample">The sample containing the data.</param>
        /// <param name="r">The order of the moment to compute.</param>
        /// <returns>The <paramref name="r"/>th raw moment of the sample.</returns>
        public static double RawMoment (this IReadOnlyCollection<double> sample, int r) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) throw new InsufficientDataException();
            if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (Mean(sample));
            } else {
                int n = 0;
                double M = 0.0;
                foreach (double value in sample) {
                    n++;
                    M += MoreMath.Pow(value, r);
                }
                return (M / n);
            }
        }

        /// <summary>
        /// Computes the given sample central moment.
        /// </summary>
        /// <param name="sample">The sample containing the data.</param>
        /// <param name="r">The order of the moment to compute.</param>
        /// <returns>The <paramref name="r"/>th central moment of the sample.</returns>
        /// <remarks>
        /// <para>This method computes the central momements of the sample data, not the estimated
        /// central moments of the underlying population; to obtain the latter, use <see cref="PopulationCentralMoment(IReadOnlyCollection{double},int)"/>.</para>
        /// </remarks>
        public static double CentralMoment(this IReadOnlyCollection<double> sample, int r) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) throw new InsufficientDataException();
            if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                double mean = Mean(sample);
                return (CentralMoment(sample, r, mean));
            }
        }

        private static double CentralMoment(IEnumerable<double> sample, int r, double mean) {
            Debug.Assert(sample != null);
            int n = 0;
            double C = 0.0;
            foreach (double value in sample) {
                n++;
                C += MoreMath.Pow(value - mean, r);
            }
            return (C / n);
        }

    }
}
