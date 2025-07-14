using Meta.Numerics.Data;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Contains methods for analyzing univariate samples.
    /// </summary>
    /// <remarks>
    /// <para>This is the central class for the analysis of univariate samples of independent, identically
    /// distributed values.</para>
    /// <para>To compute moments of the sample data, you can use methods such as
    /// <see cref="Mean(IReadOnlyCollection{double})">Mean</see>,
    /// <see cref="Variance(IReadOnlyCollection{double})">Variance</see>,
    /// <see cref="RawMoment(IReadOnlyCollection{double}, int)">RawMoment</see>, and
    /// <see cref="CentralMoment(IReadOnlyCollection{double}, int)">CentralMoment</see>.
    /// Note that these are moments of the sample data, not estimates of the moments of the underlying
    /// population from which the sample was drawn.</para>
    /// <para>To obtain estimates of the moments of the underlying population from which the sample was
    /// drawn, you can use methods such as
    /// <see cref="PopulationMean(IReadOnlyCollection{double})">PopulationMean</see>,
    /// <see cref="PopulationStandardDeviation(IReadOnlyCollection{double})">PopulationStandardDeviation</see>,
    /// <see cref="PopulationRawMoment(IReadOnlyCollection{double}, int)">PopulationRawMoment</see>,
    /// and <see cref="PopulationCentralMoment(IReadOnlyCollection{double}, int)">PopulationCentralMoment</see>.
    /// These estimates all come with associated error bars, so they return <see cref="UncertainValue"/> structures.</para>
    /// <para>You can fit a sample to any number of distributions using methods such as
    /// <see cref="FitToExponential(IReadOnlyList{double})"/>,
    /// <see cref="FitToLognormal(IReadOnlyList{double})"/>,
    /// <see cref="FitToNormal(IReadOnlyList{double})"/>,
    /// <see cref="FitToWeibull(IReadOnlyList{double})"/>.
    /// </para>
    /// <para>You can perform statistical tests on a single sample,
    /// such as the <see cref="StudentTTest(IReadOnlyCollection{double}, double)"/> or
    /// <see cref="SignTest(IReadOnlyCollection{double}, double)"/>
    /// to compare a sample to a reference value, or the <see cref="ShapiroFranciaTest(IReadOnlyList{double})"/> to test
    /// a sample for normality. You can also perform statistical tests comparing multiple samples,
    /// such as the two-sample <see cref="StudentTTest(IReadOnlyCollection{double}, IReadOnlyCollection{double})"/>,
    /// and <see cref="MannWhitneyTest(IReadOnlyList{double}, IReadOnlyList{double})" />, or the multi-sample
    /// <see cref="OneWayAnovaTest(IReadOnlyCollection{double}[])"/> and
    /// <see cref="KruskalWallisTest(IReadOnlyList{double}[])"/>.</para>
    /// <para>Most of the methods in this class are extension methods that accept as a sample any type
    /// that implements the appropriate collection interface. So, for example, given
    /// sample values in a <see cref="List{T}">List&lt;string&gt;</see> named <c>s</c>, you could estimate
    /// the variance of the population from which it was drawn either by <c>s.PopulationVariance()</c>
    /// or <c>Univariate.PopulationVariance(s)</c>.</para>
    /// </remarks>
    public static partial class Univariate {

        /// <summary>
        /// Computes the sample mean.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The mean of the sample.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        public static double Mean (this IReadOnlyCollection<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) return (Double.NaN);

            ComputeMomentsUpToFirst(sample, out int n, out double mean);

            return (mean);
        }

        /// <summary>
        /// Computes the sample variance.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The variance of the sample data.</returns>
        /// <remarks>
        /// <para>Note that some  authors use the term "sample variance" to mean "the population variance as estimated by the sample". We do not
        /// adopt this aberrant convention, which is contrary to the conventional meanings of all other sample moments. This means that
        /// what is returned by this method is the sum of squared deviations divided by n, not divided by (n-1). If you want an estimate
        /// of the population variance, use the <see cref="PopulationVariance"/> method.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        public static double Variance (this IReadOnlyCollection<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) return (Double.NaN);

            ComputeMomentsUpToSecond(sample, out int n, out double mean, out double sumOfSquares);
            Debug.Assert(sample.Count == n);

            return (sumOfSquares / n);
        }

        /// <summary>
        /// Computes the sample skewness.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The skewness of the sample data.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Skewness"/>
        public static double Skewness (this IReadOnlyCollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) return (Double.NaN);

            ComputeMomentsUpToThird(sample, out int n, out double m1, out double c2, out double c3);
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

        // This method of computing means & variances online is called the Welford algorithm.


        internal static void ComputeMomentsUpToFirst(IEnumerable<double> sample, out int n, out double mean) {

            Debug.Assert(sample != null);

            n = 0;
            mean = 0.0;

            foreach (double value in sample) {
                n++;
                mean += (value - mean) / n;
            }

        }

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

        // Two-pass calculation are the most accurate. The second pass not only computes the moments, but
        // by tracking how far off the given mean is from the real mean, produces a corrected mean and corrections for the higher moments.

        private static void FirstPassMean(IEnumerable<double> values, out int n, out double m) {
            Debug.Assert(values != null);
            n = 0;
            m = 0.0;
            foreach (double x in values) {
                n++;
                m += (x - m) / n;
            }
        }

        // Suppose we have computed central moments based on a value of the mean m (which will we suppose to be slightly wrong).
        //   c_r = n^{-1} \sum_i (x_i - m)^r
        // We will then find
        //   c_1 = ( n^{-1} \sum_i x_i - m ) \ne 0 = \delta m
        // is a correction to the mean. The true mean is m + \delta m.

        // We can then use Taylor expansion to get corrections to all the other moments.
        //   \frac{\partial c_r}{\partial m} = -r n^{-1} \sum_i (x_i - m)^{r-1} = -r c_{r-1}
        //   \frac{\partial^2 c_r}{\partial m^2} = r(r-1) n^{-1} \sum_i (x_i - m)^{r-2} = r(r-1) c_{r-2}
        // So for any r > 3, to 2nd order in (\delta m),
        //   \delta c_r = -r c_{r-1}  (\delta m) + 1/2 r(r-1) c_{r-2} (\delta m)^2 = (\delta m) [ 1/2 r (r-1) c_{r-2} (\delta m) - r c_{r-1} ]

        // For r = 2 or 3, c_1 = (\delta m) appears, which messes with the relative order of terms
        //   \delta c_2 = -2 c_1 (\delta m) + c_0 (\delta m)^2 = -2 (\delta m)^2 + (\delta m)^2 = - (\delta m)^2
        //   \delta c_3 = -3 c_2 (\delta m) + 3 c_1 (\delta m)^2 - 1 c_0 (\delta m)^3 = -3 c_2 (\delta m) + 3 (\delta m)^3 - (\delta m)^3 = (\delta m) [ 2 (\delta m)^2 - 3 c_2 ]

        private static void CorrectionPassMoments2 (IEnumerable<double> values, ref double m, out double c2) {

            int n = 0;

            // Compute c2 assuming the original m
            c2 = 0.0;
            double delta_m = 0.0;
            foreach (double x in values) {
                n++;
                double z = x - m;
                delta_m += z;
                c2 += z * z;
            }
            delta_m /= n;
            c2 /= n;

            // Correct c2 and m
            c2 -= MoreMath.Sqr(delta_m);
            m += delta_m;
        }


        private static void CorrectionPassMoments4 (IEnumerable<double> xs, ref double m, out double c2, out double c3, out double c4) {

            int n = 0;

            // Compute moments assuming the original m
            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;
            double delta_m = 0.0;
            foreach (double x in xs) {
                n++;
                double z = x - m;
                double z2 = z * z;
                delta_m += z;
                c2 += z2;
                c3 += z * z2;
                c4 += z2 * z2;
            }
            delta_m /= n;
            c2 /= n;
            c3 /= n;
            c4 /= n;

            // Add corrections based on computed delta_m
            // They must be ordered higher to lower moments, since original lower moments appear in corrections to higher
            c4 += delta_m * (6.0 * c2 * delta_m - 4.0 * c3);
            c3 += delta_m * (2.0 * MoreMath.Sqr(delta_m) - 3.0 * c2);
            c2 -= MoreMath.Sqr(delta_m);
            m += delta_m;

        }

        private static void CorrectionPassMoments6(IEnumerable<double> xs, ref double m, out double c2, out double c3, out double c4, out double c5, out double c6) {

            int n = 0;

            // Compute moments assuming the original m
            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;
            c5 = 0.0;
            c6 = 0.0;
            double delta_m = 0.0;
            foreach (double x in xs) {
                n++;
                double z = x - m;
                double z2 = z * z;
                double z3 = z * z2;
                delta_m += z;
                c2 += z2;
                c3 += z3;
                c4 += z2 * z2;
                c5 += z2 * z3;
                c6 += z3 * z3;
            }
            delta_m /= n;
            c2 /= n;
            c3 /= n;
            c4 /= n;
            c5 /= n;
            c6 /= n;

            // Add corrections based on computed delta_m
            // They must be ordered higher to lower moments, since original lower moments appear in corrections to higher
            c6 += delta_m * (15.0 * c4 * delta_m - 6.0 * c5);
            c5 += delta_m * (10.0 * c3 * delta_m - 5.0 * c4);
            c4 += delta_m * (6.0 * c2 * delta_m - 4.0 * c3);
            c3 += delta_m * (2.0 * MoreMath.Sqr(delta_m) - 3.0 * c2);
            c2 -= MoreMath.Sqr(delta_m);
            m += delta_m;

        }

        /// <summary>
        /// Computes the sample standard deviation.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The standard deviation of the sample data.</returns>
        /// <remarks>
        /// <para>Note that a few authors use the term "sample standard deviation" to mean "the square root of the population variance as estimated by the sample".
        /// We do not adopt this aberrant convention, which is contrary to the conventional meaning of "sample moment" for all other moments.
        /// In other words, what is returned by this methods is the square root of the sum of the squared deviations divided by n, not divided by (n-1).
        /// If you want the square root of the sum of the squared deviations divided by n, which is also called the Bessel-corrected standard deviation,
        /// use the <see cref="CorrectedStandardDeviation(IReadOnlyCollection{double})"/> method. Note, however, that, unlike the variance, the standard deviation
        /// corrected in this way is not unbiased.
        /// If you want a bias-corrected estimate of the standard deviation of the populace that was sampled, use the <see cref="PopulationStandardDeviation"/> method.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        public static double StandardDeviation (this IReadOnlyCollection<double> sample) {
            return (Math.Sqrt(Variance(sample)));
        }

        /// <summary>
        /// Computes the Bessel-corrected standard deviation.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The Bessel-corrected standard deviation.</returns>
        /// <remarks>
        /// <para>This probably isn't the quantity you want, even though it's what many software
        /// packages refer to as the "standard deviation", and some authors even misleadingly call the "sample standard deviation".
        /// It is the square root of the sum of the squared deviations from the mean, divided by n-1 instead of n.</para>
        /// <para>Using n-1 instead of n in the formula for <em>variance</em> produces an
        /// unbiased estimator of the <see cref="PopulationVariance"/>. But using n-1 instead of n
        /// in the formula for standard deviation, which is how this property is
        /// computed, does <em>not</em> produce an unbiased estimator of the population's standard deviation.
        /// Our implementation of <see cref="PopulationStandardDeviation"/> does a better job of reducing
        /// bias, so you should use it instead if that is what you are trying to estimate. The main
        /// reason this property exists at all is to satisfy users who want to compute the exact
        /// same value that another software package did.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <seealso cref="StandardDeviation"/>
        /// <seealso cref="PopulationStandardDeviation"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Bessel%27s_correction"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation"/>
        /// <seealso href="https://en.wikipedia.org/wiki/Bias_of_an_estimator"/>
        public static double CorrectedStandardDeviation (this IReadOnlyCollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) return (Double.NaN);

            ComputeMomentsUpToSecond(sample, out int n, out double mean, out double sumOfSquares);
            Debug.Assert(sample.Count == n);

            return (Math.Sqrt(sumOfSquares / (n-1)));
        }

        /// <summary>
        /// Computes the given sample raw moment.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="r">The order of the moment to compute.</param>
        /// <returns>The <paramref name="r"/>th raw moment of the sample data.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        public static double RawMoment (this IReadOnlyCollection<double> sample, int r) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) throw new InsufficientDataException();
            if (r == 0) {
                return 1.0;
            } else if (r == 1) {
                return Mean(sample);
            } else {
                int n = 0;
                double M = 0.0;
                foreach (double value in sample) {
                    n++;
                    M += MoreMath.Pow(value, r);
                }
                return M / n;
            }
        }

        /// <summary>
        /// Computes the given sample central moment.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="r">The order of the moment to compute.</param>
        /// <returns>The <paramref name="r"/>th central moment of the sample.</returns>
        /// <remarks>
        /// <para>This method computes the central moments of the sample data, not the estimated
        /// central moments of the underlying population; to obtain the latter, use <see cref="PopulationCentralMoment(IReadOnlyCollection{double},int)"/>.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
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
