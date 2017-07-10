using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Statistics {


    /// <summary>
    /// Contains methods for analyzing univariate samples.
    /// </summary>
    public static class Univariate {

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
            if (sample.Count < 1) throw new InsufficientDataException();
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
            double mean, sumOfSquares;
            ComputeMomentsUpToSecond(sample, out mean, out sumOfSquares);
            return (sumOfSquares / sample.Count);
        }

        private static void ComputeMomentsUpToSecond (IReadOnlyCollection<double> sample, out double mean, out double sumOfSquares) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();
            int n = 0;
            mean = 0.0;
            sumOfSquares = 0.0;
            foreach (double value in sample) {
                n++;
                double delta = value - mean;
                mean += delta / n;
                sumOfSquares += (value - mean) * delta;
            }
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
        /// <para>Note this method can be called either as a static method of the <see cref="Univariate"/> class, or as an
        /// extension method of any collection that implements <see cref="IReadOnlyCollection{Double}"/>.</para>
        /// </remarks>
        public static double StandardDeviation (this IReadOnlyCollection<double> sample) {
            return (Variance(sample));
        }

        /// <summary>
        /// Estimates the population mean.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>An estimate, with uncertainty, of the popoulation mean.</returns>
        /// <remarks>
        /// <para>Note this method can be called either as a static method of the <see cref="Univariate"/> class, or as an
        /// extension method of any collection that implements <see cref="IReadOnlyCollection{Double}"/>.</para>
        /// </remarks>
        public static UncertainValue PopulationMean (this IReadOnlyCollection<double> sample) {
            double mean, sumOfSquares;
            ComputeMomentsUpToSecond(sample, out mean, out sumOfSquares);
            return (new UncertainValue(mean, Math.Sqrt(sumOfSquares / sample.Count)));
        }

        public static UncertainValue PopulationVariance (this IReadOnlyCollection<double> sample) {
            throw new NotImplementedException();
        }

        public static UncertainValue PopulationStandardDeviation (this IReadOnlyCollection<double> sample) {
            throw new NotImplementedException();
        }

        public static double Minimum(this ICollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            double minimum = Double.NaN;
            foreach (double value in sample) {
                if (!(value > minimum)) minimum = value;
            }
            return (minimum);
        }

        public static double Maximum(this ICollection<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            double maximum = Double.NaN;
            foreach (double value in sample) {
                if (!(value < maximum)) maximum = value;
            }
            return (maximum);
        }

        /// <summary>
        /// Computes the median of the sample.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <returns>The median of the sample.</returns>
        public static double Median(this IReadOnlyList<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();

            List<int> order = new List<int>(sample.Count);
            for (int i = 0; i < sample.Count; i++) {
                order.Add(i);
            }

            order.Sort((int i, int j) => sample[i].CompareTo(sample[j]));

            int m = sample.Count / 2;
            if (sample.Count % 2 == 0) {
                return (sample[order[m]]);
            } else {
                return (0.5 * (sample[m] + sample[m + 1]));
            }

        }

    }
}
