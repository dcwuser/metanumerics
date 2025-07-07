using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;

namespace Meta.Numerics.Statistics {

    


    /// <summary>
    /// Tracks summary statistics for a stream of data points.
    /// </summary>
    /// <remarks>
    /// <para>Unlike the methods of the <see cref="Univariate"/> class, which by and large expect
    /// to work with sample data that is stored entirely within memory, this class is designed
    /// to compute the most important summary statistics on streamed data, i.e. individual
    /// data points which need only be presented once and may afterwards be discarded.</para>
    /// <para>The <see cref="Combine(SummaryStatistics, SummaryStatistics)"/> method allows you
    /// to combine two SummaryStatistics objects formed by processing disjoint data streams to
    /// obtain a combined SummaryStatistics object which has accurate summary statistics
    /// for the combined data set.</para>
    /// </remarks>
    public sealed class SummaryStatistics {

        /// <summary>
        /// Initializes a new summary for an empty sample.
        /// </summary>
        public SummaryStatistics () { }


        /// <summary>
        /// Initializes a new summary of the given sample values.
        /// </summary>
        /// <param name="values">An enumerator which can iterate over the values to summarize.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is <see langword="null"/>.</exception>
        public SummaryStatistics (IEnumerable<double> values) : this() {
            Add(values);
        }

        // n tracks the number of elements
        // m tracks the mean
        // s2 = n c_2 tracks the sum associated with the sample variance

        private int n;
        private double m, s2, s3, s4;
        private double m_correction, c2_correction;
        private double min = Double.MaxValue;
        private double max = Double.MinValue;

        /// <summary>
        /// Adds a new value to the sample.
        /// </summary>
        /// <param name="value">The value to add.</param>
        public void Add (double value) {

            n++;

            double d = value - m;
            double e = d / n;
            m = KahanSum(m, e, ref m_correction);
            double ddme = d * (d - e);
            s2 = KahanSum(s2, ddme, ref c2_correction);

            s3 += -3.0 * e * s2 + ddme * (d + e);
            s4 += -e * (4.0 * s3 + 6.0 * e * s2) + ddme * (e * e + d * e + d * d);

            if (value < min) min = value;
            if (value > max) max = value;
        }

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
        // Kahan suggested a way to do a running sum more accurately by keeping track of a correction term
        // This is essentially FastTwoSum with a re-used error term.

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double KahanSum(double augend, double addend, ref double correction) {
            double y = addend - correction;
            double sum = augend + y;
            correction = (sum - augend) - y;
            return sum;
        }

        /// <summary>
        /// Adds new values to the sample.
        /// </summary>
        /// <param name="values">An enumerator which can iterate over the values to add.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is <see langword="null"/>.</exception>
        public void Add (IEnumerable<double> values) {
            if (values == null) throw new ArgumentNullException(nameof(values));
            foreach (double value in values) Add(value);
        }

        /// <summary>
        /// Gets the number of values observed.
        /// </summary>
        public int Count {
            get {
                return n;
            }
        }

        /// <summary>
        /// Gets the mean of the observed data.
        /// </summary>
        public double Mean {
            get {
                if (n < 1) throw new InsufficientDataException();
                return m;
            }
        }

        /// <summary>
        /// Get the variance of the observed data.
        /// </summary>
        /// <remarks><para>This is not the estimated population variance that is obtained via Bessel's correction (dividing by n-1). It
        /// is the actual sample variance. If you want the estimated population variance, use <see cref="PopulationVariance"/>.</para></remarks>
        public double Variance {
            get {
                return s2 / n;
            }
        }

        /// <summary>
        /// Gets the standard deviation of the observed data.
        /// </summary>
        public double StandardDeviation {
            get {
                return Math.Sqrt(Variance);
            }
        }

        /// <summary>
        /// Gets the skewness of the observed data.
        /// </summary>
        public double Skewness {
            get {
                return s3 / n / Math.Pow(s2 / n, 3.0 / 2.0);
            }
        }

        /// <summary>
        /// Gets the smallest value observed.
        /// </summary>
        public double Minimum {
            get {
                if (n < 1) throw new InsufficientDataException();
                return min;
            }
        }

        /// <summary>
        /// Gets the largest value observed.
        /// </summary>
        public double Maximum {
            get {
                if (n < 1) throw new InsufficientDataException();
                return max;
            }
        }

        /// <summary>
        /// Estimates the mean of the underlying population.
        /// </summary>
        /// <value>As estimate, with uncertainty, of the mean of the distribution from which the sample was drawn.</value>
        public UncertainValue PopulationMean {
            get {
                if (n < 2) throw new InsufficientDataException();
                return Univariate.EstimateFirstCumulant(n, m, s2 / n);
            }
        }

        /// <summary>
        /// Estimates the variance of the underlying population.
        /// </summary>
        /// <value>As estimate, with uncertainty, of the variance of the distribution from which the sample was drawn.</value>

        public UncertainValue PopulationVariance {
            get {
                if (n < 4) throw new InsufficientDataException();
                return Univariate.EstimateSecondCumulant(n, s2 / n, s4 / n);
                //return (new UncertainValue(s2 / (n - 1), Math.Sqrt((s4 / n - MoreMath.Sqr(s2 / n)) / (n - 1))));
            }
        }

        /// <summary>
        /// Estimates the standard deviation of the underlying population.
        /// </summary>
        /// <value>An estimate, with uncertainty, of the standard deviation of the
        /// distribution from which the sample was drawn.</value>
        public UncertainValue PopulationStandardDeviation {
            get {
                if (n < 4) throw new InsufficientDataException();
                return Univariate.EstimateStandardDeviation(n, s2 / n, s4 / n);
                //double c2 = s2 / n;
                //double c4 = s4 / n;
                //double v = (c4 - c2 * c2) / (n - 1);
                //double s = Math.Sqrt(s2 / (n - 1)) * (1.0 + v / (8.0 * c2 * c2));
                //double ds = 0.5 * Math.Sqrt(v) / s;
                //
                //return (new UncertainValue(s, ds));
            }
        }


        /// <summary>
        /// Combines two sample summaries.
        /// </summary>
        /// <param name="a">One sample summary.</param>
        /// <param name="b">Another sample summary.</param>
        /// <returns>Summary statistics of the combined sample.</returns>
        public static SummaryStatistics Combine (SummaryStatistics a, SummaryStatistics b) {

            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            double d = b.m - a.m;
            double d2 = d * d;

            SummaryStatistics ab = new SummaryStatistics();
            ab.n = a.n + b.n;
            ab.m = (a.n * a.m + b.n * b.m) / ab.n;
            ab.s2 = a.s2 + b.s2 + d2 * a.n * b.n / ab.n;
            ab.s3 = a.s3 + b.s3 + (d / ab.n) * (d2 * a.n * b.n / ab.n * (a.n - b.n) + 3.0 * (a.n * b.s2 - b.n * a.s2));
            ab.s4 = a.s4 + b.s4 + (d / ab.n) * (d2 * a.n * b.n / ab.n * (d / ab.n) * (a.n * a.n - a.n * b.n + b.n * b.n) + 6.0 * (d / ab.n) * (a.n * a.n * b.s2 + b.n * b.n * a.s2) + 4.0 * (a.n * b.s3 - b.n * a.s3));

            ab.min = Math.Min(a.min, b.min);
            ab.max = Math.Max(a.max, b.max);

            return ab;
        }

    }


    // Online algorithms for the computation of moments are discussed in  https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    // The basic idea of incremental update is called Welford's algorithm after Welford 1962 paper.
    // There are different variants for moments beyond the second.
    // Pébaÿ 2008 gave formulas for arbitrary order and Meng 2015 (https://arxiv.org/pdf/1510.04923) gave the form we use here.

    // That form is:    
    //   d = x_n - M_{n-1}
    //   M_n = M_{n-1} + d/n
    //   C_{2,n} = C_{2,n-1} + d [d  -  d/n]
    //   C_{3,n} = C_{3,n-1} - 3 (d/n) C_{2,n} + d [ d^2 - (d/n)^2) ]
    //   C_{4,n} = C_{4,n-1} - 4 (d/n) C_{3,n} - 6 (d/n)^2 C_{2,n} + d [ d^3 - (d/n)^3 ]
    //   C_{r,n} = C_{r,n-1}  - \sum_{k=1}^{r-2} \binom{r, k} (d/n)^k C_{r-k,n} + d [ d^{r-1} - (d/n)^{r-1} ]
    // To get actual central moments, divide by n.

    // To speed calculation we use factored forms and record some and re-use some intermediate results
    // such as e = d/n and d(d-e).

    // Before deciding how to compute moments from both streaming and in-memory data, I compared a bunch of algorithms, as summarized below.
    // My basic takeaway was that, for in-memory data, doing a correction pass after an initial computation was a good approach, and
    // using Kahan summation in that case wasn't necessary. But for single-pass calculations, Kahan summation was useful.

    /* 
    10,000 variates from Uniform(0,1)

    m                    c2                    c3                      c4
    0.497796628542336    0.083525558341007958  3.3604243569687955E-06  0.012510385532363581 DoubleDouble two-pass            0        0        0        0        Reference standard
    0.497796628542336    0.0835255583410078    3.360424356969006E-06   0.01251038553236361  Then correction pass             0        1.8E-15  6.2E-14  2.4E-15  This is noise level in single precision
    0.497796628542336    0.0835255583410078    3.3604243569681827E-06  0.01251038553236361  Instead Kahan correction pass    0        1.8E-15  1.8E-13  2.4E-15  Actually slightly worse c_3
    0.49779662854233536  0.083525558341008457  3.3604243569657479E-06  0.012510385532363636 Naive Meng online                1.0E-15  6.0E-15  9.1E-13  4.4E-15  Not bad but worse errors than noise
    0.497796628542336    0.083525558341007861  3.3604243569693723E-06  0.012510385532363607 Then correction pass             0        1.2E-15  1.7E-13  2.1E-15  Very good, basically noise level
    0.497796628542336    0.083525558341007861  3.36042435696903E-06    0.012510385532363607 Instead Kahan correction pass    0        1.2E-15  6.9E-14  2.1E-15  This time slightly better c_3
    0.497796628542336    0.083525558341007958  3.3604243569683224E-06  0.012510385532363633 Kahan Meng online                0        0        1.4E-13  4.2E-15  Slightly better than Naive Meng
    0.49779662854233875  0.083525558341007791  3.3604243562786297E-06  0.012510385532363602 Naive two-pass                   5.6E-15  2.0E-15  2.0E-10  1.7E-15  Worse than naive Meng; do not recommend
    0.497796628542336    0.083525558341007791  3.360424356968336E-06   0.012510385532363602 Then correction pass             0        2.0E-15  1.4E-13  1.7E-15  Correction pass heals most wounds

    10,000 variates from Logistic(1,4)

    0.91893369512340162  53.246389219061378    -8.3354983404628111     12237.808724847224   DoubleDouble two-pass            0        0        0        0        Reference standard
    0.91893369512340162  53.2463892190614      -8.335498340462852      12237.808724847257   Then correction pass             0        4.0E-16  4.9E-15  2.7E-15  Gives noise level in single prevision
    0.91893369512340162  53.2463892190614      -8.3354983404628573     12237.808724847257   Instead Kahan correction pass    0        4.0E-16  5.5E-15  2.7E-15  Actually slightly worse c3
    0.91893369512340273  53.2463892190613      -8.3354983404628271     12237.808724847206   Naive Meng online                1.2E-15  1.5E-15  1.9E-15  1.5E-15  Very good
    0.91893369512340162  53.2463892190614      -8.3354983404628218     12237.808724847258   Then correction pass             0        4.0E-16  1.3E-15  2.8E-15  Makes most slightly better, but c4 slightly worse
    0.91893369512340151  53.2463892190614      -8.33549834046281       12237.808724847258   Instead Kahan correction pass    1.2E-16  4.0E-16  2.1E-16  2.8E-15  Tiny bit better than naive for c3, worse for m
    0.91893369512340162  53.246389219061378    -8.335498340462836      12237.808724847204   Kahan Meng online                0        0        3.0E-15  1.6E-15  Better than naive Meng for m, c2, not for c3, c4
    0.91893369512340217  53.246389219061392    -8.3354983404628857     12237.808724847257   Naive two-pass                   6.0E-16  2.7E-16  9.0E-15  2.7E-15  Slightly better than Meng; consistent with noise
    0.91893369512340162  53.246389219061392    -8.3354983404627934     12237.808724847257   Then correction pass             0        2.7E-16  2.1E-15  2.7E-15  Moves in right direction
     */

}
