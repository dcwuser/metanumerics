using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Functions;
using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {
    public static partial class Univariate {

        /// <summary>
        /// Performs a z-test.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="referenceMean">The mean of the comparison population.</param>
        /// <param name="referenceStandardDeviation">The standard deviation of the comparison population.</param>
        /// <returns>A test result indicating whether the sample mean is significantly different from that of the comparison population.</returns>
        /// <remarks>
        /// <para>A z-test determines whether the sample is compatible with a normal population with known mean and standard deviation.
        /// In most cases, Student's t-test (<see cref="StudentTTest(IReadOnlyCollection{double},System.Double)"/>), which does not assume a known population standard deviation,
        /// is more appropriate.</para>
        /// </remarks>
        /// <example>
        /// <para>Suppose a standardized test exists, for which it is known that the mean score is 100 and the standard deviation is 15
        /// across the entire population. The test is administered to a small sample of a subpopulation, who obtain a mean sample score of 95.
        /// You can use the z-test to determine how likely it is that the subpopulation mean really is lower than the population mean,
        /// that is that their slightly lower mean score in your sample is not merely a fluke.</para>
        /// </example>
        /// <exception cref="InsufficientDataException"><see cref="Sample.Count"/> is zero.</exception>
        /// <seealso cref="StudentTTest(IReadOnlyCollection{double},double)"/>
        public static TestResult ZTest (this IReadOnlyCollection<double> sample, double referenceMean, double referenceStandardDeviation) {
            return (ZTest(sample, referenceMean, referenceStandardDeviation, TestType.TwoTailed));
        }

        /// <summary>
        /// Performs a z-test with the given sidedness.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="referenceMean">The mean of the comparison population.</param>
        /// <param name="referenceStandardDeviation">he standard deviation of the comparison population.</param>
        /// <param name="type">The sidedness of the test to perform.</param>
        /// <returns>A test result indicating whether the sample mean is significantly different from that of the comparison population
        /// in the direction indicated by <paramref name="type"/>.</returns>
        /// <seealso cref="ZTest(IReadOnlyCollection{double},double, double)"/>
        public static TestResult ZTest (this IReadOnlyCollection<double> sample, double referenceMean, double referenceStandardDeviation, TestType type) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) throw new InsufficientDataException();
            double z = (sample.Mean() - referenceMean) / (referenceStandardDeviation / Math.Sqrt(sample.Count));
            return (new TestResult("z", z, type, new NormalDistribution()));
        }

        /// <summary>
        /// Tests whether the sample median is compatible with the given reference value.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="referenceMedian">The reference median.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>The sign test is a non-parametric alternative to the Student t-test (<see cref="StudentTTest(IReadOnlyCollection{double},double)"/>).
        /// It tests whether the sample is consistent with the given reference median.</para>
        /// <para>The null hypothesis for the test is that the median of the underlying population from which the sample is
        /// drawn is the reference median. The test statistic is simply number of sample values that lie above the median. Since
        /// each sample value is equally likely to be below or above the population median, each draw is an independent Bernoulli
        /// trial, and the total number of values above the population median is distributed according to a binomial distribution
        /// (<see cref="BinomialDistribution"/>).</para>
        /// <para>The left probability of the test result is the chance of the sample median being so low, assuming the sample to have been
        /// drawn from a population with the reference median. The right probability of the test result is the chance of the sample median
        /// being so high, assuming the sample to have been drawn from a population with the reference median.</para>
        /// </remarks>
        /// <seealso cref="StudentTTest(IReadOnlyCollection{double},double)"/>
        public static TestResult SignTest (this IReadOnlyCollection<double> sample, double referenceMedian) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) throw new InsufficientDataException();

            // count the number of entries that exceed the reference median
            int W = 0;
            foreach (double value in sample) {
                if (value > referenceMedian) W++;
            }

            // W should be distributed binomially
            return (new TestResult("W", W, TestType.TwoTailed, new DiscreteAsContinuousDistribution(new BinomialDistribution(0.5, sample.Count))));

        }

        /// <summary>
        /// Tests whether the sample mean is compatible with the reference mean.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="referenceMean">The reference mean.</param>
        /// <returns>The result of the test. The test statistic is a t-value. If t &gt; 0, the one-sided likelyhood
        /// to obtain a greater value under the null hypothesis is the (right) propability of that value. If t &lt; 0, the
        /// corresponding one-sided likelyhood is the (left) probability of that value. The two-sided likelyhood to obtain
        /// a t-value as far or farther from zero as the value obtained is just twice the one-sided likelyhood.</returns>
        /// <remarks>
        /// <para>The test statistic of Student's t-test is the difference between the sample mean and the reference mean,
        /// measured in units of the sample mean uncertainty. Under the null hypothesis that the sample was drawn from a normally
        /// distributed population with the given reference mean, this statistic can be shown to follow a Student distribution
        /// (<see cref="StudentDistribution"/>). If t is far from zero, with correspondingly small left or right tail probability,
        /// then the sample is unlikely to have been drawn from a population with the given reference mean.</para>
        /// <para>Because the distribution of a t-statistic assumes a normally distributed population, this
        /// test should only be used only on sample data compatible with a normal distribution. The sign test (<see cref="SignTest"/>)
        /// is a non-parametric alternative that can be used to test the compatibility of the sample median with an assumed population median.</para>
        /// </remarks>
        /// <example>
        /// <para>In some country, the legal limit blood alcohol limit for drivers is 80 on some scale. Because they
        /// have noticed that the results given by their measuring device fluctuate, the police perform three
        /// seperate measurements on a suspected drunk driver.
        /// They obtain the results 81, 84, and 93. They argue that, because all three results exceed the
        /// limit, the court should be very confident that the driver's blood alcohol level did, in fact, exceed
        /// the legal limit. You are the driver's lawyer. Can you make an argument to that the court shouldn't be so
        /// sure?</para>
        /// <para>Here is some code that computes the probability of obtaining such high measured values,
        /// assuming that the true level is exactly 80.</para>
        /// <code lang="c#">
        /// Sample values = new Sample();
        /// values.Add(81, 84, 93);
        /// TestResult result = values.StudentTTest(80);
        /// return(result.RightProbability);
        /// </code>
        /// <para>What level of statistical confidence do you think should a court require in order to pronounce a defendant guilty?</para>
        /// </example>
        /// <exception cref="InsufficientDataException">There are fewer than two data points in the sample.</exception>
        /// <seealso cref="StudentDistribution" />
        /// <seealso href="https://en.wikipedia.org/wiki/Student%27s_t-test"/>
        public static TestResult StudentTTest (this IReadOnlyCollection<double> sample, double referenceMean) {
            return (StudentTTest(sample, referenceMean, TestType.TwoTailed));
        }

        /// <summary>
        /// Tests whether the sample mean differs from the reference mean in the specified direction.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="referenceMean">The reference mean.</param>
        /// <param name="type">The sidedness of the test to perform.</param>
        /// <returns>A test result indicating whether the sample mean is significantly different from the reference mean
        /// in the direction indicated by <paramref name="type"/>.</returns>
        public static TestResult StudentTTest (this IReadOnlyCollection<double> sample, double referenceMean, TestType type) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 2) throw new InsufficientDataException();

            int n;
            double mean, sumOfSquaredDeviations;
            ComputeMomentsUpToSecond(sample, out n, out mean, out sumOfSquaredDeviations);

            double sigma = Math.Sqrt(sumOfSquaredDeviations/ (n - 1));
            double se = sigma / Math.Sqrt(n);

            double t = (mean - referenceMean) / se;
            int dof = n - 1;

            return (new TestResult("t", t, type, new StudentDistribution(dof)));
        }

        /// <summary>
        /// Tests whether the variances of two samples are compatible.
        /// </summary>
        /// <param name="a">The first sample.</param>
        /// <param name="b">The second sample.</param>
        /// <returns>The result of the test.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is <see langword="null"/>.</exception>
        public static TestResult FisherFTest (IReadOnlyCollection<double> a, IReadOnlyCollection<double> b) {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if ((a.Count < 2) || (b.Count < 2)) throw new InsufficientDataException();

            int aCount;
            double aMean, aSumOfSquaredDeviations;
            ComputeMomentsUpToSecond(a, out aCount, out aMean, out aSumOfSquaredDeviations);
            Debug.Assert(aCount == a.Count);
            Debug.Assert(aSumOfSquaredDeviations >= 0.0);

            int bCount;
            double bMean, bSumOfSquaredDeviations;
            ComputeMomentsUpToSecond(b, out bCount, out bMean, out bSumOfSquaredDeviations);
            Debug.Assert(bCount == b.Count);
            Debug.Assert(bSumOfSquaredDeviations >= 0.0);

            // compute population variances
            double va = aSumOfSquaredDeviations / (aCount - 1);
            double vb = bSumOfSquaredDeviations / (bCount - 1);

            // compute the ratio
            double F = va / vb;

            // right-tailed is only right if F > 1.
            return (new TestResult("F", F, TestType.RightTailed, new FisherDistribution(a.Count - 1, b.Count - 1)));
        }

        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="distribution">The distribution.</param>
        /// <returns>The test result. The test statistic is the D statistic and the probability is the chance of
        /// obtaining such a large value of D under the assumption that the sample is drawn from the given distribution.</returns>
        /// <remarks>
        /// <para>The null hypothesis of the Kolmogorov-Smirnov (KS) test is that the sample is drawn from the given continuous distribution.
        /// The test statistic D is the maximum deviation of the sample's empirical distribution function (EDF) from
        /// the distribution's cumulative distribution function (CDF). A high value of the test statistic, corresponding
        /// to a low right tail probability, indicates that the sample distribution disagrees with the given distribution
        /// to a degree unlikely to arise from statistical fluctuations.</para>
        /// <para>For small sample sizes, we compute the null distribution of D exactly. For large sample sizes, we use an accurate
        /// asymptotic approximation. Therefore it is safe to use this method for all sample sizes.</para>
        /// <para>A variant of this test, <see cref="Sample.KolmogorovSmirnovTest(Sample, Sample)"/>, allows you to non-parametrically
        /// test whether two samples are drawn from the same underlying distribution, without having to specify that distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">There is no data in the sample.</exception>
        /// <seealso cref="KolmogorovDistribution"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public static TestResult KolmogorovSmirnovTest (this IReadOnlyList<double> sample, ContinuousDistribution distribution) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (distribution == null) throw new ArgumentNullException(nameof(distribution));

            int n = sample.Count;
            if (n < 1) throw new InsufficientDataException();

            // Compute the D statistic, which is the maximum of the D+ and D- statistics.
            double DP, DM;
            ComputeDStatistics(sample, distribution, out DP, out DM);
            double D = Math.Max(DP, DM);

            ContinuousDistribution DDistribution;
            if (n < 32) {
                DDistribution = new TransformedDistribution(new KolmogorovExactDistribution(n), 0.0, 1.0 / n);
            } else {
                DDistribution = new TransformedDistribution(new KolmogorovAsymptoticDistribution(n), 0.0, 1.0 / Math.Sqrt(n));
            }
            return (new TestResult("D", D, TestType.RightTailed, DDistribution));

        }


        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="sample">The sample.</param>
        /// <param name="distribution">The distribution.</param>
        /// <returns>The test result. The test statistic is the V statistic and the chance to obtain such a large
        /// value of V under the assumption that the sample is drawn from the given distribution.</returns>
        /// <remarks>
        /// <para>Like the Kolmogorov-Smirnov test (<see cref="Univariate.KolmogorovSmirnovTest(IReadOnlyList{Double},ContinuousDistribution)"/>),
        /// Kuiper's test compares the EDF of the sample to the CDF of the given distribution.</para>
        /// <para>For small sample sizes, we compute the null distribution of V exactly. For large sample sizes, we use an accurate
        /// asymptotic approximation. Therefore it is safe to use this method for all sample sizes.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">There is no data in the sample.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Kuiper%27s_test"/>
        public static TestResult KuiperTest (this IReadOnlyList<double> sample, ContinuousDistribution distribution) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (distribution == null) throw new ArgumentNullException(nameof(distribution));

            int n = sample.Count;
            if (n < 1) throw new InsufficientDataException();

            // Compute the V statistic, which is the sum of the D+ and D- statistics.
            double DP, DM;
            ComputeDStatistics(sample, distribution, out DP, out DM);
            double V = DP + DM;

            ContinuousDistribution VDistribution;
            if (n < 32) {
                VDistribution = new TransformedDistribution(new KuiperExactDistribution(n), 0.0, 1.0 / n);
            } else {
                VDistribution = new TransformedDistribution(new KuiperAsymptoticDistribution(n), 0.0, 1.0 / Math.Sqrt(n));
            }
            return (new TestResult("V", V, TestType.RightTailed, VDistribution));

        }

        private static void ComputeDStatistics (IReadOnlyList<double> sample, ContinuousDistribution distribution, out double D1, out double D2) {

            int[] order = GetSortOrder(sample);

            D1 = 0.0;
            D2 = 0.0;

            double P1 = 0.0;
            for (int i = 0; i < sample.Count; i++) {

                // Compute the theoretical CDF value at the data-point.
                double x = sample[order[i]];
                double P = distribution.LeftProbability(x);

                // compute the experimental CDF value at x+dx
                double P2 = (i + 1.0) / sample.Count;

                // Compare the the theoretical and experimental CDF values at x-dx and x+dx;
                // Update D if it is the largest separation yet observed.
                double DU = P2 - P;
                if (DU > D1) D1 = DU;
                double DL = P - P1;
                if (DL > D2) D2 = DL;

                // Remember the experimental CDF value at x+dx;
                // it will be the experimental CDF value at x-dx for the next step
                P1 = P2;

            }

        }

        /// <summary>
        /// Performs a Shapiro-Francia test of normality on the sample.
        /// </summary>
        /// <param name="sample">The sample to test for normality.</param>
        /// <returns>The result of the test.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">There are fewer than 16 values in the sample.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Shapiro%E2%80%93Francia_test"/>
        public static TestResult ShapiroFranciaTest (this IReadOnlyList<double> sample) {

            if (sample == null) throw new ArgumentNullException(nameof(sample));

            int n = sample.Count;
            if (n < 16) throw new InsufficientDataException();

            // Determine V = 1 - r^2
            double mx = Mean(sample);
            double cxx = 0.0;
            double cmm = 0.0;
            double cmx = 0.0;
            int[] order = GetSortOrder(sample);
            for (int i = 0; i < n; i++) {
                double zx = sample[order[i]] - mx;
                // We could use symmetry to make half as many calls to NormalOrderStatistic,
                // but the code would be less straightforward, because we can't proceed linearly
                // through the data.
                double m = NormalOrderStatistic(i + 1, n);
                cxx += zx * zx;
                cmm += m * m;
                cmx += zx * m;
            }
            // W = r^2 = cmx / sqrt(cxx * cmm). But if we compute r^2 like this and then 1 - r^2, we will
            // get significant cancellation errors. So re-form V = 1 - r^2 analytically.
            double q = cxx * cmm;
            double p = Math.Sqrt(cxx * cmm);
            double lnV = Math.Log((p - cmx) * (p + cmx) / q);

            // The observed mean and variance of log(1 - W) can be fit to six digits over n=16-8192
            // by simple polynomials in log(n). To get these results, I ran simulations for n = 16, 24, 32, 48, 64, ..., 8192.
            // The number of simulated samples was 10^8 for n=10-100, 10^7 for n=100-1000, and 10^6 for n=1000-10000.
            // Following Royston, I transformed to log(1 - W) and did polynomial fits to K1 and log(K2) as
            // functions of log(n), increasing the order until I was able to reproduce both cumulants to
            // within errors for all n=16-8192. Note that K3 is detectably non-zero, so we really should
            // introduce a correction for the remaining skewness, perhaps via an Edgeworth-style expansion
            // or by fitting to a skew-normal distribution.
            double t = Math.Log(n);
            double mu = -2.41348986 + t * (0.34785948 + t * (-0.31462118 + t * (0.04298010 + t * (-0.00311513 + t * 0.00009230))));
            double sigma = Math.Exp((0.56767288 + t * (-1.18529022 + t * (0.31916472 + t * (-0.04912452 + t * (0.00383692 + t * -0.00011891))))) / 2.0);
            NormalDistribution nullDistribution = new NormalDistribution(mu, sigma);

            // Need to handle 3-15 separately.

            // N = 3: W' = W, 3/4 <= W <= 1, P = (6 / Pi)(Asin(Sqrt(W)) - Asin(Sqrt(3/4)))

            return (new TestResult("ln(1-W')", lnV, TestType.RightTailed, nullDistribution));

        }

        private static double NormalOrderStatistic (int i, int n) {

            // This uses the David and Johnson asymptotic expansion, as presented
            // in Arnold, Balakrishnan, and Nagaraja, "A First Course in Order Statistics",
            // Section 5.5, specialized to the normal distribution.

            // It is more accurate than simpler approximations, but still looses
            // accuracy in the tails. I have verified that dropping the last term
            // doesn't change the moments in my simulations to within their accuracy,
            // so it's not necessary to add more terms. If we ever do want to add
            // another term, see Childs and Balakrisnan, "Series approximations
            // for moments of order statistics using MAPLE", Computations Statistics
            // and Data Analysis 38 (2002) 331, which gives the lengthy next term.

            double p = (double) i / (n + 1);
            double q = (double) (n - i + 1) / (n + 1);

            double f0 = AdvancedMath.Probit(p, q);
            double p0 = Math.Exp(-f0 * f0 / 2.0) / Global.SqrtTwoPI;

            double f2 = f0 / (p0 * p0);
            double f3 = (1.0 + 2.0 * f0 * f0) / (p0 * p0 * p0);
            double f4 = f0 * (7.0 + 6.0 * f0 * f0) / (p0 * p0 * p0 * p0);
            double f5 = (7.0 + 46.0 * f0 * f0 + 24.0 * f0 * f0 * f0 * f0) / (p0 * p0 * p0 * p0 * p0);
            double f6 = f0 * (127.0 + 326.0 * f0 * f0 + 120.0 * f0 * f0 * f0 * f0) / (p0 * p0 * p0 * p0 * p0 * p0);

            double qmp = q - p;

            double m = f0;

            m += p * q / (n + 2) * f2 / 2.0;

            m += p * q / (n + 2) / (n + 2) * (qmp / 3.0 * f3 + p * q / 8.0 * f4);

            m += p * q / (n + 2) / (n + 2) / (n + 2) * (-qmp / 3.0 * f3 + (qmp * qmp - p * q) / 4.0 * f4 + p * q * qmp / 6.0 * f5 + p * p * q * q / 48.0 * f6);

            return (m);

        }

        /// <summary>
        /// Tests whether the sample is compatible with another sample.
        /// </summary>
        /// <param name="a">The first sample.</param>
        /// <param name="b">The second sample.</param>
        /// <returns>The test result. The test statistic is the D statistic and the likelihood is the right probability
        /// to obtain a value of D as large or larger than the one obtained.</returns>
        /// <remarks>
        /// <para>The two-sample Kolmogorov-Smirnov test is a variation of the single-sample test (<see cref="KolmogorovSmirnovTest(IReadOnlyList{double}, ContinuousDistribution)"/>)
        /// that tests whether two independent samples are drawn from the same underlying distribution. The null hypothesis of the test is that both samples are drawn from
        /// the same population.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is <see langword="null"/>.</exception>
        /// <exception cref="InsufficientDataException">One or both of the samples is empty.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public static TestResult KolmogorovSmirnovTest (IReadOnlyList<double> a, IReadOnlyList<double> b) {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            // we must have data to do the test
            if (a.Count < 1) throw new InsufficientDataException();
            if (b.Count < 1) throw new InsufficientDataException();

            // the samples must be sorted
            int[] aOrder = GetSortOrder(a);
            Debug.Assert(aOrder.Length == a.Count);
            int[] bOrder = GetSortOrder(b);

            // keep track of where we are in each sample
            int ia = 0;
            int ib = 0;

            // keep track of the cdf of each sample
            double da = 0.0;
            double db = 0.0;

            // find the maximum cdf seperation
            // this is a variant on the standard merge sort algorithm
            double d = 0.0;
            while ((ia < a.Count) && (ib < b.Count)) {
                if (a[aOrder[ia]] < b[bOrder[ib]]) {
                    // the next data point is in sample a
                    da = 1.0 * (ia + 1) / a.Count;
                    ia++;
                } else {
                    // the next data point is in sample b
                    db = 1.0 * (ib + 1) / b.Count;
                    ib++;
                }
                double dd = Math.Abs(db - da);
                if (dd > d) d = dd;
            }

            // return the result
            ContinuousDistribution nullDistribution;
            if (AdvancedIntegerMath.BinomialCoefficient(a.Count + b.Count, a.Count) < Int64.MaxValue) {
                nullDistribution = new DiscreteAsContinuousDistribution(new KolmogorovTwoSampleExactDistribution(a.Count, b.Count), Interval.FromEndpoints(0.0, 1.0));
            } else {
                nullDistribution = new TransformedDistribution(new KolmogorovDistribution(), 0.0, Math.Sqrt(1.0 / a.Count + 1.0 / b.Count));
            }

            return (new TestResult("D", d, TestType.RightTailed, nullDistribution));

        }

        /// <summary>
        /// Tests whether one sample mean is compatible with another sample mean.
        /// </summary>
        /// <param name="a">The first sample, which must contain at least two entries.</param>
        /// <param name="b">The second sample, which must contain at least two entries.</param>
        /// <returns>The result of the test. The statistic is the Student's t and the probability
        /// is the chance of obtaining such an extreme value of t if the two samples are drawn
        /// from the same distribution.</returns>
        /// <remarks>
        /// <para>Given two samples, a back-of-the-envelope way to determine whether their means differ in a statistically
        /// significant way is to compare their <see cref="PopulationMean"/> values. If their error bars overlap, they
        /// are probably statistically compatible; if they do not, the difference in means is probably statistically
        /// significant. Student's t-test is a way to refine this back-of-the-envelope procedure into a statistical test
        /// that can determine exactly how likely a given separation of means is under the null hypothesis that the
        /// two samples are drawn from the same distribution.</para>
        /// <para>The t-statistic is proportional to the mean of <paramref name="a"/> minus the mean of <paramref name="b"/>,
        /// so t > 0 indicates that <paramref name="a"/> has a greater mean.</para>
        /// <para>Student's t-test was one of the first statistical tests. It was described by William Sealy Gosset,
        /// a chemist who worked for the Guinness brewing company. Since Guinness was concerned that other breweries might take
        /// advantage of a technique published by one of its chemists, Gosset published his work under the pseudonym Student.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="a"/> or <paramref name="b"/> contains less than two values.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Student's_t-test"/>
        public static TestResult StudentTTest (IReadOnlyCollection<double> a, IReadOnlyCollection<double> b) {

            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (a.Count < 2) throw new InsufficientDataException();
            if (b.Count < 2) throw new InsufficientDataException();

            int na;
            double ma, aSumOfSquaredDeviations;
            ComputeMomentsUpToSecond(a, out na, out ma, out aSumOfSquaredDeviations);
            Debug.Assert(na == a.Count);

            int nb;
            double mb, bSumOfSquaredDeviations;
            ComputeMomentsUpToSecond(b, out nb, out mb, out bSumOfSquaredDeviations);
            Debug.Assert(nb == b.Count);

            // pool variances and counts
            double v = (aSumOfSquaredDeviations + bSumOfSquaredDeviations) / (na + nb - 2);
            double n = 1.0 / (1.0 / na + 1.0 / nb);

            // evaluate t
            double t = (ma - mb) / Math.Sqrt(v / n);

            return (new TestResult("t", t, TestType.TwoTailed, new StudentDistribution(na + nb - 2)));

        }

        /// <summary>
        /// Tests whether one sample median is compatible with another sample median.
        /// </summary>
        /// <param name="a">The first sample.</param>
        /// <param name="b">The second sample.</param>
        /// <returns>The result of the test. The statistic is the Mann-Whitney U value and the probability
        /// is the chance of obtaining such an extreme value of U if the two samples are drawn from the
        /// same distribution.</returns>
        /// <remarks>
        /// <para>The Mann-Whitney test is a non-parametric alternative to Student's t-test
        /// (<see cref="StudentTTest(IReadOnlyCollection{Double}, IReadOnlyCollection{Double})"/>).
        /// Essentially, it supposes that the medians of the two samples are equal and tests
        /// the likelihood of this null hypothesis. Unlike the t-test, it does not assume that the sample distributions are normal.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Mann-Whitney_U_test"/>
        public static TestResult MannWhitneyTest (IReadOnlyList<double> a, IReadOnlyList<double> b) {

            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (a.Count < 2) throw new InsufficientDataException();
            if (b.Count < 2) throw new InsufficientDataException();

            // Essentially, we want to order the entries from both samples together and find out how
            // many times "a's beat b's". In the ordering ababb, a beats b 5 times.

            // Implementing this naively would be O(N^2), so instead we use a formula that
            // relates this quantity to the sum of ranks of a's and b's; to get those ranks
            // we do separate sorts O(N ln N) and a merge sort O(N).

            // Sort the two samples.
            int[] aOrder = GetSortOrder(a);
            int[] bOrder = GetSortOrder(b);

            // Now we essentially do a merge sort, but instead of actually forming the merged list,
            // we just keep track of the ranks that elements from each sample would have in the merged list

            // Variables to track the sum of ranks for each sample
            int ra = 0;
            int rb = 0;

            // Pointers to the current position in each list
            int ia = 0;
            int ib = 0;

            // The current rank
            int r = 1;

            while (true) {

                if (a[aOrder[ia]] < b[bOrder[ib]]) {
                    ra += r;
                    ia++;
                    r++;
                    if (ia >= a.Count) {
                        while (ib < b.Count) {
                            rb += r;
                            ib++;
                            r++;
                        }
                        break;
                    }

                } else {
                    rb += r;
                    ib++;
                    r++;
                    if (ib >= b.Count) {
                        while (ia < a.Count) {
                            ra += r;
                            ia++;
                            r++;
                        }
                        break;
                    }

                }
            }

            // Relate u's to r's,
            int ua = ra - a.Count * (a.Count + 1) / 2;
            int ub = rb - b.Count * (b.Count + 1) / 2;
            Debug.Assert(ua + ub == a.Count * b.Count);

            // If possible, we want to use the exact distribution of U.
            // To compute it, we need to do exact integer arithmetic on numbers of order the total number of possible orderings,
            // which is (a.Count + b.Count!).
            // Since decimal is the built-in type that can hold the largest exact integers, we use it for the computation.
            // Therefore, to generate the exact distribution, the total number of possible orderings must be less than the capacity of a decimal. 
            ContinuousDistribution uDistribution;
            double lnTotal = AdvancedIntegerMath.LogFactorial(a.Count + b.Count) - AdvancedIntegerMath.LogFactorial(a.Count) - AdvancedIntegerMath.LogFactorial(b.Count);
            if (lnTotal > Math.Log((double) Decimal.MaxValue)) {
                double mu = a.Count * b.Count / 2.0;
                double sigma = Math.Sqrt(mu * (a.Count + b.Count + 1) / 6.0);
                uDistribution = new NormalDistribution(mu, sigma);
            } else {
                uDistribution = new DiscreteAsContinuousDistribution(new MannWhitneyExactDistribution(a.Count, b.Count));
            }

            return (new TestResult("U", ua, TestType.TwoTailed, uDistribution));

        }

        /// <summary>
        /// Performs a one-way analysis of variance (ANOVA).
        /// </summary>
        /// <param name="samples">The samples to compare.</param>
        /// <returns>ANOVA data, including an F-test comparing the between-group variance to
        /// the within-group variance.</returns>
        /// <remarks>
        /// <para>The one-way ANOVA is an extension of the Student t-test
        /// (<see cref="StudentTTest(IReadOnlyCollection{double},IReadOnlyCollection{double})"/>)
        /// to more than two groups. The test's null hypothesis is that all the groups' data are drawn
        /// from the same distribution. If the null hypothesis is rejected, it indicates
        /// that at least one of the groups differs significantly from the others.</para>
        /// <para>Given more than two groups, you should use an ANOVA to test for differences
        /// in the means of the groups rather than perform multiple t-tests. The reason
        /// is that each t-test incurs a small risk of a false positive, so multiple t-tests increase
        /// the total risk of a false positive. For example, given a 95% confidence requirement,
        /// there is only a 5% chance that an individual t-test will incorrectly diagnose a significant
        /// difference. But given 5 samples, there are 5 * 4 / 2 = 10 t-tests to be
        /// performed, giving about a 40% chance that at least one of them will incorrectly
        /// diagnose a significant difference! The ANOVA avoids the accumulation of risk
        /// by performing a single test at the required confidence level to test for
        /// any significant differences between the groups.</para>
        /// <para>A one-way ANOVA performed on just two samples is equivalent to a
        /// t-test (<see cref="Sample.StudentTTest(Sample,Sample)" />).</para>
        /// <para>ANOVA is an acronym for "Analysis of Variance". Do not be confused
        /// by the name and by the use of a ratio-of-variances test statistic: an
        /// ANOVA is primarily (although not exclusively) sensitive to changes in the
        /// <i>mean</i> between samples. The variances being compared by the test are not the
        /// variances of the individual samples; instead the test compares the variance of
        /// all samples considered together as one single, large sample to the variances of the samples
        /// considered individually. If the means of some groups differ significantly,
        /// then the variance of the unified sample will be much larger than the variances of the
        /// individual samples, and the test will signal a significant difference. Thus the
        /// test uses variance as a tool to detect shifts in mean, not because we
        /// are interested in the individual sample variances per se.</para>
        /// <para>ANOVA is most appropriate when the sample data are continuous and approximately normal,
        /// and the samples are distinguished by a nominal variable. For example, given
        /// a random sampling of the heights of members of five different political parties,
        /// a one-way ANOVA would be an appropriate test of the whether the different
        /// parties tend to attract people with different heights.</para>
        /// <para>Given a continuous independent variable, binning in order to define
        /// groups and perform an ANOVA is generally not appropriate.
        /// For example, given the incomes and heights of a large number of people,
        /// dividing these people into low-height, medium-height, and high-height groups
        /// and performing an ANOVA of the income of people in each group is not a
        /// good way to test whether height influences income, first because the result
        /// will be sensitive to the arbitrary boundaries you have chosen for the bins,
        /// and second because the ANOVA has no notion of bin ordering.
        /// In a case like this, it would be better to put the data into a <see cref="BivariateSample"/> and
        /// perform a test of association, such as a <see cref="BivariateSample.PearsonRTest" />,
        /// <see cref="BivariateSample.SpearmanRhoTest" />, or <see cref="BivariateSample.KendallTauTest" />
        /// between the two variables. If you have measurements
        /// of additional variables for each individual, a <see cref="MultivariateSample.LinearRegression(int)" />
        /// analysis would allow you to adjust for the confounding effects of the other variables. If you
        /// define arbitrary bins of continuously variable data in order to form groups, then your
        /// ANOVA results will depend on your choice of bins.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="samples"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="samples"/> contains fewer than two samples.</exception>
        /// <seealso href="https://en.wikipedia.org/wiki/Analysis_of_variance"/>
        /// <seealso href="https://en.wikipedia.org/wiki/One-way_analysis_of_variance"/>
        public static OneWayAnovaResult OneWayAnovaTest (params IReadOnlyCollection<double>[] samples) {
            return (OneWayAnovaTest((IReadOnlyCollection<IReadOnlyCollection<double>>) samples));
        }

        /// <summary>
        /// Performs a one-way analysis of variance (ANOVA).
        /// </summary>
        /// <param name="samples">The samples to compare.</param>
        /// <returns>ANOVA data, including an F-test comparing the between-group variance to
        /// the within-group variance.</returns>
        /// <remarks>
        /// <para>For detailed information, see <see cref="OneWayAnovaTest(IReadOnlyCollection{double}[])"/>.</para>
        /// </remarks>
        public static OneWayAnovaResult OneWayAnovaTest (IReadOnlyCollection<IReadOnlyCollection<double>> samples) {

            if (samples == null) throw new ArgumentNullException(nameof(samples));
            if (samples.Count < 2) throw new ArgumentException("There must be at least two samples in the sample list.", nameof(samples));

            // determine total count, mean, and within-group sum-of-squares
            int n = 0;
            double mean = 0.0;
            double SSW = 0.0;
            foreach (ICollection<double> sample in samples) {
                if (sample == null) throw new ArgumentNullException("sample");

                int n_sample;
                double mean_sample, sumOfSquaredDeviations_sample;

                ComputeMomentsUpToSecond(sample, out n_sample, out mean_sample, out sumOfSquaredDeviations_sample);
                n += n_sample;
                mean += n_sample * mean_sample;
                SSW += sumOfSquaredDeviations_sample;
            }
            mean = mean / n;

            // determine between-group sum-of-squares
            double SSB = 0.0;
            foreach (Sample sample in samples) {
                SSB += sample.Count * MoreMath.Sqr(sample.Mean - mean);
            }

            // determine degrees of freedom associated with each sum-of-squares
            int dB = samples.Count - 1;
            int dW = n - 1 - dB;

            // determine F statistic
            double F = (SSB / dB) / (SSW / dW);

            AnovaRow factor = new AnovaRow(SSB, dB);
            AnovaRow residual = new AnovaRow(SSW, dW);
            AnovaRow total = new AnovaRow(SSB + SSW, n - 1);
            return (new OneWayAnovaResult(factor, residual, total));

        }

        /// <summary>
        /// Performs a Kruskal-Wallis test on the given samples.
        /// </summary>
        /// <param name="samples">The set of samples to compare.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>Kruskal-Wallis tests for differences between the samples. It is a non-parametric alternative to the
        /// one-way ANOVA (<see cref="OneWayAnovaTest(IReadOnlyCollection{double}[])"/>) which is more appropriate
        /// when the data far from normally distributed.</para>
        /// <para>The test is essentially a one-way ANOVA performed on the <i>ranks</i> of sample values instead of the sample
        /// values themselves.</para>
        /// <para>A Kruskal-Wallis test on two samples is equivalent to a Mann-Whitney test
        /// (see <see cref="MannWhitneyTest(IReadOnlyList{double}, IReadOnlyList{double})"/>).</para>
        /// <para>As with a normal ANOVA, it is not appropriate to bin a continuous independent variable in order to form
        /// groups for a Kruskal-Wallis test. Kruskal-Wallis addresses the non-normality of the dependent variable, not
        /// the non-discreteness of the independent variable.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="samples"/> is <see langword="null"/>.</exception>
        /// <see cref="OneWayAnovaTest(IReadOnlyCollection{double}[])"/>
        /// <see href="https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance"/>
        public static TestResult KruskalWallisTest (params IReadOnlyList<double>[] samples) {
            return (KruskalWallisTest((IReadOnlyList<IReadOnlyList<double>>) samples));
        }

        /// <summary>
        /// Performs a Kruskal-Wallis test on the given samples.
        /// </summary>
        /// <param name="samples">The set of samples to compare.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>For detailed information, see <see cref="KruskalWallisTest(IReadOnlyList{double}[])"/>.</para>
        /// </remarks>
        public static TestResult KruskalWallisTest (IReadOnlyList<IReadOnlyList<double>> samples) {
            if (samples == null) throw new ArgumentNullException(nameof(samples));
            if (samples.Count < 2) throw new ArgumentException("There must be at least two samples in the sample list.", "samples");

            // sort each sample individually and compute count total from all samples
            int N = 0;
            int[][] orders = new int[samples.Count][];
            for (int i = 0; i < samples.Count; i++) {
                N += samples[i].Count;
                orders[i] = GetSortOrder(samples[i]);
            }

            // do a multi-merge sort to determine ranks sums

            // initialize a pointer to the current active index in each ordered list
            int[] p = new int[samples.Count];

            // keep track the current rank to be assigned
            int r = 0;

            // keep track of rank sums
            // this is all that we need for KW
            // the ranks of individual entries can be made to disappear from the final formula using sum identities
            int[] rs = new int[samples.Count];

            while (true) {

                // increment the rank
                // (programmers may think ranks start from 0, but in the definition of the KW test they start from 1)
                r++;

                // determine the smallest of current entries
                int j = -1;
                double f = Double.PositiveInfinity;
                for (int k = 0; k < samples.Count; k++) {
                    if ((p[k] < orders[k].Length) && (samples[k][orders[k][p[k]]] < f)) {
                        j = k;
                        f = samples[k][orders[k][p[k]]];
                    }
                }

                // test for all lists complete
                if (j < 0) break;

                // increment the pointer and the rank sum for that column
                p[j]++;
                rs[j] += r;

            }

            // compute the KW statistic
            double H = 0.0;
            for (int i = 0; i < samples.Count; i++) {
                double z = ((double) rs[i]) / orders[i].Length - (N + 1) / 2.0;
                H += orders[i].Length * (z * z);
            }
            H = 12.0 / N / (N + 1) * H;

            // use the chi-squared approximation to the null distribution
            return (new TestResult("H", H, TestType.RightTailed, new ChiSquaredDistribution(samples.Count - 1)));

        }

    }
}
