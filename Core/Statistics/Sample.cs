using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a set of independent draws of real numbers.
    /// </summary>
    /// <remarks>
    /// <para>A univariate sample is a data set which records one number for each independent
    /// observation. For example, data from a study which measured the weight of each subject could be
    /// stored in the Sample class. The class offers descriptive statistics for the sample, estimates
    /// of descriptive statistics of the underlying population distribution, and statistical
    /// tests to comare the sample distribution to other sample distributions or theoretical models.</para>
    /// </remarks>
    public sealed class Sample : ICollection<double>, IEnumerable<double>, IEnumerable {

        /// <summary>
        /// Initializes a new, empty sample.
        /// </summary>
        public Sample () {
        }

        /// <summary>
        /// Initializes a new sample from a list of values.
        /// </summary>
        /// <param name="values">Values to add to the sample.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        public Sample (IEnumerable<double> values) {
            if (values == null) throw new ArgumentNullException("values");
            Add(values);
        }

        // the sample data
        private List<double> data = new List<double>();

		// stuff to remember about the data
        private bool sorted = true;            // whether we have sorted the data
        private double mean = 0.0;           // the sample data mean
		private double variance;        // the sample data variance

        private bool HaveDistribution {
            get {
                return((variance < 0));
            }
        }

        private void EnsureSorted () {
            if (!sorted) {
                data.Sort();
            }
        }

        private void ComputeDistribution () {
            double m = 0.0;
            double s = 0.0;
            for (int n = 0; n < data.Count; n++) {
                double mm = m;
                m += (data[n] - m) / (n + 1);
                s += (data[n] - m) * (data[n] - mm);
            }
            mean = m;
            variance = s / data.Count;
        }

        // basic characterization

        /// <summary>
        /// Gets the number of measurements in the sample.
        /// </summary>
        public int Count {
            get {
                return(data.Count);
            }
        }

        // low moments

        /// <summary>
        /// Gets the sample mean.
        /// </summary>
        /// <remarks>
        /// <para>The mean is the average of all values in the sample.</para>
        /// </remarks>
        /// <seealso cref="PopulationMean"/>
        public double Mean {
            get {
                return(mean);
            }
        }

        /// <summary>
        /// Gets the sample variance.
        /// </summary>
        /// <remarks>
        /// <para>Note this is the actual variance of the sample values, not the infered variance
        /// of the underlying population; to obtain the latter use <see cref="PopulationVariance" />.</para>
        /// </remarks>
        public double Variance {
            get {
                if (!HaveDistribution) ComputeDistribution();
                return(variance);
            }
        }

        /// <summary>
        /// Gets the sample standard deviation.
        /// </summary>
        /// <remarks>
        /// <para>Note this is the actual standard deviation of the sample values, not the infered standard
        /// deviation of the underlying population; to obtain the latter use
        /// <see cref="PopulationStandardDeviation" />.</para>
        /// </remarks>        
        public double StandardDeviation {
            get {
                return(Math.Sqrt(Variance));
            }
        }

        /// <summary>
        /// Computes the given sample moment.
        /// </summary>
        /// <param name="n">The order of the moment to compute.</param>
        /// <returns>The <paramref name="n"/>th moment of the sample.</returns>
        public double Moment (int n) {
            if (n == 0) return (1.0);
            if (n == 1) return (mean);
            double m = 0.0;
            for (int i = 0; i < data.Count; i++) {
                m += Math.Pow(data[i], n);
            }
            return (m / data.Count);
        }

        /// <summary>
        /// Computes the given sample moment about its mean.
        /// </summary>
        /// <param name="n">The order of the moment to compute.</param>
        /// <returns>The <paramref name="n"/>th moment about its mean of the sample.</returns>
        public double MomentAboutMean (int n) {
            if (n == 0) return (1.0);
            if (n == 1) return (0.0);
            double m = 0.0;
            for (int i = 0; i < data.Count; i++) {
                m += Math.Pow(data[i] - mean, n);
            }
            return (m / data.Count);
        }

        /// <summary>
        /// Gets the sample median.
        /// </summary>
        public double Median {
            get {
                EnsureSorted();
				int m = data.Count/2;
				if ((data.Count % 2) == 0) {
					return( data[m] );
				} else {
					return( (data[m] + data[m+1]) / 2.0 );
				}
            }
        }

        /// <summary>
        /// Gets the interquartile range of sample measurmements.
        /// </summary>
        /// <remarks>The interquartile range is the interval between the 25th and the 75th percentile.</remarks>
        /// <seealso cref="InverseLeftProbability"/>
        public Interval InterquartileRange {
			get {
                return (Interval.FromEndpoints(InverseLeftProbability(0.25), InverseLeftProbability(0.75)));
			}
		}

        /// <summary>
        /// Gets the sample value corresponding to a given percentile score.
        /// </summary>
        /// <param name="P">The percentile, which must lie between zero and one.</param>
        /// <returns>The corresponding value.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="P"/> lies outside [0,1].</exception>
		public double InverseLeftProbability (double P) {
			if ((P<0.0) || (P>1.0)) throw new ArgumentOutOfRangeException("P", P, "Percentiles must be between 0 and 1.");
            EnsureSorted();
			double n = P * data.Count;
			double n1 = Math.Floor(n);
			double n2 = Math.Ceiling(n);
			// this could cause loss of accuracy for very large n; i would prefer to have a "fractional part" function
			double p1 = n2-n;
			double p2 = n-n1;
			return( p1 * data[(int) n1] + p2 * data[(int) n2] );
		}


        // population characteristics as derived from the sample

        /// <summary>
        /// Gets an estimate of the population mean from the sample.
        /// </summary>
        public UncertainValue PopulationMean {
            get {
                // population mean is estimated by sample mean
                // population standard deviation is estimated by Sqrt(N/(N-1)) * sample standard deviation
                // standard deviation of sample mean is population standard deviation / Sqrt(N)
                return (new UncertainValue(Mean, Math.Sqrt(Variance / (Count - 1))));
            }
        }

        /// <summary>
        /// Gets an estimate of the population variance from the sample.
        /// </summary>
        public UncertainValue PopulationVariance {
            get {
                return (PopulationMomentAboutMean(2));
            }
        }

        /// <summary>
        /// Gets an estimate of the population standard deviation from the sample.
        /// </summary>
        public UncertainValue PopulationStandardDeviation {
            get {
                return (UncertainMath.Sqrt(PopulationVariance));
            }
        }

        /// <summary>
        /// Estimates the given population moment using the sample.
        /// </summary>
        /// <param name="n">The order of the moment.</param>
        /// <returns>An estimate of the <paramref name="n"/>th moment of the population.</returns>
        public UncertainValue PopulationMoment (int n) {
            if (n < 0) {
                // we don't do negative moments
                throw new ArgumentOutOfRangeException("n", n, "Moment orders must be non-negative.");
            } else if (n == 0) {
                // the zeroth moment is exactly one for any distribution
                return (new UncertainValue(1.0, 0.0));
            } else if (n == 1) {
                // the first momemt is just the mean
                return (PopulationMean);
            } else {
                // moments of order two and higher
                double M_n = Moment(n);
                double M_2n = Moment(2 * n);
                return (new UncertainValue(M_n, Math.Sqrt((M_2n - M_n * M_n) / data.Count)));
            }
        }

        /// <summary>
        /// Estimates the given population moment about the mean using the sample.
        /// </summary>
        /// <param name="n">The order of the moment.</param>
        /// <returns>An estimate, with uncertainty, of the <paramref name="n"/>th moment about the mean
        /// of the underlying population.</returns>
        public UncertainValue PopulationMomentAboutMean (int n) {
            if (n < 0) {
                // we don't do negative moments
                throw new ArgumentOutOfRangeException("n", n, "Moment orders must be non-negative.");
            } else if (n == 0) {
                // the zeroth moment is exactly one for any distribution
                return (new UncertainValue(1.0, 0.0));
            } else if (n == 1) {
                // the first moment about the mean is exactly zero by the defintion of the mean
                return (new UncertainValue(0.0, 0.0));
            } else if (n == 2) {
                // the best estimate of the population variance is the sample variance C2 increased by the standard (N-1)/N correction factor
                // the uncertainty in the estimate can be computed exactly in terms of the C2 and C4 sample moments 
                int N = data.Count;
                double C2 = MomentAboutMean(2);
                double C4 = MomentAboutMean(4);
                return (new UncertainValue(C2 * N / (N - 1.0), Math.Sqrt((C4 - C2 * C2) / N)));
                //return (new UncertainValue(C2 * N / (N - 1.0), Math.Sqrt((C4 + C2 * (4.0 - C2)) / N)));
                //return (new UncertainValue(C2 / (1.0 - 1.0 / N), Math.Sqrt(((N - 1) * C4 - (N * N - 3) * C2 * C2 / (N - 1)) / (N * N - 3 * N + 3))));
            } else {
                // moments of order three and higher
                // in these cases we use formulas that include only the leading and first subleading term in a 1/N expansion

                // the formula for the best estimate of the n'th population moment involves C_n, C_{n-2}, and C_2
                // the formula for the uncertainty in this estimate involves C_{2n}, C_{n+1}, C_{n-1}, and C_2
                // we need all these sample moments
                int N = data.Count;
                double C_2 = MomentAboutMean(2);
                double C_nm2 = MomentAboutMean(n - 2);
                double C_nm1 = MomentAboutMean(n - 1);
                double C_n = MomentAboutMean(n);
                double C_np1 = MomentAboutMean(n + 1);
                double C_2n = MomentAboutMean(2 * n);

                // now that we have the moments, compute the estimate and its uncertainty.
                double C = C_n - (1.0 * n / N) * (0.5 * (n - 1) * C_2 * C_nm2 - C_n);
                double V = (C_2n - 2 * n * C_nm1 * C_np1 - C_n * C_n + C_2 * C_nm1 * C_nm1 * n * n) / N;
                return (new UncertainValue(C, Math.Sqrt(V)));
            }

        }

        // fitting to distributions

        /// <summary>
        /// Fits the sample to a normal distribution.
        /// </summary>
        /// <returns>A <see cref="FitResult"/> containg the mu and sigma parameters of the normal distribution that best fits the sample data,
        /// and a Kolmogorov-Smirnov test of the quality of the fit.</returns>
        /// <seealso cref="NormalDistribution" />
        /// <seealso cref="KolmogorovSmirnovTest(Meta.Numerics.Statistics.Distribution)"/>
        public FitResult FitToNormalDistribution () {

            // maximum likelyhood estimates are guaranteed to be asymptotically unbiased, but not necessarily unbiased
            // this hits home for the maximum likelyhood estimate of the variance of a normal distribution, which fails
            // to include the N/(N-1) correction factor. since we know the bias, there is no reason for us not to correct
            // it, and we do so here

            UncertainValue mu = PopulationMean;
            UncertainValue sigma = PopulationStandardDeviation;

            Distribution distribution = new NormalDistribution(mu.Value, sigma.Value);
            TestResult test = KolmogorovSmirnovTest(distribution, 2);

            return (new FitResult(mu.Value, mu.Uncertainty, sigma.Value, sigma.Uncertainty, 0.0, test));

        }

        /// <summary>
        /// Fits the sample to an exponential distribution.
        /// </summary>
        /// <returns>A <see cref="FitResult"/> containing the lambda parameter of the exponential distribution that best fits the sample data,
        /// and a Kolmogorov-Smirnov test of the quality of the fit.</returns>
        /// <seealso cref="ExponentialDistribution" />
        /// <seealso cref="KolmogorovSmirnovTest(Meta.Numerics.Statistics.Distribution)" />
        /// <exception cref="InvalidOperationException">One or more values in the sample is negative.</exception>
        public FitResult FitToExponentialDistribution () {

            // none of the data is allowed to be negative
            for (int i = 0; i < data.Count; i++) {
                if (data[i] < 0.0) throw new InvalidOperationException();
            }

            // the best-fit exponential's mean sample mean, with corresponding uncertainly

            double lambda = Mean;
            double dLambda = lambda / Math.Sqrt(Count);

            Distribution distribution = new ExponentialDistribution(lambda);
            TestResult test = KolmogorovSmirnovTest(distribution, 1);

            return (new FitResult(lambda, dLambda, test));
        }

        // tests on the sample

        /// <summary>
        /// Tests whether the sample mean is compatible with the reference mean.
        /// </summary>
        /// <param name="referenceMean">The reference mean.</param>
        /// <returns>The result of the test. The test statistic is a t-value. If t &gt; 0, the one-sided likelyhood
        /// to obtain a greater value under the null hypothesis is the (right) propability of that value. If t &lt; 0, the
        /// corresponding one-sided likelyhood is the (left) probability of that value. The two-sided likelyhood to obtain
        /// a t-value as far or farther from zero as the value obtained is just twice the one-sided likelyhood.</returns>
        /// <remarks><para>The test statistic of the student t-test is the difference between the
        /// sample and reference means, measured in units of the sample mean uncertainty. For normally
        /// distributed samples, this is known to follow a Student distribution. If t is
        /// far from zero, then the sample is unlikely to have been drawn from a population with the reference
        /// mean.</para>
        /// <para>Because the distribution of a t-statistic assumes a normally distributed population, this
        /// test should only be used on sample data compatible with a normal distribution. The Mann-Whitney
        /// U test is a less powerful non-parametric alternative that can be used to test mean compatibility
        /// on arbitrarly distributed data sets.</para></remarks>
        /// <seealso cref="StudentDistribution" />
        public TestResult StudentTTest (double referenceMean) {
            if (this.Count < 2) throw new InvalidOperationException();
            double sigma = Math.Sqrt(this.Count / (this.Count - 1.0)) * this.StandardDeviation;
            double se = sigma / Math.Sqrt(this.Count);

            double t = (mean - referenceMean) / se;
            int dof = this.Count - 1;

            return (new TestResult(t, new StudentDistribution(dof)));
        }

        /// <summary>
        /// Tests whether the sample mean is compatible with the mean of another sample.
        /// </summary>
        /// <param name="sample">The other sample.</param>
        /// <returns>The result of the test. The test statistic is a t-value. If t &gt; 0, the one-sided likelyhood
        /// to obtain a greater value under the null hypothesis is the (right) propability of that value. If t &lt; 0, the
        /// corresponding one-sided likelyhood is the (left) probability of that value. The two-sided likelyhood to obtain
        /// a t-value as far or farther from zero as the value obtained is just twice the one-sided likelyhood.</returns>
        /// <remarks><para></para></remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <seealso cref="StudentDistribution" />
        public TestResult StudentTTest (Sample sample) {
            if (sample == null) throw new ArgumentNullException("sample");

            double m1 = this.Mean;
            double s1 = Math.Sqrt(this.Count / (this.Count - 1.0)) * this.StandardDeviation;
            double m2 = sample.Mean;
            double s2 = Math.Sqrt(sample.Count / (sample.Count - 1.0)) * sample.StandardDeviation;

            // compute degrees of freedom and pooled variance
            int dof = this.Count + sample.Count - 2;
            double s = Math.Sqrt((s1 + s2) / dof);

            // compute t statistic and probability
            double t = (m1 - m2) / (s * Math.Sqrt(1.0 / this.Count + 1.0 / sample.Count));
            return (new TestResult(t, new StudentDistribution(dof)));
        }

        /*
        /// <summary>
        /// Tests whether the sample variance is compatible with a reference variance.
        /// </summary>
        /// <param name="variance">The reference variance.</param>
        /// <returns></returns>
        /// <remarks><para>The F test assumes that the population is normally distributed.</para></remarks>
        public TestResult FisherFTest (double variance) {
            if (variance <= 0.0) throw new ArgumentOutOfRangeException("variance");

            // compute population variance
            double v = this.Count / (this.Count - 1.0) * this.Variance;

            // compute ratio to reference variance
            double F = this.Count * v / variance;

            return (new TestResult(F, new ChiSquaredDistribution(this.Count)));
        }
        */

        /// <summary>
        /// Tests whether the sample variance is compatible with the variance of another sample.
        /// </summary>
        /// <param name="sample"></param>
        /// <returns></returns>
        public TestResult FisherFTest (Sample sample) {
            if (sample == null) throw new ArgumentNullException("sample");

            // compute population variances
            double v1 = this.Count / (this.Count - 1.0) * this.Variance;
            double v2 = sample.Count / (sample.Count - 1.0) * sample.Variance;

            // compute the ratio
            double F = v1 / v2;

            return (new TestResult(F, new FisherDistribution(this.Count - 1, sample.Count - 1)));

        }

        /// <summary>
        /// Tests whether the sample median is compatible with the mean of another sample.
        /// </summary>
        /// <param name="sample">The other sample.</param>
        /// <returns></returns>
        /// <remarks>
        /// <para>The Mann-Whitney test is a non-parametric alternative to Student's t-test.
        /// Essentially, it supposes that the medians of the two samples are equal and tests
        /// the likelihood of this null hypothesis.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Mann-Whitney_U_test"/>
        public TestResult MannWhitneyTest (Sample sample) {

            if (sample == null) throw new ArgumentNullException("sample");

            // essentially, we want to order the entries from the two samples and find out how
            // many times "a's beat b's". In the ordering ababb, a beats b 5 times.

            // implementing this naively would be O(N^2), so instead we use a formula that
            // relates this quantity to the sum of ranks of a's and b's; to get those ranks
            // we do seperate sorts O(N ln N) and a merge sort O(N).

            // sort the two samples
            this.EnsureSorted();
            sample.EnsureSorted();

            // now we essentially do a merge sort, but instead of actually forming the merged list,
            // we just keep track of the ranks that elements from each sample would have in the merged list

            // variables to track the sum of ranks for each sample
            int r1 = 0;
            int r2 = 0;

            // pointers to the current position in each list
            int p1 = 0;
            int p2 = 0;

            // the current rank
            int r = 1;

            while (true) {

                //Console.WriteLine("a[{0}] = {1} <? b[{2}] = {3}", p1, this.data[p1], p2, sample.data[p2]);
                if (this.data[p1] < sample.data[p2]) {
                    //Console.WriteLine("{0} a", r);
                    r1 += r;
                    p1++;
                    r++;
                    if (p1 >= this.Count) {
                        while (p2 < sample.Count) {
                            //Console.WriteLine("{0} b", r);
                            r2 += r;
                            p2++;
                            r++;
                        }
                        break;
                    }

                } else {
                    //Console.WriteLine("{0} b", r);
                    r2 += r;
                    p2++;
                    r++;
                    if (p2 >= sample.Count) {
                        while (p1 < this.Count) {
                            //Console.WriteLine("{0} a", r);
                            r1 += r;
                            p1++;
                            r++;
                        }
                        break;
                    }

                }
            }

            // relate u's to r's
            int u1 = r1 - this.Count * (this.Count + 1) / 2;
            int u2 = r2 - sample.Count * (sample.Count + 1) / 2;
            Debug.Assert(u1 + u2 == this.Count * sample.Count);

            // return the result

            return (new TestResult(u1, new MannWhitneyDistribution(this.Count, sample.Count)));

            //return (new TestResult(u1, new DistributionAdapter(new MannWhitneyDistribution2(this.Count, sample.Count))));

        }
 
        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="distribution">The test distribution.</param>
        /// <returns>The test result. The test statistic is the D statistic and the likelyhood is the right probability
        /// to obtain a value of D as large or larger than the one obtained.</returns>
        /// <remarks><para>The Kolmogorov-Smirnov test measures the departure of a sample from a hypothesized population
        /// distribution by comparing the cumulative probability function of the data to the cumulative probability function
        /// of the hypothesized population distribution. The test statistic D is the maximum seperation between the two curves.</para>
        /// <para>Under the null hypothesis N<sup>1/2</sup>D is known to be distributed according to the Kolomogorov distribution
        /// in the large-N limit. Because of the large-N assumption, this test should not be used with small (less than ~50)
        /// data sets.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="distrubution"/> is null.</exception>
        /// <seealso cref="KolmogorovDistribution"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public TestResult KolmogorovSmirnovTest (Distribution distribution) {
            if (distribution == null) throw new ArgumentNullException("distribution");
            if (this.Count < 1) throw new InvalidOperationException();

            return (KolmogorovSmirnovTest(distribution, 0));
        }

        private TestResult KolmogorovSmirnovTest (Distribution distribution, int count) {

            // compute the D statistic, which is the maximum of the D+ and D- statistics
            double D, DP, DM;
            ComputeDStatistics(distribution, out DP, out DM);
            if (DP > DM) {
                D = DP;
            } else {
                D = DM;
            }

            // reduce n by the number of fit parameters, if any
            // this is heuristic and needs to be changed to deal with specific fits
            int n = data.Count - count;

            //return (new TestResult(D, new FiniteKolmogorovDistribution(n)));
            return (new TestResult(D, new KolmogorovDistribution(1.0 / Math.Sqrt(n))));
        }

        /// <summary>
        /// Tests whether the sample is compatible with the given distribution.
        /// </summary>
        /// <param name="distribution">The test distribution.</param>
        /// <returns>The test result. The test statistic is the V statistic and the likelyhood is the right probability
        /// to obtain a value of V as large or larger than the one obtained.</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Kuiper%27s_test"/>
        public TestResult KuiperTest (Distribution distribution) {

            // compute the V statistic, which is the sum of the D+ and D- statistics
            double DP, DM;
            ComputeDStatistics(distribution, out DP, out DM);
            double V = DP + DM;

            return (new TestResult(V, new KuiperDistribution(1.0 / Math.Sqrt(data.Count))));

        }


        private void ComputeDStatistics (Distribution distribution, out double D1, out double D2) {

            EnsureSorted();

            D1 = 0.0;
            D2 = 0.0;

            double P1 = 0.0;
            for (int i = 0; i < data.Count; i++) {

                // compute the theoretical CDF value at the data-point
                double x = data[i];
                double P = distribution.LeftProbability(x);

                // compute the experimental CDF value at x+dx
                double P2 = (i + 1.0) / data.Count;

                // compare the the theoretical and experimental CDF values at x-dx and x+dx
                // update D if this is the largest seperation yet observed

                double DU = P2 - P;
                if (DU > D1) D1 = DU;
                double DL = P - P1;
                if (DL > D2) D2 = DL;

                //double D1 = Math.Abs(P - P1);
                //if (D1 > D) D = D1;
                //double D2 = Math.Abs(P - P2);
                //if (D2 > D) D = D2;

                // remember the experimental CDF value at x+dx;
                // this will be the experimental CDF value at x-dx for the next step
                P1 = P2;

            }

        }

        /// <summary>
        /// Tests whether the sample is compatible with another sample.
        /// </summary>
        /// <param name="sample">The other sample.</param>
        /// <returns>The test result. The test statistic is the D statistic and the likelyhood is the right probability
        /// to obtain a value of D as large or larger than the one obtained.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <seealso href="http://en.wikipedia.org/wiki/Kolmogorov-Smirnov_test"/>
        public TestResult KolmogorovSmirnovTest (Sample sample) {
            if (sample == null) throw new ArgumentNullException("sample");

            // we must have data to do the test
            if (this.Count < 1) throw new InvalidOperationException();
            if (sample.Count < 1) throw new InvalidOperationException();

            // the samples must be sorted
            this.EnsureSorted();
            sample.EnsureSorted();

            // keep track of where we are in each sample
            int c0 = 0;
            int c1 = 0;

            // keep track of the cdf of each sample
            double d0 = 0.0;
            double d1 = 0.0;

            // find the maximum cdf seperation
            double d = 0.0;
            while ((c0 < this.Count) && (c1 < sample.Count)) {
                //Console.WriteLine("s0[{0}]={1} vs s1[{2}]={3}", c0, this.data[c0], c1, sample.data[c1]);
                if (this.data[c0] < sample.data[c1]) {
                    // the next data point is in sample 0
                    //Console.WriteLine("x={0}", this.data[c0]);
                    d0 = 1.0 * (c0 + 1) / this.Count;
                    c0++;
                } else {
                    // the next data point is in sample 1
                    //Console.WriteLine("x={0}", sample.data[c1]);
                    d1 = 1.0 * (c1 + 1) / sample.Count;
                    c1++;
                }
                double dd = Math.Abs(d1 - d0);
                //Console.WriteLine("d0={0} d1={1} dd={2}", d0, d1, dd);
                if (dd > d) d = dd;
            }

            // compute the effective degrees of freedom;
            double ne = 1.0 / (1.0 / this.Count + 1.0 / sample.Count);
            //Console.WriteLine("ne={0}", ne);

            // return the result
            return (new TestResult(d, new KolmogorovDistribution(1.0 / Math.Sqrt(ne))));

        }

        // more miscelaneous operations

        // IEnumerable methods

        IEnumerator<double> IEnumerable<double>.GetEnumerator () {
            EnsureSorted();
            return (data.GetEnumerator());
        }

        IEnumerator IEnumerable.GetEnumerator () {
            EnsureSorted();
            return (data.GetEnumerator());
        }

        // Add, Clear, Contains, CopyTo, Remove, Count, IsReadOnly

        /// <summary>
        /// Adds a value to the sample.
        /// </summary>
        /// <param name="value">The value to add.</param>
        public void Add (double value) {
            sorted = false;
            mean += (value - mean) / (data.Count + 1);
            data.Add(value);
        }

        /// <summary>
        /// Adds a series of values to the sample.
        /// </summary>
        /// <param name="values">An enumerator of values to be added.</param>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> is null.</exception>
        public void Add (IEnumerable<double> values) {
            if (values == null) throw new ArgumentNullException("values");
            foreach (double value in values) {
                Add(value);
            }
        }

        /// <summary>
        /// Remove all values from the sample.
        /// </summary>
        public void Clear () {
            sorted = true;
            mean = 0.0;
            data.Clear();
        }

        /// <summary>
        /// Determines whether the sample contains a given value.
        /// </summary>
        /// <param name="value">The value to check for.</param>
        /// <returns>Whether the sample contains <paramref name="datum"/>.</returns>
        public bool Contains (double value) {
            return (data.Contains(value));
        }

        /// <summary>
        /// Removes a given value from the sample.
        /// </summary>
        /// <param name="value">The value to remove.</param>
        /// <returns>Whether the value was found and removed.</returns>
        public bool Remove (double value) {
            // sort order is not affected
            if (data.Remove(value)) {
                mean -= value / (data.Count + 1); // fix this
                mean *= ((double) (data.Count + 1)) / data.Count;
                return (true);
            } else {
                return (false);
            }
        }

        void ICollection<double>.CopyTo (double[] array, int start) {
            data.CopyTo(array, start);
        }

        bool ICollection<double>.IsReadOnly {
            get {
                return (false);
            }
        }

        /// <summary>
        /// Performs a maximum likelihood fit.
        /// </summary>
        /// <param name="distribution">The distribution to fit to the data.</param>
        /// <returns>The result of the fit, containg the parameters that result in the best fit,
        /// covariance matrix among those parameters, and a test of the goodness of fit.</returns>
        /// <seealso href="http://en.wikipedia.org/wiki/Maximum_likelihood"/>
        public FitResult MaximumLikelihoodFit (IParameterizedDistribution distribution) {

            // define a log likeyhood function
            Function<double[],double> L = delegate (double[] parameters) {
                distribution.SetParameters(parameters);
                double lnP = 0.0;
                foreach (double value in data) {
                    double P = distribution.Likelihood(value);
                    if (P == 0.0) throw new InvalidOperationException();
                    lnP += Math.Log(P);
                }
                return(-lnP);
            };

            // maximize it
            double[] v0 = distribution.GetParameters();
            SpaceExtremum min = FunctionMath.FindMinimum(L, v0);

            // turn this into a fit result
            FitResult result = new FitResult(min.Location(), min.Curvature().CholeskyDecomposition().Inverse(), null);

            return (result);

        }


    }

}
