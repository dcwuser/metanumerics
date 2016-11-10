using System;
using System.Collections;
using System.Collections.Generic;

using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;
using Meta.Numerics.SignalProcessing;
using Meta.Numerics.Statistics.Distributions;


namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a series a data points.
    /// </summary>
    public sealed class TimeSeries : IEnumerable, IEnumerable<double>, ICollection<double>, IList<double> {

        /// <summary>
        /// Initializes a new, empty time series.
        /// </summary>
        public TimeSeries () {

        }

        /// <summary>
        /// Initializes a new time series with the given values.
        /// </summary>
        /// <param name="values"></param>
        public TimeSeries (params double[] values) {
            Add(values);
        }

        private SampleStorage data = new SampleStorage();

        /// <summary>
        /// Gets the number of points in the time series.
        /// </summary>
        public int Count {
            get {
                return (data.Count);
            }
        }

        IEnumerator<double> IEnumerable<double>.GetEnumerator () {
            return (data.GetEnumerator());
        }

        IEnumerator IEnumerable.GetEnumerator () {
            return (data.GetEnumerator());
        }

        /// <summary>
        /// Adds a point to the time series.
        /// </summary>
        /// <param name="value">The value to add.</param>
        public void Add (double value) {
            data.Add(value);
        }

        /// <summary>
        /// Adds multiple points to the time series.
        /// </summary>
        /// <param name="values">The values to be added.</param>
        public void Add(params double[] values) {
            if (values == null) throw new ArgumentNullException("values");
            foreach (double value in values) data.Add(value);
        }

        bool ICollection<double>.Remove (double value) {
            throw new InvalidOperationException();
        }

        void ICollection<double>.CopyTo(double[] array, int arrayIndex) {
            data.CopyTo(array, arrayIndex);
        }

        /// <summary>
        /// Removes all points from the time series.
        /// </summary>
        public void Clear () {
            data.Clear();
        }

        /// <summary>
        /// Determines whether the time series contains a value.
        /// </summary>
        /// <param name="value">The value to look up.</param>
        /// <returns>True if the value appears in the time series, otherwise false.</returns>
        public bool Contains (double value) {
            return (data.Contains(value));
        }

        bool ICollection<double>.IsReadOnly {
            get {
                return (false);
            }
        }

        void IList<double>.RemoveAt(int index) {
            if ((index < 0) || (index >= data.Count)) throw new ArgumentOutOfRangeException("index");
            data.RemoveAt(index);
        }

        void IList<double>.Insert(int index, double value) {
            throw new InvalidOperationException();
        }

        /// <summary>
        /// Gets or sets the value of the time series at a given index.
        /// </summary>
        /// <param name="index">The index of the value.</param>
        /// <returns>The value.</returns>
        public double this [int index] {
            get {
                return (data[index]);
            }
            set {
                data[index] = value;
            }
        }

        /// <summary>
        /// Finds the index at which a value occurs.
        /// </summary>
        /// <param name="value">The value to find.</param>
        /// <returns>An index at which <paramref name="value"/> occurs,
        /// or -1 if the value does not occur in the time series.</returns>
        public int IndexOf (double value) {
            return (data.IndexOf(value));
        }

        /// <summary>
        /// Gets the mean of the time series.
        /// </summary>
        public double Mean {
            get {
                return (data.Mean);
            }
        }

        /// <summary>
        /// Computes the autocovariance of the series at the given lag.
        /// </summary>
        /// <param name="lag">The lag at which to compute the autocovariance.</param>
        /// <returns>The value of the autocovariance at the given lag.</returns>
        /// <remarks>
        /// <para>In a length-N time series, there are N-k lag-k observations. Nonetheless,
        /// the definition of the lag-k autocovariance requires division by N, not N-k. This
        /// counterintuitive convention insures that the autocovariance has desirable
        /// positive definiteness properties and agrees with the computation via FFT.</para>
        /// <para>The computation of an autocovariance via this method is O(N). If
        /// you need to compute more than a handfull of autocovariances, it is
        /// more efficient to call the <see cref="Autocovariance()"/>, which
        /// computes all of them in O(N log N).</para>
        /// <para>While the sample autocovariance does converge to the population
        /// autocovariance in the large-N limit, this convergence is very slow. If you
        /// want the value of the population covariance, use <see cref="PopulationStatistics()"/>
        /// to obtain a much better estimate.</para>
        /// </remarks>
        /// <seealso cref="Autocovariance()"/>
        public double Autocovariance (int lag) {

            if (lag < 0) lag = -lag;
            if (lag >= data.Count) throw new ArgumentOutOfRangeException(nameof(lag));

            double sum = 0.0;
            for (int i = lag; i < data.Count; i++) {
                sum += (data[i] - data.Mean) * (data[i - lag] - data.Mean);
            }
            return (sum / data.Count);

            // It appears to be standard practice to divide by n instead of n-k, even though
            // (1) there are n-k terms, not n, (2) dividing by n gives a biased estimator
            // of the population mean, while dividing by n-k yields an unbiased estimator.
            // I really don't undertand this practice, but it does appear to be necessary 
            // for agreement with the FFT computation.

        }

        /// <summary>
        /// Computes the autocovariance for all lags.
        /// </summary>
        /// <returns>An array of autocovariance values, with the array index equal to the lag index.</returns>
        /// <remarks>
        /// <para>The computation of the autocovariance for a given lag is an O(N) operation.
        /// Naively, the computation of the autocovariance for all N possible lags is an
        /// O(N^2) operation; this is in fact the cost of N invoations of <see cref="Autocovariance(int)"/>.
        /// However, it is possible using Fourier techniques to simultaneously compute
        /// the autocovariance for all possible lags in O(N log N) operations. This method
        /// uses this Fourier technique and should be called if you require the autocovariance
        /// for more than a handfull of lag values.</para>
        /// </remarks>
        /// <seealso cref="Autocovariance(int)"/>
        public double[] Autocovariance () {

            // To avoid aliasing, we must zero-pad with n zeros.
            Complex[] s = new Complex[2 * data.Count];
            for (int i = 0; i < data.Count; i++) {
                    s[i] = data[i] - data.Mean;
            }

            FourierTransformer fft = new FourierTransformer(s.Length);
            Complex[] t = fft.Transform(s);

            for (int i = 0; i < t.Length; i++) {
                t[i] = MoreMath.Sqr(t[i].Re) + MoreMath.Sqr(t[i].Im);
            }

            Complex[] u = fft.InverseTransform(t);

            double[] c = new double[data.Count];
            for (int i = 0; i < c.Length; i++) {
                c[i] = u[i].Re / data.Count;
            }

            return (c); 
        }

        /// <summary>
        /// Computes the power spectrum of the time series.
        /// </summary>
        /// <returns>An array giving the power in each frequency bin.</returns>
        /// <remarks>
        /// <para>For a time series of length n, the index k gives the
        /// power near period n / k, or frequency k / n. (For example,
        /// for a year-long monthly time series, the value at index 1
        /// is proportional to the yearly variation, the value at the index 4
        /// is proportional to the quarterly variation, and the value at index 12
        /// is proportional to the monthly variation.) The zeroth
        /// entry is proportional to the unfluctuating component of
        /// the signal, i.e. the mean.</para>
        /// </remarks>
        public double[] PowerSpectrum () {

            Complex[] s = new Complex[data.Count];
            for (int i = 0; i < data.Count; i++) {
                s[i] = data[i];
            }

            FourierTransformer fft = new FourierTransformer(s.Length);
            Complex[] t = fft.Transform(s);

            double[] u = new double[data.Count / 2];
            u[0] = MoreMath.Sqr(t[0].Re) + MoreMath.Sqr(t[0].Im);
            for (int i = 1; i < u.Length; i++) {
                u[i] = MoreMath.Sqr(t[i].Re) + MoreMath.Sqr(t[i].Im)
                    + MoreMath.Sqr(t[data.Count - i].Re) + MoreMath.Sqr(t[data.Count - i].Im);
            }

            return (u);

        }

        /// <summary>
        /// Computes estimates for the moments of the population from which the time series is drawn.
        /// </summary>
        /// <returns></returns>
        /// <remarks>
        /// <para>Just as is the case for the variance of a simple <see cref="Sample"/>,
        /// the sample autocovariances of a time series
        /// are not unbiased estimates of the autocovariances of the population from which
        /// the series is drawn.</para>
        /// </remarks>
        public TimeSeriesPopulationStatistics PopulationStatistics () {

            double m = data.Mean;

            double v;
            double[] g;
            PopulationStatistics(out g, out v);

            return (new TimeSeriesPopulationStatistics(m, v, g));

        }

        private void PopulationStatistics (out double[] g, out double v) {

            if (data.Count < 4) throw new InsufficientDataException();

            // For an ergodic series, the naive mean
            //   m = \frac{1}{n} \sum_{i=1}{n} y_i
            // is an unbiased estimator of the population mean E(m) = \mu. It's not
            // hard to show that the estimator has variance
            //   v = \frac{1}{n} \left[ g_0 + 2 \sum_{k=1}{n-1} \frac{n - k}{n} g_k \right]
            // where g_k = V(g_t, g_{t-k}).

            // For data with a known mean \bar{y},
            //   c_{k} = \frac{1}{n - k} \sum_{i = k}{n} ( y_i - \bar{y} ) ( y_{i-k} - \bar{y} )
            // is a provably unbiased estimator of the population covariance E(c_k) = g_k.

            // But almost never is the mean known independent of the sample. If the mean is
            // computed from the sample, this estimator is biased, and the bias is computable.
            //  E(c_k) = g_k - v
            // Here v is the same quantity as the variace of the mean.

            // For these results, see Fuller, "Introduction To Statistical Time Series",
            // Chapter 6.

            // You might think it should now be straightforward to recover unbiased estimates of
            // the g_k given the c_k, but it turns out it is very much not.

            // One simple approach would be to note that the deviation of c_k from g_k looks
            // to be supressed by 1/n, so a good first approximation should be to compute v
            // from the c_k and then use it to correct them. But this doesn't work because
            // v = 0 when computed using the c_k. This is a sum identity of the FFT, which
            // tends to make the higher c_k artificially negative so as to produce the required
            // sum.

            // Next you might note that our formula effectively defines a matrix A that
            // relates c = A g, if c and g are understood as vectors. If we just invert the
            // matrix we should be able to recover g = A^{-1} c. But A turns out to be
            // rank n-1 and so is not invertable. (Presumably because the same sum identity
            // effectively defines a linear relationship among the c's.) I also tried
            // doing an SVD on A and inverting using the pseudo-inverse, but didn't get
            // very good g's. (Presumably there is a lot of freedom to move power
            // around between the different g's, and the pseudo-inverse doesn't pick
            // what we want, which is one that supresses higher g's.)

            // What finally worked was to introduce a window function into v that supresses
            // the contribution of higher g's. The simplest is a hard window which just
            // does the sum up to some maximum j. This makes A invertible but it is
            // simpler and faster just do a few v -> g -> v -> cycles. Such an
            // approach is suggested by Okui and Ryo, Econometric Theory (2010) 26 : 1263,
            // "Asymptotically Unbiased Estimation of Autocovariances and Autocorrelations
            // with Long Panel Data"
            // (http://repository.kulib.kyoto-u.ac.jp/dspace/bitstream/2433/130692/1/S0266466609990582.pdf)

            // For my test series, I recovered all g's within population uncertainties
            // using j_max = n / 3. They typically stabilized to 3 digits within 8
            // cycles, 4 digits within 12 cycles.

            // One problem with this approach is that makes g_{ij} = g_{i - j} no longer positive
            // definite, so v is occasionally negative. 

            int n = data.Count;
            //int jmax = n / 3;
            int jmax = (int) Math.Floor(Math.Sqrt(n));

            double[] c = Autocovariance();
            for (int j = 0; j < c.Length; j++) {
                c[j] = n * c[j] / (n - j);
            }

            g = new double[n];
            Array.Copy(c, g, n);

            v = 0.0;
            for (int k = 0; k < 12; k++) {

                v = g[0];
                for (int j = 1; j < jmax; j++) {
                    v += 2.0 * (n - j) / n * g[j];
                }
                v = v / n;

                for (int j = 0; j < g.Length; j++) {
                    g[j] = c[j] + v;
                }

            }

            // To deal with v negative (or unrealisticaly small), don't let it
            // get smaller than the first truncated term. This is a very
            // ad hoc solution; we should try to do better.

            double vMin = 2.0 / n * (n - jmax) / n * Math.Abs(g[jmax]);
            if (v < vMin) {
                v = vMin;
            }

        } 

        /// <summary>
        /// Gets a sample containing the time-series values.
        /// </summary>
        /// <returns>A read-only sample containing the time-series value.</returns>
        public Sample AsSample () {
            return (new Sample(data, true));
        } 


        /// <summary>
        /// Fits an MA(1) model to the time series.
        /// </summary>
        /// <returns>The fit with parameters lag-1 coefficient, mean, and standard deviation.</returns>
        public FitResult FitToMA1() {

            if (data.Count < 4) throw new InsufficientDataException();

            // MA(1) model is
            //   y_t - \mu = u_t + \beta u_{t-1}
            // where u_t ~ N(0, \sigma) are IID.

            // It's easy to show that
            //   m = E(y_t) = \mu
            //   c_0 = V(y_t) = E((y_t - m)^2) = (1 + \beta^2) \sigma^2
            //   c_1 = V(y_t, y_{t-1}) = \beta \sigma^2
            // So the method of moments gives
            //   \beta = \frac{1 - \sqrt{1 - (2 g_1)^2}}{2}
            //   \mu = m
            //   \sigma^2 = \frac{c_0}{1 + \beta^2}
            // It turns out these are very poor (high bias, high variance) estimators,
            // but they illustrate a basic requirement that g_1 < 1/2.

            // The MLE estimator 

            int n = data.Count;

            double m = data.Mean;

            Func<double,double> fnc = (double theta) => {

                double s = 0.0;
                double uPrevious = 0.0;
                for (int i = 0; i < data.Count; i++) {
                    double u = data[i] - m - theta * uPrevious;
                    s += u * u;
                    uPrevious = u;
                };
                return (s);
            };

            Extremum minimum = FunctionMath.FindMinimum(fnc, Interval.FromEndpoints(-1.0, 1.0));

            double beta = minimum.Location;
            double sigma2 = minimum.Value / (data.Count - 3);

            // While there is significant evidence that the MLE value for \beta is biased
            // for small-n, I know of no anlytic correction.

            double[] parameters = new double[] { beta, m, Math.Sqrt(sigma2) };

            // The calculation of the variance for \mu can be improved over the MLE
            // result by plugging the values for \gamma_0 and \gamma_1 into the
            // exact formula.

            SymmetricMatrix covariances = new SymmetricMatrix(3);
            if (minimum.Curvature > 0.0) {
                covariances[0, 0] = sigma2 / minimum.Curvature;
            } else {
                covariances[0, 0] = MoreMath.Sqr(1.0 - beta * beta) / n;
            }
            covariances[1, 1] = sigma2 * (MoreMath.Sqr(1.0 + beta) - 2.0 * beta / n) / n;
            covariances[2, 2] = sigma2 / 2.0 / n;

            TimeSeries residuals = new TimeSeries();
            double u1 = 0.0;
            for (int i = 0; i < data.Count; i++) {
                double u0 = data[i] - m - beta * u1;
                residuals.Add(u0);
                u1 = u0;
            };
            TestResult test = residuals.LjungBoxTest();

            FitResult result = new FitResult(
                parameters, covariances, test
            );

            return (result);

        }

        /// <summary>
        /// Fits an AR(1) model to the time series.
        /// </summary>
        /// <returns>The fit with parameters lag-1 coefficient, mean, and standard deviation.</returns>
        public FitResult FitToAR1() {

            // AR1 model is
            //   (x_t - \mu) = \alpha (x_{t-1} - \mu) + u_{t}
            // where u_{t} \sim N(0, \sigma) are IID

            // It's easy to show
            //   m = E(x_t) = \mu
            //   c_0 = V(x_t) = E((x_t - m)^2) = \frac{\sigma^2}{1 - \alpha^2}
            //   c_1 = V(x_t, x_t-1) = E((x_t - m)(x_{t-1} - m)) = \alpha c_0
            // which gives a way to get paramters via the method of moments. In particular,
            //   \alpha = c_1 / c_0

            // For maximum likelyhood estimation (MLE), we need
            //   \log L = -\frac{1}{2} \sum_i \left[ \log (2 \pi \sigma^2)
            //            + \left[\frac{(x_i - \mu) - \alpha (x_{i-1} - \mu)}{\sigma}\right]^2

            // Treatment of the 1st value a bit subtle. We could treat it as normally distributed as
            // per m and c_0, but then it enters differently than all other values, which significantly
            // complicates the equations. Or we could regard it as given and compute log likelihood
            // conditional on it; then all values enter in the same way, but sum begins with the second
            // value. We do the latter, which is called the "conditional MLE" in the literature.

            // Differentiating once
            //   \frac{\partial L}{\partial \alpha} = \frac{1}{\sigma^2} \sum_i ( x_i - \mu - \alpha x_{i-1} ) x_{i-1}
            //   \frac{\partial L}{\partial \mu} = \frac{1}{\sigma^2} \sum_i ( x_i - \mu - \alpha x_{i-1} )
            //   \frac{\partial L}{\partial \sigma} = \sum_i \left[ \frac{(x_i - \mu - \alpha x_{i-1})^2}{\sigma^3} - \frac{1}{\sigma} \right]
            // Set equal to zero to get equations for \mu, \alpha, \sigma. First two give a 2X2 system
            // that can be solved for \mu and \alpha, then thrid for \sigma. The third equation just says
            //   \sum_i (x_i - \mu - \alpha x_{i-1})^2 = \sum_i \sigma^2
            // that \sigma is the rms of residuals. If we play a little fast and loose with index ranges
            // (e.g. ingoring difference between quantities computed over first n-1 and last n-1 values),
            // then the other two give the same results as from the method of moments.

            // Differentiating twice
            //   \frac{\partial^2 L}{\partial \alpha^2} = \frac{-1}{\sigma^2} \sum_i x_{i-1} x_{i-1} = \frac{-n (c_0 + m^2)}{\sigma^2}
            //   \frac{\partial^2 L}{\partial \mu^2} = \frac{-1}{\sigma^2} \sum_i 1 = \frac{-n}{\sigma^2}
            //   \frac{\partial^2 L}{\partial \sigma^2} = \frac{-2 n}{\sigma^2}
            //  Mixed derivatives vanish because of the first derivative conditions.

            if (data.Count < 4) throw new InsufficientDataException();

            int n = data.Count;

            // compute mean, variance, lag-1 autocorrelation

            double m = data.Mean;

            double c0 = 0.0;
            for (int i = 1; i < data.Count; i++) {
                c0 += MoreMath.Sqr(data[i] - m);
            }

            double c1 = 0.0;
            for (int i = 1; i < data.Count; i++) {
                c1 += (data[i] - m) * (data[i - 1] - m);
            }

            double alpha = c1 / c0;

            // This expression for alpha is guaranteed to be asymptotically unbiased by MLE,
            // but it is known to be biased at finite n, and in fact the bias is known.

            // See http://www.alexchinco.com/bias-in-time-series-regressions/ for simulations and explanation.
            // He cites Kendall, "Note on Bias in the Estimation of Autocorrelation", Biometrika (1954) 41 (3-4) 403-404

            // See Shaman and Stine, "The Bias of Autogregressive Coefficient Estimators",
            // Journal of the American Statistical Association (1988) Vol. 83, No. 403, pp. 842-848
            // (http://www-stat.wharton.upenn.edu/~steele/Courses/956/Resource/YWSourceFiles/ShamanStine88.pdf)
            // for derivation and formulas for AR(1)-AR(6).

            // For AR(1), MLE systematically underestimates alpha:
            //   \hat{\alpha} = \alpha - \frac{1 + 3 \alpha}{n}
            // I have confrimed the accuracy of this formula via my own simulations.

            alpha = alpha + (1.0 + 3.0 * alpha) / n;

            double sigma2 = 0.0;
            TimeSeries residuals = new TimeSeries();
            for (int i = 1; i < data.Count; i++) {
                double r = (data[i] - m) - alpha * (data[i - 1] - m);
                residuals.Add(r);
                sigma2 += MoreMath.Sqr(r);
            }
            sigma2 = sigma2 / (data.Count - 3);

            // Solution to MLE says denominator is n-1, but (i) Fuller says to use n-3,
            // (ii) simulations show n-3 is a better estimate, (iii) n-3 makes intuitive
            // sense because there are 3 parameters. I would prefer a more rigorous
            // argument, but that's good enough for now.

            // The formulas for the variances of alpha and sigma follow straightforwardly from
            // the second derivatives of the likelyhood function. For the variance of the mean
            // we use the exact formula with the \gamma_k for an AR(1) model with the fitted
            // alpha. After quite a bit of manipulation, that is
            //   v = \frac{\sigma^2}{(1-\alpha)^2} \left[ 1 -
            //         \frac{2\alpha}{n} \frac{1 - \alpha^n}{1 - \alpha^2} \right]
            // which gives a finite-n correction to the MLE result. Near \alpha \approx \pm 1,
            // we should use a series expansion to preserve accuracy.

            double[] parameters = new double[] { alpha, m, Math.Sqrt(sigma2) };

            SymmetricMatrix covariances = new SymmetricMatrix(3);
            covariances[0, 0] = (1.0 - alpha * alpha) / n;
            covariances[1, 1] = sigma2 / MoreMath.Sqr(1.0 - alpha) * (1.0 - 2.0 * alpha * (1.0 - MoreMath.Pow(alpha, n)) / (1.0 - alpha * alpha) / n) / n;
            covariances[2, 2] = sigma2 / 2.0 / n;

            TestResult test = residuals.LjungBoxTest();

             return (new FitResult(
                parameters,
                covariances,
                test
            ));   

        }

        /// <summary>
        /// Recomputes the time series as the differences between sequential values of the original series.
        /// </summary>
        /// <remarks>
        /// <para>Differencing decreases the number of values in the series by one.</para>
        /// </remarks>
        /// <seealso cref="Integrate"/>
        public void Difference () {

            if (data.Count < 2) throw new InsufficientDataException();

            SampleStorage newData = new SampleStorage();
            for (int i = 1; i < data.Count; i++) {
                newData.Add(data[i] - data[i - 1]);
            }

            data = newData;

        }

        /// <summary>
        /// Recomputes the time series as the sums of seqential values of the original series.
        /// </summary>
        /// <param name="c">The constant to be added to the first value.</param>
        /// <seealso cref="Difference"/>
        public void Integrate (double c) {

            if (data.Count < 1) throw new InsufficientDataException();

            SampleStorage newData = new SampleStorage();
            newData.Add(data[0] + c);
            for (int i = 1; i < data.Count; i++) {
                newData.Add(data[i] + newData[i - 1]);
            }

            data = newData;
        }

        /// <summary>
        /// Performs a Ljung-Box test for non-correlation.
        /// </summary>
        /// <returns>The result of the test.</returns>
        public TestResult LjungBoxTest () {
            if (data.Count < 2) throw new InsufficientDataException();
            int kMax = (int) Math.Floor(Math.Sqrt(data.Count));
            return (LjungBoxTest(kMax));
        }

        /// <summary>
        /// Performs a Ljung-Box test for non-correlation with the given number of lags.
        /// </summary>
        /// <param name="kMax">The number of lag times to test.</param>
        /// <returns>The result of the test.</returns>
        public TestResult LjungBoxTest (int kMax) {

            if (kMax < 1) throw new ArgumentOutOfRangeException("kMax");
            if (kMax >= data.Count) throw new InsufficientDataException();

            double[] c = Autocovariance();

            double Q = 0.0;
            for (int k = 1; k <= kMax; k++) {
                Q += MoreMath.Sqr(c[k] / c[0]) / (c.Length - k);
            }
            Q = c.Length * (c.Length + 2) * Q;

            return (new TestResult("Q", Q, TestType.RightTailed, new ChiSquaredDistribution(kMax)));

        }

    }

    /// <summary>
    /// Contains estimates of the moments of the population from which a time series is drawn.
    /// </summary>
    /// <remarks>
    /// <para>This class is returned by the method <see cref="TimeSeries.PopulationStatistics()"/>;
    /// see the documentation of that method for an explanation of its purpose.</para>
    /// </remarks>
    public class TimeSeriesPopulationStatistics {

        internal TimeSeriesPopulationStatistics (double m, double v, double[] g) {
            this.m = m;
            this.v = v;
            this.g = g;
        }

        private double m, v;

        private double[] g;

        /// <summary>
        /// Gets the best estimate, with uncertainty, of the mean of the time series.
        /// </summary>
        public UncertainValue Mean {
            get {
                return (new UncertainValue(m, Math.Sqrt(v)));
            }
        }

        /// <summary>
        /// Returns the autocovariance for the given lag.
        /// </summary>
        /// <param name="k">The lag index, which must lie between 0 and n-1.</param>
        /// <returns>The best estimate, with uncertainty, of the lag-<paramref name="k"/>
        /// autocovariance of the population.</returns>
        /// <remarks>
        /// <para>Note that, unlike the sample autocovariance (<see cref="TimeSeries.Autocovariance()"/>),
        /// and unlike the true population autocovariance, the estimated population autocovariance
        /// is not guaranteed to be positive definite when expressed as a matrix C<sub>i,j</sub>
        /// = c[|i-j|].</para>
        /// </remarks>
        public UncertainValue Autocovariance (int k) {
            if ((k < 0) || (k >= g.Length)) throw new ArgumentOutOfRangeException("k");

            // The task here is to estimate the variance of our estimate of g_k

            // Fuller shows that 
            //   V(c_k) = \frac{1}{n-k} \sum{j=-\infty}^{+\infty} ( g_j^2 + g_{k + j} g_{k - k} ) 
            // plus a term involving fourth cumulants that vanishes for normal errors.

            // For k = 0, this specializes to
            //   V(c_0) = \frac{2}{n} \left[ g_0^2 + 2 \sum{j=1}{n-1} g_j^2 \right]
            // For k > 0, the second term is usually smaller than the first, since
            // the signs tend to be mixed. Ignoring all cross terms we would get
            //   V(c_k) = \frac{1}{n - k} \left[ g_0^2 + g_k^2 + 2 \sum_{j=1}{n-1} g_j^2 \right]

            // It's problematic to apply this formula naively using the estimated g's.
            // The reason is that even vanishing g's will have small values due to
            // noise, and when we sum their squares they will each contribute
            // positively, and there are quite a lot of them, so we will tend
            // to significantly overestimate V(c_k).

            // I tried to deal with this by subtracting off the value we would expect if
            // all higher g's were zero, but this often gave negative values since
            // there is a lot of variance in the g's.

            // For the moment, I am applying the same cutoff I am using for the g_k
            // estimator.

            double e0 = MoreMath.Sqr(g[0]);
            for (int i = 1; i < g.Length / 3; i++) {
                e0 += 2.0 * MoreMath.Sqr(g[i]);
            }

            double e1 = MoreMath.Sqr(g[k]);
            int im = k;
            int ip = k;
            for (int j = 1; j < g.Length / 3; j++) {
                im--;
                if (im < 0) im += g.Length;
                ip++;
                if (ip >= g.Length) ip -= g.Length;
                e1 += 2.0 * g[im] * g[ip];
            }

            return (new UncertainValue(g[k], Math.Sqrt((e0 + e1) / (g.Length - k))));
        }

        /// <summary>
        /// Returns the autocorrelation for the given lag.
        /// </summary>
        /// <param name="k">The lag index, which must lie between 0 and n-1.</param>
        /// <returns>The best estimate, with uncertainty, of the lag-<paramref name="k"/>
        /// autocorrelation of the population.</returns>
        public UncertainValue Autocorrelation (int k) {
            if ((k < 0) || (k >= g.Length)) throw new ArgumentOutOfRangeException("k");
            return (new UncertainValue(g[k] / g[0], 0.0));
        }
    }
}
