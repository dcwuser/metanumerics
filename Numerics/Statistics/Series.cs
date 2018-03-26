using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;
using Meta.Numerics.SignalProcessing;
using Meta.Numerics.Statistics.Distributions;


namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Contains methods for the statistical analysis of time series.
    /// </summary>
    public static class Series {

        /// <summary>
        /// Computes the autocovariance of the series at the given lag.
        /// </summary>
        /// <param name="series">The data series from which to compute the autocovariance.</param>
        /// <param name="lag">The lag at which to compute the autocovariance.</param>
        /// <returns>The value of the autocovariance at the given lag.</returns>
        /// <remarks>
        /// <para>In a length-N time series, there are N-k lag-k observations. Nonetheless,
        /// the definition of the lag-k autocovariance requires division by N, not N-k. This
        /// counterintuitive convention ensures that the autocovariance has desirable
        /// positive definiteness properties and agrees with the computation via FFT.</para>
        /// <para>The computation of an autocovariance via this method is O(N). If
        /// you need to compute more than a handfull of autocovariances, it is
        /// more efficient to call the <see cref="Autocovariance(IReadOnlyList{double})"/>, which
        /// computes all of them in O(N log N).</para>
        /// <para>While the sample autocovariance does converge to the population
        /// autocovariance in the large-N limit, this convergence is very slow. If you
        /// want an estimate of the population autocovariance, use
        /// <see cref="SeriesPopulationStatistics(IReadOnlyList{double})"/>
        /// to obtain a much better estimate.</para>
        /// </remarks>
        /// <seealso cref="Autocovariance(IReadOnlyList{double})"/>
        public static double Autocovariance (this IReadOnlyList<double> series, int lag) {

            if (series == null) throw new ArgumentNullException(nameof(series));
            if (lag >= series.Count) throw new ArgumentOutOfRangeException(nameof(lag));
            if (lag < 0) lag = -lag;

            double mean = series.Mean();

            double sum = 0.0;
            for (int i = lag; i < series.Count; i++) {
                sum += (series[i] - mean) * (series[i - lag] - mean);
            }
            return (sum / series.Count);

            // It appears to be standard practice to divide by n instead of n-k, even though
            // (1) there are n-k terms, not n, (2) dividing by n gives a biased estimator
            // of the population mean, while dividing by n-k yields an unbiased estimator.
            // I really don't undertand this practice, but it does appear to be necessary 
            // for agreement with the FFT computation.

        }

        /// <summary>
        /// Computes the autocovariance for all lags.
        /// </summary>
        /// <param name="series">The data series.</param>
        /// <returns>An array of autocovariance values, with the array index equal to the lag index.</returns>
        /// <remarks>
        /// <para>The computation of the autocovariance for a given lag is an O(N) operation.
        /// Naively, the computation of the autocovariance for all N possible lags is therefore an
        /// O(N^2) operation; this is in fact the cost of N invocations of
        /// <see cref="Autocovariance(IReadOnlyList{double},int)"/>.
        /// However, using Fourier techniques, it is possible to simultaneously compute
        /// the autocovariance for all possible lags in O(N log N) operations. This method
        /// uses the Fourier technique and should be called if you require the autocovariance
        /// for more than a handfull of lag values.</para>
        /// </remarks>
        /// <seealso cref="Autocovariance(IReadOnlyList{double},int)"/>
        public static double[] Autocovariance (this IReadOnlyList<double> series) {

            if (series == null) throw new ArgumentNullException(nameof(series));

            // To avoid aliasing, we must zero-pad with n zeros.
            Complex[] s = new Complex[2 * series.Count];
            for (int i = 0; i < series.Count; i++) {
                s[i] = series[i] - series.Mean();
            }

            FourierTransformer fft = new FourierTransformer(s.Length);
            Complex[] t = fft.Transform(s);

            for (int i = 0; i < t.Length; i++) {
                t[i] = MoreMath.Sqr(t[i].Re) + MoreMath.Sqr(t[i].Im);
            }

            Complex[] u = fft.InverseTransform(t);

            double[] c = new double[series.Count];
            for (int i = 0; i < c.Length; i++) {
                c[i] = u[i].Re / series.Count;
            }

            return (c);
        }

        /// <summary>
        /// Computes the power spectrum of the time series.
        /// </summary>
        /// <param name="series">The data series.</param>
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
        public static double[] PowerSpectrum (this IReadOnlyList<double> series) {

            if (series == null) throw new ArgumentNullException(nameof(series));

            Complex[] s = new Complex[series.Count];
            for (int i = 0; i < series.Count; i++) {
                s[i] = series[i];
            }

            FourierTransformer fft = new FourierTransformer(s.Length);
            Complex[] t = fft.Transform(s);

            double[] u = new double[series.Count / 2];
            u[0] = MoreMath.Sqr(t[0].Re) + MoreMath.Sqr(t[0].Im);
            for (int i = 1; i < u.Length; i++) {
                u[i] = MoreMath.Sqr(t[i].Re) + MoreMath.Sqr(t[i].Im)
                    + MoreMath.Sqr(t[series.Count - i].Re) + MoreMath.Sqr(t[series.Count - i].Im);
            }

            return (u);

        }

        /// <summary>
        /// Computes estimates for the moments of the population from which the time series is drawn.
        /// </summary>
        /// <param name="series">The data series.</param>
        /// <returns>A collection of population statistics.</returns>
        /// <remarks>
        /// <para>Just as is the case for the variance of a sample (<see cref="Sample.Variance"/>),
        /// the sample autocovariances of a time series
        /// are not unbiased estimates of the autocovariances of the population from which
        /// the series is drawn. Additional computations must be performed to calculate
        /// unbiased estimates and error estimates for the time series mean and autocovariances.</para>
        /// <para>Unlike the case of a sample of uncorrelated values, these computations are complicated
        /// and the same computation is relevant for all moments. Therefore, the computation is performed
        /// once by this method, which returns a object from which all population statistics can be
        /// quickly queried.</para>
        /// </remarks>
        public static TimeSeriesPopulationStatistics SeriesPopulationStatistics (this IReadOnlyList<double> series) {

            if (series == null) throw new ArgumentNullException(nameof(series));

            double m = series.Mean();

            double v;
            double[] g;
            PopulationStatistics(series, out g, out v);

            return (new TimeSeriesPopulationStatistics(m, v, g));

        }

        private static void PopulationStatistics (IReadOnlyList<double> series, out double[] g, out double v) {

            Debug.Assert(series != null);
            if (series.Count < 4) throw new InsufficientDataException();

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
            // Here v is the same quantity as the variance of the mean.

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

            int n = series.Count;
            int jmax = (int) Math.Floor(Math.Sqrt(n));

            double[] c = Autocovariance(series);
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

            // To deal with v negative (or unrealistically small), don't let it
            // get smaller than the first truncated term. This is a very
            // ad hoc solution; we should try to do better.

            double vMin = 2.0 / n * (n - jmax) / n * Math.Abs(g[jmax]);
            if (v < vMin) {
                v = vMin;
            }

        }

        /// <summary>
        /// Fits an AR(1) model to the time series.
        /// </summary>
        /// <param name="series">The data series.</param>
        /// <returns>The fit with parameters lag-1 coefficient, mean, and standard deviation.</returns>
        public static AR1FitResult FitToAR1 (this IReadOnlyList<double> series) {

            if (series == null) throw new ArgumentNullException(nameof(series));

            // AR1 model is
            //   (x_t - \mu) = \alpha (x_{t-1} - \mu) + u_{t}
            // where u_{t} \sim N(0, \sigma) are IID

            // It's easy to show
            //   m = E(x_t) = \mu
            //   c_0 = V(x_t) = E((x_t - m)^2) = \frac{\sigma^2}{1 - \alpha^2}
            //   c_1 = V(x_t, x_t-1) = E((x_t - m)(x_{t-1} - m)) = \alpha c_0
            // which gives a way to get parameters via the method of moments. In particular,
            //   \alpha = c_1 / c_0

            // For maximum likelihood estimation (MLE), we need
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
            // Set equal to zero to get equations for \mu, \alpha, \sigma. The first two give a 2X2 system
            // that can be solved for \mu and \alpha. The third equation just says
            //   \sum_i (x_i - \mu - \alpha x_{i-1})^2 = \sum_i \sigma^2
            // that \sigma is the rms of residuals. If we play a little fast and loose with index ranges
            // (e.g. ignoring difference between quantities computed over first n-1 and last n-1 values),
            // then the other two give the same results as from the method of moments.

            // Differentiating twice
            //   \frac{\partial^2 L}{\partial \alpha^2} = \frac{-1}{\sigma^2} \sum_i x_{i-1} x_{i-1} = \frac{-n (c_0 + m^2)}{\sigma^2}
            //   \frac{\partial^2 L}{\partial \mu^2} = \frac{-1}{\sigma^2} \sum_i 1 = \frac{-n}{\sigma^2}
            //   \frac{\partial^2 L}{\partial \sigma^2} = \frac{-2 n}{\sigma^2}
            //  Mixed derivatives all vanish because of the first derivative conditions.

            int n = series.Count;
            if (n < 4) throw new InsufficientDataException();

            // compute mean, variance, lag-1 autocorrelation

            double m = series.Mean();

            double c0 = 0.0;
            for (int i = 1; i < n; i++) {
                c0 += MoreMath.Sqr(series[i] - m);
            }

            double c1 = 0.0;
            for (int i = 1; i < n; i++) {
                c1 += (series[i] - m) * (series[i - 1] - m);
            }

            double alpha = c1 / c0;
            Debug.Assert(Math.Abs(alpha) <= 1.0);

            // This expression for alpha is guaranteed to be asymptotically unbiased by MLE,
            // but it is known to be biased at finite n, and in fact the bias is known.

            // See http://www.alexchinco.com/bias-in-time-series-regressions/ for simulations and explanation.
            // He cites Kendall, "Note on Bias in the Estimation of Autocorrelation", Biometrika (1954) 41 (3-4) 403-404

            // See Shaman and Stine, "The Bias of Autogregressive Coefficient Estimators",
            // Journal of the American Statistical Association (1988) Vol. 83, No. 403, pp. 842-848
            // (http://www-stat.wharton.upenn.edu/~steele/Courses/956/Resource/YWSourceFiles/ShamanStine88.pdf)
            // for derivation and formulas for AR(1)-AR(6).

            // For AR(1), MLE systematically underestimates alpha:
            //   E(\hat{\alpha}) = \alpha - \frac{1 + 3 \alpha}{n}
            // I have confirmed the accuracy of this formula via my own simulations.
            // One problem with trying to correct for this is that it can push alpha over one,
            // which is not only impossible but produces negative variance estimates.
            // So we don't accept such large alphas.
            alpha = Math.Min(alpha + (1.0 + 3.0 * alpha) / n, 1.0);
            Debug.Assert(alpha <= 1.0);

            double sigma2 = 0.0;
            List<double> residuals = new List<double>(n);
            for (int i = 1; i < n; i++) {
                double r = (series[i] - m) - alpha * (series[i - 1] - m);
                residuals.Add(r);
                sigma2 += MoreMath.Sqr(r);
            }
            sigma2 = sigma2 / (n - 3);

            // Solution to MLE says denominator is n-1, but (i) Fuller says to use n-3,
            // (ii) simulations show n-3 is a better estimate, (iii) n-3 makes intuitive
            // sense because there are 3 parameters. I would prefer a more rigorous
            // argument, but that's good enough for now.

            // The formulas for the variances of alpha and sigma follow straightforwardly from
            // the second derivatives of the likelihood function. For the variance of the mean
            // we use the exact formula with the \gamma_k for an AR(1) model with the fitted
            // alpha. After quite a bit of manipulation, that is
            //   v = \frac{\sigma^2}{(1-\alpha)^2} \left[ 1 -
            //         \frac{2\alpha}{n} \frac{1 - \alpha^n}{1 - \alpha^2} \right]
            // which gives a finite-n correction to the MLE result. Near \alpha \approx \pm 1,
            // we should use a series expansion to preserve accuracy.

            double varAlpha = (1.0 - alpha * alpha) / n;
            double varMu = sigma2 / MoreMath.Sqr(1.0 - alpha) * (1.0 - 2.0 * alpha * (1.0 - MoreMath.Pow(alpha, n)) / (1.0 - alpha * alpha) / n) / n;
            double varSigma = sigma2 / 2.0 / n;

            return (new AR1FitResult(
                new UncertainValue(m, Math.Sqrt(varMu)),
                new UncertainValue(alpha, Math.Sqrt(varAlpha)),
                new UncertainValue(Math.Sqrt(sigma2), Math.Sqrt(varSigma)),
                residuals
            ));
        }

        /// <summary>
        /// Fits an MA(1) model to the time series.
        /// </summary>
        /// <param name="series">The data series.</param>
        /// <returns>The fit, with parameters lag-1 coefficient, mean, and standard deviation.</returns>
        public static MA1FitResult FitToMA1 (this IReadOnlyList<double> series) {

            if (series == null) throw new ArgumentNullException(nameof(series));
            if (series.Count < 4) throw new InsufficientDataException();

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
            // but they do illustrate the basic requirement that g_1 < 1/2.

            // The MLE estimator 

            int n = series.Count;

            double m = series.Mean();

            Func<double, double> fnc = (double theta) => {

                double s = 0.0;
                double uPrevious = 0.0;
                for (int i = 0; i < series.Count; i++) {
                    double u = series[i] - m - theta * uPrevious;
                    s += u * u;
                    uPrevious = u;
                };
                return (s);
            };

            Extremum minimum = FunctionMath.FindMinimum(fnc, Interval.FromEndpoints(-1.0, 1.0));

            double beta = minimum.Location;
            double sigma2 = minimum.Value / (n - 3);

            // While there is significant evidence that the MLE value for \beta is biased
            // for small-n, I know of no analytic correction.

            //double[] parameters = new double[] { beta, m, Math.Sqrt(sigma2) };

            // The calculation of the variance for \mu can be improved over the MLE
            // result by plugging the values for \gamma_0 and \gamma_1 into the
            // exact formula.

            double varBeta;
            if (minimum.Curvature > 0.0) {
                varBeta = sigma2 / minimum.Curvature;
            } else {
                varBeta = MoreMath.Sqr(1.0 - beta * beta) / n;
            }
            double varMu = sigma2 * (MoreMath.Sqr(1.0 + beta) - 2.0 * beta / n) / n;
            double varSigma = sigma2 / 2.0 / n;

            List<double> residuals = new List<double>(n);
            double u1 = 0.0;
            for (int i = 0; i < n; i++) {
                double u0 = series[i] - m - beta * u1;
                residuals.Add(u0);
                u1 = u0;
            };

            return (new MA1FitResult(
                new UncertainValue(m, Math.Sqrt(varMu)),
                new UncertainValue(beta, Math.Sqrt(varBeta)),
                new UncertainValue(Math.Sqrt(sigma2), Math.Sqrt(varSigma)),
                residuals
           ));
        }
        /*
        /// <summary>
        /// Fits an AR(1) model to the time series.
        /// </summary>
        /// <returns>The fit with parameters lag-1 coefficient, mean, and standard deviation.</returns>
        public FitResult FitToAR1 () {

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
        */

        /// <summary>
        /// Performs a Ljung-Box test for non-correlation.
        /// </summary>
        /// <param name="series">The data series.</param>
        /// <returns>The result of the test.</returns>
        public static TestResult LjungBoxTest (this IReadOnlyList<double> series) {
            if (series == null) throw new ArgumentNullException(nameof(series));
            if (series.Count < 2) throw new InsufficientDataException();
            int kMax = (int) Math.Floor(Math.Sqrt(series.Count));
            return (LjungBoxTest(series, kMax));
        }

        /// <summary>
        /// Performs a Ljung-Box test for non-correlation with the given number of lags.
        /// </summary>
        /// <param name="series">The data series.</param>
        /// <param name="kMax">The number of lag times to test.</param>
        /// <returns>The result of the test.</returns>
        /// <remarks>
        /// <para>If you are unsure how many lag times to test, use the default
        /// overload <see cref="LjungBoxTest(IReadOnlyList{double})"/>.</para>
        /// </remarks>
        public static TestResult LjungBoxTest (this IReadOnlyList<double> series, int kMax) {
            if (series == null) throw new ArgumentNullException(nameof(series));
            if (kMax < 1) throw new ArgumentOutOfRangeException(nameof(kMax));
            if (kMax >= series.Count) throw new InsufficientDataException();

            double[] c = Autocovariance(series);

            double Q = 0.0;
            for (int k = 1; k <= kMax; k++) {
                Q += MoreMath.Sqr(c[k] / c[0]) / (c.Length - k);
            }
            Q = c.Length * (c.Length + 2) * Q;

            return (new TestResult("Q", Q, TestType.RightTailed, new ChiSquaredDistribution(kMax)));
        }

    }

}
