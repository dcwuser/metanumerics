using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Contains estimates of the moments of the population from which a time series is drawn.
    /// </summary>
    /// <remarks>
    /// <para>This class is returned by the method <see cref="Series.SeriesPopulationStatistics(IReadOnlyList{double})"/>;
    /// see the documentation of that method for an explanation of its purpose.</para>
    /// </remarks>
    public sealed class TimeSeriesPopulationStatistics {

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
        /// Returns the auto-covariance for the given lag.
        /// </summary>
        /// <param name="k">The lag index, which must lie between 0 and n-1.</param>
        /// <returns>The best estimate, with uncertainty, of the lag-<paramref name="k"/>
        /// auto-covariance of the population.</returns>
        /// <remarks>
        /// <para>Note that, unlike the sample auto-covariance (<see cref="TimeSeries.Autocovariance()"/>),
        /// and unlike the true population auto-covariance, the estimated population auto-covariance
        /// is not guaranteed to be positive definite when expressed as a matrix C<sub>i,j</sub>
        /// = c[|i-j|].</para>
        /// </remarks>
        public UncertainValue Autocovariance (int k) {
            if ((k < 0) || (k >= g.Length)) throw new ArgumentOutOfRangeException(nameof(k));

            // The task here is to estimate the variance of our estimate of g_k

            // Fuller shows that 
            //   V(c_k) = \frac{1}{n-k} \sum_{j=-\infty}^{+\infty} ( g_j^2 + g_{k + j} g_{k - k} ) 
            // plus a term involving fourth cumulants that vanishes for normal errors.

            // For k = 0, this specializes to
            //   V(c_0) = \frac{2}{n} \left[ g_0^2 + 2 \sum_{j=1}^{n-1} g_j^2 \right]
            // For k > 0, the second term is usually smaller than the first, since
            // the signs tend to be mixed. Ignoring all cross terms we would get
            //   V(c_k) = \frac{1}{n - k} \left[ g_0^2 + g_k^2 + 2 \sum_{j=1}^{n-1} g_j^2 \right]

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

        /*
        /// <summary>
        /// Returns the autocorrelation for the given lag.
        /// </summary>
        /// <param name="k">The lag index, which must lie between 0 and n-1.</param>
        /// <returns>The best estimate, with uncertainty, of the lag-<paramref name="k"/>
        /// autocorrelation of the population.</returns>
        public UncertainValue Autocorrelation (int k) {
            if ((k < 0) || (k >= g.Length)) throw new ArgumentOutOfRangeException(nameof(k));
            return (new UncertainValue(g[k] / g[0], 0.0));
        }
        // fix computation of uncertainty before making this public
        */
    }
}
