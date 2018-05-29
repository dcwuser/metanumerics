using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;


namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents an ordered series of data points.
    /// </summary>
    public sealed class TimeSeries : IEnumerable, IEnumerable<double>, IReadOnlyCollection<double>, IList<double> {

        /// <summary>
        /// Initializes a new, empty time series.
        /// </summary>
        public TimeSeries () {

        }

        /// <summary>
        /// Initializes a new time series with the given values.
        /// </summary>
        /// <param name="values">The first values in the time series.</param>
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
            if (values == null) throw new ArgumentNullException(nameof(values));
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
            if ((index < 0) || (index >= data.Count)) throw new ArgumentOutOfRangeException(nameof(index));
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
            return (Series.Autocovariance(data, lag));
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
            return (Series.Autocovariance(data));
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
            return (Series.PowerSpectrum(data));
        }

        /// <summary>
        /// Computes estimates for the moments of the population from which the time series is drawn.
        /// </summary>
        /// <returns>A collection of population statistics.</returns>
        /// <remarks>
        /// <para>Just as is the case for the variance of a simple <see cref="Sample"/>,
        /// the sample autocovariances of a time series
        /// are not unbiased estimates of the autocovariances of the population from which
        /// the series is drawn. Additional computations must be performed to calculate
        /// unbiased estimates and error estimates for the time series mean and autocovariances.</para>
        /// <para>Unlike the case of a sime <see cref="Sample"/>, these computations are complicated
        /// and the same computation is relevant for all moments. Therefore, the computation is performed
        /// once by this method, which returns a object from which all population statistics can be
        /// quickly queried.</para>
        /// </remarks>
        public TimeSeriesPopulationStatistics PopulationStatistics () {
            return (Series.SeriesPopulationStatistics(data));
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
        public MA1FitResult FitToMA1() {
            return (Series.FitToMA1(data));
        }

        /// <summary>
        /// Fits an AR(1) model to the time series.
        /// </summary>
        /// <returns>The fit with parameters lag-1 coefficient, mean, and standard deviation.</returns>
        public AR1FitResult FitToAR1() {
            return (data.FitToAR1());
        }

        /// <summary>
        /// Re-computes the time series as the differences between sequential values of the original series.
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
        /// Re-computes the time series as the sums of sequential values of the original series.
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
            return (Series.LjungBoxTest(data, kMax));
        }

    }

}
