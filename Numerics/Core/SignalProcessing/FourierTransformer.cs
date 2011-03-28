using System;
using System.Collections.Generic;


namespace Meta.Numerics.SignalProcessing {

    /// <summary>
    /// Specifies the normalization convention to be used in a forward Fourier transform.
    /// </summary>
    /// <remarks>
    /// <para>The most common convention in signal processing applications is <see cref="FourierNormalization.None"/>.</para>
    /// </remarks>
    public enum FourierNormalization {

        /// <summary>
        /// The series is not normalized.
        /// </summary>
        None,
        
        /// <summary>
        /// The series is multiplied by 1/N<sup>1/2</sup>.
        /// </summary>
        Unitary,
        
        /// <summary>
        /// The series is multiplied by 1/N.
        /// </summary>
        Inverse
    }

    /// <summary>
    /// Specifies the sign convention to be used in the exponent of a forward Fourier transform.
    /// </summary>
    /// <remarks>
    /// <para>The most common convention in signal processing applications is <see cref="FourierSign.Negative"/>.</para>
    /// </remarks>
    public enum FourierSign {

        /// <summary>
        /// The exponent has positive imaginary values.
        /// </summary>
        Positive,

        /// <summary>
        /// The exponent has negative imaginary values.
        /// </summary>
        Negative
    }

    /// <summary>
    /// An engine for performing Fourier transforms on complex series.
    /// </summary>
    /// <remarks>
    /// <para>A Fourier transform decomposes a function into a sum of different frequency components. This is
    /// useful for a wide array of applications.</para>
    /// <para>Mathematically, the DFT is an N-dimensional linear transfromation
    /// with coefficients that are the Nth complex roots of unity.</para>
    /// <img src="../images/Fourier.png" />
    /// <para>An instance of the FourierTransformer class performs DFTs on series of a particular length,
    /// given by its <see cref="FourierTransformer.Length"/> property. This specialization allows certain parts of the DFT
    /// calculation, which are indepdent of the transformed series but dependent on the length of the series,
    /// to be performed only once and then used for all transforms of that length. This saves time and improves
    /// performance. If you need to perform DFTs on series with different lengths, simply create a seperate instance
    /// of the FourierTransform class for each required length.</para>
    /// <para>Many simple DFT implementations require that the series length be a power of two (2, 4, 8, 16, etc.).
    /// Meta.Numerics supports DFTs of any length. Our DFT implementation is fast -- order O(N log N) -- for any length composed
    /// of small prime factors.
    /// This includes not only all powers of two, but also all powers of ten (10, 100, 1000, etc.) and many other
    /// lengths (e.g. 3600 = 2<sup>4</sup> 3<sup>2</sup> 5<sup>2</sup>). Our DFT implementation is slow -- O(N<sup>2</sup>) -- for
    /// lengths composed of large prime factors (e.g. 3571, which is prime).</para>
    /// </remarks>
    /// <example>
    /// <para>The following code performs a simple DFT and then inverts it to re-obtain the original data.</para>
    /// <code lang="c#">
    /// // Create a Fourier transformer for length-6 series
    /// FourierTransformer ft = new FourierTransformer(6);
    /// // Create a length-6 series and transform it
    /// Complex[] x = new Complex[] { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
    /// Complex[] xt = ft.Transform(x);
    /// // Re-use the same transformer to transform a different  series
    /// Complex[] y = new Complex[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    /// Complex[] yt = ft.Transform(y);
    /// // Transform them back
    /// Complex[] xtt = ft.InverseTransform(xt);
    /// Complex[] ytt = ft.InverseTransform(yt);
    /// </code>
    /// </example>
    /// <seealso href="http://en.wikipedia.org/wiki/Discrete-time_Fourier_transform"/>
    public sealed class FourierTransformer {

        /// <summary>
        /// Initializes a new instance of the Fourier transformer.
        /// </summary>
        /// <param name="size">The series length of the transformer, which must be positive.</param>
        public FourierTransformer (int size)
            : this(size, FourierSign.Negative, FourierNormalization.None) {
        }

        /// <summary>
        /// Initializes a new instance of the Fourier transformer with the given sign and normalization conventions.
        /// </summary>
        /// <param name="size">The series length of the transformer, which must be positive.</param>
        /// <param name="signConvention">The sign convention of the transformer.</param>
        /// <param name="normalizationConvention">The normalization convention of the transformer.</param>
        /// <remarks>
        /// <para>There are multiple conventions for both the sign of the exponent and the overall normalization of
        /// Fourier transforms. The default conventions for some widely used software packages are summarized in the following
        /// table.</para>
        /// <table>
        ///     <tr><th>Software</th><th>Sign</th><th>Normalization</th></tr>
        ///     <tr><td>Meta.Numerics</td><td><see cref="FourierSign.Negative"/></td><td><see cref="FourierNormalization.None"/></td></tr>
        ///     <tr><td>Matlab</td><td><see cref="FourierSign.Negative"/></td><td><see cref="FourierNormalization.None"/></td></tr>
        ///     <tr><td>Mathmatica</td><td><see cref="FourierSign.Positive"/></td><td><see cref="FourierNormalization.Unitary"/></td></tr>
        ///     <tr><td>Numerical Recipies</td><td><see cref="FourierSign.Positive"/></td><td><see cref="FourierNormalization.None"/></td></tr>
        /// </table>
        /// </remarks>
        public FourierTransformer (int size, FourierSign signConvention, FourierNormalization normalizationConvention) {
            if (size < 1) throw new ArgumentOutOfRangeException("size");
            this.size = size;
            this.factors = Factor.Factorize(size);
            this.signConvention = signConvention;
            this.normalizationConvention = normalizationConvention;
            this.roots = FourierAlgorithms.ComputeRoots(size, +1);
        }

        private int size;
        private List<Factor> factors;
        private FourierNormalization normalizationConvention;
        private FourierSign signConvention;
        private Complex[] roots;

        /// <summary>
        /// The series length for which the transformer is specialized.
        /// </summary>
        public int Length {
            get {
                return (size);
            }
        }

        /// <summary>
        /// Gets the normalization convention used by the transformer.
        /// </summary>
        public FourierNormalization NormalizationConvention {
            get {
                return (normalizationConvention);
            }
        }

        /// <summary>
        /// Gets the normalization convention used by the transformer.
        /// </summary>
        public FourierSign SignConvention {
            get {
                return (signConvention);
            }
        }

        private int GetSign () {
            if (signConvention == FourierSign.Positive) {
                return (+1);
            } else {
                return (-1);
            }
        }

        private static void Normalize (Complex[] x, double f) {
            for (int i = 0; i < x.Length; i++) {
                x[i] = new Complex(f * x[i].Re, f * x[i].Im);
            }
        }


        /// <summary>
        /// Computes the Fourier transform of the given series.
        /// </summary>
        /// <param name="values">The series to transform.</param>
        /// <returns>The discrete Fourier transform of the series.</returns>
        public Complex[] Transform (IList<Complex> values) {
            if (values == null) throw new ArgumentNullException("values");
            if (values.Count != size) throw new DimensionMismatchException();

            // copy the original values into a new array
            Complex[] x = new Complex[size];
            values.CopyTo(x, 0);

            // normalize the copy appropriately
            if (normalizationConvention == FourierNormalization.Unitary) {
                Normalize(x, 1.0 / Math.Sqrt(size));
            } else if (normalizationConvention == FourierNormalization.Inverse) {
                Normalize(x, 1.0 / size);
            }

            // create a scratch array
            Complex[] y = new Complex[size];

            // do the FFT
            FourierAlgorithms.Fft(values.Count, factors, ref x, ref y, roots, GetSign());

            return (x);

        }

        /// <summary>
        /// Computes the inverse Fourier transform of the given series.
        /// </summary>
        /// <param name="values">The series to invert.</param>
        /// <returns>The inverse discrete Fourier transform of the series.</returns>
        public Complex[] InverseTransform (IList<Complex> values) {
            if (values == null) throw new ArgumentNullException("values");
            if (values.Count != size) throw new DimensionMismatchException();

            // copy the original values into a new array
            Complex[] x = new Complex[size];
            values.CopyTo(x, 0);

            // normalize the copy appropriately
            if (normalizationConvention == FourierNormalization.None) {
                Normalize(x, 1.0 / size);
            } else if (normalizationConvention == FourierNormalization.Unitary) {
                Normalize(x, 1.0 / Math.Sqrt(size));
            }

            // create a scratch array
            Complex[] y = new Complex[size];

            // do the FFT
            FourierAlgorithms.Fft(values.Count, factors, ref x, ref y, roots, -GetSign());

            return (x);

        }

    }

}
