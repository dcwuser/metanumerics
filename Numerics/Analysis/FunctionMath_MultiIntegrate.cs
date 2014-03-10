using System;
using System.Collections.Generic;

using Meta.Numerics;

namespace Meta.Numerics.Functions {

    public static partial class FunctionMath {

        /// <summary>
        /// Integrates the given function over the given volume.
        /// </summary>
        /// <param name="integrand">The function to integrate, which maps R<sup>d</sup> to R.</param>
        /// <param name="volume">The box defining the volume over with to integrate.</param>
        /// <returns>The value of the integral.</returns>
        /// <remarks>
        /// <para>Note that the integration region must be a hyper-rectangle. You can integrate over regions with more complex boundaries by specifying the integration
        /// volume as a bounding hyper-rectangle that encloses your desired integration region, and returing the value 0 for the integrand outside of the desired integration
        /// region. For example, to find the volume of a unit d-sphere, you can integrate</para>
        /// <para>over the volume [-1,1]<sup>d</sup>.</para>
        /// <para>Volumes with dimension greater than 8 are not currently supported. Infinite volumes are not currently supported.</para>
        /// <h2>Accuracy and Evaluation Budget</h2>
        /// <para>By default, our multidimensional integration system targets a relative accuracy of about 10<sup>-4</sup> (0.01%) for d=2, falling gradually
        /// to about 10<sup>-1</sup> (10%) for d=8. To achieve that accuracy, it allows up to about one million (1,000,000) evaluations of the integrand for d=2, rising
        /// up to about one billion (1,000,000,000) evaluations for d=8.</para>
        /// <para>You can change the accuracy demands and evaluation budget by passing an <see cref="EvaluationSettings"/> object to the integration method. By decreasing
        /// the accuracy you require or increasing the evaluation budget, you may be able to successfully complete integrals that would fail for the default settings.</para>
        /// <h2>Hints</h2>
        /// <para>Because multidimensional integrals are so computationally expensive, you should use any known symmetries to simplify the problem.</para>
        /// </remarks>
        /// <exception cref="InvalidOperationException">An integrand of dimension greater than 8 was specified.</exception>
        /// <exception cref="NonconvergenceException">The integral could not be determined to the specified accuracy within the given evaluation budget.</exception>
        public static double Integrate (Func<IList<double>, double> integrand, IList<Interval> volume) {
            return(Integrate(integrand, volume, null));
        }

        /// <summary>
        /// Integrates the given function over the given volume.
        /// </summary>
        /// <param name="integrand">The function to integrate, which maps R<sup>d</sup> to R.</param>
        /// <param name="volume">The box defining the volume over with to integrate.</param>
        /// <param name="settings">The evaluation settings to use.</param>
        /// <returns>The value of the integral.</returns>
        public static double Integrate (Func<IList<double>, double> integrand, IList<Interval> volume, EvaluationSettings settings) {

            if (integrand == null) throw new ArgumentNullException("integrand");
            if (volume == null) throw new ArgumentNullException("volume");

            // get the dimension of the problem from the volume
            int d = volume.Count;
            if ((d < 1) || (d > sobolParameters.Length)) throw new InvalidOperationException();

            // if no settings were provided, use defaults
            if (settings == null) {
                settings = new EvaluationSettings() {
                    RelativePrecision = 1.0 / (1 << (14 - 3 * d / 2)),
                    AbsolutePrecision = RelativePrecision / 128.0,
                    EvaluationBudget = 1 << (18 + 3 * d / 2)
                };
            }

            // compute the volume of the box
            double V = 1.0;
            for (int j = 0; j < volume.Count; j++) {
                V *= volume[j].Width;
            }

            // generate the appropriate number of Sobol sequences
            SobolSequence[] s = new SobolSequence[d];
            for (int j = 0; j < d; j++) {
                SobolSequenceParameters p = sobolParameters[j];
                s[j] = new SobolSequence(p.Dimension, p.Coefficients, p.Seeds);
            }

            // create one vector to store the argument
            // we don't want to recreate the vector in the
            // heap on each iteration
            double[] x = new double[d];

            // keep track of the average sampled value
            double M = 0.0; double M_old = Double.NaN;

            // set the first convergence checkpoint to ~4000 points for d=2, increasing
            // for higher d
            int i_next = 1 << (10 + d);

            for (int i = 1; i <= settings.EvaluationBudget; i++) {

                // move to the next value in each dimension's Sobol sequence
                // and construct the argument vector 
                for (int j = 0; j < d; j++) {
                    s[j].MoveNext();
                    x[j] = volume[j].LeftEndpoint + s[j].Current * volume[j].Width;
                }

                // evaluate the function
                double y = integrand(x);

                // update the mean
                M += (y - M) / i;

                // check for convergence at regular intervals
                if (i == i_next) {

                    // estimate error as change since last check
                    double dM = 2.0 * Math.Abs(M - M_old);

                    if (dM < settings.RelativePrecision * Math.Abs(M) || Math.Abs(V) * dM < settings.AbsolutePrecision) {
                        return (V * M);
                    }

                    // no convergence, so remember current value and set next checkpoint
                    M_old = M;
                    i_next = 2 * i_next;

                    // Consider how to do this better:
                    // 1. Is there any way we can do better than doubling? Say increasing by 4/3 and then 3/2, for an average increase of
                    // 40% instead of 100%? If so, how do we estimate error at each step?
                    // 2. Can we use Romberg extrapolation? I tried using both M and M_old, assuming an error term of 1/N. This gives
                    // the extrapolated value 2 M - M_old, but experimentally this is usually a worse value than M. I should try out
                    // 3-value extrapolation.

                }

            }

            throw new NonconvergenceException();

        }

        // See http://web.maths.unsw.edu.au/~fkuo/sobol/ for lists of suggested parameters

        internal static readonly SobolSequenceParameters[] sobolParameters = new SobolSequenceParameters[] {
            new SobolSequenceParameters(1, 0, new ulong[] { 1 }),
            new SobolSequenceParameters(2, 1, new ulong[] { 1, 3 }),
            new SobolSequenceParameters(3, 1, new ulong[] { 1, 3, 1 }),
            new SobolSequenceParameters(3, 2, new ulong[] { 1, 1, 1 }),
            new SobolSequenceParameters(4, 1, new ulong[] { 1, 1, 3, 3 }),
            new SobolSequenceParameters(4, 4, new ulong[] { 1, 3, 5, 13 }),
            new SobolSequenceParameters(5, 2, new ulong[] { 1, 1, 5, 5, 17 }),
            new SobolSequenceParameters(5, 4, new ulong[] { 1, 1, 5, 5, 5 }),
            new SobolSequenceParameters(5, 7, new ulong[] { 1, 1, 7, 11, 19 }), 
            new SobolSequenceParameters(5, 11, new ulong[] { 1, 1, 5, 1, 1 }),
            new SobolSequenceParameters(5, 13, new ulong[] { 1, 1, 1, 3, 11 }),
            new SobolSequenceParameters(5, 14, new ulong[] { 1, 3, 5, 5, 31 })
        };

        
    }

    // Define a Sobol sequence

    internal class SobolSequence {

         public SobolSequence (int s, ulong a, IList<ulong> m) {

            // store the binary sequence that encodes our direction-generating recurrsion
            // we could do this with a single ulong, but the literature typically specifies
            // the degree s seperately from the coefficients of x^{s-1}, \cdots, x^2 x^1, which
            // are stored as binary digits in a = (a_{s-1} \cdots a_2 a_1)
            this.s = s;
            this.a = a;

            // now we construct the "directions"

            // The ith direction number v_i = u_i / 2^i, where u_i is an whole number with fewer
            // than i bits. 

            // we will store the v's scaled up by a constant factor 2^n, where n is the maximum i, so they
            // can be manipulated as whole binary numbers; then we will divide by 2^n once at the very end
            // when we return a floating point value
            v = new ulong[n];

            // store the scale factor 1 / 2^n
            f = 1.0 / ((double) (((ulong) 1) << n));

            // the first few "directions" are given by the initialization seeds
            // the number of seed values must equal the degree of the primitive polynomial
            if (m.Count != s) throw new InvalidOperationException();
            for (int i = 0; i < s; i++) {
                // seeds must be odd and the ith must be less than 2^i
                if ((m[i] % 2 == 0) || (m[i] > (((ulong) 1) << (i + 1)))) throw new InvalidOperationException();
                // the values are m_i / 2^i, but remember we then scale up by 2^n
                v[i] = m[i] << (n - (i + 1));
            }

            // The remainder of the "directions" are constructed via the recurrance
            //   u_i = u_{i-s} + 2^1 a_1 u_{i-1} + 2^2 a_2 u_{i-2} + \cdots + 2^{s-1} a_{s-1} u_{i-(s-1)} + 2^s u_{i-s}
            // Note + in mod 2 arithmetic is XOR.
            // Use v_i = u_i / 2^(i+1) * 2^n to translate this into a recurrence for the v's
            //   v_i = v_{i-s} 2^{-s} + a_1 v_{i-1} + a_2 v_{i-2} + \cdots + a_{s-1} v_{i-(s-1)} + v_{i-s}
            // where the a's are the binary digits of a = (a_1 a_2 \cdots a_{s-1})
            for (int i = s; i < n; i++) {
                // since the first and last coefficients are known to be one, we start with them
                ulong vi = v[i - s] ^ (v[i - s] >> s);
                // proceed through the bits of a, XORing in the appropriate value if the bit is 1
                for (int j = 1; j < s; j++) {
                    if (((a >> (s - j - 1)) & 1) != 0) {
                        vi ^= v[i - j];
                    }
                }
                v[i] = vi;
            }

        }

        // degree of primitive polynomial (i.e. degree of recurrsion)
        int s;

        // polynomial coefficients, expressed
        ulong a;

        // the number of "directions"
        // the maximum number of values in the sequence is 2^n
        int n = 31;

        // the "directions"
        ulong[] v;

        // the factor to convert a given value to a floating point number on (0,1)
        double f;

        // our positition in the sequence
        // p is the count
        // q stores the "state" of the system
        uint p = 0;
        ulong q = 0;

        public double Current {
            get {
                return (f * q);
            }
        }

        public bool MoveNext () {

            // get index of first zero binary digit of p = (\cdots p_2 p_1 p_0)
            int i = 0;
            while (((p >> i) & 1) != 0) i++;

            q ^= v[i];
            p++;

            return (true);

        }

        public void Reset () {
            p = 0;
            q = 0;
        }

    }


    internal class SobolSequenceParameters {
        public SobolSequenceParameters (int dimension, ulong coefficients, ulong[] seeds) {
            this.Dimension = dimension;
            this.Coefficients = coefficients;
            this.Seeds = seeds;
        }
        public int Dimension { get; private set; }
        public ulong Coefficients { get; private set; }
        public ulong[] Seeds { get; private set; }
    }

}
