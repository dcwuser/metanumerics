using System;
using System.Collections.Generic;

using Meta.Numerics.Matrices;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Gumbel distribution.
    /// </summary>
    /// <seealso href="http://en.wikipedia.org/wiki/Gumbel_distribution"/>
    public sealed class GumbelDistribution : ContinuousDistribution {

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity));
            }
        }

        /// <summary>
        /// Initializes a new standard Gumbel distribution.
        /// </summary>
        public GumbelDistribution () : this(0.0, 1.0) {
        }

        /// <summary>
        /// Initializes a new Gumbel distribution with the given parameters.
        /// </summary>
        /// <param name="location">The location parameter.</param>
        /// <param name="scale">The scale parameter, which must be positive.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="scale"/> is negative or zero.</exception>
        public GumbelDistribution (double location, double scale) {
            if (scale <= 0.0) throw new ArgumentOutOfRangeException(nameof(scale));
            this.m = location;
            this.s = scale;
        }

        private readonly double m, s;


        /// <summary>
        /// Gets the location parameter of the distribution.
        /// </summary>
        public double Location {
            get {
                return (m);
            }
        }

        /// <summary>
        /// Gets the scale parameter of the distribution.
        /// </summary>
        public double Scale {
            get {
                return (s);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            // check for infinite e; when e = PositiveInfinity, Exp(-e) = 0, but the computer doesn't recognize that the
            // latter is even smaller than the former is big, and gives NaN rather than zero when they are multiplied
            // for some reason e ~ Infinity rather than PositiveInfinity, so we check for that; i actually thought that
            // infinities were always either PositiveInfinity or NegativeInfinity, but that appears not to be the case
            if (Double.IsInfinity(e)) {
                return (0.0);
            } else {
                return (e * Math.Exp(-e) / s);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            return (Math.Exp(-e));
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            return (-MoreMath.ExpMinusOne(-e));
        }

        /// <inheritdoc />
        public override double Hazard (double x) {
            double z = (x - m) / s;
            double e = Math.Exp(-z);
            return (e / MoreMath.ExpMinusOne(e) / s);
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (m - s * Math.Log(-Math.Log(P)));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (m - s * Math.Log(-MoreMath.LogOnePlus(-Q)));
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (m - s * Math.Log(Global.LogTwo));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (m + s * AdvancedMath.EulerGamma);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (MoreMath.Sqr(Math.PI * s) / 6.0);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                return (Math.PI / Math.Sqrt(6.0) * s);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (12.0 * Math.Sqrt(6.0) * AdvancedMath.Apery / MoreMath.Pow(Math.PI, 3));
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                return (12.0 / 5.0);
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r == 1) {
                return (0.0);
            } else {
                return (base.CentralMoment(r));
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                double[] M = RawMoments(r);
                return (M[r]);
            }
        }

        // Substituting u = e^{-z}, the raw moment integral for the standard Gumbel distribution can be re-written as 
        //   M_r = \int_{0}^{\infty} \! e^{-u} ( - \log u )^r \, du
        // Koelbig, "On the Integral \int_{0}^{\infty} \!e^{-\mu t} t^{\nu - 1} \log^m t \, dt", Mathematics of Computation 41 (1983) 171
        // derives formulas for integrals of this form. In particular, his equation (23) implies
        //   M_j = \sum_{k=0}{j-1} \frac{(j-1)!}{k!} \hat{\zeta}(j-k) M_k
        // where \hat{\zeta}(n) = \zeta(n) s^n for n > 1 and m + \gamma s for n = 1. This expresses a given raw moment in terms
        // of all the lower raw moments, making computation of the rth moment O(r^2).

        internal override double[] RawMoments (int rMax) {

            // Create an array to hold the moments.
            double[] M = new double[rMax + 1];
            M[0] = 1.0;
            if (rMax == 0) return (M);

            // Pre-compute the zeta-hat values, since we will use them repeatedly.
            double[] Z = new double[rMax + 1];
            double sj = s; // tracks s^j
            Z[1] = m + s * AdvancedMath.EulerGamma;
            for (int j = 2; j < Z.Length; j++) {
                sj *= s;
                Z[j] = AdvancedMath.RiemannZeta(j) * sj;
            }

            // Compute higher moments in turn.
            for (int j = 1; j <= rMax; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += Z[j - k] * M[k] / AdvancedIntegerMath.Factorial(k);
                }
                M[j] = AdvancedIntegerMath.Factorial(j - 1) * sum;
            }

            return (M);
            
        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (0.0);
            } else if (r == 1) {
                return (Mean);
            } else {
                double C = AdvancedIntegerMath.Factorial(r - 1) * AdvancedMath.RiemannZeta(r) * MoreMath.Pow(s, r);
                return (C);
            }
        }

        /// <summary>
        /// Find the parameters of a Gumbel distribution that best fit the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The fit result.</returns>
        public static GumbelFitResult FitToSample (IReadOnlyList<double> sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            // To do a maximum likelihood fit, start from the log probability of each data point and aggregate to
            // obtain the log likelihood of the sample
            //   z_i = \frac{x_i - m}{s}
            //   -\ln p_i = \ln s + ( z_i + e^{-z_i})
            //   \ln L = \sum_i \ln p_i

            // Take derivatives wrt m and s.
            //   \frac{\partial \ln L}{\partial m} = \frac{1}{s} \sum_i ( 1 - e^{-z_i} )
            //   \frac{\partial \ln L}{\partial s} = \frac{1}{s} \sum_i ( -1 + z_i - z_i e^{-z_i} )

            // Set derivatives to zero to get a system of equations for the maximum.
            //    n = \sum_i e^{-z_i}
            //    n = \sum_i ( z_i - z_i e^{-z_i} )
            // that is, <e^z> = 1 and <z> - <z e^z> = 1.

            // To solve this system, pull e^{m/s} out of the sum in the first equation and solve for m
            //    n = e^{m / s} \sum_i e^{-x_i / s}
            //    m = -s \ln \left( \frac{1}{n} \sum_i e^{-x_i / s} \right) = -s \ln <e^{-x/s}>
            // Substituting this result into the second equation gets us to
            //    s = \bar{x} - \frac{ <x e^{-x/s}> }{ <e^{x/s}> }  
            // which involves only s. We can use a one-dimensional root-finder to determine s, then determine m
            // from the first equation.

            // To avoid exponentiating potentially large x_i, it's better to write the problem in terms
            // of d_i, where x_i = \bar{x} + d_i.
            //    m = \bar{x} - s \ln <e^{-d/s}>
            //    s = -\frac{ <d e^{-d/s}> }{ <e^{-d/s}> }

            // To get the covariance matrix, we need the curvature matrix at the minimum, so take more derivatives
            //    \frac{\partial^2 \ln L}{\partial m^2} = - \frac{1}{s} \sum_i e^{-z_i} = - \frac{n}{s^2}
            //    \frac{\partial^2 \ln L}{\partial m \partial s} = - \frac{n}{s^2} <z e^{-z}>
            //    \frac{\partial^2 \ln L}{\partial s^2} = - \frac{n}{s^2} ( <z^2 e^{-z}> + 1 )

            // Several crucial pieces of this analysis are taken from Mahdi and Cenac, "Estimating Parameters of Gumbel Distribution
            // "using the method of moments, probability weighted moments, and maximum likelihood", Revista de Mathematica:
            // Teoria y Aplicaciones 12 (2005) 151-156 (http://revistas.ucr.ac.cr/index.php/matematica/article/viewFile/259/239) 

            // We will be needed the sample mean and standard deviation
            int n;
            double mean, stdDev;
            Univariate.ComputeMomentsUpToSecond(sample, out n, out mean, out stdDev);
            stdDev = Math.Sqrt(stdDev / n);

            // Use the method of moments to get an initial estimate of s.
            double s0 = Math.Sqrt(6.0) / Math.PI * stdDev;

            // Define the function to zero
            Func<double, double> fnc = (double s) => {
                double u, v;
                MaximumLikelihoodHelper(sample, n, mean, s, out u, out v);
                return (s + v / u);
            };

            // Zero it to compute the best-fit s
            double s1 = FunctionMath.FindZero(fnc, s0);

            // Compute the corresponding best-fit m
            double u1, v1;
            MaximumLikelihoodHelper(sample, n, mean, s1, out u1, out v1);
            double m1 = mean - s1 * Math.Log(u1);

            // Compute the curvature matrix
            double w1 = 0.0;
            double w2 = 0.0;
            foreach (double x in sample) {
                double z = (x - m1) / s1;
                double e = Math.Exp(-z);
                w1 += z * e;
                w2 += z * z * e;
            }
            w1 /= sample.Count;
            w2 /= sample.Count;
            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = (n - 2) / (s1 * s1);
            C[0, 1] = (n - 2) / (s1 * s1) * w1;
            C[1, 1] = (n - 2) / (s1 * s1) * (w2 + 1.0);
            SymmetricMatrix CI = C.CholeskyDecomposition().Inverse();
            // The use of (n-2) here in place of n is a very ad hoc attempt to increase accuracy.


            // Compute goodness-of-fit
            GumbelDistribution dist = new GumbelDistribution(m1, s1);
            TestResult test = sample.KolmogorovSmirnovTest(dist);

            return (new GumbelFitResult(m1, s1, CI[0, 0], CI[1, 1], CI[0, 1], test));
        }

        // Compute u = <e^{-d/s}> and v = <d e^{-d/s}> for a given s and a given sample.

        private static void MaximumLikelihoodHelper (IEnumerable<double> sample, int n, double mean, double s, out double u, out double v) {
            u = 0.0;
            v = 0.0;
            foreach (double x in sample) {
                double d = x - mean;
                double e = Math.Exp(-d / s);
                u += e;
                v += d * e;
            }
            u /= n;
            v /= n;
        }

    }



}
