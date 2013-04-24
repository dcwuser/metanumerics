using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a beta distribution.
    /// </summary>
    /// <remarks>
    /// <para>The beta distribution is defined on the interval [0,1]. Depending on its two shape parameters, it can take on a wide
    /// variety of forms on this interval.</para>
    /// <para>If the two shape parameters are equal, the distribution is symmetric. If the first shape parameter is less than one,
    /// the distribution has a singularity at its left endpoint. If the first shape parameter is greater than one, the distribution
    /// goes to zero at its left endpoint. The second shape parameter similarly governs the distribution's behavior at its right
    /// endpoint.</para>
    /// <para>When both shape parameters are one, the beta distribution reduces to a standard uniform distribution.</para>
    /// <img src="../images/UniformFromBeta.png" />
    /// <para>Beta distributions describe the maximum and minimum values obtained from multiple, independent draws from a standard
    /// uniform distribution. For n draws, the maximum value is distributed as B(n,1).</para>
    /// <img src="../images/BetaFromUniform.png" />
    /// <para>Similiarly, the minimum value is distributed as B(1,n).</para>
    /// <para>Because of the wide variety of shapes it can take, the beta distribution is sometimes
    /// used as an ad hoc model to describe any distribution observed on a finite interval.</para>
    /// </remarks>
    /// <seealso href="http://en.wikipedia.org/wiki/Beta_distribution"/>
    /// <seealso href="http://mathworld.wolfram.com/BetaDistribution.html"/>
    public sealed class BetaDistribution : Distribution {

        /// <summary>
        /// Instantiates a new &#x3B2; distribution.
        /// </summary>
        /// <param name="alpha">The left shape parameter, which controls the form of the distribution near x=0.</param>
        /// <param name="beta">The right shape parameter, which controls the form of the distribution near x=1.</param>
        /// <remarks>
        /// <para>The <paramref name="alpha"/> shape parameter controls the form of the distribution near x=0. The
        /// <paramref name="beta"/> shape parameter controls the form of the distribution near z=1. If a shape parameter
        /// is less than one, the PDF diverges on the side of the distribution it controls. If a shape parameter
        /// is greater than one, the PDF goes to zero on the side of the distribution it controls. If the left and right
        /// shapre parameters are equal, the distribution is symmetric about x=1/2.</para>
        /// </remarks>
        /// <seealso href="http://en.wikipedia.org/wiki/Beta_distribution" />
        public BetaDistribution (double alpha, double beta) {
            if (alpha <= 0.0) throw new ArgumentOutOfRangeException("alpha");
            if (beta <= 0.0) throw new ArgumentOutOfRangeException("beta");
            this.alpha = alpha;
            this.beta = beta;
            // cache value of B(alpha, beta) to avoid having to re-calculate it whenever needed
            this.bigB = AdvancedMath.Beta(alpha, beta);
            // get a beta generator
            this.betaRng = DeviateGeneratorFactory.GetBetaGenerator(alpha, beta);
        }

        private readonly double alpha, beta;
        private readonly double bigB;
        private readonly IDeviateGenerator betaRng;

        /// <summary>
        /// Gets the left shape parameter.
        /// </summary>
        public double Alpha {
            get {
                return (alpha);
            }
        }

        /// <summary>
        /// Gets the right shape parameter.
        /// </summary>
        public double Beta {
            get {
                return (beta);
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, 1.0));
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if ((x < 0.0) || (x > 1.0)) {
                return (0.0);
            } else {
                return (Math.Pow(x, alpha - 1.0) * Math.Pow(1.0 - x, beta - 1.0) / bigB);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x >= 1.0) {
                return (1.0);
            } else {
                // the value is I_x(a,b), which is LeftRegularizedBeta
                // since we have precomputed bigB, use it here instead of calling
                // that method, which would recompute it
                return (AdvancedMath.Beta(alpha, beta, x) / bigB);
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else if (x >= 1.0) {
                return (0.0);
            } else {
                // use 1.0 - I_x(a, b) = I_{1-x}(b, a), which is essentially the symmetry of
                // the beta distribution definition, to avoid calculating 1 - small number
                return (AdvancedMath.LeftRegularizedBeta(beta, alpha, 1.0 - x));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (alpha / (alpha + beta));
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                double ab = alpha + beta;
                return (alpha * beta / (ab + 1.0) / (ab * ab));
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                double ab = alpha + beta;
                return (2.0 * (beta - alpha) / (ab + 2.0) * Math.Sqrt((ab + 1.0) / (alpha * beta)));
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else {
                // this is just a recursive development of Beta(alpha + n, beta) / Beta(alpha, beta)
                double sum = alpha + beta;
                double M = 1.0;
                for (int i = 0; i < n; i++) {
                    M = (alpha + i) / (sum + i) * M;
                }
                return (M);
            }
        }

        // According to http://mathworld.wolfram.com/BetaDistribution.html, central moments are given by
        //   C_n = \left( -\frac{\alpha}{\alpha+\beta} \right)^n 2F1(\alpha, -n; \alpha + \beta; \frac{\alpha+\beta}{\alpha})
        // where 2F1 is a hypergeometric function. The Gauss recurrence for hypergeometric functions (Abromowitz and Stegun 15.2.22)
        //   [c - 2 b + (b - a) z] 2F1(a, b; c; z) + b (1 - z) 2F1(a, b + 1; c; z) - (c - b) 2F1(a, b - 1; c; z) = 0
        // implies
        //   C_{n+1} = \frac{n}{(\alpha + \beta + n) (\alpha + \beta)} \left[ (\beta - \alpha) C_{n} + \frac{\alpha \beta}{\alpha + \beta} C_{n - 1} \right]
        // This recurrence appears to work and not suffer from the cancelation errors that computation from the raw moments does.

        /// <inheritdoc />
        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else {

                // use recurrsion

                double s = alpha + beta;
                double t = beta - alpha;
                double u = alpha * beta / s;

                double C0 = 1.0;
                double C1 = 0.0;
                for (int i = 1; i < n; i++) {
                    double C2 = i / (s + i) / s * (t * C1 + u * C0);
                    C0 = C1;
                    C1 = C2;
                }
                return (C1);

            }
        }

        // inverse CDF

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            // we invert I_x(a, b) = P by first approximating x from P, then refining the result by root polishing
            double x0 = ApproximateInverseBeta(alpha, beta, P);
            if ((x0 == 0.0) || (x0 == 1.0)) return (x0);
            double x1 = RefineInverseBeta(x0, P);
            return (x1);
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException("Q");
            // improve this to deal with small values of Q
            return (InverseLeftProbability(1.0 - Q));
        }

        private double ApproximateInverseBeta (double a, double b, double P) {

            // try series from the left
            double x0 = BetaInverseSeries(a, b, P);
            if (x0 >= 0.0) return (x0);

            // try series from the right
            double x1 = BetaInverseSeries(b, a, 1.0 - P);
            if (x1 >= 0.0) return (1.0 - x1);

            // use the normal approximation
            // the series cases above should take care of extreme values, so this will be used only sufficiently close to the peak
            // that we do not overrun the [0,1] bounds; this appears to work but I have no proof that it always does
            double z = AdvancedMath.ApproximateProbit(P);
            return (Mean + StandardDeviation * z);
            // can we apply a transformation to remove skew?

        }

        // this series is documented at http://functions.wolfram.com/GammaBetaErf/InverseBetaRegularized/06/01/02/

        private double BetaInverseSeries (double a, double b, double P) {

            double w = a * P * bigB;
            if (w < 1.0) {
                double s = Math.Pow(w, 1.0 / a);
                // in this series, the fiducial quantity is z when a ~ b or a >> b, and b z /a^2 when b >> a
                double u = (b <= a) ? s : Math.Max(b, 1.0) * s / Math.Max(a, 1.0);
                //double u = s;
                //if (b > a) u = Math.Max(b, 1.0) * s / Math.Max(a, 1.0);
                if (u < 0.25) {
                    double ap1 = a + 1.0;
                    double s2 = (b - 1.0) / ap1 * s;
                    double s3 = (a * a + (3.0 * b - 1.0) * a + (5.0 * b - 4.0)) / ap1 / (a + 2.0) / 2.0 * s2 * s;
                    double s4 = (a * a * a * a + (6.0 * b - 1.0) * a * a * a + (8.0 * b * b + 11.0 * b - 10.0) * a * a +
                        (33.0 * b * b - 30.0 * b + 4.0) * a + (31.0 * b * b - 47.0 * b + 18.0)) / (ap1 * ap1) / (a + 2.0) / (a + 3.0) / 3.0 * s2 * s * s;
                    return (s * (1.0 + s2 + s3 + s4));
                    // including the s4 term more than doubles our operations count, but it is the only way we can get into the convergence
                    // basin of our root polisher for high enough u that we can avoid using the normal approximation when it would return
                    // an x out of bounds
                }
            }
            // signal failure by returning a negative number
            return (-1.0);

        }

        // refine an approximation to the inverse beta function by root polishing

        private double RefineInverseBeta (double x, double P) {

            double a1 = alpha - 1.0; double b1 = beta - 1.0;
            for (int i = 0; i < 8; i++) {
                double x_old = x;
                // we actually find roots of B(a,b,x) - P B(a,b) rather than that I(a, b, x) - P
                // this reduces the number of calculations of and divisions by B(a,b)
                double y = AdvancedMath.Beta(alpha, beta, x) - P * bigB;
                double yp = Math.Pow(x, a1) * Math.Pow(1.0 - x, b1);
                // dx = -y / yp is Newton's method (1st derivative); improve this by using Halley's method (2nd derivative)
                double dx = -y / (yp - y / 2.0 * (a1 / x - b1 / (1 - x)));
                x += dx;
                //if (x < 0.0) {
                //    x = x_old / 2.0;
                //} else if (x > 1.0) {
                //    x = (x_old + 1.0) / 2.0;
                //}
                if (x == x_old) return (x);
            }
            return (x);

        }

        /*
        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException("rng");
            return (betaRng.GetNext(rng));
        }
        */

        /// <summary>
        /// Computes the Beta distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <remarks>
        /// <para>The returned fit parameters are the &#x3B1; (<see cref="Alpha"/>) and  &#x3B2; (<see cref="Beta"/>) parameters, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="BetaDistribution(double,double)"/> constructor to
        /// specify a new Beta distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        /// <exception cref="InvalidOperationException">Not all the entries in <paramref name="sample" /> lie between zero and one.</exception>
        public static FitResult FitToSample (Sample sample) {
            if (sample == null) throw new ArgumentNullException("sample");
            if (sample.Count < 3) throw new InsufficientDataException();

            // maximum likelyhood calculation
            //   \log L = \sum_i \left[ (\alpha-1) \log x_i + (\beta-1) \log (1-x_i) - \log B(\alpha,\beta) \right]
            // using \frac{\partial B(a,b)}{\partial a} = \psi(a) - \psi(a+b), we have
            //   \frac{\partial \log L}{\partial \alpha} = \sum_i \log x_i -     N \left[ \psi(\alpha) - \psi(\alpha+\beta) \right]
            //   \frac{\partial \log L}{\partial \beta}  = \sum_i \log (1-x_i) - N \left[ \psi(\beta)  - \psi(\alpha+\beta) \right]
            // set equal to zero to get equations for \alpha, \beta
            //   \psi(\alpha) - \psi(\alpha+\beta) = <\log x>
            //   \psi(\beta) - \psi(\alpha+\beta) = <\log (1-x)>

            // compute the mean log of x and (1-x)
            // these are the (logs of) the geometric means
            double ga = 0.0; double gb = 0.0;
            foreach (double value in sample) {
                if ((value <= 0.0) || (value >= 1.0)) throw new InvalidOperationException();
                ga += Math.Log(value); gb += Math.Log(1.0 - value);
            }
            ga /= sample.Count; gb /= sample.Count;

            // define the function to zero
            Func<double[], double[]> f = delegate(double[] x) {
                double pab = AdvancedMath.Psi(x[0] + x[1]);
                return (new double[] {
                    AdvancedMath.Psi(x[0]) - pab - ga,
                    AdvancedMath.Psi(x[1]) - pab - gb
                });
            };

            // guess initial values using the method of moments
            //   M1 = \frac{\alpha}{\alpha+\beta} C2 = \frac{\alpha\beta}{(\alpha+\beta)^2 (\alpha+\beta+1)}
            // implies
            //   \alpha = M1 \left( \frac{M1 (1-M1)}{C2} - 1 \right)
            //   \beta = (1 - M1) \left( \frac{M1 (1-M1)}{C2} -1 \right)
            double m = sample.Mean; double mm = 1.0 - m;
            double q = m * mm / sample.Variance - 1.0;
            double[] x0 = new double[] { m * q, mm * q };

            // find the parameter values that zero the two equations
            double[] x1 = FunctionMath.FindZero(f, x0);
            double a = x1[0]; double b = x1[1];

            // take more derivatives of \log L to get curvature matrix
            //   \frac{\partial^2 \log L}{\partial\alpha^2} = - N \left[ \psi'(\alpha) - \psi'(\alpha+\beta) \right]
            //   \frac{\partial^2 \log L}{\partial\beta^2}  = - N \left[ \psi'(\beta)  - \psi'(\alpha+\beta) \right]
            //   \frac{\partial^2 \log L}{\partial \alpha \partial \beta} = - N \psi'(\alpha+\beta)
            // covariance matrix is inverse of curvature matrix
            SymmetricMatrix CI = new SymmetricMatrix(2);
            CI[0, 0] = sample.Count * (AdvancedMath.Psi(1, a) - AdvancedMath.Psi(1, a + b));
            CI[1, 1] = sample.Count * (AdvancedMath.Psi(1, b) - AdvancedMath.Psi(1, a + b));
            CI[0, 1] = sample.Count * AdvancedMath.Psi(1, a + b);
            CholeskyDecomposition CD = CI.CholeskyDecomposition();
            SymmetricMatrix C = CD.Inverse();

            // do a KS test on the result
            TestResult test = sample.KolmogorovSmirnovTest(new BetaDistribution(a, b));

            // return the results
            FitResult result = new FitResult(x1, C, test);
            return (result);

        }

    }

}