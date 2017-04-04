using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a log-normal distribution.
    /// </summary>
    /// <remarks>
    /// <para>The logrithm of a log-normal distributed variable is distributed normally.</para>
    /// <img src="../images/LogNormalFromNormal.png" />
    /// <para>The log-normal distribution is commonly used in financial engineering as a model of stock prices.
    /// If the rate of return on an asset is distributed normally, then its price will be distributed log-normally.</para>
    /// </remarks>
    /// <seealso cref="NormalDistribution"/>
    /// <seealso href="http://en.wikipedia.org/wiki/Log-normal_distribution" />
    public sealed class LognormalDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a log-normal distribution.
        /// </summary>
        /// <param name="mu">The mean of the underlying normal distribution.</param>
        /// <param name="sigma">The standard deviation of the underlying normal distribution.</param>
        /// <remarks>
        /// <para>Note that the values of &#x3BC; and  &#x3C3; parameters are
        /// <em>not</em> the mean and standard deviation of the log-normal distribution.
        /// They are the mean and standard deviation of their logrithms z = ln x. This is
        /// the standard method of characterizing a log-normal distribution.</para>
        /// </remarks>
        public LognormalDistribution (double mu, double sigma) {
            if (sigma <= 0.0) throw new ArgumentOutOfRangeException(nameof(sigma));
            this.mu = mu;
            this.sigma = sigma;
            this.normal = new NormalDistribution(mu, sigma);
        }

        /// <summary>
        /// Initializes a standard log-normal distribution.
        /// </summary>
        /// <remarks>
        /// <para>A standard log-normal distribution has &#x3BC; = 0 and &#x3C3; = 1.
        /// It is the log transform of the standard normal distribution.</para>
        /// </remarks>
        public LognormalDistribution () : this(0, 1) {
        }

        private double mu = 0.0;
        private double sigma = 1.0;
        private NormalDistribution normal;

        /// <summary>
        /// Gets the value of the &#x3BC; parameter.
        /// </summary>
        /// <remarks>
        /// <para>Note that the value of this parameter is not the mean of the distribution.
        /// It is the value given to the distribution constructor (<see cref="LognormalDistribution(double,double)"/>),
        /// which is the mean of the underlying normal distribution.</para>
        /// </remarks>
        public double Mu {
            get {
                return (mu);
            }
        }

        /// <summary>
        /// Gets the value of the &#x3C3; parameter.
        /// </summary>
        /// <remarks>
        /// <para>Note that the value of this parameter is not the standard deviation of the distribution.
        /// It is the value given to the distribution constructor (<see cref="LognormalDistribution(double,double)"/>),
        /// which is the sandard deviation of the underlying normal distribution.</para>
        /// </remarks>
        public double Sigma {
            get {
                return (sigma);
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) return (0.0);
            double z = (Math.Log(x) - mu) / sigma;
            return (Math.Exp(-z * z / 2.0) / x / (Global.SqrtTwoPI * sigma));
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (Math.Exp(mu + sigma * sigma / 2.0));
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (Math.Exp(mu));
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (MoreMath.ExpMinusOne(sigma * sigma) * Math.Exp(2.0 * mu + sigma * sigma));
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (Math.Sqrt(MoreMath.ExpMinusOne(sigma * sigma)) * (Math.Exp(sigma * sigma) + 2.0));
            }
        }

        /// <inheritdoc />
        public override double ExcessKurtosis {
            get {
                double s2 = sigma * sigma;
                return (MoreMath.ExpMinusOne(4.0 * s2) + 2.0 * MoreMath.ExpMinusOne(3.0 * s2) + 3.0 * MoreMath.ExpMinusOne(2.0 * s2));
            }
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                return (normal.LeftProbability(Math.Log(x)));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else {
                return (normal.RightProbability(Math.Log(x)));
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (Math.Exp(normal.InverseLeftProbability(P)));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (Math.Exp(normal.InverseRightProbability(Q)));
        }

        /// <inheritdoc />
        public override double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException(nameof(rng));
            return (Math.Exp(normal.GetRandomValue(rng)));
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                return (Math.Exp(r * mu + MoreMath.Sqr(r * sigma) / 2.0));
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
            } else if (r == 2) {
                return (Variance);
            } else {

                double[] K = Cumulants(r);
                return (MomentMath.CumulantToCentral(K, r));

                // Begin with M_r = e^{r \mu + \frac{1}{2} r^2 \sigma^2} and use
                //     C_r = <(x-M_1)^r> = \sum_{k=0}^{r} {r \choose k} M_k (-M_1)^{r-k}
                //         = \sum_{k=0}^{r} (-1)^{r-k} {r \choose k} e^{k \mu + \frac{1}{2} k^2 \sigma^2} e^{(r-k)(\mu + \frac{1}{2} \sigma^2)}
                //         = e^{r (\mu + \frac{1}{2} \sigma^2)} \sum_{k=0}^{r} (-1)^{r-k} {r \choose k} e^{\frac{1}{2} k (k - 1) \sigma^2}
                // So far, this is just the usual raw to central moment conversion, subject to the usual significant cancellation errors,
                // particularly for small \sigma, in which case all exponentials will be nearly equal. Write
                //     e^{\frac{1}{2} k (k - 1) \sigma^2} = (e^{\frac{1}{2} k (k - 1) \sigma^2} - 1) + 1.
                // When summed over all k, the contribution of the final +1 vanishes, so
                //     C_r = e^{r (\mu + \frac{1}{2} \sigma^2)} \sum_{k=0}^{r} (-1)^{r-k} {r \choose k} (e^{\frac{1}{2} k (k - 1) \sigma^2} - 1)
                // and since ( ) is zero for k = 0 and k = 2, we need not start the sum until k = 2.
                /*
                double s = 0.0;
                double x = sigma * sigma / 2.0;
                for (int k = 2; k <= r; k++) {
                    double ds = AdvancedIntegerMath.BinomialCoefficient(r, k) * MoreMath.ExpMinusOne(k * (k - 1) * x);
                    if ((r - k) % 2 != 0) ds = -ds;
                    s += ds;
                }
                return (Math.Exp(r * (mu + x)) * s);
                */
                /*
                // This follows from a straightforward expansion of (x-m)^n and substitution of expressions for M_k.
                // It eliminates some arithmetic but is still subject to loss of significance due to cancelation.
                double s2 = sigma * sigma / 2.0;

                double C = 0.0;
                for (int i = 0; i <= r; i++) {
                    double dC = AdvancedIntegerMath.BinomialCoefficient(r, i) * Math.Exp((MoreMath.Sqr(r - i) + i) * s2);
                    if (i % 2 != 0) dC = -dC;
                    C += dC;
                }
                return (Math.Exp(r * mu) * C);
                */
                // this isn't great, but it does the job
                // expand in terms of moments about the origin
                // there is likely to be some cancelation, but the distribution is wide enough that it may not matter
                /*
                double m = -Mean;
                double C = 0.0;
                for (int k = 0; k <= n; k++) {
                    C += AdvancedIntegerMath.BinomialCoefficient(n, k) * Moment(k) * Math.Pow(m, n - k);
                }
                return (C);
                */
            }

        }

        /// <inheritdoc />
        public override double Cumulant (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else {
                double[] K = Cumulants(r);
                return (K[r]);
            }
        }

        internal double[] Cumulants (int rMax) {

            double[] K = new double[rMax + 1];

            K[0] = 0.0;
            if (rMax == 0) return (K);

            K[1] = this.Mean;
            if (rMax == 1) return (K);

            // C. L. Mallows & J. Riordan, "The Inversion Enumerator for Labeled Trees" derives a
            // recurrence for lognormal cumulants.
            //   K_r = (M_1)^r (e^x - 1)^{r-1} J_{r-1}(x)
            // Here x = e^{\sigma^2} and J_{n}(x) is a polynomial with all positive coefficients.
            //   J_{n+1} = \sum_{k=0}^{n} {n \choose k} (1 + x + \cdots x^k) J_{k} J_{n - k}
            // with J_0 = J_1 = 1.

            // Thus the first few cumulants are
            //   K_0 = 1
            //   K_1 = M_1
            //   K_2 = (M_1)^2 (x - 1)
            //   K_3 = (M_1)^3 (x - 1)^2 (2 + x)
            //   K_4 = (M_1)^4 (x - 1)^3 (6 + 6 x^2 + 3 x^2 + x^3)
            //   K_5 = (M_1)^5 (x - 1)^4 (24 + 36 x + 30 x^2 + 20 x^3 + 4 x^4 + x^5)

            double x = Math.Exp(sigma * sigma);
            //double x = 0.5;

            // Form L_k = 1 + x + \cdots + x^{k}
            double[] L = new double[rMax];
            double xk = 1.0;
            L[0] = 1.0;
            for (int i = 1; i < L.Length; i++) {
                xk *= x;
                L[i] = L[i - 1] + xk;
            }

            double y = MoreMath.ExpMinusOne(sigma * sigma) * K[1];
            double yk = K[1];
            //double y = 0.1;
            //double yk = 1.0;

            double[] J = new double[rMax];
            J[0] = 1.0;
            for (int i = 1; i < rMax; i++) {
                J[i] = 0.0;
                IEnumerator<double> B = AdvancedIntegerMath.BinomialCoefficients(i - 1).GetEnumerator();
                for (int j = 0; j < i; j++) {
                    B.MoveNext();
                    J[i] += B.Current * L[j] * J[j] * J[(i - 1) - j];
                }
                yk *= y;
                K[i + 1] = yk * J[i];
            }

            return (K);

        }

        /// <summary>
        /// Computes the log-normal distribution that best fits the given sample.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The best fit parameters.</returns>
        /// <remarks>
        /// <para>The returned fit parameters are the &#x3BC; (<see cref="Mu"/>) and &#x3C3; (<see cref="Sigma"/>) parameters, in that order.
        /// These are the same parameters, in the same order, that are required by the <see cref="LognormalDistribution(double,double)"/> constructor to
        /// specify a new log-normal distribution.</para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is null.</exception>
        /// <exception cref="InsufficientDataException"><paramref name="sample"/> contains fewer than three values.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="sample"/> contains non-positive values.</exception>
        public static FitResult FitToSample (Sample sample) {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 3) throw new InsufficientDataException();

            // Writing out the log likelyhood from p(x), taking its
            // derivatives wrt mu and sigma, and setting them equal
            // to zero to find the minimizing values, you find that
            // the results of the normal fit are reproduced exactly
            // with x -> log x, i.e.

            // \mu = < \log x >, \sigma^2 = < (\log x - \mu)^2 >
            
            // do a one-pass computation of these quantities
            SampleSummary summary = new SampleSummary();
            foreach (double value in sample) {
                if (value <= 0.0) throw new InvalidOperationException();
                summary.Add(Math.Log(value));
            }

            // the second derivatives are also just as in the normal
            // case, including the vanishing of the mixed derivative
            // this makes direct inversion trivial

            SymmetricMatrix C = new SymmetricMatrix(2);
            C[0, 0] = summary.Variance / summary.Count;
            C[1, 1] = summary.Variance / summary.Count / 2;

            // test the fit
            ContinuousDistribution d = new LognormalDistribution(summary.Mean, Math.Sqrt(summary.Variance));
            TestResult r = sample.KolmogorovSmirnovTest(d);

            return (new FitResult(new double[] { summary.Mean, Math.Sqrt(summary.Variance) }, C, r));
        }

    }

}