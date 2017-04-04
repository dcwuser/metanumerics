using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents a Pareto or power law distribution.
    /// </summary>
    /// <seealso href="http://en.wikipedia.org/wiki/Pareto_distribution"/>
    public sealed class ParetoDistribution : ContinuousDistribution {

        /// <summary>
        /// Initializes a new Pareto distribution.
        /// </summary>
        /// <param name="mu">The scale parameter, which must be positive.</param>
        /// <param name="alpha">The shape parameter, which must be positive.</param>
        public ParetoDistribution (double mu, double alpha) {
            if (mu <= 0.0) throw new ArgumentOutOfRangeException(nameof(mu));
            if (alpha <= 0.0) throw new ArgumentOutOfRangeException(nameof(alpha));
            this.m = mu;
            this.a = alpha;
        }

        private readonly double m;
        private readonly double a;

        /// <summary>
        /// Gets the scale parameter of the Pareto distribution.
        /// </summary>
        public double ScaleParameter {
            get {
                return (m);
            }
        }

        /// <summary>
        /// Gets the shape parameter of the Pareto distribution.
        /// </summary>
        /// <remarks>
        /// <para>For a given shape parameter &#x3B1;, the probability density falls off as x<sup>-(&#x3B1;+1)</sup>.</para>
        /// </remarks>
        public double ShapeParameter {
            get {
                return (a);
            }
        }

        /// <summary>
        /// Gets the Gini coefficient corresponding to the distribution.
        /// </summary>
        public double GiniCoefficient {
            get {
                return (1.0 / (2.0 * a - 1.0));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                if (a > 1.0) {
                    return (a * m / (a - 1.0));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (Math.Pow(2.0, 1.0 / a) * m);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                if (a > 2.0) {
                    return (m / (a - 1.0) * Math.Sqrt(a / (a - 2.0)));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                if (a > 2.0) {
                    double mm = m / (a - 1.0);
                    return (mm * mm * a / (a - 2.0));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                if (a > 3.0) {
                    return (2.0 * (a + 1.0) / (a - 3.0) * Math.Sqrt((a - 2.0) / a));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double RawMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r >= a) {
                return (Double.PositiveInfinity);
            } else {
                // Straightforward integration yields M_r = \frac{\alpha m^r}{\alpha - r}
                return (a / (a - r) * MoreMath.Pow(m, r));
            }
        }

        /// <inheritdoc />
        public override double CentralMoment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else if (r >= a) {
                return (Double.PositiveInfinity);
            } else {

                // Using C_r = \left( \frac{\alpha m}{1 - \alpha} \right)^2 F(-\alpha, -r; 1 - \alpha; 1 - \alpha^{-1}) and the recurrence relation
                // on the b-value of the hypergeometric function, we get the recurrence
                //   (\alpha - 1 - r) C_{r+1} = \frac{r m}{\alpha - 1} \left[ (\alpha + 1) C_r + \frac{\alpha m}{\alpha - 1} C_{r-1} \right]
                // We implement this recurrence here to compute the rth central moment.

                // Start with C_0 = 1 amd C_1 = 0
                double C0 = 1.0;
                double C1 = 0.0;

                // \alpha - 1 will be used repeatedly
                double a1 = a - 1.0;

                // Recurr upward and return C_r
                for (int s = 1; s < r; s++) {
                    double C2 = s * m / a1 / (a1 - s) * ((a + 1.0) * C1 + a * m * C0 / a1);
                    C0 = C1;
                    C1 = C2;
                }
                return (C1);

            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x < m) {
                return (0.0);
            } else {
                return (a / m * Math.Pow(m / x, a + 1.0));
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= m) {
                return (0.0);
            } else {
                double r = (x - m) / m;
                if (a * r < 0.5) {
                    // close to the left border, use a series expansion to
                    // avoid loss of accuracy due to cancelation
                    return (LeftProbabilitySeries((x - m) / m));
                } else {
                    return (1.0 - Math.Pow(m / x, a));
                }
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= m) {
                return (1.0);
            } else {
                return (Math.Pow(m / x, a));
            }
        }

        private double LeftProbabilitySeries (double r) {
            // expand 1 - (1+r)^(-alpha) for small r
            double df = a * r;
            double f = df;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                df *= -r * (a + k) / (k + 1);
                f += df;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        /// <inheritdoc />
        public override double Hazard (double x) {
            if (x < m) {
                return (0.0);
            } else {
                return (a / x);
            }
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException(nameof(P));
            return (m / Math.Pow(1.0 - P, 1.0 / a));
        }

        /// <inheritdoc />
        public override double InverseRightProbability (double Q) {
            if ((Q < 0.0) || (Q > 1.0)) throw new ArgumentOutOfRangeException(nameof(Q));
            return (m / Math.Pow(Q, 1.0 / a));
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(m, Double.PositiveInfinity));
            }
        } 


    }

}