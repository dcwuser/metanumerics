using System;
using System.Collections.Generic;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics {


    public class ParetoDistribution : Distribution {

        public ParetoDistribution (double mu, double alpha) {
            if (mu <= 0.0) throw new ArgumentOutOfRangeException("mu");
            if (alpha <= 0.0) throw new ArgumentOutOfRangeException("alpha");
            this.mu = mu;
            this.alpha = alpha;
        }

        private double mu;
        private double alpha;

        public double ScaleParameter {
            get {
                return (mu);
            }
        }

        public double ShapeParameter {
            get {
                return (alpha);
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                if (alpha > 1.0) {
                    return (alpha * mu / (alpha - 1.0));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                return (Math.Pow(2.0, 1.0 / alpha) * mu);
            }
        }

        /// <inheritdoc />
        public override double StandardDeviation {
            get {
                if (alpha > 2.0) {
                    return (mu / (alpha - 1.0) * Math.Sqrt(alpha / (alpha - 2.0)));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                if (alpha > 2.0) {
                    double mm = mu / (alpha - 1.0);
                    return (mm * mm * alpha / (alpha - 2.0));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                if (alpha > 3.0) {
                    return (2.0 * (alpha + 1.0) / (alpha - 3.0) * Math.Sqrt((alpha - 2.0) / alpha));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                if (alpha > n) {
                    return (alpha / (alpha - n) * MoreMath.Pow(mu, n));
                } else {
                    return (Double.PositiveInfinity);
                }
            }
        }

        public override double MomentAboutMean (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (0.0);
            } else {
                return (CentralMomentFromRawMoment(n));
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x < mu) {
                return (0.0);
            } else {
                return (alpha / mu * Math.Pow(mu / x, alpha + 1.0));
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x <= mu) {
                return (0.0);
            } else {
                double r = (x - mu) / mu;
                if (alpha * r < 0.5) {
                    // close to the left border, use a series expansion to
                    // avoid loss of accuracy due to cancelation
                    return (LeftProbabilitySeries((x - mu) / mu));
                } else {
                    return (1.0 - Math.Pow(mu / x, alpha));
                }
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x <= mu) {
                return (1.0);
            } else {
                return (Math.Pow(mu / x, alpha));
            }
        }

        private double LeftProbabilitySeries (double r) {
            // expand 1 - (1+r)^(-alpha) for small r
            double df = alpha * r;
            double f = df;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                df *= -r * (alpha + k) / (k + 1);
                f += df;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();
        }

        /// <inheritdoc />
        public override double InverseLeftProbability (double P) {
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            return (mu / Math.Pow(1.0 - P, 1.0 / alpha));
        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(mu, Double.PositiveInfinity));
            }
        } 


    }

}