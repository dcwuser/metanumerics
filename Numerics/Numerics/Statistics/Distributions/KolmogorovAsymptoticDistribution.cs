using System;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    internal class KolmogorovAsymptoticDistribution : Distribution {

        // The implementation of this distribution is based on an asymptotic expansion of the CDF
        // for finite-n Kolmogorov-Smirnov distribution
        //    P(x) = K_0(x) + K_1(x) / \sqrt{n} + K_2(x) / n + \cdots
        // as described in
        //   Pelz and Good, "Approximating the Lower Tail-areas of the Kolmogorov-Smirnov One-sample Statistic",
        //   Journal of the  Royal Statistical Society Series B, 38 (1976) 152-156
        // which is based on earlier work of Li-Chien (1956) and Korolyuk (1960).
        // By differentiating the expressions, we obtain p(x). By integrating p(x) x^m, we obtain expressions
        // for the moments. Thus the approximation is self-consistent.

        // The leading term is
        //   K_0 = \sum_{k=-\infty}^{+\infty} (-1)^k e^{-2 k^2 x^2} = 1 - 2 \sum_{k=1}^{\infty} (-1)^{k-1} e^{-2 k^2 x^2}
        // This form is useful for large x. To transform it into a form useful for small x, note that this is just the theta
        // function (http://en.wikipedia.org/wiki/Theta_function) \vartheta(1/2; 2 i x^2 / \pi). A Jacobi transformation relates
        // relates \vartheta(z; \tau) to \vartheta(z/\tau; -1/\tau). Apply this transformation to get
        //   K_0 = \sum_{k=-\infty}^{+\infty} (-1)^k e^{-2 k^2 x^2} =
        //   \frac{1}{x} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \exp\left\{-\frac{(2k-1)^2 \pi^2}{8 x^2} \right\}
        // a form useful to small x.

        // The first correction is given by \frac{1}{n^{1/2}} times
        //   K_1 = - 2 x / 3 \sum_{k=-\infty}{+\infty} (-1)^k k^2 e^{-2 k^2 x^2} = 4 / 3 \sum_{k=1}^{\infty} (-1)^{k-1} k^2 e^{-2 k^2 x^2}
        // To get a form of this useful for small-x, take derivatives of our previous theta relationship to obtain
        //   \sum_{k=-\infty}^{+\infty} (-1)^k k^2 e^{-2 k^2 x^2} =
        //   \frac{1}{4 x^3} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \left[ 1 - \frac{(2k-1)^2 \pi^2}{4 x^2} \right] \exp\left\{-\frac{(2k-1)^2 \pi^2}{8 x^2} \right\}
        // from which follows
        //   P_1 = - \frac{1}{6 x^2} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \left[ 1 - \frac{(2k-1)^2 \pi^2}{4 x^2} \right] \exp\left\{-\frac{(2k-1)^2 \pi^2}{8 x^2} \right\}
        // a form useful for small x.

        // To compute the PDF instead of the CDF we just take derivatives of the relevant expressions. For example
        //   p_0 = P_0' = - 4 x \sum_{k=-\infty}^{+\infty} (-1)^k k^2 e^{-2 k^2 x^2} =
        //   - \frac{1}{x^2} \sqrt{\frac{\pi}{2}} \sum_{k=-\infty}^{+\infty} \left[ 1 - \frac{(2k-1)^2 \pi^2}{4 x^2} \right] \exp\left\{-\frac{(2k-1)^2 \pi^2}{8 x^2} \right\}
        // For very small n and very extreme values of x, we can obtain negative p(x), but in practice this is extremely unlikely.
        // This is a problem well known for Edgeworth expansions, and is really unavoidable for any expansion of p(x). Since the
        // leading term must integrate to one, and the whole series must also integrate to one, higher terms must integrate to zero.
        // Hence they must have regions of negative values.

        public KolmogorovAsymptoticDistribution (int n) {
            if (n < 1) throw new ArgumentOutOfRangeException(nameof(n));
            this.n = n;
        }

        private int n;

        // The dividing line between where P-series and Q-series are used.

        const double x0 = 1.0;


        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        public override double LeftProbability (double x) {
            if (x < 0.0) {
                return (0.0);
            } else if (x < x0) {
                return (P(x));
            } else {
                return (1.0 - Q(x));
            }
        }

        public override double RightProbability (double x) {
            if (x < 0.0) {
                return (1.0);
            } else if (x < x0) {
                return (1.0 - P(x));
            } else {
                return (Q(x));
            }
        }

        private double P (double x) {

            // Pelz and good write arguments of exponents as \pi^2 k^2 / 2 x^2 and \pi^2 (k + 1/2)^2 / 2 x^2, i.e.
            // squared integer and half-integer multiples of \pi^2 / 2 x^2. We find it easier to define k = 2j
            // so that these become squared even and odd integers. That way we can loop exclusively over integers,
            // testing for convergence after each one. 

            double x2 = x * x;

            double u = Math.PI / x / 2.0;
            double u2 = u * u;
            double a = -u2 / 2.0;

            double c0 = Global.SqrtTwoPI / x;
            double c1 = c0 / 6.0 / x / Math.Sqrt(n);
            double c2 = c0 / 72.0 / x2 / n;

            double a0 = 6.0 * x2 + 2.0;
            double a2 = (2.0 * x2 - 5.0) * u2;
            double a4 = (1.0 - 2.0 * x2) * u2 * u2;

            double s = 0.0;

            for (int j = 1; j < Global.SeriesMax; j++) {

                double s_old = s;

                int j2 = j * j;

                double d0, d1, d2;
                if (j % 2 == 0) {
                    d0 = 0.0;
                    d1 = 0.0;
                    d2 = -Math.PI * Math.PI / 2.0 * j2;
                } else {
                    d0 = 1.0;
                    d1 = u2 * j2 - 1.0;
                    d2 = a0 + a2 * j2 + a4 * j2 * j2;
                }
                
                double ds = (c0 * d0 + c1 * d1 + c2 * d2) * Math.Exp(a * j2);

                s += ds;

                if (s == s_old) return (s);

            }

            throw new NonconvergenceException();
        }

        private double Q (double x) {

            if (x > Global.SqrtMax) return (0.0);

            // In the first version of this, we computed K0, K1, and K2 in three seperate methods.
            // This allows us to distinguish contribution of each and to do checks like higher
            // corrections integrating to zero to preserve unitarity, but it means we unnecessarily
            // calculate identical exponentials repeatedly. Now we do it all together.

            double x2 = x * x;
            double a = -2.0 * x2;

            // Pre-factors for each term
            double c0 = 2.0;
            double c1 = -4.0 / 3.0 * x / Math.Sqrt(n);
            double c2 = -1.0 / 9.0 / n;

            // Some coefficients used in computing O(1/n) term
            double a0 = 1.0 - 12.0 * x2;
            double a2 = 4.0 * x2 * (2.0 * x2 - 1.0);
            double b0 = -1.0;
            double b2 = 4.0 * x2;

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;

                int k2 = k * k;

                // O(1) term
                double d0 = 1.0;

                // O(1/\sqrt{n}) term
                double d1 = k2;

                // O(1/n) term
                double d2 = (a0 + a2 * k2) * k2;
                if (k % 2 != 0) {
                    d2 += b0 + b2 * k2;
                }

                // Add 'em up
                double ds = (c0 * d0 + c1 * d1 + c2 * d2) * Math.Exp(a * k2);
                if (k % 2 == 0) ds = -ds;
                s += ds;

                if (s == s_old) return (s);
            }

            throw new NonconvergenceException();

        }

        public override double ProbabilityDensity (double x) {
            if (x < 0.0) {
                return (0.0);
            } else if (x < x0) {
                return (PPrime(x));
            } else {
                return (-QPrime(x));
            }
        }

        private double PPrime (double x) {

            double x2 = x * x;

            double u = Math.PI / x / 2.0;
            double u2 = u * u;
            double a = -u2 / 2.0;

            double c0 = Global.SqrtTwoPI / x2;
            double c1 = c0 / 6.0 / x / Math.Sqrt(n);
            double c2 = c0 / 72.0 / x2 / n;

            double u4 = u2 * u2;
            double a0 = -6.0 * (1.0 + x2);
            double a2 = 27.0 * u2;
            double a4 = 12.0 * (x2 - 1.0) * u4;
            double a6 = (1.0 - 2.0 * x2) * u4 * u2;

            double s = 0.0;

            for (int j = 1; j < Global.SeriesMax; j++) {

                double s_old = s;

                int j2 = j * j;

                double d0, d1, d2;
                if (j % 2 == 0) {
                    d0 = 0.0;
                    d1 = 0.0;
                    d2 = -Math.PI * Math.PI / 2.0 * (u2 * j2 - 3.0) * j2;
                } else {
                    int j4 = j2 * j2;
                    d0 = u2 * j2 - 1.0;
                    d1 = u4 * j4 - 5.0 * u2 * j2 + 2.0;
                    d2 = a0 + a2 * j2 + a4 * j4 + a6 * j4 * j2;
                }

                double ds = (c0 * d0 + c1 * d1 + c2 * d2) * Math.Exp(a * j2);

                s += ds;

                if (s == s_old) return (s);

            }

            throw new NonconvergenceException();

        }

        private double QPrime (double x) {

            if (x > Global.SqrtMax) return (0.0);

            double x2 = x * x;
            double a = -2.0 * x2;

            // Pre-factors for each term
            double c0 = -8.0 * x;
            double c1 = 4.0 / 3.0 / Math.Sqrt(n);
            double c2 = 4.0 / 9.0 * x / n;

            // Calculate some coefficients that appear so as not to re-calculate them in the loop.
            double a10 = -1.0;
            double a12 = 4.0 * x2;
            double a20 = 6.0;
            double a22 = 3.0 - 20.0 * x2;
            double a24 = 4.0 * (2.0 * x2 - 1.0) * x2;
            double b20 = -3.0;
            double b22 = 4.0 * x2;

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;

                int k2 = k * k;

                // O(1) term
                double d0 = 1.0;

                // O(1/\sqrt{n}) term
                double d1 = a10 + a12 * k2;

                // O(1/n) term
                double d2 = a20 + a22 * k2 + a24 * k2 * k2;
                if (k % 2 != 0) {
                    d2 += b20 + b22 * k2;
                }

                // Add 'em up
                double ds = k2 * (c0 * d0 + c1 * d1 + c2 * d2) * Math.Exp(a * k2);
                if (k % 2 == 0) ds = -ds;
                s += ds;

                if (s == s_old) return (s);
            }

            throw new NonconvergenceException();

        }

        // The P-series expression for p(x) can be multiplied by x^r and integrated term-by term to give expressions
        // for <x^r> involving powers of 1/k, i.e. zeta series. 

        public override double Moment (int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return (1.0);
            } else {

                double q = Math.Pow(2.0, r / 2.0);
                double d0 = AdvancedMath.DirichletEta(r);
                double d1 = AdvancedMath.DirichletEta(r - 1);
                double g0 = AdvancedMath.Gamma(r / 2.0 + 1.0);

                double m0 = 2.0 / q * g0 * d0;
                double m1 = -Global.SqrtTwo / 3.0 / q * AdvancedMath.Gamma(r / 2.0 + 1.0 / 2.0) * d1  * r;

                // The O(1/N) expression has some canceling infinities for r = 1 so we form the limit by hand for that cases
                double m2;
                if (r == 1) {
                    m2 = Global.SqrtHalfPI / 36.0 * (3.0 * Global.LogTwo - 1.0);
                } else {
                    m2 = -1.0 / 18.0 / q * g0 * (
                        r * (r - 4) * d0 -
                        2 * (r - 1) * AdvancedMath.DirichletEta(r - 2) +
                        2 * (r - 1) * (1.0 - MoreMath.Pow(2, -r)) / (1.0 - MoreMath.Pow(2, 1-r)) * d0
                    );
                }

                return (m0 + m1 / Math.Sqrt(n) + m2 / n);
            }
        }

        public override double Mean {
            get {
                double m0 = Global.SqrtHalfPI * Global.LogTwo;
                double m1 = -1.0 / 6.0;
                double m2 = Global.SqrtHalfPI / 36.0 * (3.0 * Global.LogTwo - 1.0);
                return (m0 + m1 / Math.Sqrt(n) + m2 / n);
            }
        }

    }
}
