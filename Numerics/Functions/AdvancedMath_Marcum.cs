using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Text;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        // Marcum Q function
        //   Q_m(a, b) = a^{1-m} \int_{b}^{\infty} \! dx \, x^m e^{-(a^2 + x^2)/2} I_{m-1}(a x)
        // where I_{\nu}(z) is modified Bessel function.

        // To get series development, expand I_{m-1}(a x) using series definition and integrate term-by-term.
        //   Q_m(a, b) = e^{-a^2/2} \sum_{k=0}^{\infty} \frac{1}{k!} \left( \frac{a^2}{2} \right)^k Q(k+m, b^2/2)
        // where Q(x, y) = \Gamma(x, y) / \Gamma(y) is right regularized incomplete Gamma function.

        // Writing in terms of left regularized incomplete Gamma functions
        // P(x, y) = 1 - Q(x, y) = \gamma(x, y) / \Gamma(y) instead gives a complementary function
        //   P_m(a, b) = e^{-a^2/2} \sum_{k=0}^{\infty} \frac{1}{k!} \left( \frac{a^2}{2} \right)^k P(k+m, b^2/2)
        // with
        //   P_m(a, b) + Q_m(a, b) = 1

        // From limits of regularized incomplete Gamma functions it follows that
        //                  Q_{\infty}(a, b) = 1  Q increases with m
        //   Q_m(a, 0) = 1  Q_m(a, \infty) = 0    Q decreases with a
        //                  Q_m(\infty, b) = 1    Q increases with b

        // The transition from P small to Q small occurs around
        //   b^2/2 \approx a^2/2 + m

        // A lot of good results in Gil, Segura, and Temme, "Computation of the Marcum Q-function"
        // https://arxiv.org/abs/1311.0681
        // Note they use transformed arguments in a way that as far as I can tell no one else does.

        // See https://arxiv.org/abs/1404.0302 for another asymptotic expansion.

        // See also http://www.phaselockedsystems.com/NoncentralChiSquared.pdf for more.

        internal static void Marcum (double m, double a, double b, out double P, out double Q) {

            if (m < 0.0) throw new ArgumentOutOfRangeException(nameof(m));
            if (a < 0.0) throw new ArgumentOutOfRangeException(nameof(a));
            if (b < 0.0) throw new ArgumentOutOfRangeException(nameof(b));

            double x = 0.5 * a * a;
            double y = 0.5 * b * b;

            // These series work well for x < ~30-60, i.e. a < ~10.
            // I tried to use asymptotic expansion of I to derive a series for large a but surprisingly little luck.
            // Still need to figure out what to do for large a. Temme's asymptotic expansions are painful.
            if (y < x + m) {
                P = MarcumP_Series(m, a, b);
                Q = 1.0 - P;
            } else {
                Q = MarcumQ_Series(m, a, b);
                P = 1.0 - Q;
            }

        }

        // Use the recurrsion Q(\mu + 1, y) = Q(\mu, y) + \frac{y^{\mu} e^{-y}}{\Gamma(\mu + 1)}

        private static double MarcumQ_Series (double m, double a, double b) {
            double x = 0.5 * a * a;
            double y = 0.5 * b * b;
            double q = AdvancedMath.RightRegularizedGamma(m, y);
            double dq = AdvancedMath.PowerOverFactorial(y, m) * Math.Exp(-y); // Use Poisson probability instead
            double t = 1.0;
            double s = q;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                t *= x / k;
                q += dq;
                Debug.Assert(0.0 <= q && q <= 1.0);
                s += t * q;
                if (s == s_old) return Math.Exp(-x) * s;
                dq *= y / (m + k);
            }
            throw new NonconvergenceException();
        }

        // Recurrsion P(\mu + 1, y) = P(\mu, y) - \frac{y^{\mu} e^{-y}}{\Gamma(\mu + 1)} is upward unstable.
        // Temme guesses required k and recurrs down, but that guess is quite complicated.

        private static double MarcumP_Series(double m, double a, double b) {
            double x = 0.5 * a * a;
            double y = 0.5 * b * b;
            double p = AdvancedMath.LeftRegularizedGamma(m, y);
            double t = 1.0;
            double s = p;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double s_old = s;
                t *= x / k;
                p = AdvancedMath.LeftRegularizedGamma(m + k, y);
                s += t * p;
                if (s == s_old) {
                    return Math.Exp(-x) * s;
                }
            }
            throw new NonconvergenceException();
        }

#if FUTURE

        public static void TemmeMarcumLargeMu (double mu, double x, double y, out double P , out double Q) {

            x = 0.5 * x * x;
            y = 0.5 * y * y;

            mu -= 1.0;
            x /= mu;
            y /= mu;

            double z = TemmeZeta(x, y);
            double ep = 0.5 * AdvancedMath.Erfc(z * Math.Sqrt(0.5 * mu));
            double em = 0.5 * AdvancedMath.Erfc(-z * Math.Sqrt(0.5 * mu));

            Q = 0.5 * AdvancedMath.Erfc(-z * Math.Sqrt(0.5 * mu));
            P = 0.5 * AdvancedMath.Erfc(z * Math.Sqrt(0.5 * mu));

            double uSquared = 1.0 / (2.0 * x + 1.0);
            double u = Math.Sqrt(uSquared);
            double f01 = 1.0 / 24.0 * uSquared * (3.0 - 5.0 * uSquared * uSquared);
            double f10 = 1.0 / 6.0 * u * (3.0 + uSquared);

        }

        public static double TemmeZeta (double x, double y) {
            //double s = y - x - 1.0;
            double t = Math.Sqrt(1.0 + 4.0 * x * y);
            double z = x + y - t + Math.Log((1.0 + t) / (2.0 * y));
            z = Math.Sqrt(2.0 * z);
            //if (s < 0.0) z = -z;
            return z;
        }

#endif

    }
}
