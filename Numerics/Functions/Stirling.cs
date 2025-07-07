using System;
using System.Diagnostics;
using System.Runtime.Serialization;

namespace Meta.Numerics.Functions
{
    // The Stirling approximation is an asymptotic series for large |z| (A&S 6.1.40, DLMF 5.11.1):
    //   \ln \Gamma(z) = z \ln z - z - \frac{1}{2} \ln(\frac{z}{2\pi}) + S(z)
    //     S(z) = \sum_{k=1}^{\infty} \frac{B_{2k}}{2k(2k-1) z^{2k-1}}
    // Exponentiate to get \Gamma(z). Differentiate to get \psi(z):
    //   \psi(z) = \ln z - \frac{1}{2z} + S'(z)
    //     S'(z) = -\sum_{k=1}{\infty} \frac{B_{2k}}{2k z^{2k}}

    // One way to get \ln k! = \ln \Gamma(k+1) is just to plug in z = k + 1.
    // But one often finds the following form instead:
    //   \ln k! = \ln \Gamma(k + 1) = \ln (k \Gamma(k)) = \ln k + \ln \Gamma(k)
    //          = k \ln k - k + \frac{1}{2} \ln(2 \pi k ) + S(k)
    // I spent a few minutes trying to get this by manipulating the series with z = k + 1
    // before realizing I could prove it using the \Gamma recursion relation.


    internal static class Stirling {

        // The Stirling sum converges to full double precision within 12 terms for all x > 16.

        private static double Sum(double x) {
            double rxPower = 1.0 / x;
            double rxSquared = rxPower * rxPower;
            double f = 0.5 * AdvancedIntegerMath.Bernoulli[1] * rxPower;
            for (int k = 2; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                double f_old = f;
                rxPower *= rxSquared;
                f += AdvancedIntegerMath.Bernoulli[k] / (2 * k * (2 * k - 1)) * rxPower;
                if (f == f_old) {
                    return f;
                }
            }
            throw new NonconvergenceException();
        }

        private static Complex Sum (Complex z) {
            Complex rzPower = 1.0 / z;
            Complex rzSquared = ComplexMath.Sqr(rzPower);
            Complex f = 0.5 * AdvancedIntegerMath.Bernoulli[1] * rzPower;
            for (int k = 2; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                Complex f_old = f;
                rzPower *= rzSquared;
                f += AdvancedIntegerMath.Bernoulli[k] / (2 * k * (2 * k - 1)) * rzPower;
                if (f == f_old) return f;
            }
            throw new NonconvergenceException();
        }

        private static double SumPrime(double x) {
            double rxSquared = MoreMath.Sqr(1.0 / x);
            double rxPower = rxSquared;
            double f = -0.5 * AdvancedIntegerMath.Bernoulli[1] * rxPower;
            for (int k = 2; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                double f_old = f;
                rxPower *= rxSquared;
                f -= AdvancedIntegerMath.Bernoulli[k] / (2 * k) * rxPower;
                if (f == f_old) return f;
            }
            throw new NonconvergenceException();
        }

        public static double LogGamma(double x) {
            // Sum from smallest to largest terms to minimize error.
            return Sum(x) + halfLogTwoPi - x + (x - 0.5) * Math.Log(x);
        }

        public static Complex LogGamma (Complex z) {
            return (z - 0.5) * ComplexMath.Log(z) - z + halfLogTwoPi + Sum(z);
        }

        private static readonly double halfLogTwoPi = 0.5 * Math.Log(2.0 * Math.PI);

        public static double Gamma(double x) {
            return (Math.Sqrt(2.0 * Math.PI / x) * Math.Pow(x / Math.E, x) * Math.Exp(Sum(x)));
        }

        // Let's try to apply Stirling series to computation of Pochammer (x)_y = \frac{\Gamma(x + y)}{\Gamma(x)}
        // Note Stirling applies as long as x >> 1 and (x + y) >> 1; we don't require y >> 1 and in fact it can be tiny.
        // Plugging in naively gets you to
        //   \ln (x)_y = \ln \Gamma(x + y) - \ln \Gamma(x)
        //             = (x - 1/2) \ln (1 + y/x) + y \ln (x + y) - y + S(x + y) - S(x)
        // Already at this point a few terms have cancelled, which is good for avoiding overflow and maintaining precision.
        // But we can do better. Note that every term here is, for small y, proportional to y.
        // Some are explicitly so. \ln (1 + y/x) is because \ln (1 + e) ~ e for small e. And we expect [S(x + y) - S(x)] ~ y
        // because linear approximation should be good in the small limit. Let us define
        //   \ln (1 + y/x) = log1p(y/x) = \ell
        // If we simply evaluate the sums, the behavior of [S(x + y) - S(x)] for small y will depend on cancellations.
        // To get around this, write
        //   S(x + y) - S(x) = \sum_{k} \frac{B_{2k}}{2k(2k-1)} \left[ \frac{1}{(x + y)^{2k-1}} - \frac{1}{x^{2k-1}} \right]
        //                   = \sum_{k} \frac{B_{2k}}{2k(2k-1)} \frac{1}{x^{2k-1}} \left[ \left( 1 + y/x \right)^{-(2k-1)} - 1 \right]
        // Now use a combination of log1p and expm1 to compute term in brackets (this techinique is used for compound
        // interest calculations, too).
        //   \left( 1 + y/x \right) = e^{ \ln (1 + y/x) } = e^{\ell}
        //   \left( 1 + y/x \right)^{-(2k-1)} = e^{-(2k-1) \ell}
        //   \left( 1 + y/x \right)^{-(2k-1)} - 1 = e^{-(2k-1) \ell} - 1 = expm1(-(2k-1) \ell)
        // So
        //   S(x + y) - S(x) = \sum_{k} \frac{B_{2k}}{2k(2k-1)} \frac{expm1(-(2k-1) \ell)}{x^{2k-1}}
        // Now each term ~y because \ell ~ y and expm1(e) ~ e.

        // As written, the proportionality of \ln(x)_y ~ y is still implicit. It's possible to make it explicit by using
        // "reduced" versions of log1p and expm1 which factor out the proportionality to the argument, e.g.
        //   \ln(1 + e) = e rlog1p(e)
        // I did the work and carried it through, and doing so even makes the result prettier and saves a few flops.
        // But it requires coding these reduced functions and doesn't improve accuracy in any way I can figure.

        public static double LogPochhammer (double x, double y) {
            double z = x + y;
            double ell = MoreMath.LogOnePlus(y / x);
            double sp = -ReducedLogPochhammerSum2(x, z) * y;
            return sp + (x - 0.5) * ell + y * Math.Log(z / Math.E); 
        }

        public static double Pochhammer (double x, double y) {
            // Invoking Pow instead of Exp(LogPochhammer) appears to have
            // higher accuracy for large y and not have lower accuracy for small y.
            double z = x + y;
            double sp = -ReducedLogPochhammerSum2(x, z) * y;
            return Math.Pow(z / x, x - 0.5) * Math.Pow(z / Math.E, y) * Math.Exp(sp);
            //return Math.Sqrt(x / z) * Math.Pow(z / x, x) * Math.Pow(z / Math.E, y) * Math.Exp(s);
        }

        private static double ReducedLogPochhammerSum2 (double x, double z) {
            // First term is 1/2 B_2 R_1
            double R = 1.0;
            double S = 0.5 * AdvancedIntegerMath.Bernoulli[1];
            // We will need some recripricol powers of x and z
            double rx = 1.0 / x;
            double rz = 1.0 / z;
            double rxz = rx + rz;
            double rx2 = MoreMath.Sqr(rx);
            double rz2 = MoreMath.Sqr(rz);
            double rzPower = rz;
            for (int k = 2; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                // Higher terms are B_{2k} R_{2k-1} / (2k * (2k-1))
                double S_old = S;
                // We need to run R recurrsion twice
                // R_{k+1} = R_k / x + 1 / z^k requires 3 flops + 1 to advance power of z so 4 flops, run twice would be 8 flops
                // but R_{k+2} = R_k / x^2 + 1 / z^k * (1/x + 1/z) is still only 3 flops + 1 to advance power of z, so faster this way
                R = R * rx2 + rzPower * rxz;
                int two_k = 2 * k;
                S += AdvancedIntegerMath.Bernoulli[k] / (two_k * (two_k - 1)) * R;
                if (S == S_old) return S * rx * rz;
                rzPower *= rz2;
            }
            throw new NonconvergenceException();
        }

        // For \ln\Beta, perform analytic cancelations of \ln\Gamma(x) - \ln\Gamma(y) - \ln\Gamma(x+y)
        // so that they don't cancel numerically.
        //   B(x, y) = \frac{\Gamma(x) \Gamma(y)}{\Gamma(x + y)}
        //   \ln B(x, y) = \ln\Gamma(x) + \ln\Gamma(y) - \ln\Gamma(x + y)
        //               =   (x - 1/2) \ln x - x + 1/2 \ln(2\pi) + S(x) 
        //                 + (y - 1/2) \ln y - y + 1/2 \ln(2\pi) + S(y)
        //                 - (x + y - 1/2) \ln(x + y) - (x + y) - 1/2 \ln(2\pi) - S(x + y)
        //               = x \ln(\frac{x}{x + y}) + y \ln(\frac{y}{x + y}) + 1/2 \ln(\frac{2\pi(x+y)}{xy}) + S(x) + S(y) - S(x + y)
        // There appears to be no case where there is significant cancellation among sums. For x ~ y, leading term is 1/x + 1/x - 1/(2x) ~ 3/2 1/x.
        // For y << x, leading term is 1/x + 1/y - 1/(x + y) ~ 1/y (1/x nearly cancels with 1/(x+y), but are much smaller than 1/y).

        public static double Beta(double x, double y) {
            Debug.Assert(x > 0.0);
            Debug.Assert(y > 0.0);
            if (x < y) Global.Swap(ref x, ref y);
            Debug.Assert(y <= x);
            double z = x + y;
            double r = y / x;
            Debug.Assert(r <= 1.0);
            if (r > 0.25) {
                // For not-so-different x and y, this form is more accurate.
                return Global.SqrtTwoPI * Math.Sqrt(z / x / y) * Math.Pow(x / z, x) * Math.Pow(y / z, y) * Math.Exp(Sum(x) + Sum(y) - Sum(z));
            } else {
                // For very different x and y, this form is more accurate.
                return Global.SqrtTwoPI / Math.Sqrt(y) * Math.Pow(r, y) * Math.Exp((0.5 - z) * MoreMath.LogOnePlus(r) + Sum(x) + Sum(y) - Sum(z));
            }
        }

        // Applying the Stirling series to the Beta function allow us to cancel the linear terms analytically,
        // and to combine the log terms in a way that involves ratios.

        public static double LogBeta(double x, double y) {
            // Note that all three leading terms in this expression have the same sign (negative), so cancellation
            // isn't an issue.
            // Note that Sum(x) + Sum(y) - Sum(x+y) doesn't suffer any significant calcelation. The inidividual sums
            // are also supressed relative to the leading terms.
            double xy = x + y;
            double s = Sum(x) + Sum(y) - Sum(xy);
            return (x - 0.5) * Math.Log(x / xy) + (y - 0.5) * Math.Log(y / xy) + 0.5 * Math.Log(2.0 * Math.PI / xy) + s;
        }

        public static double Psi(double x) {
            return Math.Log(x) - 0.5 / x + SumPrime(x);
        }

        // The factor x^n / n! appears in multiple places, e.g. factors in front of Bessel function series.
        // Naive evaluation can go wrong because x^n or n! can overflow even when ratio is representable.
        // For example, x = 15, n = 200, which is below Bessel series limit: 16^256 ~ 1.0E308 and 256! ~ 8.5E506,
        // but ratio 2.1E-199 is representable.

        // x^n / n! is just e^x * PoissonProbability(x, n), but we can't re-use the logic there naively.
        // In that case, maximum was 1 and if we underflowed that was fine. Here, maximum may overflow
        // but we still want representable values on the sides, which may occur for underflows of PoissonProbability.
        // So we take a different approach.
        //   P = \frac{x^n}{n!}
        //   \ln P = n \ln x - \ln n!
        //   \ln P = n \ln x - n \ln n + n - \frac{1}{2} \ln(2 \pi n) - S(n)
        //   \ln P = n \ln (e x / n) - \frac{1}{2} \ln(2 \pi n) - S(n)
        //   P = \left( \frac{e x}{n} \right)^n \frac{e^{-S(n)}}{\sqrt{2 \pi n}
        // It's easy to show that maximum is at n ~ x. Only factor that can overflow
        // or underflow is power, and will only do so when result overflows or underflows.

        public static double PowerOverFactorial(double x, double nu) {
            return Math.Exp(-Sum(nu)) / Math.Sqrt(2.0 * Math.PI * nu) * Math.Pow(Math.E * x / nu, nu);
        }

        // Compute \frac{x^a (1-x)^b}{B(a,b)}
        // For large a, b, B(a,b) is often unrepresentably small but this ratio is still representable because x^a (1-x)^b is also small.
        // To compute this ratio, we bring x and (1-x) inside the a and b powers that we need to take anyway to compute B(a,b).

        public static double PowOverBeta(double a, double b, double x) {
            return(a / (a + b) * b * BinomialProbability(x, a, 1.0 - x, b));
        }

        // Our approach to binomial and Poisson probabilities is taken from
        //   Catherine Loader, Fast and Accurate Computation of Binomial Probabilities, 2000
        //   http://octave.1599824.n4.nabble.com/attachment/3829107/0/loader2000Fast.pdf
        // This great paper introduces the D-function we use here; I have not seen it anywhere else.

        // Poisson probability is
        //   P = \frac{\lambda^k e^{-\lambda}}{k!}
        //   \ln P = k \ln \lambda - \lambda - \ln(k!)
        //   \ln P = k \ln \lambda - \lambda - k \ln k + k - \frac{1}{2} \ln(2\pi k) - S(k)
        // This form is not bad, but there is still cancellation between terms (which you would expect
        // in order to produce a maximum at k ~ \lambda). Isolate this cancelation by writing as
        //   \ln P = -\frac{1}{2} \ln(2\pi k) - S(k)
        //     -\lambda \left[ \frac{k}{\lambda} \ln\left(\frac{k}{\lambda}\right) + 1 - \frac{k}{\lambda} \right]
        //   \ln P = -\frac{1}{2} \ln(2\pi k) - S(k) - D(k/\lambda)
        // where
        //   D(x) = x \ln x - x + 1
        // Since D(x) >= 0, all terms are negative.

        public static double PoissonProbability(double lambda, double k) {
            Debug.Assert(lambda > 0.0);
            Debug.Assert(k >= 16.0);
            double kappa = (k - lambda) / lambda;
            return Math.Exp(-Sum(k) - lambda * Dee(kappa)) / Math.Sqrt(2.0 * Math.PI * k); 
            //return (Math.Exp(-Sum(k) - lambda * D(k / lambda)) / Math.Sqrt(2.0 * Math.PI * k));
        }

        // D(x) = (1+x) ln(1+x) - x
        // D >= 0, D(0) = 0, D ~ 1/2 x^2 + O(x^3)

        private static double Dee (double x) {
            Debug.Assert(x >= -1.0);
            if (Math.Abs(x) < 0.5) {
                double xk = x * x;
                double f = 0.5 * xk;
                for (int k = 3; k < Global.SeriesMax; k++) {
                    double f_old = f;
                    xk *= -x;
                    double df = xk / (k * (k - 1));
                    f += df;
                    if (f == f_old) return f;
                }
                throw new NonconvergenceException();
            } else {
                double y = 1.0 + x;
                return y * Math.Log(y) - x;
            }
        }
        

        // Binomial probability is
        //   P = \binom{n}{k} p^k q^{n-k} = \frac{n! p^k q^{n-k}}{k! (n-k)!}
        // where q = 1 - p. Taking logs and using the Stirling series
        //   \ln P = k \ln p + (n-k) \ln q + \ln n! - \ln k! - \ln(n-k)!
        //         = \frac{1}{2} \ln\left(\frac{n}{2 \pi n (n-k)}\right) + S(n) - S(k) - S(n-k) +
        //           + k \ln\left(\frac{k}{np}\right) + (n-k) \ln\left(\frac{n-k}{nq}\right)
        // There are cancellations among the S, but they tend not to be too big since S(x) ~ \frac{1}{12x}.
        // But the last two terms can exhibit large cancellations. Loader's insight was that they could
        // be written in terms of D(x):
        //   np D(k/np) + nq D((n-k)/nq) =
        //     np \left[ \frac{k}{np} \ln\left(\frac{k}{np}\right) - \frac{k}{np} + 1 \right] +
        //     nq \left[ \frac{n-k}{nq} \ln\left(\frac{n-k}{nq}\right) - \frac{n-k}{nq} + 1 \right]
        //     = k \ln\left(\frac{k}{np}\right) + (n-k) \ln\left(\frac{n-k}{nq}\right)
        //       - k - (n - k) + n(p + q)
        // Since p + q = 1, the last three terms cancel, leaving the last two terms above. So
        //   P = \sqrt{\frac{n}{2 \pi n (n-k)}} e^{[S(n) - S(k) - S(n-k)] + np D(k/np) + nq D((n-k)/nq)}
        // and since D(x) >= 0, there is no more cancelation among these terms.

        public static double BinomialProbability(double p, double k, double q, double nmk) {
            Debug.Assert(Math.Abs(p + q - 1.0) < 1.0E-15);
            Debug.Assert(k > 0);
            Debug.Assert(nmk > 0);
            // Zero probabilities result in 0 * Infinity = NaN in our formulas, so handle them here.
            if ((p == 0.0) || (q == 0.0)) return (0.0);
            Debug.Assert(p > 0.0);
            Debug.Assert(q > 0.0);
            double n = k + nmk;
            double np = n * p;
            double nq = n * q;
            double S1 = Sum(n) - Sum(k) - Sum(nmk);
            double D1 = np * D(k / np) + nq * D(nmk / nq);
            return (Math.Sqrt(n / (2.0 * Math.PI * k * (n - k))) * Math.Exp(S1 - D1));
        }

        // We formed D(x) = x \ln x - x + 1 to isolate cancellations, but there are still cancelations hiding inside it.
        // By graphing and linearizing, you can discover that the region where we need to worry about
        // cancelations is x ~ 1. Carefully doing a full expansion of (1 + y) \ln (1 + y) gives
        //   D(1 + y) = \sum_{k=2}^{\infty} \frac{(-1)^k y^k}{k (k-1)} = 1/2 y^2 + \cdots

        private static double D (double x) {
            Debug.Assert(x >= 0.0);
            if (x == 0.0) {
                return (1.0);
            } else {
                double y = x - 1.0;
                // Ideally, we would get y = x - 1 as a direct input,
                // but that's not how it arises in our problems, so
                // we accept the loss from this one subtraction.
                if (Math.Abs(y) < 0.25) {
                    double yk = y * y;
                    double f = 0.5 * yk;
                    for (int k = 3; k < Global.SeriesMax; k++) {
                        double f_old = f;
                        yk *= -y;
                        double df = yk / (k * (k - 1));
                        f += df;
                        if (f == f_old) return (f);
                    }
                    throw new NonconvergenceException();
                }
                else {
                    return (x * Math.Log(x) - x + 1.0);
                }
            }
         }
    }

}
