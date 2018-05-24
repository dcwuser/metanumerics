using System;
using System.Diagnostics;

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

        // The Sterling sum converges to full double precision within 12 terms for all x > 16.

        private static double Sum(double x) {

            double xx = x * x; // x^2 
            double xk = x; // tracks x^{2k - 1}
            double f = AdvancedIntegerMath.Bernoulli[1] / (2.0 * xk); // k = 1 term
            for (int k = 2; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                double f_old = f;
                xk *= xx;
                f += AdvancedIntegerMath.Bernoulli[k] / ((2 * k) * (2 * k - 1)) / xk;
                if (f == f_old) return (f);
            }

            throw new NonconvergenceException();

        }

        private static double SumPrime(double x) {

            double xx = x * x;
            double xk = xx; // tracks x^{2k}
            double f = -AdvancedIntegerMath.Bernoulli[1] / (2.0 * xk); // k=1 term
            for (int k = 2; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                double f_old = f;
                xk *= xx;
                f -= AdvancedIntegerMath.Bernoulli[k] / (2 * k) / xk;
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();

        }

        // Computes (S(x + y) - S(x)) / y, i.e. the rate at which S(x + y) - S(x) changes with y

        private static double ReducedPochhammerSum(double x, double y) {
            double L = MoreMath.ReducedLogOnePlus(1.0 / x, y);
            double xx = x * x;
            double xk = x;
            double f = AdvancedIntegerMath.Bernoulli[1] / (2.0 * xk) * MoreMath.ReducedExpMinusOne(-L, y);
            for (int k = 2; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                double f_old = f;
                xk *= xx;
                f += AdvancedIntegerMath.Bernoulli[k] / ((2 * k) * (2 * k - 1)) / xk * MoreMath.ReducedExpMinusOne((1 - 2 * k) * L, y);
                if (f == f_old) return (f);
            }
            throw new NonconvergenceException();

        }

        public static double LogGamma(double x) {
            // we-write to use (x-0.5) form to eliminate one one by storing log(2\pi)?
            return (x * Math.Log(x) - x - 0.5 * Math.Log(x / (2.0 * Math.PI)) + Sum(x));
        }

        public static double Gamma(double x) {
            // return (Math.Sqrt(2.0 * Math.PI / x) * Math.Pow(x / Math.E, x) * Math.Exp(Sum(x)));
            return (Math.Exp(LogGamma(x)));
        }

        public static double ReducedLogPochhammer(double x, double y) {
            double L = MoreMath.ReducedLogOnePlus(1.0 / x, y);
            return ((x - 0.5) * L + Math.Log(x + y) - 1.0 + ReducedPochhammerSum(x, y));
        }

        public static double Beta(double x, double y) {
            double xy = x + y;
            return (
                Math.Sqrt(2.0 * Math.PI * xy / x / y) *
                Math.Pow(x / xy, x) * Math.Pow(y / xy, y) *
                Math.Exp(Sum(x) + Sum(y) - Sum(xy))
            );
        }

        public static double LogBeta(double x, double y) {
            double xy = x + y;
            return (x * Math.Log(x / xy) + y * Math.Log(y / xy) - Math.Log(x / xy * y / (2.0 * Math.PI)) / 2.0 + Sum(x) + Sum(y) - Sum(xy));
        }

        public static double Psi(double x) {
            return (Math.Log(x) - 0.5 / x + SumPrime(x));
        }

        // Compute \frac{x^{\nu}}{\Gamma(\nu + 1)}

        public static double PowOverGammaPlusOne(double x, double nu) {
            // This should probably be routed to PoissonProbability, but need to deal with e^{-x} factor.
            return (Math.Pow(x * Math.E / nu, nu) / Math.Sqrt(Global.TwoPI * nu) / Math.Exp(Sum(nu)));
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
            return (Math.Exp(-Sum(k) - lambda * D(k / lambda)) / Math.Sqrt(2.0 * Math.PI * k));
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
