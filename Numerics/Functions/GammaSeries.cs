using System;
using System.Diagnostics;

namespace Meta.Numerics.Functions {

    // This class evaluates Gamma and related function on the basis of the series
    // DLMF 5.7.3 (https://dlmf.nist.gov/5.7):
    //   \ln\Gamma(1 + x) = -\ln(1 + x) + x (1 - \gamma) +
    //     \sum_{k=2}^{\infty} \frac{\zeta(k) - 1}{k} (-z)^k

    // This is more accurate than Lanczos in the region -1/2 < x < 5/2, where
    // \ln\Gamma has zeros that Lanczos relies on delicate cancelations to reproduce. 
    // Since the series can take up to ~30 terms near half-integers, it's slower
    // than Lanczos for those values. Very close to the integers, though, it's faster
    // than Lanczos.

    internal static class GammaSeries {

        // This is the largest radius we are willing to compute, determined by number of
        // zeta values we record and how accuracy drops off as we get further away.

        public static readonly double Limit = 0.75;

        // We record (-1)^k [\zeta(k) - 1] explicitly for two reasons. First, it's faster than
        // invoking AdvancedMath.Zeta. Second, \zeta(k) gets very close to 1, and we
        // want to know the difference to high accuracy. An alternative to storing
        // externally computed values is to sum the zeta series without the first unit
        // term; this converges quickly for large k.

        private static double[] altZetaMinusOne = new double[]
        {
            -1.5,
            Double.NaN,
            6.44934066848226436472415E-1,
            -2.02056903159594285399738E-1,
            8.23232337111381915160037E-2,
            -3.69277551433699263313655E-2,
            1.73430619844491397145179E-2,
            -8.34927738192282683979755E-3,
            4.07735619794433937868524E-3,
            -2.00839282608221441785277E-3,
            9.94575127818085337145959E-4,
            -4.94188604119464558702282E-4,
            2.46086553308048298637998E-4,
            -1.22713347578489146751836E-4,
            6.12481350587048292585451E-5,
            -3.05882363070204935517285E-5,
            1.52822594086518717325715E-5,
            -7.63719763789976227360029E-6,
            3.81729326499983985646164E-6,
            -1.90821271655393892565696E-6,
            9.53962033872796113152039E-7,
            -4.76932986787806463116720E-7,
            2.38450502727732990003648E-7,
            -1.19219925965311073067789E-7,
            5.96081890512594796124402E-8,
            -2.98035035146522801860637E-8,
            1.49015548283650412346585E-8,
            -7.45071178983542949198100E-9,
            3.72533402478845705481920E-9,
            -1.86265972351304900640391E-9,
            9.31327432419668182871765E-10,
            -4.65662906503378407298923E-10,
            2.3283118336765054920E-10,
            -1.1641550172700519776E-10,
            5.8207720879027008892E-11,
            -2.9103850444970996869E-11,
            1.4551921891041984236E-11,
            -7.2759598350574810145E-12,
            3.6379795473786511902E-12,
            -1.8189896503070659476E-12,
            9.0949478402638892825E-13,
            -4.5474737830421540268E-13,
            2.2737368458246525152E-13,
            -1.1368684076802278493E-13,
            5.6843419876275856093E-14,
            -2.8421709768893018555E-14,
            1.4210854828031606770E-14,
            -7.1054273952108527129E-15,
            3.5527136913371136733E-15,
            -1.7763568435791203275E-15,
            8.8817842109308159031E-16,
            -4.4408921031438133642E-16,
            2.2204460507980419840E-16
        };

        // Since
        //  \Gamma(2 + x) = (1 + x) \Gamma(1 + x),
        //  \ln\Gamma(2 + x) = \ln(1 + x) + \Gamma(1 + x)
        // The log1p term precisely cancels this term in the series, so it's
        // actually more accurate and efficient to first define LogGammaTwoPlus
        // and then evaluate LogGammaOnePlus from it.

        private static double ReducedSeries(double x) {
            Debug.Assert(Math.Abs(x) <= Limit);
            double s = 1.0 - AdvancedMath.EulerGamma;
            double xPower = x;
            for (int k = 2; k < altZetaMinusOne.Length; k++) {
                double s_old = s;
                s += altZetaMinusOne[k] * xPower / k;
                if (s == s_old)  return s;
                xPower *= x;
            }
            throw new NonconvergenceException();
        }

        public static double ReducedLogGammaTwoPlus(double x) {
            return ReducedSeries(x);
        }

        public static double ReducedLogGammaOnePlus(double x) {
            return ReducedLogGammaTwoPlus(x) - MoreMath.ReducedLogOnePlus(x);
        }

        private static double ReducedBetaSeries(double x, double y) {
            Debug.Assert(Math.Abs(x) <= Limit && Math.Abs(y) <= Limit && Math.Abs(x + y) <= Limit);
            double z = x + y;
            double xPower = 1.0;
            double yPower = 1.0;
            double q = 2.0;
            double s = 0.0;
            for (int k = 2; k < altZetaMinusOne.Length; k++) {
                double s_old = s;
                double ds = altZetaMinusOne[k] / k * q;
                s += ds;
                if ((k % 2 == 0) && (s == s_old)) return x * y * s;
                xPower *= x;
                yPower *= y;
                q = z * q + (xPower + yPower);
            }
            throw new NonconvergenceException();
        }

        public static double LogGammaTwoPlus(double x) {
            return x * ReducedSeries(x);
        }

        public static double LogGammaOnePlus(double x) {
            return LogGammaTwoPlus(x) - MoreMath.LogOnePlus(x);
        }

        public static double GammaOnePlus(double x) {
            return GammaTwoPlus(x) / (1.0 + x);
        }

        public static double GammaTwoPlus(double x) {
            return Math.Exp(LogGammaTwoPlus(x));
        }

        public static double PsiTwoPlus(double x) {
            double s = 1.0 - AdvancedMath.EulerGamma;
            double xPower = x;
            for (int k = 2; k < altZetaMinusOne.Length; k++) {
                double s_old = s;
                s += altZetaMinusOne[k] * xPower;
                if (s == s_old) {
                    return (s);
                }
                xPower *= x;
            }
            throw new NonconvergenceException();
        }

        public static double PsiOnePlus(double x) {
            return (PsiTwoPlus(x) - 1.0 / (1.0 + x));
        }

        // For Pochhammer symbol, we need S(x + y) - S(x).
        // Each term needs (x + y)^n - x^n, and since x^n cancels, every term has at least one power of y.
        // Let (x + y)^n - x^n = y R_n(x, y). For example
        //   R_1 = 1
        //   R_2 = 2 x - y
        //   R_3 = 3 x^2 + 3 x y + y^2
        // We can derive a recurrence for R's using the binomial recurrence.
        //   R_{n + 1} = (x + y) R_n + x^n
        // We can use this to derive the whole series with a factor y extracted.

        public static double ReducedPochammerSeries(double x, double y) {
            Debug.Assert(Math.Abs(x) <= Limit && Math.Abs(x + y) <= Limit);
            double s = 1.0 - AdvancedMath.EulerGamma;
            double z = x + y;
            double xPower = x;
            double r = 1.0;
            for (int k = 2; k < altZetaMinusOne.Length; k++) {
                double s_old = s;
                r = z * r + xPower;
                double ds = altZetaMinusOne[k] / k * r;
                s += ds;
                if (s == s_old) return s;
                xPower *= x;
            }
            throw new NonconvergenceException();
        }


        public static double LogPochhammerTwoPlus (double x, double y) {
            return y * ReducedPochammerSeries(x, y);
        }

        public static double PochhammerTwoPlus (double x, double y) {
            return Math.Exp(LogPochhammerTwoPlus(x, y));
        }

        public static double LogBetaOnePlus(double x, double y) {
            double s = - MoreMath.LogOnePlus(x) - MoreMath.LogOnePlus(y);
            return s - ReducedBetaSeries(x, y);
        }

        public static double BetaOnePlus (double x, double y) {
            return 1.0 / ((1.0 + x) * (1.0 + y) * Math.Exp(ReducedBetaSeries(x, y)));
        }

        public static double BetaNearOne (double x, double y) {
            Debug.Assert(0.5 <= x && x <= 1.5);
            Debug.Assert(0.5 <= y && y <= 1.5);
            return 1.0 / x / y / Math.Exp(ReducedBetaSeries(x - 1.0, y - 1.0));
        }

        // To compute series for log Beta, we need P_k = (x + y)^k - x^k - y^k for each term.
        // Naive computation would have cancelation error.
        //   P_2 = (x + y)^2 - x^2 - y^2 = 2 x y
        //   P_3 = (x + y)^3 - x^3 - y^3 = 3 x^2 y + 3 x y^2 = 3 x y (x + y)
        //   P_4 = (x + y)^4 - x^4 - y^4 = 4 x^3 y + 6 x^2 y^2 + 4 x y^3 = x y (4 x^2 + 6 x y + 4 y^2)
        // Notice we just get binomial series without first and last terms, so can always extract a factor x y.
        // Define P_k = x y Q_k and derive a recursion for Q_k using the binomial recursion.
        //   Q_{k+1} = (x + y) Q_k + x^{k-1} + y^{k-1}
        // This implies for Qs
        //   Q_2 = 2
        //   Q_3 = 2 (x + y) + x + y = 3 (x + y)
        //   Q_4 = 3 (x + y) (x + y) + x^2 + y^2 = 3 x^2 + 6 x y + 3 y^2 + x^2 + y^2 = 4 x^2 + 6 x y + 4 y^2
        // which agrees with our results for Ps above.

    }
}

