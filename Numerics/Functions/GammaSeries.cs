using System;

namespace Meta.Numerics.Functions {

    // This class evaluates Gamma and related function on the basis of the series
    // DLMF 5.7.3 (https://dlmf.nist.gov/5.7):
    //   \ln\Gamma(1 + x) = -\ln(1 + x) + x (1 - \gamma) +
    //     \sum_{k=2}^{\infty} \frac{\zeta(k) - 1}{k} (-z)^k

    // This is more accurate than Lanczos in the region 1/2 < x < 5/2, where
    // \ln\Gamma has zeros that Lanczos relies on delicate cancelations to reproduce. 
    // Since the series can take up to ~30 terms near half-integers, it's slower
    // than Lanczos for those values. Very close to the integers, though, it's faster
    // than Lanczos.

    internal static class GammaSeries {

        // We record \zeta(k) - 1 explicitly for two reasons. First, it's faster than
        // invoking AdvancedMath.Zeta. Second, \zeta(k) gets very close to 1, and we
        // want to know the difference to high accuracy. An alternative to storing
        // externally computed values is to sum the zeta series without the first unit
        // term; this converges quickly for large k.

        private static double[] zetaMinusOne = new double[]
        {
            -1.5,
            Double.NaN,
            6.44934066848226436472415E-1,
            2.02056903159594285399738E-1,
            8.23232337111381915160037E-2,
            3.69277551433699263313655E-2,
            1.73430619844491397145179E-2,
            8.34927738192282683979755E-3,
            4.07735619794433937868524E-3,
            2.00839282608221441785277E-3,
            9.94575127818085337145959E-4,
            4.94188604119464558702282E-4,
            2.46086553308048298637998E-4,
            1.22713347578489146751836E-4,
            6.12481350587048292585451E-5,
            3.05882363070204935517285E-5,
            1.52822594086518717325715E-5,
            7.63719763789976227360029E-6,
            3.81729326499983985646164E-6,
            1.90821271655393892565696E-6,
            9.53962033872796113152039E-7,
            4.76932986787806463116720E-7,
            2.38450502727732990003648E-7,
            1.19219925965311073067789E-7,
            5.96081890512594796124402E-8,
            2.98035035146522801860637E-8,
            1.49015548283650412346585E-8,
            7.45071178983542949198100E-9,
            3.72533402478845705481920E-9,
            1.86265972351304900640391E-9,
            9.31327432419668182871765E-10,
            4.65662906503378407298923E-10
        };

        // Since
        //  \Gamma(2 + x) = (1 + x) \Gamma(1 + x),
        //  \ln\Gamma(2 + x) = \ln(1 + x) + \Gamma(1 + x)
        // The log1p term precisely cancels this term in the series, so it's
        // actually more accurate and efficient to first define LogGammaTwoPlus
        // and then evaluate LogGammaOnePlus from it.

        private static double Series (double x) {
            double s = 0.0;
            double xMinus = -x;
            double xPower = xMinus;
            for (int k = 2; k < zetaMinusOne.Length; k++) {
                double s_old = s;
                xPower *= xMinus;
                s += zetaMinusOne[k] * xPower / k;
                if (s == s_old) {
                    return (s);
                }
            }
            throw new NonconvergenceException();
        }

        public static double LogGammaTwoPlus (double x) {
            return ((1.0 - AdvancedMath.EulerGamma) * x + Series(x));
        }

        public static double LogGammaOnePlus (double x) {
            return ((1.0 - AdvancedMath.EulerGamma) * x - MoreMath.LogOnePlus(x) + Series(x));
        }

        public static double GammaOnePlus (double x) {
            return (GammaTwoPlus(x) / (1.0 + x));
        }

        public static double GammaTwoPlus (double x) {
            return (Math.Exp(LogGammaTwoPlus(x)));
        }

        public static double PsiTwoPlus (double x) {
            double s = 1.0 - AdvancedMath.EulerGamma;
            double xMinus = -x;
            double xPower = 1.0;
            for (int k = 2; k < zetaMinusOne.Length; k++) {
                double s_old = s;
                xPower *= xMinus;
                s -= zetaMinusOne[k] * xPower;
                if (s == s_old) {
                    return (s);
                }
            }
            throw new NonconvergenceException();
        }

        public static double PsiOnePlus (double x) {
            return (PsiTwoPlus(x) - 1.0 / (1.0 + x));
        }

    }
}

