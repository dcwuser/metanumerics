using System;

namespace Meta.Numerics.Functions
{

    // This class handles the Lanczos approximation to the \Gamma function and the corresponding approximations to associated functions.
    // For basic background to the Lanczos approximation, see http://en.wikipedia.org/wiki/Lanczos_approximation and
    // http://mathworld.wolfram.com/LanczosApproximation.html and http://www.boost.org/doc/libs/1_53_0/libs/math/doc/sf_and_dist/html/math_toolkit/backgrounders/lanczos.html.
    // The basic Lanczos formula is:
    //   \Gamma(z+1) = \sqrt{2 \pi} (z + g + 1/2)^(z+1/2) e^{-(z + g + 1/2)} \left[ c_0 + \frac{c_1}{z+1} + \frac{c_2}{z+2} + \cdots + \frac{c_N}{z+N} \right]
    // Given a value of g, the c-values can be computed using a complicated set of matrix equations that require high precision.
    // We write this as:
    //   \Gamma(z) = \sqrt{2 \pi} (z + g - 1/2)^(z-1/2) e^{-(z + g - 1/2)} \left[ c_0 + \frac{c_1}{z} + \frac{c_2}{z+1} + \cdots + \frac{c_N}{z+N-1} \right]
    //             = \sqrt{2 \pi} (\frac{z + g - 1/2}{e})^{z-1/2} e^{-g} \left[ c_0 + \frac{c_1}{z} + \frac{c_2}{z+1} + \cdots + \frac{c_N}{z+N-1} \right]

    internal static class Lanczos {

        // These are listed at http://www.mrob.com/pub/ries/lanczos-gamma.html as the values used by GSL, although I don't know if they still are.
        // Measured deviations at integers 2, 2, 4, 11, 1, 17, 22, 21 X 10^(-16) so this is clearly worse than Godfrey's coefficients, although it
        // does manage with slightly fewer terms.
        /*
        private const double LanczosG = 7.0;
        private static readonly double[] LanczosC = new double[] {
            0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
        };
        */

        // Godfrey's coefficients, claimed relative error < 10^(-15), documented at http://my.fit.edu/~gabdo/gamma.txt and in NR 3rd edition section 6.1.
        // Measured relative deviation at integers 1, 1, 4, 1, 4, 5, 6, 3 X 10^(-16) so this appears about right.
        // These improves to 1, 1, 2, 1, 3, 3, 3 X 10^(-16) when we pull the 1/e into Math.Pow(t/e, z+1/2) instead of calling Math.Exp(-t) separately.

        private const double LanczosG = 607.0 / 128.0;

        private static readonly double[] LanczosC = new double[] {
            0.99999999999999709182,
            57.156235665862923517,
            -59.597960355475491248,
            14.136097974741747174,
            -0.49191381609762019978,
            3.3994649984811888699e-5,
            4.6523628927048575665e-5,
            -9.8374475304879564677e-5,
            1.5808870322491248884e-4,
            -2.1026444172410488319e-4,
            2.1743961811521264320e-4,
            -1.6431810653676389022e-4,
            8.4418223983852743293e-5,
            -2.6190838401581408670e-5,
            3.6899182659531622704e-6
        };


        // These coefficients are given by Pugh in his thesis (http://web.viu.ca/pughg/phdThesis/phdThesis.pdf) table 8.5, p. 116
        // (We record them in his form; to get the usual form, multiply by Exp(-\gamma+1/2) \sqrt{2} / \pi.)
        // He claims that they "guarantee 16 digit floating point accuracy in the right-half plane", and since he uses fewer coefficients
        // than Godfrey that would be fantastic. But we measure relative deviations at the integers of 2, 0, 35, 23, 45, 42, 6 X 10^(-16),
        // making this relatively bad.
        // Unfortunately, we didn't do these measurements ourselves at first, so we actually used these coefficients until version 3.
        // Perhaps this really does give 53 bits of accuracy if you do the calculation with more bits. The fact that G is relatively large
        // makes the initial coefficients relatively large, which probably leads to cancellation errors.
        /*
        private const double LanczosG = 10.900511;

        private static readonly double[] LanczosC = new double[] {
            +2.48574089138753565546e-5,
            1.05142378581721974210,
            -3.45687097222016235469,
            +4.51227709466894823700,
            -2.98285225323576655721,
            +1.05639711577126713077,
            -1.95428773191645869583e-1,
            +1.70970543404441224307e-2,
            -5.71926117404305781283e-4,
            +4.63399473359905636708e-6,
            -2.71994908488607703910e-9,
        };
        */

        // From LanczosG, we derive several values that we need only compute once.

        private const double LanczosGP = LanczosG - 0.5;
        private static readonly double LanczosExpG = Math.Exp(-LanczosG);
        private static readonly double LanczosExpGP = Math.Exp(-LanczosGP);

        private static double Sum(double x) {
            double s = LanczosC[0] + LanczosC[1] / x;
            for (int i = 2; i < LanczosC.Length; i++) {
                x += 1.0;
                s += LanczosC[i] / x;
            }
            return (s);
        }

        private static Complex Sum(Complex z) {
            Complex s = LanczosC[0] + LanczosC[1] / z;
            for (int i = 2; i < LanczosC.Length; i++) {
                z += 1.0;
                s += LanczosC[i] / z;
            }
            return (s);
        }

        private static double LogSumPrime(double x) {
            double q = LanczosC[0] + LanczosC[1] / x;
            double p = LanczosC[1] / (x * x);
            for (int i = 2; i < LanczosC.Length; i++) {
                x += 1.0;
                q += LanczosC[i] / x;
                p += LanczosC[i] / (x * x);
            }
            return (-p / q);
        }

        private static Complex LogSumPrime(Complex z) {
            Complex q = LanczosC[0] + LanczosC[1] / z;
            Complex p = LanczosC[1] / (z * z);
            for (int i = 2; i < LanczosC.Length; i++) {
                z += 1.0;
                q += LanczosC[i] / z;
                p += LanczosC[i] / (z * z);
            }
            return (-p / q);
        }

        // To calculate \frac{\Gamma(x + y)}{\Gamma(x)}, we need \frac{S(x+y)}{S(x)}. A little bit of algebra shows
        //   \frac{S(x+y)}{S(x)} = 1 - y \left[ \sum_{k=1} \frac{c_k}{(x + y + (k-1))(x + (k-1))} \right]
        //                               \left[ c_0 + \sum_{k=1} \frac{c_k}{x + (k-1)} \right]
        // Expressing the result as the rate at which the ratio moves away from 1 as y increases allows us to
        // derive accurate results for small y.

        private static double RatioOfSumsResidual(double x, double y) {
            double z = x + y;
            double p = LanczosC[1] / (x * z);
            double q = LanczosC[0] + LanczosC[1] / x;
            for (int i = 2; i < LanczosC.Length; i++) {
                x += 1.0;
                z += 1.0;
                p += LanczosC[i] / (x * z);
                q += LanczosC[i] / x;
            }
            return (p / q);
        }

        public static double Gamma(double x) {
            double t = x + LanczosGP;
            return (
                Global.SqrtTwoPI *
                Math.Pow(t / Math.E, x - 0.5) * LanczosExpG *
                Sum(x)
            );
        }

        public static Complex Gamma (Complex z) {
            Complex t = z + LanczosGP;
            return (
                Global.SqrtTwoPI *
                ComplexMath.Pow(t / Math.E, z - 0.5) * LanczosExpG *
                Sum(z)
            );
        }

        public static double LogGamma(double x) {
            double t = x + LanczosGP;
            return (
                Math.Log(Global.SqrtTwoPI * Sum(x)) +
                (x - 0.5) * Math.Log(t) - t
            );
        }

        public static Complex LogGamma(Complex z) {
            Complex t = z + LanczosGP;
            return (
                Math.Log(Global.SqrtTwoPI) +
                (z - 0.5) * ComplexMath.Log(t) - t +
                ComplexMath.Log(Sum(z))
            );
        }

        public static double Psi(double x) {
            double t = x + LanczosGP;
            return (Math.Log(t) - LanczosG / t + LogSumPrime(x));
        }

        public static Complex Psi(Complex z) {
            Complex t = z + LanczosGP;
            return (ComplexMath.Log(t) - LanczosG / t + LogSumPrime(z));
        }

        public static double ReducedLogPochhammer(double x, double y) {
            double r = MoreMath.ReducedLogOnePlus(-RatioOfSumsResidual(x, y), y) +
                (x - 0.5) * MoreMath.ReducedLogOnePlus(1.0 / (x + LanczosGP), y) +
                Math.Log(x + y + LanczosGP) - 1.0;
            return (r);
        }

        // If we just compute Exp( LogGamma(x) + LogGamma(y) - LogGamma(x+y) ) then several leading terms in the sum cancel,
        // potentially introducing cancelation error. So we write out the ratios explicitly and take the opportunity
        // to write the result in terms of some naturally occurring ratios.

        public static double Beta(double x, double y) {
            double tx = x + LanczosGP;
            double ty = y + LanczosGP;
            double txy = x + y + LanczosGP;
            return (
                Global.SqrtTwoPI * LanczosExpGP *
                Math.Pow(tx / txy, x) * Math.Pow(ty / txy, y) * Math.Sqrt(txy / tx / ty) *
                Sum(x) * Sum(y) / Sum(x + y)
            );
        }

        public static double LogBeta(double x, double y) {
            double tx = x + LanczosGP;
            double ty = y + LanczosGP;
            double txy = x + y + LanczosGP;
            // For very small x, Sum(x) explodes as 1/x. So for x or y < ~1.0E-150, Sum(x) * Sum(y) / Sum(x + y)
            // would suffer intermediate overflow, but / Sum(x + y) * Sum(x) * Sum(y) does not. There is no
            // corresponding problem for very large x, for sum Sum(x) -> ~1.0.
            return (
                0.5 * Math.Log(2.0 * Math.PI / txy) + (x - 0.5) * Math.Log(tx / txy) + (y - 0.5) * Math.Log(ty / txy) +
                Math.Log(LanczosExpGP / Sum(x + y) * Sum(x) * Sum(y))
            );
        }

    }

}
