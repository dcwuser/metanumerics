using System;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Represents the distribution of the Kolmogorov-Smirnov D statistic.
    /// </summary>
    /// <remarks><para>The D statistic in a Kolmogorov-Smirnov test is distributed (under the null hypothesis) according to a Kolmogorov disribution.</para></remarks>
    /// <seealse cref="Sample.KolmogorovSmirnovTest(Meta.Numerics.Statistics.Distributions.Distribution)" />
    public class KolmogorovDistribution : Distribution {

        /// <summary>
        /// Instantiates a new asymptotic Kolmogorov distribution.
        /// </summary>
        public KolmogorovDistribution () { }

        // the sample size; when N=0 we will report the asymptotic distribution

        internal KolmogorovDistribution (double scale) {
            this.scale = scale;
        }

        private double scale = 1.0;

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {

            if (x < scale) {
                return (AsymptoticPPrime(x/scale)/scale);
            } else {
                return (AsymptoticQPrime(x/scale)/scale);
            }

        }

        // the asymptotic PDF for x <~ 1

        private static double AsymptoticPPrime (double x) {

            if (x <= 0.0) return (0.0);

            double p = 0.0;
            for (int k = 1; k < Global.SeriesMax; k += 2) {
                double p_old = p;
                double z = k * Math.PI / x / 2.0;
                double dp = Math.Exp(-z * z / 2.0) * (z * z - 1.0);
                p += dp;
                if (p == p_old) return (Global.SqrtTwoPI / (x * x) * p);
            }

            throw new NonconvergenceException();
        }

        // the asymptotic PDF for x >~ 1

        private static double AsymptoticQPrime (double x) {

            double p = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double p_old = p;
                double kx = k * x;
                double dp = k * k * Math.Exp(-2.0 * kx * kx);
                if (k % 2 == 0) dp = -dp;
                p += dp;
                if (p == p_old) return (8.0 * p * x);
            }

            throw new NonconvergenceException();
        }


        /// <inheritdoc />
        public override double LeftProbability (double x) {

            if (x < scale) {
                return (AsymptoticP(x / scale));
            } else {
                return (1.0 - AsymptoticQ(x / scale));
            }

        }

        /// <inheritdoc />
        public override double RightProbability (double x) {

            if (x < scale) {
                return (1.0 - AsymptoticP(x / scale));
            } else {
                return (AsymptoticQ(x / scale));
            }

        }

        // implements \frac{\sqrt{2\pi}}{x1} \sum{k=0}{\infty} e^{ \frac{(2k+1)^2 \pi^2}{8 x1^2} }
        // convergence is rapid; 4 terms at x~1 and still just 10 terms at x~3

        private static double AsymptoticP (double x) {

            if (x <= 0.0) {
                return (0.0);
            } else {

                double p = 0.0;
                for (int k = 1; k < Global.SeriesMax; k += 2) {
                    double p_old = p;
                    double z = k * Math.PI / x / 2.0;
                    double dp = Math.Exp(-z * z / 2.0);
                    p += dp;
                    if (p == p_old) return (Global.SqrtTwoPI / x * p);
                }

                throw new NonconvergenceException();
            }
        }

        // implements \sum_{k=-\infty}^{\infty} (-1)^k e^{-2 k^2 x1^2}
        // convergence is very rapid; 5 terms at x~1 and just 2 terms at x~3

        private static double AsymptoticQ (double x) {
            double xx = x * x;
            double f = 0.0;
            int sign = -1;
            for (int k = 1; k < Global.SeriesMax; k++) {
                double f_old = f;
                sign = -sign;
                double df = sign * Math.Exp(-(2 * k * k) * xx);
                f = f_old + df;
                if (f == f_old) {
                    return (2.0 * f);
                }
            }
            throw new NonconvergenceException();
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                return (Math.Sqrt(Global.HalfPI) * Global.LogTwo * scale);
            }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                double ln2 = Global.LogTwo;
                return (Math.PI / 2.0 * (Math.PI / 6.0 - ln2 * ln2) * scale * scale);

            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
                // this constant was determined empiricaly
                return (0.82757355518991 * scale);
            }
        }

        /// <inheritdoc />
        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else if (n == 1) {
                return (Mean);
            } else if (n == 2) {
                return (Math.PI * Math.PI / 12.0 * scale * scale);
            } else {
                return (AdvancedMath.Gamma(n / 2.0 + 1.0) * AdvancedMath.DirichletEta(n) / Math.Pow(2.0, n / 2.0 - 1.0) * Math.Pow(scale, n));
            }
        }

        /// <inheritdoc />
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
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /*

        // an exact formula for 1/2 <= t <= 1

        private static double Smallest_P (int n, double t) {

            Debug.Assert((0.5 <= t) && (t <= 1.0));

            return (Math.Exp(n * Math.Log((2.0 * t - 1.0) / n) + AdvancedIntegerMath.LogFactorial(n)));

        }

        // an exact formula for n-1 <= t <= n

        private static double Smallest_Q (int n, double t) {

            Debug.Assert((n - 1) <= t);

            return (2.0 * Math.Pow(1.0 - t / n, n));

        }
        */

    }

    // we will build this assuming n <~ 50, because it gets very ugly with increasing n

    public class PolynomialFiniteKolmogorovDistribution : Distribution {

        public PolynomialFiniteKolmogorovDistribution (int n) {
            if (n < 1) throw new ArgumentOutOfRangeException("n");
            this.n = n;
        }

        private int n;

        // extreme left value for 1/2 < t < 1

        private static double P0 (int n, double t) {
            return (AdvancedIntegerMath.Factorial(n) * MoreMath.Pow((2.0 * t - 1.0) / n, n));
        }

        // extreme right value for n / 2 < t < n

        private static double Q0 (int n, double t) {

            double Q = 2.0 * MoreMath.Pow((n - t) / n, n);

            if (t > n - 1) return (Q);

            Q += 2.0 * t * MoreMath.Pow((n - 1 - t) / n, n - 1);

            if (t > n - 2) return (Q);

            for (int j = 2; n - j - t > 0.0; j++) {
                Q += 2.0 * AdvancedIntegerMath.BinomialCoefficient(n, j) * (t / n) * MoreMath.Pow((t + j) / n, j - 1) * MoreMath.Pow((n - j - t) / n, n - j);
            }

            return (Q);

        }

        public override double ProbabilityDensity (double x) {
            double t = x;
            if (t < 0.5) {
                return (0.0);
            } else if (t < 1.0) {
                throw new NotImplementedException();
                // very easy, derivate of P0
            } else if (t < n / 2) {
                // hard
                throw new NotImplementedException();
            } else if (t < n) {
                // easy, derivate of -Q0
                throw new NotImplementedException();
            } else {
                return (0.0);
            }
        }

        public override double LeftProbability (double x) {
            double t = x;
            if (t < 0.5) {
                return (0.0);
            } else if (t < 1.0) {
                return (P0(n, t));
            } else if (t < n / 2) {
                // hard
                throw new NotImplementedException();
            } else if (t < n) {
                return (1.0 - Q0(n, t));
            } else {
                return (1.0);
            }
        }

        public override double RightProbability (double x) {
            double t = x;
            if (t < 0.5) {
                return (1.0);
            } else if (t < 1.0) {
                return (1.0 - P0(n, t));
            } else if (t < n / 2) {
                // hard
                throw new NotImplementedException();
            } else if (t < n) {
                return (Q0(n, t));
            } else {
                return (0.0);
            }
        }

    }


    public class AsymptoticFiniteKolmogorovDistribution : Distribution {

        // This class uses an expansion of P(x) in powers of 1/\sqrt{n} that is explored in
        //   Pelz and Good, Journal of the Royal Statistical Society. Series B (Methodological), Vol. 38, No. 2 (1976), pp. 152-156
        //   Simard and L'Ecuyer, Journal of Statistical Software
        // By differentiating the expressions, we obtain p(x). By integrating p(x) x^m, we obtain expressions for the moments.
        // Thus the approximation is self-consistent.

        // For very small n and very extreme values of x, it can give negative p(x), but in practice this is extremely unlikely.
        // This is a problem well known for Edgeworth expansions, and is really unavoidable for any expansion of p(x). Since the
        // leading term must integrate to one, and the whole series must also integrate to one, higher terms must integrate to zero.
        // Hence they must have regions of negative values.

        public AsymptoticFiniteKolmogorovDistribution (int n) {
            if (n < 1) throw new ArgumentOutOfRangeException("n");
            this.n = n;
        }

        int n;

        // Asymptotic result as series useful for small x is
        //   P(x) = \frac{\sqrt{2 \pi}}{x} \sum_{j=1}^{\infty} e^{-\frac{(2j-1)^2 \pi^2}{8x^2}}
        // This can derived from the Q-result below by writing the series as a Jacobi theta function
        // and using the \tau \rightarrow 1/\tau transformation for that function (see http://en.wikipedia.org/wiki/Theta_function)

        private double P0 (double x) {
            double a = MoreMath.Pow2(Math.PI / x) / 8.0;
            double f = 0.0;
            for (int k = 1; k < 100; k += 2) {
                double f_old = f;
                double df = Math.Exp(-a * (k * k));
                f += df;
                if (f == f_old) { return (Math.Sqrt(2.0 * Math.PI) / x * f); }
            }
            throw new NonconvergenceException();
        }

        // Asymptotic result as series useful for large x is
        //  P = \sum_{k=-\infty}^{\infty} (-1)^k e^{-2 k^2 x^2}
        //    = 1 - 2 \sum_{k=1}^{\infty} (-1)^{k-1} e^{-2 k^2 x^2}
        // and since P = 1 - Q the second term is Q.

        private double Q0 (double x) {
            double a = 2.0 * x * x;
            double f = 0.0;
            for (int k = 1; k < 50; k++) {
                double f_old = f;
                double df = Math.Exp(-a * (k * k));
                if (k % 2 == 0) df = -df;
                f += df;
                if (f == f_old) { return (2.0 * f); }
            }
            throw new NonconvergenceException();
        }

        // Derivative of P_0, which is leading term in p(x) = P'(x) useful for small x
        //   P' = \frac{\sqrt{2\pi}}{x^2} \sum_{k=1}^{\infty} \left( \frac{(2k-1)^2 \pi^2}{4 x^2} - 1 \right) e^{-\frac{(2k-1)^2 \pi^2}{8 x^2}}
        // This is derived straightforwardly by term-by-term differentiation of the P0 series

        private double P0Prime (double x) {
            double a = MoreMath.Pow2(Math.PI / x) / 8.0;
            double f = 0.0;
            for (int k = 1; k < 100; k += 2) {
                double f_old = f;
                double z = a * (k * k);
                double df = (2.0 * z - 1.0) * Math.Exp(-z);
                f += df;
                if (f == f_old) return (Global.SqrtTwoPI / (x * x) * f);
            }
            throw new NonconvergenceException();
        }

        // Derivative of Q_0, which is the leading term in p(x) = -Q'(x) useful for large x
        //   Q' = 8 x \sum_{k=1}^{\infty} (-1)^k k^2 e^{-2 k^2 x^2}
        // This is derived straightforwardly by term-by-term differentiation of the Q_0 series

        private double Q0Prime (double x) {
            double a = 2.0 * x * x;
            double f = 0.0;
            for (int k = 1; k < 50; k++) {
                double f_old = f;
                int k2 = k * k;
                double df = k2 * Math.Exp(-a * k2);
                if (k % 2 != 0) df = -df;
                f += df;
                if (f == f_old) return (8.0 * x * f);
            }
            throw new NonconvergenceException();
        }

        // The next terms in the expansion give the 1/\sqrt{n} corrections
        // It happens that P_1 = P_0' / 6, so use this to reduce the code
        // we have to write and maintain

        private double P1 (double x) {
            return (P0Prime(x) / 6.0);
        }

        private double Q1 (double x) {
            return (Q0Prime(x) / 6.0);
        }

        // Derivative of Q_1, which is the first correction term to p(x) = -Q'(x)
        //    Q_1' = \frac{4}{3} \sum_{k=1}^{\infty} (-1)^k k^2 ( 1 - 4 k^2 x^2 ) e^{-2 k^2 x^2}
        // This form, useful for large x, is obtained by straightforward differentiation of the expression
        // for Q_1 = Q_0' / 6

        private double Q1Prime (double x) {
            double a = 2.0 * x * x;
            double f = 0.0;
            for (int k = 1; k < 50; k++) {
                double f_old = f;
                int k2 = k * k;
                double z = a * k2;
                double df = k2 * (1.0 - 2.0 * z) * Math.Exp(-z);
                if (k % 2 != 0) df = - df;
                f += df;
                if (f == f_old) return (4.0 / 3.0 * f);
            }
            throw new NonconvergenceException();
        }

        private double P1Prime (double x) {
            double f = 0.0;
            for (int k = 1; k < 100; k += 2) {
                double f_old = f;
                double z = (k * k) * (Math.PI * Math.PI) / (x * x) / 8.0;
                double df = (4.0 * z * z - 10 * z + 2.0) * Math.Exp(-z);
                f += df;
                if (f == f_old) return (Global.SqrtTwoPI / 6.0 / (x * x * x) * f);
            }
            throw new NotImplementedException();
        }

        public override double LeftProbability (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x < 1.0) {
                return (P0(x) + P1(x) / Math.Sqrt(n));
            } else {
                return (1.0 - Q0(x) - Q1(x) / Math.Sqrt(n));
            }
        }

        public override double RightProbability (double x) {
            if (x <= 0.0) {
                return (1.0);
            } else if (x < 1.0) {
                return (1.0 - P0(x) - P1(x) / Math.Sqrt(n));
            } else {
                //return (KQ0(x));
                return (Q0(x) + Q1(x) / Math.Sqrt(n));
            }
        }

        public override double ProbabilityDensity (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else if (x < 1.0) {
                return (P0Prime(x) + P1Prime(x) / Math.Sqrt(n));
            } else {
                return (-Q0Prime(x) - Q1Prime(x) / Math.Sqrt(n));
            }
        }

        public override double Mean {
            get {
                return (Math.Sqrt(Global.HalfPI) * Global.LogTwo - 1.0 / 6.0 / Math.Sqrt(n));
            }
        }

        public override double Variance {
            get {
                return (Math.PI / 2.0 * (Math.PI / 6.0 - Global.LogTwo * Global.LogTwo) - Math.Sqrt(Global.HalfPI) * Global.LogTwo / 12.0 / Math.Sqrt(n));
            }
        }

        public override double Moment (int m) {
            if (m < 0) {
                throw new ArgumentOutOfRangeException("m");
            } else if (m == 0) {
                return (1.0);
            } else {
                // we can get these expressions just by integrating Q_0' and Q_1' term by term
                double M0 = AdvancedMath.DirichletEta(m) * AdvancedMath.Gamma(m / 2.0 + 1.0) / Math.Pow(2.0, m / 2.0 - 1.0);
                double M1 = -2.0 / 3.0 * AdvancedMath.DirichletEta(m - 1) * AdvancedMath.Gamma((m + 1) / 2.0) / Math.Pow(2.0, (m + 1) / 2.0) * m;
                return (M0 + M1 / Math.Sqrt(n));
            }
        }

    }

#if FUTURE
    public class KolmogorovExactDistribution : Distribution {

        public KolmogorovExactDistribution (int size) {
            if (size < 1) throw new ArgumentOutOfRangeException("size");

            N = size;
            sqrtN = Math.Sqrt(N);
        }

        int N;
        double sqrtN;

        int maxN = 256;

        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.5 / N, 1.0));
            }
        }

        public override double ProbabilityDensity (double d) {

            double t = d * N;
            if (2.0 * t < N) {
                return (DurbinPPrime(N, t) * N);
            } else {
                return (DurbinQPrime(N, t) * N);
            }

        }

        /// <inheritdoc />
        public override double LeftProbability (double d) {

            if (d <= 0.5 / N) {
                return (0.0);
            } else if (d >= 1.0) {
                return (1.0);
            } else {

                // use Durbin formulas for small N

                double t = d * N;
                if (2.0 * t < N) {
                    // the Durbin series formula is faster than the durbin matrix formula,
                    // but it has alternating sign terms and and can suffer a catastrophic
                    // loss of accuracy
                    return (MatrixP(N, t));
                    //return (DurbinP(N, t));
                } else {
                    return (1.0 - DurbinQ(N, t));
                }


            }

        }

        /// <inheritdoc />
        public override double RightProbability (double d) {

            if (d <= 0.5 / N) {
                return (1.0);
            } else if (d >= 1.0) {
                return (0.0);
            } else {

                // use Durbin formulas for small N

                double t = d * N;
                if (2.0 * t < N) {
                    return (1.0 - MatrixP(N, t));
                    //return (1.0 - DurbinP(N, t));
                } else {
                    return (DurbinQ(N, t));
                }

            }

        }

        // Durbin's formula for exact P_n(t) that holds for n > Truncate(2 t), i.e. small t
        // we holding an array of P_m(t) for m < n; this keeps us having to reevaluate repeatedly

        // the formula unfortunately involves canceling terms of increasing size;
        // it breaks down when n gets large

        private static double DurbinP (int n, double t) {

            if (t <= 0.5) return (0.0);

            int t2 = (int) Math.Truncate(2.0 * t);

            //double s = 0.0;
            //for (int j = 1; j <= t2; j++) {
            //    double ds = Math.Exp(
            //        j * Math.Log(2.0 * t - j) - AdvancedIntegerMath.LogFactorial(j) +
            //        (n - j) * Math.Log(n - j) - AdvancedIntegerMath.LogFactorial(n - j)
            //    ) * DurbinP1(n - j, t);
            //    if (j % 2 == 0) ds = - ds;
            //    s += ds;
            //}
            //s = Math.Exp(AdvancedIntegerMath.LogFactorial(n) - n * Math.Log(n)) * s;
            //return (s);


            // an array for P_m(t)
            double[] P = new double[n+1];
            P[0] = 1.0;

            // populate up to m = t2 using the Durbin Q formula
            for (int m = 1; m <= t2; m++) {
                P[m] = 1.0 - DurbinQ(m, t);
            }

            // populate higher m using the recurrsion relation
            for (int m = t2 + 1; m <= n; m++) {
                double B = 1.0; // binomial coefficient (n j)
                double s = 0.0;
                for (int j = 1; j <= t2; j++) {
                    B = B * (m - (j - 1)) / j;
                    double C = Math.Pow((2.0 * t - j) / m, j) * Math.Pow(1.0 * (m - j) / m, m - j);
                    double ds = B * C * P[m - j];
                    if (j % 2 == 0) ds = -ds;
                    s += ds;
                }
                P[m] = s;
            }

            /*
            for (int i = 0; i <= n; i++) {
                Console.WriteLine("p[{0}]={1}", i, P[i]);
            }
            */

            return (P[n]);

        }

        // Durbin's formula for exact Q_n(t) that holds for 2t > n

        private static double DurbinQ (int n, double t) {

            if (t >= n) return (0.0);

            double s = Math.Pow(1.0 - t / n, n); // j = 0 term
            int jmax = (int) Math.Truncate(n - t);
            double B = 1.0; // binomial coefficient ( n j )
            for (int j = 1; j <= jmax; j++) {
                B = B * (n - (j-1)) / j;
                double C = Math.Pow((t + j) / n, j - 1) * Math.Pow((n - j - t) / n, n - j) * t / n;
                s += B * C;
            }
            return (2.0 * s);

        }

        private static double DurbinQPrime (int n, double t) {

            if (t >= n) return (0.0);

            double s = 2.0 * Math.Pow(1.0 - t / n, n - 1);
            int jmax = (int) Math.Truncate(n - t);
            if (jmax == (n - t)) jmax--;
            for (int j = 1; j <= jmax; j++) {
                double B = Math.Exp(AdvancedIntegerMath.LogFactorial(n-1) - AdvancedIntegerMath.LogFactorial(j) - AdvancedIntegerMath.LogFactorial(n - j));
                double C = Math.Pow((t + j) / n, j - 1) * Math.Pow((n - j - t) / n, n - j);
                double D = 1.0 + (j-1) * t / (t + j) - (n-j) * t / (n - j - t);
                s += -2.0 * B * C * D;
            }
            return (s);
        }

        private static double DurbinPPrime (int n, double t) {

            int t2 = (int) Math.Truncate(2.0 * t);

            // arrays for PDF and CDF for m <= n
            double[] p = new double[n+1];
            double[] P = new double[n+1];
            p[0] = 0.0;
            P[0] = 1.0;

            // populate up to m = t2 using the Durbin Q formula
            for (int m = 1; m <= t2; m++) {
                P[m] = 1.0 - DurbinQ(m, t);
                p[m] = DurbinQPrime(m, t);
            }

            // compute higher p and P using the recursion formula
            for (int m = t2 + 1; m <= n; m++) {
                double B = 1.0; // binomial coefficient (m j)
                double s = 0.0;
                double sp = 0.0;
                double jmax = t2;
                if (t2 == 2.0 * t) jmax--;
                for (int j = 1; j <= jmax; j++) {
                    B = B * (m - (j - 1)) / j;
                    double C = Math.Pow((2.0 * t - j) / m, j) * Math.Pow(1.0 * (m - j) / m, m - j);
                    double ds = B * C * P[m - j];
                    double dsp = B * C * (p[m - j] + 2.0 * j / (2.0 * t - j) * P[m - j]);
                    if (j % 2 == 0) {
                        ds = - ds;
                        dsp = - dsp;
                    }
                    s += ds;
                    sp += dsp;
                }
                P[m] = s;
                p[m] = sp;
            }

            /*
            for (int i = 0; i <= n; i++) {
                Console.WriteLine("P[{0}]={1} p[{0}]={2}", i, P[i], p[i]);
            }
            */

            // return the desired PDF value
            return (p[n]);

        }

        // Durbin's matrix form, also programed by Marsaglia
        // all number are positive, so this does not suffer from the cancelation probmlems of Durbin's recursion

        private static double MatrixP (int n, double t) {

            // compute stuff used in matrix entries
            int tp = (int) Math.Truncate(t) + 1;
            double h = tp - t;
            int p = 2 * tp - 1;


            // construct the matrix
            SquareMatrix H = new SquareMatrix(p);

            // superdiagonal
            for (int j = 1; j < p; j++) {
                H[j-1,j] = 1.0;
            }

            // diagonal and subdiagonals
            double F = 1.0; // factorial
            double hh = h; // power of h
            for (int i = 1; i < p; i++) {
                H[i - 1, 0] = (1.0 - hh) / F;
                H[p-1, p-i] = H[i-1,0];
                for (int j = i+1; j < p; j++) {
                    H[j - 1, j - i] = 1.0 / F;
                }
                hh = hh * h;
                F = F * (i+1);
            }

            // lower-left element
            double g = 1.0 - 2.0 * hh;
            if (h > 0.5) g = g + Math.Pow(2.0 * h - 1.0, p);
            g = g / F;
            H[p-1,0] = g;

            // raise the matrix to the nth power
            SquareMatrix HN = MatrixPower(H, n);

            // return the appropriate element
            double hf = Math.Exp(AdvancedIntegerMath.LogFactorial(n) - n * Math.Log(n));
            return (hf * HN[tp-1, tp-1]);


        }

        private static SquareMatrix MatrixPower (SquareMatrix A, int n) {

            SquareMatrix B = null;

            SquareMatrix D = A.Clone();

            while (true) {
                if (n % 2 != 0) {
                    if (B == null) {
                        B = D.Clone();
                    } else {
                        B = B * D;
                    }
                }
                n = n / 2;
                if (n == 0) break;
                D = D * D;
            }

            return (B);


        }

        ///<inheritdoc />
        public override double Mean {
            get {
                if (N < maxN) {

                    if (N == 1) {
                        return ((3.0 / 4.0) / N);
                    } else if (N == 2) {
                        return ((13.0 / 12.0) / N);
                    } else if (N == 3) {
                        return ((293.0 / 216.0) / N);
                    }

                    throw new NotImplementedException();
                } else {
                    return (base.Mean / sqrtN);
                }
            }
        }

        ///<inheritdoc />
        public override double Variance {
            get {
                if (N < maxN) {
                    if (N == 1) {
                        return (1.0 / 48.0 / N);
                    } else if (N == 2) {
                        return ((7.0 / 72.0) / N);
                    }
                    throw new NotImplementedException();
                } else {
                    return (base.Variance / N);
                }
            }
        }

        ///<inheritdoc />
        public override double Moment (int n) {
            throw new NotImplementedException();
        }

        ///<inheritdoc />
        public override double MomentAboutMean (int n) {
            throw new NotImplementedException();
        }

    }

#endif
}
