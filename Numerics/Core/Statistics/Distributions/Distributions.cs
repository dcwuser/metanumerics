using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// Represents a probability distribution.
    /// </summary>
	public abstract class Distribution {

        /// <summary>
        /// Returns the probability density at the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The probability density p(x1).</returns>
		public abstract double ProbabilityDensity (double x);

        /// <summary>
        /// Returns the cumulative probability to the left of (below) the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The integrated probability P(x1) to obtain a result below the reference point.</returns>
		public abstract double LeftProbability (double x);

        /// <summary>
        /// Return the cumulative probability to the right of (above) the given point.
        /// </summary>
        /// <param name="x">The reference point.</param>
        /// <returns>The integrated probability 1-P(x1) to obtain a result above the reference point.</returns>
		public virtual double RightProbability (double x) {
			return( 1.0 - LeftProbability(x) );
		}

        /// <summary>
        /// Returns the point at which the cumulative distribution function attains a given value. 
        /// </summary>
        /// <param name="P">The left cumulative probability P, which must lie between 0 and 1.</param>
        /// <returns>The point x1 at which the left cumulative probability attains the value P.</returns>
		public virtual double InverseLeftProbability (double P) {
            // find x where LeftProbability(x) = P 
            if ((P < 0.0) || (P > 1.0)) throw new ArgumentOutOfRangeException("P");
            Function<double,double> f = delegate(double x) {
                return(LeftProbability(x) - P);
            };
            double y = FunctionMath.FindZero(f, Mean);
            return(y);
            // since we have the PDF = CDF', change to a method using Newton's method
		}

        /// <summary>
        /// Returns the given moment of the distribution.
        /// </summary>
        /// <param name="n">The order of the moment to determine.</param>
        /// <returns>The moment M<sub>n</sub> about the origin.</returns>
        /// <seealso cref="MomentAboutMean"/>
		public virtual double Moment (int n) {
			// inheritors may implement
			throw new NotImplementedException();
		}

        /// <summary>
        /// Returns the given moment of the distribution, about the mean. 
        /// </summary>
        /// <param name="n">The order of the moment to determine.</param>
        /// <returns>The moment of order n about the mean.</returns>
        /// <seealso cref="Moment" />
		public virtual double MomentAboutMean (int n) {
			// inheritors may implement
			throw new NotImplementedException();
		}

        internal virtual double Cumulant (int n) {
            // inheritors may implement
            throw new NotImplementedException();
        }

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
		public virtual double Mean {
			get {
				return(Moment(1));
			}
		}

        /// <summary>
        /// Gets the standard deviation of the distribution.
        /// </summary>
		public virtual double StandardDeviation {
			get {
				return( Math.Sqrt(Variance) );
			}
		}

        /// <summary>
        /// Gets the variance of the distribution.
        /// </summary>
        public virtual double Variance {
            get {
                return (MomentAboutMean(2));
            }
        }

        /// <summary>
        /// Ges the skewness of the distribution.
        /// </summary>
        /// <remarks>
        /// <para>The skewness of a distribution is a measurement of its asymmetry about its mean.
        /// It is the third moment about the mean, measured in units of the cubed standard deviation.</para>
        /// </remarks>
        public virtual double Skewness {
            get {
                return (MomentAboutMean(3) / Math.Pow(MomentAboutMean(2), 3.0 / 2.0));
            }
        }

        /// <summary>
        /// Gets the median of the distribution.
        /// </summary>
        /// <remarks>The median is the point with equal integrated probability above and below, i.e. with P(x1) = 0.5.</remarks>
        public virtual double Median {
            get {
                return (InverseLeftProbability(0.5));
            }
        }

        /// <summary>
        /// Gets the interval over which the distribution is nonvanishing.
        /// </summary>
        public virtual Interval Support {
            get {
                return (Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity));
            }
        }

        /// <summary>
        /// Computes the expectation value of the given function.
        /// </summary>
        /// <param name="f">The function.</param>
        /// <returns>The expectation value of the function.</returns>
        public virtual double ExpectationValue (Function<double, double> f) {
            return (FunctionMath.Integrate(f, Support));
        }

        /// <summary>
        /// Returns a random value.
        /// </summary>
        /// <param name="rng">A random number generator.</param>
        /// <returns>A number distributed according to the distribution.</returns>
        /// <remarks>
        /// <para>Note that the random number generator <paramref name="rng"/> will be advanced by this method. The next call to its
        /// generator methods will not give the same value as it would had it not been passed to this method.</para>
        /// </remarks>
        public double GetRandomValue (Random rng) {
            if (rng == null) throw new ArgumentNullException("rng");
            return (InverseLeftProbability(rng.NextDouble()));
        }

        // compute central moments from raw moments
        // this is subject to loss of precision from cancelation, so be careful

        internal virtual double CentralMomentFromRawMoment (int n) {

            double m = Mean;

            double mm = 1.0;
            double C = Moment(n);
            for (int k = 1; k <= n; k++) {
                mm = mm * (-m);
                C += AdvancedIntegerMath.BinomialCoefficient(n,k) * mm * Moment(n - k);
            }

            return (C);

        }

        // compute raw moments from central moments
        // this doesn't suffer from the cancelation problem of the reverse calculation

        internal virtual double RawMomentFromCentralMoments (int n) {

            double m = Mean;

            double mm = 1.0;
            double M = MomentAboutMean(n);
            int B = 1; // binomial coefficient; use recurrence B(n,k) = (n-k+1/k) B(n,k-1)
            for (int k = 1; k <= n; k++) {
                B = B * (n - k + 1) / k;
                mm = mm * m;
                M += B * mm * MomentAboutMean(n - k);
            }

            return (M);

        }

    }

    /// <summary>
    /// Represents the distribution of the Kolmogorov-Smirnov D statistic.
    /// </summary>
    /// <remarks><para>The D statistic in a Kolmogorov-Smirnov test is distributed (under the null hypothesis) according to a Kolmogorov disribution.</para></remarks>
    /// <seealse cref="Sample.KolmogorovSmirnovTest(Meta.Numerics.Statistics.Distribution)" />
    public class KolmogorovDistribution : Distribution {

        /// <summary>
        /// Instantiates a new asymptotic Kolmogorov distribution.
        /// </summary>
        public KolmogorovDistribution () {}

        // the sample size; when N=0 we will report the asymptotic distribution

        internal KolmogorovDistribution (double scale) {
            this.scale = scale;
        }

        private double scale = 1.0;

        /// <inheritdoc />
        public override double ProbabilityDensity (double d) {

            if (d < 1.0) {
                return (AsymptoticPPrime(d));
            } else {
                return (AsymptoticQPrime(d));
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
        public override double LeftProbability (double d) {

            if (d < scale) {
                return (AsymptoticP(d/scale));
            } else {
                return (1.0 - AsymptoticQ(d/scale));
            }

        }

        /// <inheritdoc />
        public override double RightProbability (double d) {

            if (d < scale) {
                return (1.0 - AsymptoticP(d/scale));
            } else {
                return (AsymptoticQ(d/scale));
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
                    //Console.WriteLine("{0} {1} {2} {3} {4}", k, z, Math.Exp(-z * z / 2.0), dp, p);
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
                sign = - sign;
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
                return (Math.Sqrt(Math.PI / 2.0) * Global.LogTwo * scale);
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
        public override double StandardDeviation {
            get {
                return (Math.Sqrt(Variance));
            }
        }

        /// <inheritdoc />
        public override double Median {
            get {
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

#if FUTURE

    public class FiniteKolmogorovDistribution : KolmogorovDistribution {

        public FiniteKolmogorovDistribution (int size) {
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

            if (N < maxN) {

                double t = d * N;
                if (2.0 * t < N) {
                    return (DurbinPPrime(N, t) * N);
                } else {
                    return (DurbinQPrime(N, t) * N);
                }

            } else {

                return (base.ProbabilityDensity(d * sqrtN) * sqrtN);

            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double d) {

            if (d <= 1.0 / (2.0 * N)) {
                return (0.0);
            } else if (d >= 1.0) {
                return (1.0);
            } else {

                if (N < maxN) {

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

                } else {
                    return (base.LeftProbability(d * sqrtN));
                }

            }

        }

        /// <inheritdoc />
        public override double RightProbability (double d) {

            if (d <= 1.0 / (2.0 * N)) {
                return (1.0);
            } else if (d >= 1.0) {
                return (0.0);
            } else {

                if (N < maxN) {

                    // use Durbin formulas for small N

                    double t = d * N;
                    if (2.0 * t < N) {
                        return (1.0 - MatrixP(N, t));
                        //return (1.0 - DurbinP(N, t));
                    } else {
                        return (DurbinQ(N, t));
                    }

                } else {

                    return (base.RightProbability(d * sqrtN));

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

    /// <summary>
    /// Represents the asymptotic distribution of Kuiper's V statistic.
    /// </summary>
    public class KuiperDistribution : Distribution {

        /// <summary>
        /// Instantiates a new Kuiper distribution.
        /// </summary>
        public KuiperDistribution () {
        }

        internal KuiperDistribution (double scale) {
            this.scale = scale;
        }

        private double scale = 1.0;

        /// <inheritdoc />
        public override double LeftProbability (double x) {
            if (x < 1.0) {
                return (AsymptoticP(x/scale));
            } else {
                return (1.0 - AsymptoticQ(x/scale));
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double x) {
            if (x < 1.0) {
                return (1.0 - AsymptoticP(x/scale));
            } else {
                return (AsymptoticQ(x/scale));
            }
        }

        // series \sqrt{2\pi}{x^3} \sum_{k=1}^{\infty} k^2 \pi^2 e^{-k^2 \pi^2 / 2 x^2}
        // useful for small x

        private static double AsymptoticP (double x) {

            if (x <= 0.0) return (0.0);

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;
                double z = k * Math.PI / x;
                double z2 = z * z;
                double ds = z2 * Math.Exp(-z2 / 2.0);
                s += ds;

                if (s == s_old) return (Global.SqrtTwoPI / x * s);

            }

            throw new NonconvergenceException();

        }

        // series \sum_{k=1}^{\infty} (4 k^2 x^2 - 1) e^{-2 k^2 x^2}
        // useful for large x

        private static double AsymptoticQ (double x) {

            double s = 0.0;
            for (int k = 1; k < Global.SeriesMax; k++) {

                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (4.0 * z2 - 1.0) * Math.Exp(-2.0 * z2);
                s += ds;

                if (s == s_old) return (2.0 * s);
            }

            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double x) {
            if (x < 1.0) {
                return (AsymptoticPPrime(x/scale)/scale);
            } else {
                return (AsymptoticQPrime(x/scale)/scale);
            }
        }

        private static double AsymptoticQPrime (double x) {

            double s = 0.0;
            for (int k = 1; k < 20; k++) {

                double s_old = s;
                double z = k * x;
                double z2 = z * z;
                double ds = (k * k) * (4.0 * z2 - 3.0) * Math.Exp(-2.0 * z2);
                s += ds;

                if (s == s_old) return (8.0 * x * s);
            }

            throw new NonconvergenceException();

        }

        private static double AsymptoticPPrime (double x) {

            if (x <= 0.0) return (0.0);

            double s = 0.0;
            for (int k = 1; k < 20; k++) {

                double s_old = s;
                double z = Math.PI * k / x;
                double z2 = z * z;
                double ds = z2 * (z2 - 3.0) * Math.Exp(-z2 / 2.0);
                s += ds;

                if (s == s_old) return (Global.SqrtTwoPI / (x * x) * s);

            }

            throw new NonconvergenceException();

        }

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        /// <inheritdoc />
        public override double Mean {
            get {
                // from numerical integral; would be nice to get an analytic result
                return (1.2533141373155 * scale);
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
                return (Math.PI * Math.PI / 6.0 * scale * scale);
            } else {
                return (AdvancedMath.RiemannZeta(n) * AdvancedMath.Gamma(1 + n / 2.0) * (n - 1) / Math.Pow(2.0, n / 2.0 - 1) * Math.Pow(scale, n));
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

    }



    /// <summary>
    /// Represents the distribution of the Mann-Whitney statistic.
    /// </summary>
    /// <remarks>
    /// <para>The expected use of this class is as the TestResult.Distribution of a
    /// Mann-Whitney U-test. You probably don't want to use this distribution directly
    /// unless you are studying the Mann-Whitney U-test.</para>
    /// </remarks>
    public class MannWhitneyDistribution : Distribution {

        /// <summary>
        /// Instantiates a new Mann-Whitney distribution.
        /// </summary>
        /// <param name="m">The number of elements in the first sample.</param>
        /// <param name="n">The number of elements in the second sample.</param>
        public MannWhitneyDistribution (int m, int n) {

            if (m < 1) throw new ArgumentOutOfRangeException("m");
            if (n < 1) throw new ArgumentOutOfRangeException("n");

            mc = m;
            nc = n;

            total = Math.Exp(AdvancedIntegerMath.LogFactorial(n + m) - AdvancedIntegerMath.LogFactorial(n) - AdvancedIntegerMath.LogFactorial(m));

            if (0.99 * total < ((double) Decimal.MaxValue)) {

                // for a total small enough to fit into a decimal, compute the exact count of partitions
                // contributing to each u; this operation is m^2 n, but it's integer math and we are
                // only doing it up to m~n~50, so it's still pretty fast

                counts = GaussianBinomialCoefficients(m + n, m);

                decimal sum = 0M;
                for (int u = 0; u <= m * n; u++) {
                    sum += counts[u];
                }
                total = (double) sum;

            } else {

                // for larger values, we will use a normal approximation

                normal = new NormalDistribution(Mean, StandardDeviation);

            }

            
            // i'd like to use an edgeworth approximation, and I have formulas for the required higher
            // cumulants, but edgeworth breaks down for extreme values (giving negative probabilities
            // and probabilities larger than one); i could transition back to a normal approximation
            // at those extremes, but i don't know how to do this and keep all my quantities self-consistent,
            // i.e. p is derivative of P, moment functions return actual moment integrals, etc.

            /*
            double M1 = m * n / 2.0;
            double C2 = 0.0;
            double C4 = 0.0;
            double C6 = 0.0;
            for (int i = 0; i <= m * n; i++) {
                double z = (i - M1);
                double z2 = z * z;
                double z4 = z2 * z2;
                double z6 = z4 * z2;
                C2 += ((double) counts[i]) * z2;
                C4 += ((double) counts[i]) * z4;
                C6 += ((double) counts[i]) * z6;
            }
            C2 = C2 / ((double) total);
            C4 = C4 / ((double) total);
            C6 = C6 / ((double) total);

            double m2 = m * m;
            double n2 = n * n;
            double m3 = m2 * m;
            double n3 = n2 * n;
            double m4 = m2 * m2;
            double n4 = n2 * n2;

            double S1 = m * n / 2.0;
            double S2 = m * n * (m + 1.0 + n) / 12.0;
            double S3 = 0.0;
            double S4 = - m * n * (m + n + 1.0) * (m2 + m + m * n + n + n2) / 120.0;            
            double S6 = m * n * (m + n + 1.0) * (
                2.0 * m4 + 4.0 * m3 * n + 4.0 * m3 + 7.0 * m2 * n + m2 - m +
                2.0 * n4 + 4.0 * m * n3 + 4.0 * n3 + 7.0 * m * n2 + n2 - n +
                6.0 * m2 * n2 + 2.0 * m * n ) / 504.0;

            Console.WriteLine("M1 {0} {1}", M1, S1);
            Console.WriteLine("C2 {0} {1}", C2, S2);
            Console.WriteLine("C4 {0} {1}", C4, S4 + 3.0 * S2 * S2);
            Console.WriteLine("C6 {0} {1}", C6, S6 + 15.0 * S4 * S2 + 10.0 * S3 * S3 + 15.0 * S2 * S2 * S2);
            */

        }

        private int mc, nc;
        private double total;
        private decimal[] counts;
        private NormalDistribution normal;

        public override double  ProbabilityDensity(double x) {
 	        throw new NotImplementedException();
        }

        private double Probability (int u) {
            if ((u < 0) || (u > mc*nc)) return (0.0);

            if (counts != null) {
                return ((double) counts[u] / total);
            } else {
                return (normal.ProbabilityDensity(u));
            }

        }

        /// <inheritdoc />
        public override double  LeftProbability(double x) {
            return( LeftInclusiveProbability((int) Math.Truncate(x)) );
        }   

        private double LeftInclusiveProbability (int u) {

            if (u < 0) return (0.0);
            if (u > mc*nc) return (1.0);

            if (counts != null) {
                double P = 0;
                for (int i = 0; i <= u; i++) {
                    P += (double) counts[i];
                }
                return (P / total);
            } else {
                return( normal.LeftProbability(u) );
            }


        }

        /// <inheritdoc />
        public override double RightProbability (double u) {
            return( RightExclusiveProbability((int) Math.Truncate(u)) );
        }

        private double RightExclusiveProbability (int u) {

            if (u < 0) return (1.0);
            if (u > mc*nc) return (0.0);

            if (counts != null) {
                double Q = 0;
                for (int i = u + 1; i <= mc * nc; i++) {
                    Q += (double) counts[i];
                }
                return (Q / total);
            } else {
                return( normal.RightProbability(u) );
            }
        }

        
        // this routine is based on the recurrsion
        // [ m n ] = ( 1 - q^m ) / ( 1 - q^(m-n) ) [ m-1 n ]
        // and the starting point [ n n ] = 1

        // the coefficients are integers and get large quickly as m and n increase
        // we use decimal because it handles larger integers than long
        // we can't use double because the calculation requires delicate cancelations
        // among large intermediate values, thus necessicating exact integer arithmetic
        // look into using an arbitrary-sized integer structure in the future

        private decimal[] GaussianBinomialCoefficients (int m, int n) {

            if (m < 0) throw new ArgumentOutOfRangeException("m");
            if (n < 0) throw new ArgumentOutOfRangeException("n");

            Debug.Assert(m >= n);

            // create  an array to hold our coefficients
            decimal[] c = new decimal[(m - n) * n + 1];

            // start with [n n] = 1 * q^0
            c[0] = 1;

            // keep track of current degree of our polynomial
            int d = 0;

            // create a scratch array for intermediate use
            // it needs to be larger than the previous array by (m-n) to hold intermediate polynomials
            decimal[] b = new decimal[c.Length + (m-n)];

            // interate from [n n] up to [m n]
            for (int k = n + 1; k <= m; k++) {

                // multiply by (1-q^k)
                for (int i = 0; i <= d; i++) {
                    b[i] = c[i];
                }
                d = d + k;
                for (int i = k; i <= d; i++) {
                    b[i] = b[i] - c[i - k];
                }

                // divide by (1-q^(k-n))
                for (int i = d - (k - n); i >= 0; i--) {
                    c[i] = -b[k - n + i];
                    b[k - n + i] = b[k - n + i] + c[i];
                    b[i] = b[i] - c[i];
                }
                d = d - (k - n);

            }

            // we're done
            return (c);


        }

        // these are some other ways we tried to generate counts; they are not as efficient as
        // the gaussian binomial coefficient technique, but I keep them around for the record

        /*
        public static int[] ComputeProbabilities (int n, int m) {

            int[] counts = new int[n * m + 1];

            int nm = n + m;

            // initialze to the lexographically first sequence
            // we could save memory by doing this witha BitArray, but we loose
            // about a factor of 2 in performance
            bool[] sequence = new bool[nm];
            for (int i = 0; i < n; i++) {
                sequence[i] = false;
            }
            for (int i = n; i < nm; i++) {
                sequence[i] = true;
            }
            int W = 0;

            while (true) {

                counts[W] = counts[W] + 1;

                // find the descending-order string at the end (the first 0<-1 transition)
                int j = nm - 1;
                int n0 = 0;
                int n1 = 0;
                while (j > 0) {

                    if (!sequence[j]) {
                        n0++;
                    } else {
                        n1++;
                        if (!sequence[j - 1]) break;
                    }

                    j--;

                }

                // check whether we're done
                if (j == 0) break;

                // flip the 0 to a 1
                sequence[j - 1] = true;

                // put the remainder of the string in ascending-order (0s, then 1s)
                for (int i = 0; i < (n0 + 1); i++) {
                    sequence[j] = false;
                    j++;
                }
                for (int i = 0; i < (n1 - 1); i++) {
                    sequence[j] = true;
                    j++;
                }

                // change W
                int dW = (n0 + 1) * (n1 - 1) - n1;
                W -= dW;

            }

            int sum = 0;
            for (int i = 0; i < counts.Length; i++) {
                sum += counts[i];
                Console.Write(counts[i] + " ");
            }
            Console.WriteLine();
            Console.WriteLine(sum);


            return(counts);

        }

        public static int CountOrderings (int m, int n, int u) {
            
            if ((m < 0) || (n < 0) || (u < 0) || (u > m*n)) return(0);

            if (u == 0) return (1);

            return (CountOrderings(m - 1, n, u - n) + CountOrderings(m, n - 1, u));

        }
        */

        /// <inheritdoc />
        public override double Mean {
            get {
                return (mc * nc / 2.0);
            }
        }

        /// <inheritdoc />
        public override double Median {
	        get { 
		        return (Mean);
	        }
        }

        /// <inheritdoc />
        public override double Variance {
            get {
                return (mc * nc * (mc + nc + 1) / 12.0);
            }
        }

        /// <inheritdoc />
        public override double Skewness {
            get {
                return (0.0);
            }
        }

        // add moments

    }


    /// <summary>
    /// Represents an parameterized likelihood distribution.
    /// </summary>
    public interface IParameterizedDistribution {

        /// <summary>
        /// Gets the parameter values of the distribution.
        /// </summary>
        /// <returns></returns>
        double[] GetParameters ();

        /// <summary>
        /// Sets the parameter values of the distribution.
        /// </summary>
        /// <param name="parameters">A list of parameter values.</param>
        void SetParameters (IList<double> parameters);

        /// <summary>
        /// Gets the likelihood of a value, given the current parameters.
        /// </summary>
        /// <param name="x">The value.</param>
        /// <returns>The likelihood of the value.</returns>
        double Likelihood (double x);

    }

    /*
    public class WaldDistribution : Distribution {


        public WaldDistribution (double mu, double lambda) {
            if (mu <= 0.0) throw new ArgumentOutOfRangeException("mu");
            if (lambda <= 0.0) throw new ArgumentOutOfRangeException("lambda");
            this.mu = mu;
            this.lambda = lambda;
        }

        double mu;
        double lambda;

        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.0, Double.PositiveInfinity));
            }
        }

        public override double ProbabilityDensity (double x) {
            return (Math.Sqrt(lambda / (2.0 * Math.PI * x * x * x)) * Math.Exp(-lambda * (x - mu) / (2.0 * x * mu * mu)));
        }

        public override double Mean {
            get {
                return (mu);
            }
        }


        public override double Moment (int n) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException("n");
            } else if (n == 0) {
                return (1.0);
            } else {
                // use recursion relation M_{n+1} = ((2n-1) (M_n / lambda) + M_{n-1}) mu^2
                double MM = mu;
                if (n == 1) return (MM);
                double mu2 = mu * mu;
                double M0 = m2 * (lambda + mu) / lambda;
                if (n == 2) return (M0);
                for (int k = 2; k < n; k++) {
                    double MP = ((2 * k - 1) * (M0 / lambda) + MM) * mu2;
                    MM = M0;
                    M0 = MP;
                }
                return (M0);
            }
        }

    }
    */
#if FUTURE
    internal static class DistributionMath {

        public static double ComputeCentralMomentFromRawMoments (double[] rawMoments, int n) {

            throw new NotImplementedException();
        }

        public static double ComputeRawMomentFromCentralMoments (double mean, double[] centralMoments, int n) {

            if (centralMoments == null) throw new ArgumentNullException("C");
            if (centralMoments.Length <= n) throw new InvalidOperationException();


            throw new NotImplementedException();
        }

        public static double ComputeRawMomentFromCulumants (double[] cumulants, int n) {

            if (cumulants == null) throw new ArgumentNullException("cumulants");
            if (n < 0) throw new ArgumentNullException("n");
            if (cumulants.Length <= n) throw new InvalidOperationException();

            double K = 0.0;
            IntegerPartitionEnumerator e = new IntegerPartitionEnumerator(n);
            while (e.MoveNext()) {
                double K1 = 1.0;
                int[] fs = e.Current;
                foreach (int f in fs) {
                    K1 *= cumulants[f];
                }
                K += K1;
            }

            return (K);

        }

        public static double ComputeCentralMomentFromCumulants (double[] cumulants, int n) {
            throw new NotImplementedException();
        }

        public static double ComputeCumulantFromRawMoments (double[] rawMoments, int n) {

            if (rawMoments == null) throw new ArgumentNullException("rawMoments");
            if (n < 0) throw new ArgumentOutOfRangeException("n");
            if (rawMoments.Length <= n) throw new InvalidOperationException();

            double K = rawMoments[n];
            for (int k = n - 1; k >= 0; k++) {
            }

            return (K);
        }

        public static double CumputeCumulantFromCentralMoments (int mean, double[] centralMoments, int n) {
            throw new NotImplementedException();
        }

    }
#endif

    // Deviates
    // Maximum likelyhood estimation
    // Cumulants

	// Gamma Distribution

}

