using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.ObjectModel;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.SignalProcessing;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace FutureTest {

    [TestClass]
    public class FutureTest {


        [TestMethod]
        public void KelvinTest () {
            Console.WriteLine(AdvancedMath.KelvinBer(0.0, 1.0));
        }

        [TestMethod]
        public void TestFinitePQ () {

            int n = 32;
            double t = 16.25;

            Console.WriteLine(FiniteQ(n, t));
            Console.WriteLine(1.0 - MatrixP(n, t));

        }

        private double FiniteQ (int n, double t) {

            if (t < n / 2.0) throw new ArgumentOutOfRangeException();

            double Q = 2.0 * MoreMath.Pow((n - t) / n, n);

            if (t > n - 1) return (Q);

            Q += 2.0 * t * MoreMath.Pow((n - 1 - t) / n, n - 1);

            if (t > n - 2) return (Q);

            for (int j = 2; n - j - t > 0.0; j++) {
                Q += 2.0 * AdvancedIntegerMath.BinomialCoefficient(n, j) * (t / n) * MoreMath.Pow((t + j) / n, j - 1) * MoreMath.Pow((n - j - t) / n, n - j);
            }

            return (Q);

        }

        private double K0Q (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {
                double p = 1.0;
                for (int k = 1; k < 50; k++) {
                    double p_old = p;
                    double z = k * x;
                    double dp = 2.0 * Math.Exp(-2.0 * z * z);
                    if (k % 2 != 0) dp = -dp;
                    p += dp;
                    if (p == p_old) { Console.WriteLine("KQ {0}", k); return (p); }
                }
                throw new NonconvergenceException();
            }
        }

        private double K0 (double x) {
            if (x <= 0.0) {
                return (0.0);
            } else {

                double p = 0.0;
                for (int k = 1; k < 100; k += 2) {
                    double p_old = p;
                    double z = k * Math.PI / x / 2.0;
                    double dp = Math.Exp(-z * z / 2.0);
                    p += dp;
                    if (p == p_old) { Console.WriteLine("KP {0}", k); return (Math.Sqrt(2.0 * Math.PI) / x * p); }
                }

                throw new NonconvergenceException();
            }
        }

        private double K1Q (double x) {
            double x2 = x * x;
            double p = 0.0;
            for (int k = 1; k < 100; k++) {
                double p_old = p;
                int k2 = k * k;
                double dp = k2 * Math.Exp(-2.0 * k2 * x2);
                if (k % 2 == 0) dp = -dp;
                p += dp;
                if (p == p_old) return (4.0 / 3.0 * x * p);
            }
            throw new NonconvergenceException();
        }

        private double K1 (double x) {
            double p = 0.0;
            for (int j = 0; j < 100; j++) {
                double p_old = p;
                int k = 2 * j + 1;
                double z = k * Math.PI / x / 2.0;
                double dp = Math.Exp(-z * z / 2.0) * (z * z - 1.0);
                p += dp;
                if (p == p_old) return (Math.Sqrt(2.0 * Math.PI) / 6.0 / (x * x) * p);
            }
            throw new NonconvergenceException();
        }

        public double K2 (double x) {
            double p = 0.0;
            for (int j = 0; j < 100; j++) {
                double p_old = p;
                int k = 2 * j + 1;
                double z = k * Math.PI / x / 2.0;
                double dp = Math.Exp(-z * z / 2.0) * (z * z * z * z * (1.0 - 2.0 * x * x) + z * z * (2.0 * x * x - 5.0) + (6.0 * x * x + 2.0));
                p += dp;
                if (p == p_old) return(Math.Sqrt(2.0 * Math.PI) * 72.0 / (x * x * x) * p );
            }
            throw new NonconvergenceException();
        }



        private double KP1Prime (double x) {
            double f = 0.0;
            for (int k = 1; k < 100; k += 2) {
                double f_old = f;
                double z = (k * k) * (Math.PI * Math.PI) / (x * x) / 8.0;
                double df = (4.0 * z * z - 10 * z + 2.0) * Math.Exp(-z);
                f += df;
                if (f == f_old) return (Math.Sqrt(2.0 * Math.PI) / 6.0 / (x * x * x) * f);
            }
            throw new NotImplementedException();
        }

        private double KQ1Prime (double x) {
            double a = 2.0 * x * x;
            double f = 0.0;
            for (int k = 1; k < 50; k++) {
                double f_old = f;
                int k2 = k * k;
                double z = a * k2;
                double df = k2 * (1.0 - 2.0 * z) * Math.Exp(-z);
                if (k % 2 != 0) df = -df;
                f += df;
                if (f == f_old) return (4.0 / 3.0 * f);
            }
            throw new NonconvergenceException();
        }

        [TestMethod]
        public void ComparePrimes () {
            for (double x = 0.4; x < 2.0; x += 0.2) {
                Console.WriteLine("{0} {1} {2}", x, KP1Prime(x), KQ1Prime(x));
            }
        }

        [TestMethod]
        public void TestKSSeries2 () {

            for (double x = 0.2; x <= 2.0; x += 0.2) {

                Console.WriteLine("x={0} {1} {2}", x, K1(x), K1Q(x));

            }

        }

        [TestMethod]
        public void TestFKDClass () {

            int n = 10;

            int m = 10;

            AsymptoticFiniteKolmogorovDistribution D = new AsymptoticFiniteKolmogorovDistribution(n);
            Interval i = Interval.FromEndpoints(0.0, Double.PositiveInfinity);

            //Console.WriteLine(D.ProbabilityDensity(1.1));

            double Q = FunctionMath.Integrate(x => D.ProbabilityDensity(x) * MoreMath.Pow(x, m), i);
            Console.WriteLine("{0} {1}", Q, D.Moment(m));


        }

        [TestMethod]
        public void TestFKDClass2 () {
            int n = 100;
            AsymptoticFiniteKolmogorovDistribution D = new AsymptoticFiniteKolmogorovDistribution(n);
            Console.WriteLine(D.LeftProbability(0.0874483967333 * Math.Sqrt(n)));
        }

        [TestMethod]
        public void TestKSSeries () {

            int n = 512;

            KolmogorovDistribution KD = new KolmogorovDistribution();
            AsymptoticFiniteKolmogorovDistribution FKD = new AsymptoticFiniteKolmogorovDistribution(n);

            for (double x = 0.2; x < 2.0; x += 0.2) {

                double t = Math.Sqrt(n) * x;

                double P0 = K0(x);
                double P1 = K0(x) + K1(x) / Math.Sqrt(n);
                double P2 = K0(x) + K1(x) / Math.Sqrt(n) + K2(x) / n;
                double PM = MatrixP(n, t);
                Console.WriteLine("x={0} PK={1} P0={2} P1={3} P1'={4} PM={5}", x, KD.LeftProbability(x), P0, P1, FKD.LeftProbability(x), PM);

            }

        }

        private static double MatrixP (int n, double t) {

            // compute stuff used in matrix entries
            int tp = (int) Math.Truncate(t) + 1;
            double h = tp - t;
            int p = 2 * tp - 1;


            // construct the matrix
            SquareMatrix H = new SquareMatrix(p);

            // superdiagonal
            for (int j = 1; j < p; j++) {
                H[j - 1, j] = 1.0;
            }

            // diagonal and subdiagonals
            double F = 1.0; // factorial
            double hh = h; // power of h
            for (int i = 1; i < p; i++) {
                H[i - 1, 0] = (1.0 - hh) / F;
                H[p - 1, p - i] = H[i - 1, 0];
                for (int j = i + 1; j < p; j++) {
                    H[j - 1, j - i] = 1.0 / F;
                }
                hh = hh * h;
                F = F * (i + 1);
            }

            // lower-left element
            double g = 1.0 - 2.0 * hh;
            if (h > 0.5) g = g + Math.Pow(2.0 * h - 1.0, p);
            g = g / F;
            H[p - 1, 0] = g;

            // raise the matrix to the nth power
            SquareMatrix HN = MatrixPower(H, n);

            // return the appropriate element
            double hf = Math.Exp(AdvancedIntegerMath.LogFactorial(n) - n * Math.Log(n));
            return (hf * HN[tp - 1, tp - 1]);


        }

        private static SquareMatrix MatrixPower (SquareMatrix A, int n) {

            SquareMatrix B = null;

            SquareMatrix D = A.Copy();

            while (true) {
                if (n % 2 != 0) {
                    if (B == null) {
                        B = D.Copy();
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

        [TestMethod]
        public void FftTimer2 () {

            int n = 65536;

            FourierTransformer t = new FourierTransformer(n);

            Complex[] x = new Complex[n];
            for (int i = 0; i < n; i++) {
                x[i] = new Complex(2.0, -2.0);
            }

            Stopwatch s = Stopwatch.StartNew();
            Complex[] y = t.Transform(x);
            s.Stop();

            Console.WriteLine(s.ElapsedMilliseconds);

        }

        [TestMethod]
        public void FftTimer () {

            int nmax = 4000;

            FourierTransformer[] t = new FourierTransformer[nmax];
            for (int n = 2; n < nmax; n++) {
                t[n] = new FourierTransformer(n);
            }

            Stopwatch s = Stopwatch.StartNew();

            for (int n = 2; n < nmax; n++) {
                if (n % 4 != 0) continue;
                Complex[] y = new Complex[n];
                for (int i = 0; i < n; i++) y[i] = i;
                Complex[] yt = t[n].Transform(y);

            }

            s.Stop();
            Console.WriteLine(s.ElapsedMilliseconds);
        }

        private void FermatFactor (int n) {

            int c = 0;
            int s = (int)Math.Floor(Math.Sqrt(n));
            int x = 2 * s + 1; int y = 1; int r = s * s - n;
            while (r != 0) {
                r += x; x += 2;
                c++;
                do {
                    r -= y; y += 2;
                    c++;
                } while (r > 0);

            }
            int p = (x - y) / 2;
            int q = (x + y - 2) / 2;
            if (p == 1) {
                //Console.WriteLine(q);
            } else {
                FermatFactor(p);
                FermatFactor(q);
            }
        }

        private void PollardRhoFactor (int n) {

            int x = 5; int y = 2; int k = 1; int l = 1;

            for (int c = 0; c < 10000; c++) {
            //while (true) {
                int g = (int) AdvancedIntegerMath.GCF(Math.Abs(y - x), n);
                if (g == n) {
                    Console.WriteLine("n={0}", n);
                    return;
                } else if (g == 1) {
                    k--;
                    if (k == 0) {
                        y = x;
                        l = 2 * l;
                        k = l;
                    }
                    //Console.WriteLine("{0}^2 mod {1}", x, n);
                    x = AdvancedIntegerMath.PowMod(x, 2, n) + 1;
                    if (x == n) x = 0;
                } else {
                    //Console.WriteLine("g={0}", g);
                    n = n / g;
                    x = x % n;
                    y = y % n;
                }
            }



        }

#if FUTURE
        [TestMethod]
        public void TestFactor () {
            /*
            Stopwatch s1 = Stopwatch.StartNew();
            FermatFactor(1157625);
            s1.Stop(); Console.WriteLine(s1.ElapsedMilliseconds);
            
            Stopwatch s2 = Stopwatch.StartNew();
            PollardRhoFactor(37);
            s2.Stop(); Console.WriteLine(s2.ElapsedMilliseconds);
            */
            // for (3*5*7)^3 = 1157625, Pollard's Rho fails to factor 125 = 5^3
            
            int c = 0;
            Stopwatch s1 = Stopwatch.StartNew();

            //for (int i = 1; i < 1000000; i+=2) {
            int i = 220211;
                List<Factor> factors = AdvancedIntegerMath.Factor(i);
                //FermatFactor(i);
                
                int m = 1;
                foreach (Factor factor in factors) {
                    Console.WriteLine(factor);
                    if (!AdvancedIntegerMath.IsPrime(factor.Value)) {
                        c++;
                        Console.WriteLine("for {0}, factor {1} is not prime", i, factor.Value);
                    }
                    //Console.WriteLine("{0} {1}", factor.Value, factor.Multiplicity);
                    m *= (int)MoreMath.Pow(factor.Value, factor.Multiplicity);
                }
                if (m != i) {
                    Console.WriteLine("for {0}, factors do not multiply to number", i);
                }
                //Assert.IsTrue(m == i);
                //Console.WriteLine(m);
                
            // }
            s1.Stop(); Console.WriteLine(s1.ElapsedMilliseconds);

            Console.WriteLine(c);
            
        }
#endif
#if FUTURE
        [TestMethod]
        public void BigAdd () {

            BigFloat a = new BigFloat(new byte[] { 2, 2, 5 }, 0);
            Console.WriteLine(a);
            BigFloat b = new BigFloat(new byte[] { 3, 3, 3 }, 0);
            Console.WriteLine(b);
            BigFloat c = a + b;
            Console.WriteLine(c);

        }
#endif
        [TestMethod]
        public void EinTest () {

            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 3.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 2.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 2.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 1.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 1.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 0.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, 0.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -0.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -1.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -1.5)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -2.0)));
            Console.WriteLine(AdvancedComplexMath.Ein(new Complex(-50.0, -2.5)));

        }

        [TestMethod]
        public void ErfTest () {
            Console.WriteLine(AdvancedComplexMath.Erf(new Complex(1.0E-10, 10.0)));

        }

        [TestMethod]
        public void MultiRootTest () {

            // here is a system with a zeros at (3,4) and (1,0)
            Func<double[], double[]> f = delegate (double[] u) {
                double x = u[0]; double y = u[1];
                double a = 2.0 * y - x;
                double b = (x * x + x * (y * y - 2.0) - 4.0 * y) / (x + 4.0);
                double c = Math.Sqrt(x * x + y * y);
                return (new double[] { b - a, c - a });
            };

            FunctionMath.FindZero(f, new double[] { 1.0, 1.0 });

        }

        [TestMethod]
        public void KruskalWallis () {

            // start with some data lists
            double[][] data = new double[5][];
            data[0] = new double[] { 13, 4, 12, 1 };
            data[1] = new double[] { 19, 7, 9, 17 };
            data[2] = new double[] { 8, 20, 18, 5 };
            data[3] = new double[] { 15, 14, 16, 2 };
            data[4] = new double[] { 6, 10, 3, 11 };

            // sort each sample individually and compute count total from all samples
            int N = 0;
            for (int i = 0; i < data.Length; i++) {
                N += data[i].Length;
                Array.Sort<double>(data[i]);
            }
            Console.WriteLine("N={0}", N);

            // do a multi-merge sort to determine ranks sums

            // initialize a pointer to the current active index in each ordered list
            int[] p = new int[data.Length];

            // keep track the current rank to be assigned
            int r = 0;

            // keep track of rank sums
            // this is all that we need for KW
            // the ranks of individual entries can be made to disappear from the final formula using sum identities
            int[] rs = new int[data.Length];

            while (true) {
                
                // increment the rank
                // (programmers may think ranks start from 0, but in the definition of the KW test they start from 1)
                r++;

                // determine the smallest of current entries
                int j = -1;
                double f = Double.PositiveInfinity;
                for (int k = 0; k < data.Length; k++) {
                    if ((p[k] < data[k].Length) && (data[k][p[k]] < f)) {
                        j = k;
                        f = data[k][p[k]];
                    }
                }

                // test for all lists complete
                if (j < 0) break;

                // increment the pointer and the rank sum for that column
                p[j]++;
                rs[j] += r;

            }

            double H0 = 0.0;
            for (int i = 0; i < data.Length; i++) {
                double z = ((double) rs[i]) / data[i].Length - (N + 1) / 2.0;
                Console.WriteLine("{0} {1}", rs[i], z);
                H0 += data[i].Length * (z * z);
            }
            H0 = 12.0 / N / (N + 1) * H0;
            Console.WriteLine("H={0}", H0);

            Sample s1 = new Sample(data[0]);
            Sample s2 = new Sample(data[1]);
            Sample s3 = new Sample(data[2]);
            Sample s4 = new Sample(data[3]);
            Sample s5 = new Sample(data[4]);

            // total ranks
            Console.WriteLine("{0} {1} {2} {3} {4}", s1.Mean * s1.Count, s2.Mean * s2.Count, s3.Mean * s3.Count, s4.Mean * s4.Count, s5.Mean * s5.Count);

            // statistic
            N = 20;
            double H = 12.0 / N / (N + 1) * (s1.Count * MoreMath.Pow(s1.Mean, 2) + s2.Count * MoreMath.Pow(s2.Mean, 2) +
                s3.Count * MoreMath.Pow(s3.Mean, 2) + s4.Count * MoreMath.Pow(s4.Mean, 2) + s5.Count * MoreMath.Pow(s5.Mean, 2))
                - 3.0 * (N + 1);
            Console.WriteLine("H={0}", H);

            Distribution DH = new ChiSquaredDistribution(4);
            //Console.WriteLine(DH.InverseRightProbability(0.05));
            Console.WriteLine(DH.RightProbability(H));

            Console.WriteLine("--");
            TestResult ar = Sample.OneWayAnovaTest(s1, s2, s3, s4, s5).Result;
            Console.WriteLine(ar.Statistic);
            Console.WriteLine(ar.RightProbability);

            Console.WriteLine("---");
            TestResult kw = Sample.KruskalWallisTest(s1, s2, s3, s4, s5);
            Console.WriteLine(kw.Statistic);
            Console.WriteLine(kw.RightProbability);

        }


        [TestMethod]
        public void RealFourier () {
            double[] x = new double[6];
            for (int i = 0; i < x.Length; i++) {
                x[i] = i;
            }

            Complex[] xc = new Complex[x.Length];
            for (int i = 0; i < x.Length; i++) {
                xc[i] = x[i];
            }
            FourierTransformer xft = new FourierTransformer(x.Length);
            Complex[] xct = xft.Transform(xc);
            for (int i = 0; i < x.Length; i++) {
                Console.WriteLine("{0} {1}", i, xct[i]);
            }

            Assert.IsTrue(x.Length % 2 == 0);
            Complex[] z = new Complex[x.Length / 2];
            for (int i = 0; i < z.Length; i++) {
                z[i] = new Complex(x[2 * i], x[2 * i + 1]);
                Console.WriteLine("  {0} {1}", i, z[i]);
            }
            FourierTransformer zft = new FourierTransformer(z.Length);
            Complex[] zt = zft.Transform(z);
            for (int i = 0; i < z.Length; i++) {
                Console.WriteLine("  {0} {1}", i, zt[i]);
            }
            for (int i = 1; i < z.Length; i++) {
                Complex p = (zt[i] + zt[z.Length-i].Conjugate) / 2.0;
                Complex q = -ComplexMath.I * (zt[i] - zt[z.Length-i].Conjugate) / 2.0;
                double t = -2.0 * Math.PI * i / x.Length;
                Console.WriteLine("{0} {1}", i, p + q * new Complex(Math.Cos(t), Math.Sin(t)));
            }

        }

        public static double ApproximateInverseBetaSeries (double a, double b, double P) {

            double bigB = AdvancedMath.Beta(a, b);

            double z = Math.Pow(a * P * bigB, 1.0 / a);

            double z2 = (b - 1.0) / (a + 1.0) * z;

            double z3 = z2 * (a * a + 3.0 * b * a - a + 5.0 * b - 4.0) / (a + 1.0) / (a + 2.0) / 2.0 * z;

            double z4 = z2 * (a * a * a * a + (6.0 * b - 1.0) * a * a * a + (b + 2.0) * (8.0 * b - 5.0) * a * a +
                (33.0 * b * b - 30.0 * b + 4.0) * a + b * (31.0 * b - 47.0) + 18.0) / MoreMath.Pow(a + 1.0, 2) / (a + 2.0) / (a + 3.0) / 3.0 * z * z;

            Console.WriteLine("z={0} z2={1} z3={2} z4={3}", z, z2, z3, z4);

            return (z * (1.0 + z2 + z3 + z4));

        }

        public static double ApproximateInverseBeta (double a, double b, double P) {

            if (P > 0.5) {
                return (1.0 - ApproximateInverseBeta(b, a, 1.0 - P));
            } else {
                double bigB = AdvancedMath.Beta(a, b);
                double z = Math.Pow(a * P * bigB, 1.0 / a);
                double z2 = (b - 1.0) / (a + 1.0) * z;
                if (z2 < 0.25) {
                    double z3 = z2 * (a * a + 3.0 * b * a - a + 5.0 * b - 4.0) / (a + 1.0) / (a + 2.0) / 2.0 * z;
                    return (z * (1.0 + z2 + z3));
                } else {
                    throw new NotImplementedException();
                }
            }

        }

        public static double RefineInverseBeta (double a, double b, double P, double x) {

            for (int i = 0; i < 8; i++) {
                double x_old = x;
                double y = AdvancedMath.LeftRegularizedBeta(a, b, x) - P;
                double yp = Math.Pow(x, a - 1.0) * Math.Pow(1.0 - x, b - 1.0) / AdvancedMath.Beta(a, b);
                double dx = -y / yp;
                x += dx;
                if (x == x_old) return (x);
            }
            return (x);
        }

        [TestMethod]
        public void TestBeta () {

            double a = 200.0; double b = 200.0; double P = 1.0E-5;
            double x1 = ApproximateInverseBetaSeries(a, b, P);
            if ((0.0 < x1) && (x1 < 1.0)) {
                Console.WriteLine("x1 {0} {1}", x1, AdvancedMath.LeftRegularizedBeta(a, b, x1));
            }

            double x2 = 1.0 - ApproximateInverseBetaSeries(b, a, 1.0 - P);
            if ((0.0 < x2) && (x2 < 1.0)) {
                Console.WriteLine("x2 {0} {1}", x2, AdvancedMath.LeftRegularizedBeta(a, b, x2));
            }

            //x1 = RefineInverseBeta(a, b, P, x1);
            //Console.WriteLine("{0} {1}", x1, AdvancedMath.LeftRegularizedBeta(a, b, x1));

            NormalDistribution N = new NormalDistribution();
            double m = a / (a + b); double s = Math.Sqrt(a * b / (a + b + 1.0)) / (a + b);
            double x3 = m + s * N.InverseLeftProbability(P);
            if ((0.0 < x3) && (x3 < 1.0)) {
                Console.WriteLine("x3 {0} {1}", x3, AdvancedMath.LeftRegularizedBeta(a, b, x3));
            }

            //Console.WriteLine(AdvancedMath.Beta(a, b, 0.35) / AdvancedMath.Beta(a, b));
            //Console.WriteLine(AdvancedMath.Beta(a, b, 0.40) / AdvancedMath.Beta(a, b));
            //Console.WriteLine(AdvancedMath.Beta(a, b, 0.45) / AdvancedMath.Beta(a, b));

        }

    }

}
