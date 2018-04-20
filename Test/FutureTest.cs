using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.ObjectModel;
using System.Text;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Analysis;
using Meta.Numerics.Matrices;
using Meta.Numerics.SignalProcessing;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {



    [TestClass]
    public class FutureTest {

        //[TestMethod]
        public void FranciaShapiro () {

            Random rng = new Random(1);
            NormalDistribution d = new NormalDistribution();
            Sample w = new Sample();

            int n = 2048;
            for (int i = 0; i < 1000000; i++) {

                Sample s = new Sample();
                for (int j = 0; j < n; j++) {
                    s.Add(d.GetRandomValue(rng));
                }
                TestResult r = s.ShapiroFranciaTest();
                w.Add(r.Statistic.Value);
            }

        }

        //[TestMethod]
        public void StandardDeviation () {

            Random rng = new Random(21212);
            ContinuousDistribution d = new ExponentialDistribution();
            //ContinuousDistribution d = new NormalDistribution();
            //ContinuousDistribution d = new ChiSquaredDistribution(2.0);
            Sample a = new Sample();
            Sample b = new Sample();

            int n = 100;
            for (int i = 0; i < 10000; i++) {
                Sample s = new Sample();
                for (int j = 0; j < n; j++) {
                    s.Add(d.GetRandomValue(rng));
                }
                a.Add(s.PopulationStandardDeviation.Value);
                double c2 = s.CentralMoment(2);
                double c4 = s.CentralMoment(4);
                double u = Math.Sqrt(c2 * n / (n - 1)) * (1.0 + (c4 - c2 * c2) / (8.0 * c2 * c2 * n));
                b.Add(u);
            }

        }

        [TestMethod]
        public void SkewnessVariance () {

            int n = 100;
            Random rng = new Random(11111);
            //ContinuousDistribution d = new BetaDistribution(0.75, 0.25);
            ContinuousDistribution d = new NormalDistribution();
            double C2 = d.CentralMoment(2);
            double C3 = d.CentralMoment(3);
            double C4 = d.CentralMoment(4);
            double C5 = d.CentralMoment(5);
            double C6 = d.CentralMoment(6);
            double V2 = (C4 - C2 * C2) / n;
            double V3 = (C6 - C3 * C3 + 9.0 * C2 * C2 * C2 - 6.0 * C4 * C2) / n;
            double V23 = (C5 - 4.0 * C3 * C2) / n;
            double VG = (V3 + 9.0 / 4.0 * MoreMath.Sqr(C3 / C2) * V2 - 3.0 * (C3 / C2) * V23) / (C2 * C2 * C2);

            Sample us = new Sample();
            Sample k2s = new Sample();
            Sample k3s = new Sample();
            Sample g0s = new Sample();
            Sample g1s = new Sample();

            for (int i = 0; i < 10000; i++) {
                Sample s = new Sample();
                for (int j = 0; j < n; j++) {
                    s.Add(d.GetRandomValue(rng));
                }
                double u = s.Skewness;
                us.Add(u);

                double c2 = s.CentralMoment(2);
                double c3 = s.CentralMoment(3);
                double c4 = s.CentralMoment(4);
                double c5 = s.CentralMoment(5);
                double c6 = s.CentralMoment(6);
                double v2 = (c4 - c2 * c2) / (n - 1);
                double v3 = (c6 - c3 * c3 + 9.0 * c2 * c2 * c2 - 6.0 * c4 * c2) / (n - 1);
                double v32 = (c5 - 4.0 * c3 * c2) / (n - 1);

                double k2 = c2 * n / (n - 1);
                double k3 = c3 * n * n / (n - 1) / (n - 2);
                double g0 = k3 / Math.Pow(k2, 1.5);
                double r = c3 / c2;
                double g1 = g0 * (1.0 - 15.0 / 8.0 * v2 / (c2 * c2) + 3.0 / 2.0 * v32 / (c3 * c2));
                double vg = g0 * g0 * (v3 / (c3 * c3) - 3.0 * v32 / (c3 * c2) + 9.0 / 4.0 * v2 / (c2 * c2));

                k2s.Add(k2);
                k3s.Add(k3);
                g0s.Add(g0);
                g1s.Add(g1);

                //double v = (u * u / n) * (v3 / (c3 * c3) + 9.0 / 4.0 * v2 / (c2 * c2) - 3.0 * v32 / (c3 * c2));
            }

        }

        private int PoissonDeviate (double lambda, Random rng) {
            double t0 = Math.Exp(-lambda);
            int k = 0;
            double t = rng.NextDouble();
            while (t > t0) {
                k++;
                t *= rng.NextDouble();
            }
            return (k);
        }

        private int PoissonDeviateNR (double lambda, Random rng) {

            double sqlam = Math.Sqrt(lambda);
            double loglam = Math.Log(lambda);

            while (true) {
                double u = 0.64 * rng.NextDouble();
                double v = -0.68 + 1.28 * rng.NextDouble();
                int k = (int) Math.Floor(sqlam * (v / u) + lambda + 0.5);
                if (k < 0) continue;
                double u2 = u * u;
                double lfac = AdvancedIntegerMath.LogFactorial(k);
                double p = sqlam * Math.Exp(-lambda + k * loglam - lfac);
                if (u2 < p) return (k);
            }

        }

        private int PoissonDeviatePRUA (double lambda, Random rng) {
            double c = 1.715527770; // 2.0 * Math.Sqrt(2.0 / Math.E); // 1.7155
            double d = 0.898916162; //  3.0 - 2.0 * Math.Sqrt(3.0 / Math.E); // 0.8989
            double a = lambda + 0.5;
            double r = Math.Sqrt(a);
            double h = c * r + d;
            double logLambda = Math.Log(lambda);

            while (true) {
                double u = rng.NextDouble();
                double v = rng.NextDouble();
                double x = a + h * (v - 0.5) / u;
                if (x < 0.0) continue;
                // test x >= b
                int k = (int) Math.Floor(x);
                // fast accept
                // fast reject
                double t = k * logLambda - AdvancedIntegerMath.LogFactorial(k) - lambda;
                if (u * u < Math.Exp(t)) return (k);
            }

        }

        [TestMethod]
        public void ComparePoissonGenerators () {

            int n = 10000000;
            double lambda = 15.0;
            Random rng = new Random(1);

            PoissonDistribution p = new PoissonDistribution(lambda);

            int s0 = 0;
            Stopwatch t0 = Stopwatch.StartNew();
            for (int i = 0; i < n; i++) {
                s0 += p.GetRandomValue(rng);
            }
            t0.Stop();

            int s1 = 0;
            Stopwatch t1 = Stopwatch.StartNew();
            for (int i = 0; i < n; i++) {
                s1 += PoissonDeviateNR(lambda, rng);
            }
            t1.Stop();
        }

        [TestMethod]
        public void AiryZeros () {

            for (int k = 1; k < 4; k++) {

                double t = 3.0 / 8.0 * Math.PI * (4 * k - 3);
                double x0 = -AiryZeroT(t);
                double x1 = FunctionMath.FindZero(AdvancedMath.AiryBi, x0);
                Debug.WriteLine($"{k} {x0} {x1}");

            }

        }

        public static double AiryAiZero (int k) {

            if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

            //double t = 3.0 / 8.0 * Math.PI * (4 * k - 1);
            //double x0 = -AiryZeroT(t);
            //double x1 = FunctionMath.FindZero(AdvancedMath.AiryAi, x0);

            double a = AiryZeroT( 3.0 / 8.0 * Math.PI * (4 * k + 0) );
            double b = AiryZeroT( 3.0 / 8.0 * Math.PI * (4 * k - 2) );
            Interval ab = Interval.FromEndpoints(-a, -b);
            double x = FunctionMath.FindZero(AdvancedMath.AiryAi, ab);

            return (x);

        }

        public static double AiryBiZero (int k) {

            if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

            //double t = 3.0 / 8.0 * Math.PI * (4 * k - 3);
            //double x0 = -AiryZeroT(t);
            //double x1 = FunctionMath.FindZero(AdvancedMath.AiryBi, x0);

            // Problem for k = 0, makes b NaN
            double a = AiryZeroT( 3.0 / 8.0 * Math.PI * (4 * k - 2) );
            double b = AiryZeroT( 3.0 / 8.0 * Math.PI * (4 * k - 4) );
            Interval ab = Interval.FromEndpoints(-a, -b);
            double x = FunctionMath.FindZero(AdvancedMath.AiryBi, ab);

            return (x);

        }

        private static double AiryZeroT (double t) {

            double t2 = 1.0 / (t * t);
            double s = 1.0 + t2 * (5.0 / 48.0 + t2 *(-5.0 / 36.0 + t2 * 77125.0 / 82944.0));
            return (Math.Pow(t, 2.0 / 3.0) * s);

        }

        [TestMethod]
        public void NonlinearRegressionLinearRegressionAgreementTest () {

            BivariateSample s = new BivariateSample();
            s.Add(1.0, 2.0);
            s.Add(3.0, 2.0);
            s.Add(3.0, 4.0);
            s.Add(5.0, 4.0);

            LinearRegressionResult linearFit = s.LinearRegression();
            NonlinearRegressionResult nonlinearFit = s.NonlinearRegression(
                (IReadOnlyList<double> c, double x) => c[0] + c[1] * x,
                new double[] { 1.0, 1.0 } 
            );

            Assert.IsTrue(TestUtilities.IsNearlyEqual(linearFit.Parameters.ValuesVector, nonlinearFit.Parameters.ValuesVector, Math.Sqrt(TestUtilities.TargetPrecision)));

        }

        public void JacobiPeriod () {

            double k = 0.875;

            int n = 0;
            double a = 1.0;
            double b = Math.Sqrt(1.0 - k * k);
            double c2 = a * a - b * b;
            double s = c2;
            while (a != b) {

                double ap = (a + b) / 2.0;
                double bp = Math.Sqrt(a * b);
                double c2p = MoreMath.Sqr(c2 / 4.0 / ap);

                n++;
                a = ap;
                b = bp;
                c2 = c2p;
                s += MoreMath.Pow(2.0, n) * c2;

            }

            // This works, but s~1 for k~1, so don't use to compute E for k~1.

        }

        [TestMethod]
        public void TwoWayAnova () {

            Sample[,] samples = new Sample[,] {
                { new Sample(54, 49, 59, 39, 55), new Sample(25, 29, 47, 26, 28) },
                { new Sample(53, 72, 43, 56, 52), new Sample(46, 51, 33, 47, 41) },
                { new Sample(33, 30, 26, 25, 29), new Sample(18, 21, 34, 40, 24) }
            };

            /*
            Sample[,] samples = new Sample[,] {
                { new Sample(4, 5, 6, 5), new Sample(7, 9, 8, 12), new Sample(10, 12, 11, 7) },
                { new Sample(6, 6, 4, 4), new Sample(13, 15, 12, 12), new Sample(12, 13, 10, 13) }
            };
            */

            Sample.TwoWayAnovaTest(samples);

        }


        [TestMethod]
        public void BinomialInverseCdf () {

            foreach (int n in new int[] { 2, 8, 32, 128 }) {
                foreach (double p in new double[] { 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99 }) {
                    BinomialDistribution d = new BinomialDistribution(p, n);
                    foreach (double P in new double[] { 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99 }) {

                        double Q = 1.0 - P;

                        double kMin = 0.0;
                        string sMin = "Minimum";
                        double kMax = n;
                        string sMax = "Maximum";

                        double mu = n * p;

                        double kMarkov = mu / Q;
                        if (kMarkov < kMax) {
                            kMax = kMarkov;
                            sMax = "Markov";
                        }

                        double kMedian;
                        if (P < 0.5) {
                            kMedian = Math.Ceiling(mu);
                            if (kMedian < kMax) {
                                kMax = kMedian;
                                sMax = "Median";
                            }
                        } else {
                            kMedian = Math.Floor(mu);
                            if (kMedian > kMin) {
                                kMin = kMedian;
                                sMin = "Median";
                            }
                        }

                        double kChernovLower = mu - Math.Sqrt(-2.0 * mu * Math.Log(P));
                        if (kChernovLower > kMin) {
                            kMin = kChernovLower;
                            sMin = "Chernov";
                        }

                        
                        double kChernovUpper = 1.0 + mu + Math.Sqrt(-2.0 * mu * Math.Log(Q));
                        if (kChernovUpper < kMax) {
                            kMax = kChernovUpper;
                            sMax = "Chernov";
                        }
                        

                        double kCantelli = mu + Math.Sqrt(n * p * (1.0 - p) * P / Q);
                        if (kCantelli < kMax) {
                            kMax = kCantelli;
                            sMax = "Cantelli";
                        }

                        int k = d.InverseLeftProbability(P);
                        Assert.IsTrue(kMin <= k);
                        Assert.IsTrue(k <= kMax);

                        Console.WriteLine("n={0} p={1}: P={2} {3} ({4}) <= {5} <= {6} ({7})",
                            n, p, P, kMin, sMin, k, kMax, sMax);
                    }
                }
            }

        }



        [TestMethod]
        public void NegativeBinomialInverseCdf () {

            foreach (double r in new double[] { 0.1, 2, 34 }) {
                foreach (double p in new double[] { 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99 }) {
                    double q = 1.0 - p;
                    NegativeBinomialDistribution d = new NegativeBinomialDistribution(r, p);
                    foreach (double P in new double[] { 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99 }) {

                        double Q = 1.0 - P;

                        double kMin = 0.0;
                        string sMin = "Minimum";
                        double kMax = Double.PositiveInfinity;
                        string sMax = "Maximum";

                        double mu = r * p / q;

                        double kMarkov = mu / Q;
                        if (kMarkov < kMax) {
                            kMax = kMarkov;
                            sMax = "Markov";
                        }

                        /*
                        double kMedian;
                        if (P < 0.5) {
                            kMedian = Math.Ceiling(mu);
                            if (kMedian < kMax) {
                                kMax = kMedian;
                                sMax = "Median";
                            }
                        } else {
                            kMedian = Math.Floor(mu);
                            if (kMedian > kMin) {
                                kMin = kMedian;
                                sMin = "Median";
                            }
                        }
                        */

                        /*
                        double kChernovLower = mu - Math.Sqrt(-2.0 * mu * Math.Log(P));
                        if (kChernovLower > kMin) {
                            kMin = kChernovLower;
                            sMin = "Chernov";
                        }


                        double kChernovUpper = 1.0 + mu + Math.Sqrt(-2.0 * mu * Math.Log(Q));
                        if (kChernovUpper < kMax) {
                            kMax = kChernovUpper;
                            sMax = "Chernov";
                        }
                        */

                        double kCantelliLower = mu - Math.Sqrt(p * r * Q / P) / q;
                        if (kCantelliLower > kMin) {
                            kMin = kCantelliLower;
                            sMin = "Cantelli";
                        }

                        double kCantelliUpper = mu + Math.Sqrt(p * r * P / Q) / q;
                        if (kCantelliUpper < kMax) {
                            kMax = kCantelliUpper;
                            sMax = "Cantelli";
                        }

                        int k = d.InverseLeftProbability(P);
                        Assert.IsTrue(kMin <= k);
                        Assert.IsTrue(k <= kMax);

                        Console.WriteLine("r={0} p={1}: P={2} {3} ({4}) <= {5} <= {6} ({7})",
                            r, p, P, kMin, sMin, k, kMax, sMax);
                    }
                }
            }


        }

        [TestMethod]
        public void PoissonInverseCdf () {

            foreach (double mu in new double[] { 0.1, 2.0, 30.0 }) {
                PoissonDistribution d = new PoissonDistribution(mu);
                foreach (double P in new double[] { 0.01, 0.1, 0.3 , 0.5, 0.7, 0.9, 0.99 }) {

                    double Q = 1.0 - P;

                    double kMin = 0.0;
                    string sMin = "Zero";
                    double kMax = Double.PositiveInfinity;
                    string sMax = "Infinity";

                    double kMarkov = mu / Q;
                    if (kMarkov < kMax) {
                        kMax = kMarkov;
                        sMax = "Markov";
                    }

                    double kMedian;
                    if (P < 0.5) {
                        kMedian = mu + 1.0 / 3.0;
                        if (kMedian < kMax) {
                            kMax = kMedian;
                            sMax = "Median";
                        }
                    } else {
                        kMedian = mu - Math.Log(2.0);
                        if (kMedian > kMin) {
                            kMin = kMedian;
                            sMin = "Median";
                        }
                    }

                    double kChernov = mu - Math.Sqrt(-2.0 * mu * Math.Log(P));
                    if (kChernov > kMin) {
                        kMin = kChernov;
                        sMin = "Chernov";
                    }

                    double kCantelli = mu + Math.Sqrt(mu * P / Q); 
                    if (kCantelli < kMax) {
                        kMax = kCantelli;
                        sMax = "Cantelli";
                    }

                    int k = d.InverseLeftProbability(P);
                    Assert.IsTrue(kMin <= k);
                    Assert.IsTrue(k <= kMax);

                    Console.WriteLine("mu={0} P={1} {2} ({3}) <= {4}<= {5} ({6})", 
                        mu, P, kMin, sMin, k, kMax, sMax);
                }
            }

        }

        [TestMethod]
        public void PredictionVariance () {

            double a0 = -2.0;
            double b0 = 3.0;

            double[] xValues = TestUtilities.GenerateUniformRealValues(1.0, 10.0, 8);

            ContinuousDistribution eDistribuiton = new NormalDistribution();

            Sample zSample = new Sample();
            Sample vSample = new Sample();
            Sample uSample = new Sample();

            BivariateSample pSample = new BivariateSample();
            Sample sSample = new Sample();

            Sample cbbSample = new Sample();
            Sample cabSample = new Sample();
            Sample caaSample = new Sample();

            for (int i = 0; i < 10000; i++) {
                Random rng = new Random(i);
                BivariateSample sample = new BivariateSample();
                foreach (double x in xValues) {
                    double y = a0 + b0 * x + eDistribuiton.GetRandomValue(rng);
                    sample.Add(x, y);
                }

                int n = sample.Count;
                double mx = sample.X.Mean;
                double my = sample.Y.Mean;
                double cxx = sample.X.Variance;
                double cyy = sample.Y.Variance;
                double cxy = sample.Covariance;

                double b = cxy / cxx;
                double a = my - b * mx;
                pSample.Add(a, b);

                double r = cxy / Math.Sqrt(cxx * cyy);
                ContinuousDistribution rDistribution = new PearsonRDistribution(n);
                double pr = rDistribution.RightProbability(r) * 2.0;

                double F = r * r / (1.0 - r * r) * (n - 2);
                ContinuousDistribution fDistribution = new FisherDistribution(1, n - 2);
                double pF = fDistribution.RightProbability(F);

                BivariateSample residuals = new BivariateSample();

                double s2 = 0.0;
                foreach(XY point in sample) {
                    double rs = point.Y - (a + b * point.X);
                    residuals.Add(point.X, rs);
                    s2 += rs * rs;
                }
                s2 = s2 / (n - 2);
                sSample.Add(s2);

                double cbb = s2 / cxx / n;
                double cab = -mx * cbb;
                double caa = (cxx + mx * mx) * cbb;

                cbbSample.Add(cbb);
                cabSample.Add(cab);
                caaSample.Add(caa);

                LinearRegressionResult fit = sample.LinearRegression();

                double x0 = 0.0;
                double yp = fit.Intercept.Value + x0 * fit.Slope.Value;
                double y0 = a0 + b0 * x0 + eDistribuiton.GetRandomValue(rng);
                zSample.Add(yp - y0);

                ColumnVector c = new ColumnVector(1.0, x0);
                double v = c.Transpose * (fit.Parameters.CovarianceMatrix) * c;
                vSample.Add(v);

                double u = s2 * (1.0 + 1.0 / n + MoreMath.Sqr(x0 - mx) / (cxx * n));
                uSample.Add(u);
               
            }

        }


        [TestMethod]
        public void Seperate () {

            int m = 3;

            int n = 100;
            double[] y = new double[n];

            Random rng = new Random(1);
            for (int i = 0; i < y.Length; i++) {
                y[i] = 1.0 + 2.0 * i / n - 1.0 * i * i / n / n + (i % 3 - 1.0) + rng.NextDouble();
            }

            // define moving average filter
            double[] f;
            if (m % 2 == 0) {
                f = new double[m + 1];
                for (int i = 1; i < m; i++) f[i] = 1.0 / m;
                f[0] = 1.0 / m / 2.0;
                f[m] = 1.0 / m / 2.0;
            } else {
                f = new double[m];
                for (int i = 0; i < m; i++) f[i] = 1.0 / m;
            }

            // extract the trend
            double[] t = new double[y.Length - f.Length];
            for (int i = 0; i < t.Length; i++) {
                for (int j = 0; j < f.Length; j++) {
                    t[i] += y[i + j] * f[j];
                }
            }

            // extract the seasonal cycle
            double[] s = new double[m];
            int[] c = new int[m];
            int si = 0;
            for (int i = 0; i < t.Length; i++) {
                s[si % m] += y[i + m / 2] - t[i];
                c[si % m]++;
                si++;
            }
            for (int i = 0; i < s.Length; i++) s[i] = s[i] / c[i];

            // extract the remainder
            double[] z = new double[t.Length];
            for (int i = 0; i < z.Length; i++) {
                z[i] = y[i + m / 2] - t[i] - s[i % m];
            }

        }

        [TestMethod]
        public void Filter () {


            double[] y = new double[] { 3.0, 1.5, 2.0, 4.5, 3.0, 3.5, 6.0, 4.5, 5.0, 7.5, 6.0, 6.5, 9.0, 7.5, 8.0 };
            //double[] y = new double[] { 1.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 5.0, 3.0, 4.0, 5.0, 6.0, 4.0, 5.0, 6.0, 7.0 };

            double[] f = new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
            //double[] f = new double[] { 0.25, 0.5, 0.25 };
            //double[] f = new double[] { 0.25, 0.25, 0.25, 0.25 };

            double[] z = new double[y.Length - f.Length];
            for (int i = 0; i < z.Length; i++) {
                for (int j = 0; j < f.Length; j++) {
                    z[i] += y[i + j] * f[j];
                }
            }

            Complex[] yp = new Complex[y.Length + f.Length];
            for (int i = 0; i < y.Length; i++) yp[i] = y[i];

            Complex[] fp = new Complex[yp.Length];
            for (int i = 0; i < f.Length; i++) fp[i] = f[i];

            FourierTransformer fft = new FourierTransformer(yp.Length);
            Complex[] ypt = fft.Transform(yp);
            Complex[] fpt = fft.Transform(fp);

            Complex[] zpt = new Complex[fft.Length];
            for (int i = 0; i < zpt.Length; i++) {
                zpt[i] = ypt[i] * fpt[i];
            }
            Complex[] zp = fft.InverseTransform(zpt);

        }

        [TestMethod]
        public void PowerSpectrum () {

            TimeSeries series = new TimeSeries(1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0);
            //TimeSeries series = new TimeSeries(0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0);
            double[] ps = series.PowerSpectrum();

        }

        [TestMethod]
        public void CoulombAsymptotic () {

            
            double F1 = AdvancedMath.CoulombF(4, -30.0, 1.0);
            double G1 = AdvancedMath.CoulombG(4, -30.0, 1.0);
            Console.WriteLine("{0} {1}", F1, G1);
            

            for (int L = 0; L < 8; L++) {

                foreach(double eta in new double[] { -30.1, -5.1, -0.31, 0.0, 0.51, 3.1, 50.1}) {

                    foreach (double rho in new double[] { 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0 }) {
                        //double rho = 32.0 + (L * L + eta * eta) / 2.0;
                        //for (int i = 0; i < 8; i++) {

                        SolutionPair r = AdvancedMath.Coulomb(L, eta, rho);
                        //SolutionPair r = AdvancedMath.Coulomb_Asymptotic(L, eta, rho);
                        Console.WriteLine("{0} {1} {2} : {3} {4} {5} {6} : {7}",
                            L, eta, rho,
                            r.FirstSolutionValue, r.FirstSolutionDerivative,
                            r.SecondSolutionValue, r.SecondSolutionDerivative,
                            r.FirstSolutionDerivative * r.SecondSolutionValue - r.FirstSolutionValue * r.SecondSolutionDerivative
                        );

                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(
                            r.FirstSolutionDerivative * r.SecondSolutionValue, - r.FirstSolutionValue * r.SecondSolutionDerivative, 1.0,
                            TestUtilities.TargetPrecision * 8.0
                        ));

                        //rho = 1.25 * rho;
                    }
                }

            }

        }

    }

#if FUTURE

    // Uniform asymptotic exansion for Bessel functions of large order is summarized in A&S 9.3 and DLMF 
    
    // The basic form is:
    //    J_{\nu}(\nu z) = \left(\frac{4 \zeta}{1 - z^2}\right)^{1/4}

    // Here \zeta is determined from z via a complicated expression that maps 0->\infinity, 1->0, and \infinity->-\infinity.
    
    // What makes this difficult is that the expressions for \zeta and a_k(\zeta) and b_k(\zeta) have canceling singularities
    // at z = 1 (\zeta = 0). While the functions are ultimimately analytic, the singularities that must cancel are quite
    // strong. (E.g. for b_0(\zeta) singularities of order \zeta^{-3} and \zeta^{-3/2} cancel near \zeta = 0, leaving a constant
    // value. The degree of cancelation gets higher with increasing k.) This makes the practical evaluation of these functions
    // very much more difficult than a cursory inspection leads one to expect.

    // We follow N. M. Temme, "Numerical algorithms for Airy-type asymptotic expansions", Numerical Algorithms 15 (1997) 207
    // (available at http://oai.cwi.nl/oai/asset/2178/2178A.pdf)

    internal static class BesselUae {

        private static readonly double CbrtTwo = Math.Pow(2.0, 1.0 / 3.0);

        private static readonly double[] c = new double[] {
            0.0,
            -1.0,
            3.0 / 10.0,
            -32.0 / 175.0,
            1037.0 / 7875.0,
            -103727.0 / 1010625.0,
            33060241.0 / 394143750.0,
            -4393499056.0 / 62077640625.0,
            15356175508.0 / 251266640625.0,
            -296160295945538.0 / 5514046428515625.0,
            2095680660298957.0 / 43788015755859375.0,
            -1178774779896086224.0 / 27333708572091796875.0,
            305517380952020655742.0 / 7790106943046162109375.0,
            -970178984722577536907354.0 / 26992720557654951708984375.0,
            43257635097401811452791703.0 / 1304648160286655999267578125.0,
            -1413145527493419079807132911728.0 / 45940812876930433457845458984375.0,
            8477921788590209817825424263485828.0 / 295629130863047339301235528564453125.0,
            -19355768761912897160895124453273481957.0 / 720842364087730428996179297149658203125.0
        };

        // Zeta is given by A&S 9.3.38 and 9.3.39:
        
        //

        public static void EvaluateZeta (double z, out double zeta, out double phi) {

            double e = z - 1.0;
            if (Math.Abs(e) < 0.125) {
                double ek = e;
                zeta = c[1] * e;
                phi = c[1];
                for (int k = 2; k < c.Length; k++) {
                    phi += c[k] * ek;
                    ek *= e;
                    zeta += c[k] * ek;
                }
                zeta = CbrtTwo * zeta;
                phi = CbrtTwo * Math.Pow(-phi / (1.0 + e / 2.0), 1.0 / 4.0);
            } else {
                if (z < 1.0) {
                    double s = Math.Sqrt(1.0 - z * z);
                    zeta = Math.Pow(3.0 / 2.0 * (Math.Log((1.0 + s) / z) - s), 2.0 / 3.0);
                } else {
                    zeta = -Math.Pow(3.0 / 2.0 * (Math.Sqrt(z * z - 1.0) - Math.Acos(1.0 / z)), 2.0 / 3.0);
                }
                phi = Math.Pow(4.0 * zeta / (1.0 - z * z), 1.0 / 4.0);
            }

        }

        // A & S give expressions for a_k(\zeta) and b_k(\zeta). They start out

        //   a_0 = 1

        //   b_0 = 

        // For general k,

        //   a_k = \sum_{j=0}^{2k} \mu_j \zeta^{-3j/2} u_{2k-j}
        //   b_k = - \sum_{j=0}^{2k+1} \lambda_j \zeta^{-3(j+1)/2} u_{2k+1-j}

        // Here u_j is a polynomial of order 3j, even for even j and odd for odd j. The argument of u is t=1/\sqrt{1-z^2}. For even j, this
        // means u is actually a polynomial of order 3j/2 in t^2=1/(1-z^2), which is no problem. For odd j and z > 1 (\zeta < 0), the argument of
        // u is actually pure imaginary, but the result is multiplied by \zeta to a half-power, which is also pure imaginary, giving a real
        // product. We deal with this situation by taking out one power of t and writing \zeta^{1/2} t = \phi^2 / 2.

        // Finally, the coefficients \mu_j and \lambda_j are pure numbers, given by A&S 9.3.41
        //   \mu_j = -\frac{6j+1}{6j-1} \lambda_j
        // So
        //   \lambda_0 = 1, \lambda_1 = 5/48, \lambda_2 = 385/4608, \lambda_3 = 85085/663552, \lambda_4 = 37182145/127401984
        //   \mu_0 = 1, \mu_1 = -7/48, \mu_2 = -455/4608, \mu_3 = -95095/663552, \mu_4 = -40415375/127401984

        // All these expressions have canceling singularities at \zeta = 0 (z = 1), and therefore cannot be used near \zeta ~ 0.

        // Temme's paper derives a systematic approach for producing the power series for these quantities near \zeta ~ 0.
        // Actually, he uses the trivially re-scaled parameter \eta where \zeta = 2^{1/3} \eta.

        // All Temme's series are derived from the series for
        //   \psi = \frac{5}{16 \zeta^2} + \frac{\zeta z^2 (z^2 + 4)}{4 (z^2 - 1)^3}
        //        = \sum_{k=0}^{\infty} \psi^{(k)} \eta^k = 2^{1/3} \left[ \frac{1}{70} + \frac{2}{75} \eta + \frac{138}{13475} \eta^2 + \cdots \right]
        // He shows that
        //   b_0 = \sum_{k=0}^{\infty} b_0^{(k)} \eta^k  b_0^{(k)} = \frac{\psi^{(k)}}{2k + 1}
        // and 

        private static readonly double[] PsiCoefficients = {
            1.0 / 70.0,
            2.0 / 75.0,
            138.0 / 13475.0,
            -296.0 / 73125.0,
            -38464.0 / 7074375.0,
            -117472.0 / 72515625.0,
            4835327648.0 / 6986122171875.0,
            196235264.0 / 251266640625.0,
            1660313255231104.0 / 7920224923939453125.0, 
            -66399236608.0 / 712340926171875.0,
            -3601874963621281792.0 / 37027051519416943359375.0,
            -31198440677844992.0 / 1259270353112255859375.0,
            339118989445333019205632.0 / 30305783606970600781787109375.0,
            1778526038024190427136.0 / 158627254804774859619140625.0,
            40948053418677371831115518818104241381.0 / 14820745152267765143433288288000000000000000.0,
            -28685174387749486839072824921671801053073.0 / 22675740082969680669452931080640000000000000000.0,
            -27122073903808716897938907648268793696783.0 / 21978025003493690495008225508928000000000000000.0
        };

        private static readonly double[] b0 = ComputeBZeroCoefficients();

        private static readonly double[] a1 = ComputeACoefficients(-1.0 / 225.0, b0);

        private static readonly double[] b1 = ComputeBCoefficients(a1);

        private static readonly double[] a2 = ComputeACoefficients(151439.0 / 218295000.0, b1);

        public static double[] ComputeBZeroCoefficients () {
            double[] bCoefficients = new double[PsiCoefficients.Length];
            for (int i = 0; i < bCoefficients.Length; i++) {
                bCoefficients[i] = PsiCoefficients[i] / (2 * i + 1);
            }
            return (bCoefficients);
        }

        public static double[] ComputeACoefficients (double aZeroCoefficient, double[] bCoefficients) {
            double[] aCoefficients = new double[bCoefficients.Length - 2];
            aCoefficients[0] = aZeroCoefficient;
            for (int i = 0; i < aCoefficients.Length - 1; i++) {
                double t = 0.0;
                for (int j = 0; j <= i; j++) {
                    t += PsiCoefficients[j] * bCoefficients[i - j];
                }
                t /= (i + 1);
                t -= (i + 2) * bCoefficients[i + 2] / 2.0;
                aCoefficients[i + 1] = t;
                Console.WriteLine("a {0} {1:R}", i + 1, t);
            }
            return (aCoefficients);
        }

        public static double[] ComputeBCoefficients (double[] aCoefficients) {
            double[] bCoefficients = new double[aCoefficients.Length - 2];
            for (int i = 0; i < bCoefficients.Length; i++) {
                double t = 0.0;
                for (int j = 0; j <= i; j++) {
                    t += PsiCoefficients[j] * aCoefficients[i - j];
                }
                t -= (i + 1) * (i + 2) * aCoefficients[i + 2] / 2.0;
                t /= (2 * i + 1);
                bCoefficients[i] = t;
                Console.WriteLine("b {0} {1:R}", i, t);
            }
            return (bCoefficients);
        }

        private static double U0 (double t2) { return(1.0); }

        private static double U1 (double t2) { return ((3.0 - 5.0 * t2) / 24.0); }

        private static double U2 (double t2) { return(t2 * (81.0 - 462.0 * t2 + 385.0 * t2 * t2) / 1152.0); }

        private static double U3 (double t2) { return (t2 * (30375.0 - 369603.0 * t2 + 765765.0 * t2 * t2 - 425425.0 * t2 * t2 * t2) / 414720.0); }

        private static double U4 (double t2) { return (t2 * t2 * (4465125.0 - 94121676.0 * t2 + 349922430.0 * t2 * t2 - 446185740.0 * t2 * t2 * t2 + 185910725.0 * t2 * t2 * t2 * t2) / 39813120.0); }

        public static double EvaluateB0 (double zeta, double z, double phi) {
            double eta = zeta / CbrtTwo;
            if (Math.Abs(eta) < 0.25) {
                double etak = eta;
                double f = b0[0] + b0[1] * eta;
                for (int k = 2; k < b0.Length; k++) {
                    etak *= eta;
                    f += b0[k] * etak;
                }
                return (CbrtTwo * f);
            } else {
                double t2 = 1.0 / (1.0 - z * z);
                return (-((phi * phi / 2.0) / zeta * U1(t2) + 5.0 / 48.0 / (zeta * zeta)));
            }

        }

        public static double EvaluateA1 (double zeta, double z, double phi) {
            double eta = zeta / CbrtTwo;
            if (Math.Abs(eta) < 0.25) {
                double etak = eta;
                double f = a1[0] + a1[1] * eta;
                for (int k = 2; k < a1.Length; k++) {
                    etak *= eta;
                    f += a1[k] * etak;
                }
                return (f);
            } else {
                double t2 = 1.0 / (1.0 - z * z);
                double result = U2(t2) - 7.0 / 48.0 * (phi * phi / (zeta * zeta) / 2.0) * U1(t2) - 455.0 / 4608.0 / (zeta * zeta * zeta);
                return (result);
            }
        }

        public static double EvaluateB1 (double zeta, double z, double phi) {
            double eta = zeta / CbrtTwo;
            if (Math.Abs(eta) < 0.25) {
                double etak = eta;
                double f = b1[0] + b1[1] * eta;
                for (int k = 2; k < b1.Length; k++) {
                    etak *= eta;
                    f += b1[k] * etak;
                }
                return (CbrtTwo * f);
            } else {
                double t2 = 1.0 / (1.0 - z * z);
                double result = -(
                    (phi * phi / 2.0) / zeta  * U3(t2) +
                    5.0 / 48.0 / (zeta * zeta) * U2(t2) +
                    385.0 / 4608.0 * (phi * phi / 2.0) / (zeta * zeta * zeta * zeta) * U1(t2) +
                    85085.0 / 663552.0 / (zeta * zeta * zeta * zeta * zeta)
                );
                return (result);
            }
        }

        public static double EvaluateA2 (double zeta, double z, double phi) {
            double eta = zeta / CbrtTwo;
            if (Math.Abs(eta) < 0.25) {
                double etak = eta;
                double f = a2[0] + a2[1] * eta;
                for (int k = 2; k < a2.Length; k++) {
                    etak *= eta;
                    f += a2[k] * etak;
                }
                return (f);
            } else {
                return (0.0);
            }
        }

        public static double Evaluate (double nu, double x) {

            // Write argument as x = \nu z.
            double z = x / nu;

            // Compute \zeta and \phi.
            double zeta, phi;
            EvaluateZeta(z, out zeta, out phi);
            Console.WriteLine("z'={0:R} zeta'={1:R} phi'={2:R}", z, zeta, phi);

            // Compute Airy functions.
            double nu23 = Math.Pow(nu, 2.0 / 3.0);
            SolutionPair airy = AdvancedMath.Airy(nu23 * zeta);
            Console.WriteLine("arg'={0:R} Ai={1:R} Ai'={2:R}", nu23 * zeta, airy.FirstSolutionValue, airy.FirstSolutionDerivative);

            // Compute a and b coefficients.
            const double a0 = 1.0;
            double b0 = EvaluateB0(zeta, z, phi);
            double a1 = EvaluateA1(zeta, z, phi);
            double b1 = EvaluateB1(zeta, z, phi);
            double a2 = EvaluateA2(zeta, z, phi);

            // Combine to compute Bessel functions.
            double nu2 = nu * nu;
            double nu13 = Math.Pow(nu, 1.0 / 3.0);
            double nu43 = nu23 * nu23;
            double a = a0 + a1 / nu2 + a2 / (nu2 * nu2);
            double b = b0 + b1 / nu2;
            double J = (phi / nu13) * (airy.FirstSolutionValue * a + airy.FirstSolutionDerivative / nu43 * b);
            double Y = (phi / nu13) * (airy.SecondSolutionValue * a + airy.SecondSolutionDerivative / nu43 * b);

            return (J);

        }

    }


    

    [TestClass]
    public class FutureTest {

        // See functions at http://www.sfu.ca/~ssurjano/optimization.html

        [TestMethod]
        public void STA () {

            // Has local minimum when any coordinate is at one of two values, global when all coordinates at one of them.
            Func<IList<double>, double> fStyblinskiTang = (IList<double> x) => {
                double fst = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    double x1 = x[i];
                    double x2 = MoreMath.Sqr(x1);
                    fst += x2 * (x2 - 16.0) + 5.0 * x1;
                }
                return (fst / 2.0);
            };

            // Asymptotic "minima" at infinity, and (3,1/2)->0
            Func<IList<double>, double> fBeale = (IList<double> x) =>
                MoreMath.Sqr(1.5 - x[0] + x[0] * x[1]) +
                MoreMath.Sqr(2.25 - x[0] + x[0] * x[1] * x[1]) +
                MoreMath.Sqr(2.625 - x[0] + x[0] * x[1] * x[1] * x[1]);


            Func<IList<double>, double> fEggholder = (IList<double> x) => {
                double y47 = x[1] + 47.0;
                return(-y47 * Math.Sin(Math.Sqrt(Math.Abs(y47 + x[0] / 2.0))) - x[0] * Math.Sin(Math.Sqrt(Math.Abs(x[0] - y47))));
            };

            // Many local minima, global minimum at (0,0)->0
            Func<IList<double>, double> fAckley = (IList<double> x) => {
                double s = 0.0;
                double c = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    s += x[i] * x[i];
                    c += Math.Cos(2.0 * Math.PI * x[i]);
                }
                return (-20.0 * Math.Exp(-0.2 * Math.Sqrt(s / x.Count)) - Math.Exp(c / x.Count) + 20.0 + Math.E);
            };

            // Burkin has a narrow valley, not aligned with any axis, punctuated with many tiny "wells" along its bottom.
            // The deepest well is at (-10,1)-> 0.
            Func<IList<double>, double> fBurkin = (IList<double> x) => 100.0 * Math.Sqrt(Math.Abs(x[1] - 0.01 * x[0] * x[0])) + 0.01 * Math.Abs(x[0] + 10.0);

            Func<IList<double>, double> threeHumpCamel = (IList<double> x) => {
                double x2 = x[0] * x[0];
                return (2.0 * x2 - 1.05 * x2 * x2 + x2 * x2 * x2 / 6.0 + (x[0] + x[1]) * x[1]);
            };

            // Easom has many lobal extra ripples and then a deep but narrow global minimum at (\pi, \pi) -> -1
            Func<IList<double>, double> fEasom = (IList<double> x) => -Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Exp(-(MoreMath.Sqr(x[0] - Math.PI) + MoreMath.Sqr(x[1] - Math.PI)));

            // Test over [-500,500], minimum at (420.969,...) -> -418.983*d, many local minima
            Func<IList<double>, double> fSchwefel = (IList<double> x) => {
                double s = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    s += x[i] * Math.Sin(Math.Sqrt(Math.Abs(x[i])));
                }
                return (-s);
            };

            Func<IList<double>, double> function = fSchwefel;

            int n = 4;

            ColumnVector start = new ColumnVector(n);
            //for (int i = 0; i < start.Dimension; i++) start[i] = 1.0;

            Interval[] box = new Interval[n];
            for (int i = 0; i < box.Length; i++) {
                box[i] = Interval.FromEndpoints(-500.0, 500.0);
            }
            //box[0] = Interval.FromEndpoints(-15.0, -5.0);
            //box[1] = Interval.FromEndpoints(-3.0, 3.0);

            MultiFunctionMath.FindGlobalMinimum(function, box);
            //MultiFunctionMath.FindExtremum_Amobea(fStyblinskiTang, start, new EvaluationSettings() { RelativePrecision = 1.0E-10, AbsolutePrecision = 1.0E-12, EvaluationBudget = 1000 });

        }

        /*
        [TestMethod]
        public void Rosenbock () {

            Func<IList<double>, double> fRosenbrock = (IList<double> x) =>
                MoreMath.Sqr(1.0 - x[0]) + 100.0 * MoreMath.Sqr(x[1] - x[0] * x[0]);
            //Func<IList<double>, double> f = (IList<double> x) => 7.0 + 6.0 * MoreMath.Sqr(x[0] - 1.0) + 5.0 * MoreMath.Sqr(x[1] - 2.0) + 4.0 * (x[0] - 1.0) * (x[1] - 2.0);
            Func<IList<double>, double> fGoldsteinPrice = (IList<double> p) => {
                double x = p[0]; double yy = p[1]; return (
                    (1 + MoreMath.Pow(x + yy + 1, 2) * (19 - 14 * x + 3 * x * x - 14 * yy + 6 * x * yy + 6 * yy * yy)) *
                    (30 + MoreMath.Pow(2 * x - 3 * yy, 2) * (18 - 32 * x + 12 * x * x + 48 * yy - 36 * x * yy + 27 * yy * yy))
                );
            };
            Func<IList<double>, double> fBeale = (IList<double> x) =>
                MoreMath.Sqr(1.5 - x[0] + x[0] * x[1]) +
                MoreMath.Sqr(2.25 - x[0] + x[0] * x[1] * x[1]) +
                MoreMath.Sqr(2.625 - x[0] + x[0] * x[1] * x[1] * x[1]);
            Func<IList<double>, double> fStyblinskiTang = (IList<double> x) => {
                double fst = 0.0;
                for (int i = 0; i < x.Count; i++) {
                    double x1 = x[i];
                    double x2 = MoreMath.Sqr(x1);
                    fst += x2 * (x2 - 16.0) + 5.0 * x1;
                }
                return (fst / 2.0);
            };
            Func<IList<double>, double> fMcCormick = (IList<double> x) => Math.Sin(x[0] + x[1]) + MoreMath.Sqr(x[0] - x[1]) - 1.5 * x[0] + 2.5 * x[1] + 1.0;
            // Levi causes some problems we will need to work out.
            Func<IList<double>, double> fLevi = (IList<double> x) =>
                MoreMath.Sqr(Math.Sin(3.0 * Math.PI * x[0])) +
                MoreMath.Sqr(x[0] - 1.0) * (1.0 + MoreMath.Sqr(Math.Sin(3.0 * Math.PI * x[1]))) +
                MoreMath.Sqr(x[1] - 1.0) * (1.0 + MoreMath.Sqr(Math.Sin(2.0 * Math.PI * x[1])));

            var f = fStyblinskiTang;

            double[] y = new double[] { 1.0, 1.0 };
            //double[][] points; double[] values;
            QuadraticInterpolationModel model = QuadraticInterpolationModel.Construct(f, y, 1.0);
            int evaluationcount = 6;

            double D = 1.0;
            while (true) {

                // find maximum and minimum points
                int iMax = 0; double fMax = model.values[0];
                for (int i = 1; i < model.values.Length; i++) {
                    if (model.values[i] > fMax) { iMax = i; fMax = model.values[i]; }
                }

                // find the minimum of the model
                double[] z = model.FindMinimum(D);
                double expectedValue = model.Evaluate(z);
                double[] point = new double[z.Length];
                for (int i = 0; i < point.Length; i++) point[i] = model.origin[i] + z[i];
                double value = f(point);
                evaluationcount++;

                // Adjust the trust region radius based on the ratio of the actual reduction to the expected reduction.
                double r = (model.MinimumValue - value) / (model.MinimumValue - expectedValue);
                if (r < 0.25) {
                    // If we achieved less than 25% of the expected reduction, reduce the trust region.
                    D = D / 2.0;
                } else if (r > 0.75) {
                    // If we achieved at least 75% of the expected reduction, increase the trust region.
                    D = 2.0 * D;
                }

                // if the new value is less than the maximum, accept it
                if (value < fMax) {
                    model.ReplacePoint(iMax, point, z, value);
                } else {
                    // Otherwise, reduce the trust region and try again
                }


            }
        }
        */

        [TestMethod]
        public void Bell1 () {

            RectangularMatrix A = new RectangularMatrix(new double[,] {
                {1.0, 0.0, 0.0, 0.0, 2.0},
                {0.0, 0.0, 3.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 4.0, 0.0, 0.0, 0.0}
            });

            //A = A.Transpose();

            SingularValueDecomposition svd = A.SingularValueDecomposition();

            Console.WriteLine("values");
            for (int i = 0; i < svd.Dimension; i++) {
                Console.WriteLine(svd.SingularValue(i));
            }


            Console.WriteLine("left");
            WriteMatrix(svd.LeftTransformMatrix());


            Console.WriteLine("right");
            WriteMatrix(svd.RightTransformMatrix());

            Console.WriteLine("product");
            WriteMatrix(svd.LeftTransformMatrix().Transpose() * A * svd.RightTransformMatrix());

        }


        [TestMethod]
        public void TestCumulant2 () {

            int n = 8;

            //Distribution[] set = new Distribution[] { new UniformDistribution(), new ExponentialDistribution(2.0), new GammaDistribution(3.0), new LogisticDistribution(-4.0, 3.0), new NormalDistribution(-1.0, 2.0), new WaldDistribution(1.0, 2.0) };
            DiscreteDistribution[] set = new DiscreteDistribution[] { new PoissonDistribution(0.25), new DiscreteUniformDistribution(0, 10) };

            foreach (DiscreteDistribution d in set) {
            //foreach (Distribution d in set) {
                Console.WriteLine(d.GetType().Name);

                // From cumulants to central and raw moments

                double[] inK = new double[n];
                for (int r = 0; r < n; r++) inK[r] = d.Cumulant(r);

                double[] outC = MomentMath.CumulantToCentral(inK);
                for (int r = 0; r < n; r++) Console.WriteLine("r={0} K={1} -> C={2} v C={3}", r, inK[r], d.MomentAboutMean(r), outC[r]);
                for (int r = 0; r < n; r++) Assert.IsTrue(Math.Abs(outC[r] - d.MomentAboutMean(r)) <= 1.0E-14 * Math.Abs(outC[r]));

                double[] outM = MomentMath.CumulantToRaw(inK);
                for (int r = 0; r < n; r++) Assert.IsTrue(Math.Abs(outM[r] - d.Moment(r)) <= 1.0E-14 * Math.Abs(outM[r]));

                // From central moments to cumulants and raw moments

                double[] inC = new double[n];
                for (int r = 0; r < n; r++) inC[r] = d.MomentAboutMean(r);

                double[] outK = MomentMath.CentralToCumulant(d.Mean, inC);
                for (int r = 0; r < n; r++) Console.WriteLine("r={0} C={1:R} -> K={2:R} v K={3:R}", r, inC[r], outK[r], d.Cumulant(r));
                for (int r = 0; r < n; r++) Assert.IsTrue(Math.Abs(outK[r] - d.Cumulant(r)) <= 1.0E-13 * Math.Abs(outK[r]));
                // moved to 10E-13 due to K_4 of Logistic; why is there even that much disagreement?

                double[] outM2 = MomentMath.CentralToRaw(d.Mean, inC);
                for (int r = 0; r < n; r++) Assert.IsTrue(Math.Abs(outM2[r] - d.Moment(r)) <= 1.0E-14 * Math.Abs(outM2[r]));

                // From raw moments to central moments and cumulants
                // This is unstable.

            }

        }

        [TestMethod]
        public void TestNormal () {
            NormalDistribution n = new NormalDistribution();
            Console.WriteLine(n.InverseLeftProbability(1.0E-10));
            Console.WriteLine(n.InverseRightProbability(1.0E-10));
            Console.WriteLine(n.InverseLeftProbability(1.0E-300));
            Console.WriteLine(n.InverseRightProbability(1.0E-300));
            Console.WriteLine(n.InverseLeftProbability(1.0));
            //Console.WriteLine(n.InverseLeftProbability(0.26));
        }
      

        [TestMethod]
        public void Repro () {
            double I = 0.9998063099306;
            double a = 0.00034509911609819255;
            double b = 6.8453983996634218;
            var bd = new BetaDistribution(a, b);
            double x = bd.InverseLeftProbability(I);
            Console.WriteLine(x);
            Console.WriteLine(bd.LeftProbability(x));
            Console.WriteLine(I);
        }

        [TestMethod]
        public void BetaExercise () {

            double[] c = new double[] { 1.0 / 32768.0, 1.0 / 1024.0, 1.0 / 64.0, 1.0 / 8.0, 1.0 / 2.0, 1.0, 2.0, 8.0, 64.0, 1024.0, 32768.0 };

            double[] p = new double[] { 1.0 / 65536.0, 1.0 / 256.0, 1.0 / 16.0, 1.0 / 4.0, 1.0 - 1.0 / 4.0, 1.0 - 1.0 / 16.0, 1.0 - 1.0 / 256.0 };

            Stopwatch s = Stopwatch.StartNew();
            double sum = 0.0;
            foreach (double a in c) {
                foreach (double b in c) {
                    BetaDistribution BD = new BetaDistribution(a, b);
                    foreach (double x in p) {
                        //Console.WriteLine("a={0} b={1} P={2}", a, b, x);
                        //Console.WriteLine(BD.InverseLeftProbability(x));
                        sum += BD.InverseLeftProbability(x);
                    }
                }
            }
            s.Stop();
            Console.WriteLine(sum);
            Console.WriteLine(s.ElapsedMilliseconds);

        }

        /*
        [TestMethod]
        public void TestSin () {
            Console.WriteLine(Math.Sin(Math.Pow(2, 20)));
            Console.WriteLine(Math.Sin(Math.Pow(2, 30)));
            Console.WriteLine(Math.Sin(Math.Pow(2, 62)));
            Console.WriteLine(Math.Sin(Math.Pow(2, 63)));
            Console.WriteLine(Math.Sin(Math.Pow(2, 64)));
            Console.WriteLine(Math.Sin(Math.Pow(2, 65)));
        }

        [TestMethod]
        public void TestSplitting () {

            long granularity = 65536;
            //long granularity = 1048576;
            //long granularity = 16777216;

            double dsMax = 0.0; long ksMax = 0L;
            double dcMax = 0.0; long kcMax = 0L;

            double ym = Math.Pow(2.0, -6);

            for (long k = 30000L * granularity; k < 32768L * granularity; k++) {

                double x = k / ((double) granularity);

                double s0 = Math.Sin(x);
                double c0 = Math.Cos(x);

                if ((Math.Abs(s0) > 4.0 * ym) && (Math.Abs(c0) > 4.0 * ym)) continue;
                if ((Math.Abs(s0) < ym) || (Math.Abs(c0) < ym)) continue;
                //if ((Math.Abs(s0) < 0.0000152587890625) || (Math.Abs(c0) < 0.0000152587890625)) continue;

                double s1 = RangeReduction.Sin(x);
                double c1 = RangeReduction.Cos(x);

                double ds = Math.Abs((s0 - s1) / s0);
                double dc = Math.Abs((c0 - c1) / c0);

                if (ds > dsMax) { dsMax = ds; ksMax = k; }
                if (dc > dcMax) { dcMax = dc; kcMax = k; }  

            }

            double xsMax = ksMax / ((double) granularity);
            Console.WriteLine("Sin({0}/{1} = {2}) = {3:R} v {4:R} ({5})", ksMax, granularity, xsMax, Math.Sin(xsMax), RangeReduction.Sin(xsMax), dsMax);

            double xcMax = kcMax / ((double) granularity);
            Console.WriteLine("Cos({0}/{1} = {2}) = {3:R} v {4:R} ({5})", kcMax, granularity, xcMax, Math.Cos(xcMax), RangeReduction.Cos(xcMax), dcMax);

        }
        */

        [TestMethod]
        public void BesselExactTest () {

            double y;

            y = 2.0 / Math.PI;
            Console.WriteLine("{0:R}", (AdvancedMath.BesselJ(0.5, Math.PI / 2.0) - y) / y);

            y = MoreMath.Sqr(2.0 / Math.PI);
            Console.WriteLine("{0:R}", (AdvancedMath.BesselJ(1.5, Math.PI / 2.0) - y) / y);

            y = (2.0 / Math.PI) * (3.0 * MoreMath.Sqr(2.0 / Math.PI) - 1.0);
            Console.WriteLine("{0:R}", (AdvancedMath.BesselJ(2.5, Math.PI / 2.0) - y) / y);

            y = 0.0;
            Console.WriteLine("{0:R}", AdvancedMath.BesselJ(0.5, Math.PI) - y);

            y = Math.Sqrt(2.0) / Math.PI;
            Console.WriteLine("{0:R}", (AdvancedMath.BesselJ(1.5, Math.PI) - y) / y);

            y = 3.0 * Math.Sqrt(2.0) / MoreMath.Sqr(Math.PI);
            Console.WriteLine("{0:R}", (AdvancedMath.BesselJ(2.5, Math.PI) - y) / y);
        }

        [TestMethod]
        public void TestGammaValues () {

            Console.WriteLine((AdvancedMath.Gamma(0.5) - Math.Sqrt(Math.PI)) / Math.Sqrt(Math.PI));
            Console.WriteLine(AdvancedMath.Gamma(1) - 1.0);
            Console.WriteLine(AdvancedMath.Gamma(2) - 1.0);
            Console.WriteLine((AdvancedMath.Gamma(3) - 2.0) / 2.0);
            Console.WriteLine((AdvancedMath.Gamma(4) - 6.0) / 6.0);
            Console.WriteLine((AdvancedMath.Gamma(5) - 24.0) / 24.0);
            Console.WriteLine((AdvancedMath.Gamma(6) - 120.0) / 120.0);
            Console.WriteLine((AdvancedMath.Gamma(7) - 720.0) / 720.0);
            Console.WriteLine((AdvancedMath.Gamma(8) - 5040.0) / 5040.0);
        }

        [TestMethod]
        public void TestConstants () {

            double tp0 = 6.2831853071795864769;
            double tp1 = 2.0 * Math.PI;
            Console.WriteLine("{0:R} {1:R} {2:R}", tp0, tp1, (tp1 - tp0) / tp0);

            double stp0 = 2.5066282746310005024;
            double stp1 = Math.Sqrt(2.0 * Math.PI);
            Console.WriteLine("{0:R} {1:R} {2:R}", stp0, stp1, (stp1-stp0) / stp0);

            double sp0 = 1.7724538509055160273;
            double sp1 = Math.Sqrt(Math.PI);
            Console.WriteLine("{0:R} {1:R} {2:R}", sp0, sp1, (sp1 - sp0) / sp0);

            double spt0 = 1.2533141373155002512;
            double spt1 = Math.Sqrt(Math.PI / 2.0);
            Console.WriteLine("{0:R} {1:R} {2:R}", spt0, spt1, (spt1 - spt0) / spt0);

            double st0 = 1.7320508075688772935;
            double st1 = Math.Sqrt(3.0);
            Console.WriteLine("{0:R} {1:R} {2:R}", st0, st1, (st1-st0) / st0);

            double lt0 = 0.69314718055994530942;
            double lt1 = Math.Log(2.0);
            Console.WriteLine("{0:R} {1:R} {2:R}", lt0, lt1, (lt1 - lt0) / lt0);


        }

        [TestMethod]
        public void AiryTest () {
            Console.WriteLine("{0:R}", AdvancedMath.AiryAi(5.125));
            Console.WriteLine("{0:R}", AdvancedMath.Airy(5.125).FirstSolutionValue);
        }

        [TestMethod]
        public void ComputeChiTest () {
            foreach (double z in new double[] { 1.0E-4, 1.0E-2, 0.5, 0.75, 0.90, 0.95, 0.99, 1.0, 1.01, 1.05, 1.10, 1.15, 1.25, 1.5, 2.0, 10.0, 1.0E2, 1.0E4 }) {
                double zeta, phi;
                BesselUae.EvaluateZeta(z, out zeta, out phi);
                double a0 = 1.0;
                double b0 = BesselUae.EvaluateB0(zeta, z, phi);
                double a1 = BesselUae.EvaluateA1(zeta, z, phi);
                double b1 = BesselUae.EvaluateB1(zeta, z, phi);
                Console.WriteLine("z={0:R} zeta={1:R} phi={2:R} b0={3} a1={4} b1={5}", z, zeta, phi, b0, a1, b1);

                double nu = 1024.0;
                double nu3 = Math.Pow(nu, 1.0 / 3.0);
                SolutionPair p = AdvancedMath.Airy(nu3 * nu3 * zeta);
                Console.WriteLine("arg={0:R} Ai={1:R} Ai'={2:R}", nu3 * nu3 * zeta, p.FirstSolutionValue, p.FirstSolutionDerivative);
                double J0 = phi * (p.FirstSolutionValue / nu3 * a0);
                double J1 = phi * (p.FirstSolutionValue / nu3 * a0 + p.FirstSolutionDerivative / MoreMath.Pow(nu3, 5) * b0);
                double J2 = phi * (p.FirstSolutionValue / nu3 * (a0 + a1 / (nu * nu)) + p.FirstSolutionDerivative / MoreMath.Pow(nu3, 5) * b0);
                double J3 = phi * (p.FirstSolutionValue / nu3 * (a0 + a1 / (nu * nu)) + p.FirstSolutionDerivative / MoreMath.Pow(nu3, 5) * (b0 + b1 / (nu * nu)));
                Console.WriteLine("J0={0:R} J1={1:R} J2={2:R} J3={3:R}", J0, J1, J2, J3);
                Console.WriteLine("J'={0:R}", BesselUae.Evaluate(nu, nu * z));
                //try {
                    Console.WriteLine("old J={0}", AdvancedMath.BesselJ(nu, nu * z));
                //} catch (NonconvergenceException) { }
            }
        }

        public void ComputeBesselExpansionParameter (double z, out double xi, out double phi, out double[] a, out double[] b) {
            a = new double[] { 1.0, 0.0 };
            b = new double[2];
            double e = 1.0 - z;
            if (Math.Abs(e) < 1.0 / 32.0) {
                double[] xi_c = new double[] { 0.0, 1.0, 3.0 / 10.0, 32.0 / 175.0, 1037.0 / 7875.0, 103727.0 / 1010635.0, 33060241.0 / 394143750.0 };
                double[] phi_c = new double[] { 1.0, 1.0 / 5.0, 3.0 / 35.0, 73.0 / 1575.0, 35209.0 / 1212750.0, 380069.0 / 18768750.0, 1897703867.0 / 124155281250.0 };
                double[] b0_c = new double[] { 1.0 / 70.0, 2.0 / 225.0, 953.0 / 202125.0, 17942.0 / 7882875.0, 694817.0 / 709458750.0, 16629512.0 / 50253328125.0, 9633362212.0 / 367603095234375.0 };
                double[] a1_c = new double[] { -1.0 / 225.0, -71.0 / 38500.0, 25591.0 / 45045000.0, 2982172.0 / 1773646875.0, 25025359.0 / 13400887500.0, 5379178799.0 / 3334268437500.0, 23287042920833.0 / 18905302040625000.0 };
                double ek = 1.0;
                xi = xi_c[0];
                phi = phi_c[0];
                b[0] = b0_c[0];
                a[1] = a1_c[0];
                for (int i = 1; i < xi_c.Length; i++) {
                    ek *= e;
                    xi += xi_c[i] * ek;
                    phi += phi_c[i] * ek;
                    b[0] += b0_c[i] * ek;
                    a[1] += a1_c[i] * ek;
                }
                xi = Math.Pow(2.0, 1.0 / 3.0) * xi;
                phi = Math.Pow(2.0, 1.0 / 3.0) * phi;
                b[0] = Math.Pow(2.0, 1.0 / 3.0) * b[0];
                b[1] = Math.Pow(2.0, 1.0 / 3.0) * 1213.0 / 1023750.0;
            } else {
                if (z < 1.0) {
                    double w = Math.Sqrt((1.0 - z) * (1.0 + z));
                    xi = Math.Pow(3.0 / 2.0 * (Math.Log((1.0 + w) / z) - w), 2.0 / 3.0);
                } else {
                    double w = Math.Sqrt((z - 1.0) * (z + 1.0));
                    xi = -Math.Pow(3.0 / 2.0 * (w - Math.Acos(1.0 / z)), 2.0 / 3.0);
                }
                phi = Math.Pow(4.0 * xi / (1.0 - z) / (1.0 + z), 1.0 / 4.0);
                b[0] = -5.0 / 48.0 / (xi * xi) + phi * phi / 48.0 / xi * (5.0 / (1.0 - z * z) - 3.0);
            }
        }

        [TestMethod]
        public void BesselUAE () {

            double nu = 4000.0;
            double x = 4000.0;

            double z = x / nu;
            double xi, phi;
            double[] a, b;
            ComputeBesselExpansionParameter(z, out xi, out phi, out a, out b);
            Console.WriteLine("z={0} xi={1} phi={2}", z, xi, phi);
            Console.WriteLine("a0={0} b0={1} a1={2} b1={3}", a[0], b[0], a[1], b[1]);

            double y = Math.Pow(nu, 2.0 / 3.0) * xi;
            SolutionPair p = AdvancedMath.Airy(y);
            Console.WriteLine("y={0} A={1} A'={2} B={3} B'={4}", y, p.FirstSolutionValue, p.FirstSolutionDerivative, p.SecondSolutionValue, p.SecondSolutionDerivative);

            double J0 = phi * (a[0] * p.FirstSolutionValue / Math.Pow(nu, 1.0 / 3.0));
            double J1 = phi * (a[0] * p.FirstSolutionValue / Math.Pow(nu, 1.0 / 3.0) + b[0] * p.FirstSolutionDerivative / Math.Pow(nu, 5.0 / 3.0));
            double J2 = phi * (p.FirstSolutionValue * (a[0] / Math.Pow(nu, 1.0 / 3.0) + a[1] / Math.Pow(nu, 7.0 / 3.0)) + p.FirstSolutionDerivative * b[0] / Math.Pow(nu, 5.0 / 3.0));
            double J3 = phi * (p.FirstSolutionValue * (a[0] / Math.Pow(nu, 1.0 / 3.0) + a[1] / Math.Pow(nu, 7.0 / 3.0)) + p.FirstSolutionDerivative * (b[0] / Math.Pow(nu, 5.0 / 3.0) + b[1] / Math.Pow(nu, 11.0 / 3.0)));
            double Y0 = -phi * (a[0] * p.SecondSolutionValue / Math.Pow(nu, 1.0 / 30));
            double Y1 = -phi * (a[0] * p.SecondSolutionValue / Math.Pow(nu, 1.0 / 3.0) + b[0] * p.SecondSolutionDerivative / Math.Pow(nu, 5.0 / 3.0));
            double Y2 = -phi * (p.SecondSolutionValue * (a[0] / Math.Pow(nu, 1.0 / 3.0) + a[1] / Math.Pow(nu, 7.0 / 3.0)) + p.SecondSolutionDerivative * b[0] / Math.Pow(nu, 5.0 / 3.0));
            double Y3 = -phi * (p.SecondSolutionValue * (a[0] / Math.Pow(nu, 1.0 / 3.0) + a[1] / Math.Pow(nu, 7.0 / 3.0)) + p.SecondSolutionDerivative * (b[0] / Math.Pow(nu, 5.0 / 3.0) + b[1] / Math.Pow(nu, 11.0 / 3.0)));

            Console.WriteLine("{0} {1} {2} {3}", J0, J1, J2, J3);
            Console.WriteLine("{0} {1} {2} {3}", Y0, Y1, Y2, Y3);


        }

        [TestMethod]
        public void BesselFix () {

            SolutionPair p = AdvancedMath.Airy(-3.0 / 2.0);
            Console.WriteLine("{0:R} {1:R} {2:R} {3:R}", p.FirstSolutionValue, p.FirstSolutionDerivative, p.SecondSolutionValue, p.SecondSolutionDerivative);
            Console.WriteLine(p.FirstSolutionValue * p.SecondSolutionDerivative - p.SecondSolutionValue * p.FirstSolutionDerivative);

        }

        [TestMethod]
        public void BesselYSeries () {

            Console.WriteLine("0: {0}", AdvancedMath.BesselY(16, 12.0));
            Console.WriteLine("1: {0}", BesselY_Series(16, 12.0) - BesselY_Series2(16, 12.0));
            //Console.WriteLine("2: {0}", BesselY_Series2(16, 12.0));

        }


        private static double BesselY (int n, double x) {

            double x2 = x / 2.0;
            double z = x2 * x2;

            double dy = -1.0 / (MoreMath.Pow(x2, n) / AdvancedIntegerMath.Factorial(n - 1)) / Math.PI;
            double y = dy;

            for (int k = 1; k < n; k++) {
                double y_old = y;
                dy *= z / ((n - k) * k);
                y += dy;
                Console.WriteLine("{0} {1} {2}", k, dy, y);
                if (y == y_old) return (y);
            }

            z = - z;

            double a = 2.0 / Math.PI * MoreMath.Pow(x2, n) / AdvancedIntegerMath.Factorial(n);
            double b = Math.Log(x2) + AdvancedMath.EulerGamma - AdvancedIntegerMath.HarmonicNumber(n) / 2.0;

            y += a * b;

            for (int k = 1; k < 128; k++) {

                double y_old = y;

                a *= z / (k * (n + k));
                b -= (1.0 / k + 1.0 / (n + k)) / 2.0;

                y += a * b;

                Console.WriteLine("k={0} a={1} b={2} y={3}", k, a, b, y);

                if (y == y_old) return (y);

            }

            throw new NonconvergenceException(); 

        }

        private static double BesselY_Series2 (int n, double x) {

            double z = x * x / 4.0;

            if (n == 0) return (0);

            double a = 1.0 / Math.PI * MoreMath.Pow(2.0 / x, n) * AdvancedIntegerMath.Factorial(n - 1);

            double y = a;
            Console.WriteLine("{0} {1} {2}", 0, a, y);

            for (int k = 1; k < n; k++) {
                double y_old = y;
                a *= z / ((n - k) * k);
                y += a;
                Console.WriteLine("{0} {1} {2}", k, a, y);
                if (y == y_old) return (y);
            }

            return (y);
        }

        // Start from

        // Y_n(x) = \frac{2}{\pi} \log\left(\frac{x}{2}\right) J_n(x) - \frac{1}{\pi} \left(

        // and insert

        //   \psi(k) = -\gamma + \sum_{j=1}^{k-1} \frac{1}{j} 

        private static double BesselY_Series (int n, double x) {

            double z = - x * x / 4.0;

            double a = 2.0 / Math.PI * MoreMath.Pow(x / 2.0, n) / AdvancedIntegerMath.Factorial(n);
            double b = Math.Log(x / 2.0) + AdvancedMath.EulerGamma - AdvancedIntegerMath.HarmonicNumber(n) / 2.0;

            double y = a * b;

            Console.WriteLine("k={0} a={1} b={2} y={3}", 0, a, b, y);


            for (int k = 1; k < 128; k++) {

                double y_old = y;

                a *= z / (k * (n + k));
                b -= (1.0 / k + 1.0 / (n + k)) / 2.0;

                y += a * b;

                Console.WriteLine("k={0} a={1} b={2} y={3}", k, a, b, y);

                if (y == y_old) return (y);

            }

            throw new NonconvergenceException();
        }


        [TestMethod]
        public void CumulantTest () {

            Distribution d = new GammaDistribution(2.0);

            double[] K = new double[8];
            K[0] = 0.0;
            for (int r = 1; r < K.Length; r++) {
                K[r] = 2.0 * AdvancedIntegerMath.Factorial(r - 1);
            }

            double[] M = MomentMath.CumulantToRaw(K);
            for (int r = 0; r < M.Length; r++) {
                Console.WriteLine("{0} {1}", M[r], d.Moment(r));
            }

            Console.WriteLine("---");

            double[] C = MomentMath.CumulantToCentral(K);
            for (int r = 0; r < C.Length; r++) {
                Console.WriteLine("{0} {1}", C[r], d.MomentAboutMean(r));
            }

            Console.WriteLine("---");

            double[] KP = MomentMath.RawToCumulant(M);
            for (int r = 0; r < KP.Length; r++) {
                Console.WriteLine("{0} {1}", KP[r], K[r]);
            }

        }

        [TestMethod]
        public void TwoSampleKS () {

            int n = 2 * 3 * 3 * 2; int m = 2 * 2 * 3;
            int lcm = (int) AdvancedIntegerMath.LCM(n, m);
            int gcf = (int) AdvancedIntegerMath.GCF(n, m);
            Console.WriteLine("{0} {1} {2} {3}", n, m, lcm, gcf);
            KolmogorovTwoSampleExactDistribution d0 = new KolmogorovTwoSampleExactDistribution(n, m);
            
            for (int k = d0.Minimum; k <= d0.Maximum; k++) {
                Console.WriteLine("{0} {1} {2}", k, d0.LatticePathSum(k * gcf), d0.LatticePathSum(k * gcf - 1));
            }

            //for (int c = 1; c <= n; c++) {
            //    Console.WriteLine("{0} {1}", c, EqualKS(n, c));
            //}
        }

        [TestMethod]
        public void ChiSquareDistribution () {

            int B = 8;
            Random rng = new Random(0);
            //BinomialDistribution d = new BinomialDistribution(0.25, 6);
            //DiscreteUniformDistribution d = new DiscreteUniformDistribution(0, 4);
            //BernoulliDistribution d = new BernoulliDistribution(0.25);
            DiscreteDistribution d = new PoissonDistribution(2.0);

            Sample s = new Sample();
            ChiSquaredDistribution nullDistribution = null;
            for (int i = 0; i < 512; i++) {

                Histogram h = new Histogram(B);
                for (int j = 0; j < 1024; j++) {
                    int k = d.GetRandomValue(rng);
                    h.Add(k);
                    //if (k < h.Count) h[k].Increment();
                    //h[d.GetRandomValue(rng)].Increment();
                }
                TestResult cr = h.ChiSquaredTest(d);
                nullDistribution = (ChiSquaredDistribution) cr.Distribution;
                //Console.WriteLine(((ChiSquaredDistribution) cr.Distribution).DegreesOfFreedom);
                s.Add(cr.Statistic);

            }

            Console.WriteLine(nullDistribution.DegreesOfFreedom);
            TestResult kr = s.KuiperTest(nullDistribution);
            Console.WriteLine(kr.LeftProbability);

        }

        [TestMethod]
        public void TwoSampleKS2 () {

            int n = 2 * 3 * 3; int m = 2 * 2 * 3;
            Random rng = new Random(0);
            NormalDistribution d = new NormalDistribution();

            Histogram h = new Histogram((int) AdvancedIntegerMath.LCM(n, m) + 1);

            //int[] h = new int[AdvancedIntegerMath.LCM(n, m) + 1];

            int count = 1000;
            for (int i = 0; i < count; i++) {

                Sample A = new Sample();
                for (int j = 0; j < n; j++) A.Add(d.GetRandomValue(rng));
                Sample B = new Sample();
                for (int j = 0; j < m; j++) B.Add(d.GetRandomValue(rng));

                TestResult r = Sample.KolmogorovSmirnovTest(A, B);
                int k = (int) Math.Round(r.Statistic * AdvancedIntegerMath.LCM(n, m));
                //Console.WriteLine("{0} {1}", r.Statistic, k);
                h[k].Increment();
                //h[k] = h[k] + 1;

            }

            KolmogorovTwoSampleExactDistribution ks = new KolmogorovTwoSampleExactDistribution(n, m);
            double chi2 = 0.0; int dof = 0;
            for (int i = 0; i < h.Count; i++) {
                double ne = ks.ProbabilityMass(i) * count;
                Console.WriteLine("{0} {1} {2}", i, h[i].Counts, ne);
                if (ne > 4) {
                    chi2 += MoreMath.Sqr(h[i].Counts - ne) / ne;
                    dof++;
                }
            }
            Console.WriteLine("chi^2={0} dof={1}", chi2, dof);

            TestResult r2 = h.ChiSquaredTest(ks);
            ChiSquaredDistribution rd = (ChiSquaredDistribution) r2.Distribution;
            Console.WriteLine("chi^2={0} dof={1} P={2}", r2.Statistic, rd.DegreesOfFreedom, r2.RightProbability);
        }


        private static double EqualKS (int n, int c) {

            int sign = -1;
            double sum = 0.0;
            for (int m = n - c; m >= 0; m -= c) {
                sign = -sign;
                sum += sign * AdvancedIntegerMath.BinomialCoefficient(2 * n, m);
            }
            return (2.0 * sum);
        }

        [TestMethod]
        public void Test2 () {
            double value = -1.0;
            bool sign; long mantissa; int exponent;
            Frexp(value, out sign, out mantissa, out exponent);
            Console.WriteLine("{0} = {1} {2} X 2^({3})", value, sign, mantissa, exponent); 
        }

        private static void Frexp (double value, out bool sign, out long mantissa, out int exponent) {
            long bits = BitConverter.DoubleToInt64Bits(value);
            sign = (bits < 0);
            exponent = (int) ((bits >> 52) & 0x7ffL);
            mantissa = bits & 0xfffffffffffffL;

            mantissa = mantissa | (1L << 52);
            exponent -= 1075;
        }

        [TestMethod]
        public void Test12 () {

            double R = 5734161139222659 * Math.Pow(2.0, -54.0);
            double C1 = 7074237752028440 * Math.Pow(2.0, -51.0);
            double C2 = 4967757600021504 * Math.Pow(2.0, -105.0);
            double C3 = 7744522442262976 * Math.Pow(2.0, -155.0);

            double x = 1000000000.0;
            Console.WriteLine(x);
            Console.WriteLine("S {0}", Math.Sin(x));

            double j = Math.Round(x * R);

            x = -j * C1 + x;
            x = -j * C2 + x;
            x = -j * C3 + x;

            Console.WriteLine(x);
            Console.WriteLine("S {0}", Math.Sin(x));

        }

        [TestMethod]
        public void ET () {

            SquareMatrix A = new SquareMatrix(4);
            A[0, 0] = 1.0; A[0, 1] = 1.0; A[0, 2] = 1.0; A[0, 3] = 1.0;
            A[1, 0] = 1.0; A[1, 1] = 1.0; A[1, 2] = 1.0; A[1, 3] = 1.0;
            A[2, 0] = 0.0; A[2, 1] = 0.0; A[2, 2] = 1.0; A[2, 3] = 1.0;
            A[3, 0] = 0.0; A[3, 1] = 0.0; A[3, 2] = 1.0; A[3, 3] = 1.0;
            A.Eigensystem();
        }

        // We want to find the rotation that brings a 2 X 2 matrix into triangular form.
        //   (  c  s ) ( a11  a12 ) ( c  -s )
        //   ( -s  c ) ( a21  a22 ) ( s   c )
        // Multiplying out gives
        //   ( c^2 a11 + cs a21 + cs a12 + s^2 a22  c^2 a12 - cs a11 + cs a22 - s^2 a21 )
        //   ( c^2 a21 - cs a11 + cs a22 - s^2 a12  c^2 a22 - cs a21 - cs a12 + s^2 a11 )

        // We want a21' = 0. Since the rotation must vanish when a21 = 0 and grow as a21 grows, make the Ansatz that s ~ a21, i.e.
        //   s = \frac{a_{21}}{\sqrt{a_{21}^2 + b^2}}  c = \frac{b}{\sqrt{a_{21}^2 + b^2}}
        // It is then straightforward to derive
        //   b = \frac{a_{11} - a_{22} \pm q}{2}
        // where q^2 = ( a_{11} - a_{22} )^2 + 4 a12 a21 is the same descriminant as appears in the eigenvalue problem, so this rotation
        // exists iff the eigenvalues are real.

        private void TwoByTwoSchur (ref double a11, ref double a12, ref double a21, ref double a22, out double s, out double c) {

            // compute some quantities we will use 
            double u = a11 + a22;
            double v = a11 - a22;
            double w = a11 * a22 - a12 * a21;
            double q2 = v * v + 4.0 * a12 * a21;
            // note u is the trace and w is the determinant

            if (q2 >= 0.0) {

                // the descriminant is positive so the eigenvalues are real
                double q = Math.Sqrt(q2);
 
                // find the rotation sets a21' = 0
                // in the equation for b, choose the sign so as to minimize cancelation
                double b = (v >= 0.0) ? (v + q) / 2.0 : (v - q) / 2.0;
                double rho = MoreMath.Hypot(a21, b);
                s = a21 / rho;
                c = b / rho;

                // Note that a12' - a21' = a12 - a21, and since a21' = 0, a12' = a12 - a21
                a12 = a12 - a21;
                a21 = 0.0;

                // the eigenvalues are (u \pm q) / 2
                // we avoid cancelation by computing the non-canceling one first and
                // computing the other using the fact that their product equals the determinant
                if (u >= 0.0) {
                    a11 = (u + q) / 2.0;
                    a22 = w / a11;
                } else {
                    a22 = (u - q) / 2.0;
                    a11 = w / a22;
                }
                // we have placed the u + q eigenvalue in the a11 slot. if v >= 0, that is where the rotation puts it
                // if v < 0, the u - q eigenvalue belongs these, so we need to swap them
                if (v < 0) { double t = a11; a11 = a22; a22 = t; }

            } else {

                // In the q2 < 0 case, we can't zero a21. But we can rotate so a11' = a22'

                double r = - (a12 + a21) / v;
                double t = Math.Sign(r) / (Math.Abs(r) + MoreMath.Hypot(1.0, r));
                c = 1.0 / MoreMath.Hypot(1.0, t);
                s = t * c;

                // Since rotations preserve the trace, the equal diagonal elements must equal the average of the previous elements 
                a11 = u / 2.0;
                a22 = a11;

                //double q = Math.Sqrt(-q2);
                //e1 = new Complex(u / 2.0, q / 2.0);
                //e2 = e1.Conjugate;

                //s = 0.0;
                //c = 1.0;

            }
        }

        [TestMethod]
        public void EigenTest () {
            
            Random rng = new Random(1);
            for (int i = 0; i < 100; i++) {

                double a00 = 1.0 - 2.0 * rng.NextDouble();
                double a01 = 1.0 - 2.0 * rng.NextDouble();
                double a10 = 1.0 - 2.0 * rng.NextDouble();
                double a11 = 1.0 - 2.0 * rng.NextDouble();

                SquareMatrix A = new SquareMatrix(2);
                A[0, 0] = a00;
                A[0, 1] = a01;
                A[1, 0] = a10;
                A[1, 1] = a11;

                double s, c;
                TwoByTwoSchur(ref a00, ref a01, ref a10, ref a11, out s, out c);

                if (s > 1.0) continue;

                SquareMatrix T = new SquareMatrix(2);
                T[0, 0] = c;
                T[0, 1] = s;
                T[1, 0] = -s;
                T[1, 1] = c;

                SquareMatrix S = T * A * T.Transpose();

                Console.WriteLine("{0} {1}", a00, S[0,0]);

            }
            
            //A[0, 0] = 1.0; A[0, 1] = 4.0;
            //A[1, 0] = 2.0; A[1, 1] = 3.0;
            //A.Eigenvalues();
        }

        /* ORDER STATISTICS */

        // 24-point Gauss-Hermite integration
        // We chose this because the smallest weight is 10^{-16}

        private static double[] xs = new double[] {
            0.22441454747251558515,
            0.67417110703721223600,
            1.1267608176112450721,
            1.5842500109616941485,
            2.0490035736616989118,
            2.5238810170114269742,
            3.0125461375655648257,
            3.5200068130345247113,
            4.0536644024481495039,
            4.6256627564237872650,
            5.2593829276680443674,
            6.0159255614257397173
        };

        private static double[] ws = new double[] {
            0.42693116386869924965,
            0.28617953534644301790,
            0.12773962178455916065,
            3.744547050323074601E-2,
            7.04835581007267097E-3,
            8.23692482688417458E-4,
            5.68869163640437977E-5,
            2.15824570490233363E-6,
            4.01897117494142968E-8,
            3.04625426998756390E-10,
            6.58462024307817006E-13,
            1.66436849648910887E-16
        };

        private static readonly double SqrtTwo = Math.Sqrt(2.0);
        private static readonly double SqrtPI = Math.Sqrt(Math.PI);

        private static double GaussHermiteIntegrate (Func<double, double> f) {
            double y = 0.0;
            for (int i = 0; i < xs.Length; i++) {
                double x = SqrtTwo * xs[i];
                y += (ws[i] / SqrtPI) * (f(x) + f(-x));
            }
            return (y);
        }

        private static double NormalMeanOrderStatisticExpansion (int i, int n) {

            double p = i / (n + 1.0);
            double q = 1.0 - p;

            double FI = Math.Sqrt(2.0) * AdvancedMath.InverseErf(2.0 * p - 1.0);
            double FI2 = 2.0 * Math.PI * Math.Exp(FI * FI) * FI;

            //double FI3 = Math.Pow(2.0 * Math.PI * Math.Exp(FI * FI), 3.0 / 2.0) * (1.0 + 2.0 * FI * FI);
            //double FI4 = Math.Pow(2.0 * Math.PI * Math.Exp(FI * FI), 2.0) * FI * (7.0 + 6.0 * FI * FI);

            return (FI + p * q / 2.0 * FI2 / (n + 2));

            //return (FI + p * q / 2.0 * FI2 / (n + 2) + p * q / (n + 2) / (n + 2) * ((q - p) / 3.0 * FI3 + p * q / 8.0 * FI4));

        }

        private static double NormalMeanOrderStatisticExpansion2 (int i, int n) {

            double p = i / (n + 1.0);
            double q = 1.0 - p;

            double FI = Math.Sqrt(2.0) * AdvancedMath.InverseErf(2.0 * p - 1.0);
            double FI2 = 2.0 * Math.PI * Math.Exp(FI * FI) * FI;

            double FI3 = Math.Pow(2.0 * Math.PI * Math.Exp(FI * FI), 3.0 / 2.0) * (1.0 + 2.0 * FI * FI);
            double FI4 = Math.Pow(2.0 * Math.PI * Math.Exp(FI * FI), 2.0) * FI * (7.0 + 6.0 * FI * FI);

            return (FI + p * q / 2.0 * FI2 / (n + 2) + p * q / (n + 2) / (n + 2) * ((q - p) / 3.0 * FI3 + p * q / 8.0 * FI4));


        }

        [TestMethod]
        public void TestNormalOrderStatistic () {

            int n = 100;
            //int r = 3 * n / 4;
            int r = 52;
            Distribution d = new NormalDistribution();

            double C = Math.Exp(AdvancedIntegerMath.LogFactorial(n) - AdvancedIntegerMath.LogFactorial(r - 1) - AdvancedIntegerMath.LogFactorial(n - r));

            double m = GaussHermiteIntegrate(x => C * MoreMath.Pow(d.LeftProbability(x), r - 1) * MoreMath.Pow(d.RightProbability(x), n - r) * x);
            //double m = GaussHermiteIntegrate(x => 1.0);

            double m2 = FunctionMath.Integrate(
                //x => 1.0 * Math.Exp(-x * x / 2.0) / Math.Sqrt(2.0 * Math.PI),
                x => C * MoreMath.Pow(d.LeftProbability(x), r - 1) * MoreMath.Pow(d.RightProbability(x), n - r) * x * Math.Exp(-x * x / 2.0) / Math.Sqrt(2.0 * Math.PI),
                Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity)
            );

            Console.WriteLine(m);
            Console.WriteLine(m2);
            Console.WriteLine(NormalMeanOrderStatisticExpansion(r, n));
            Console.WriteLine(NormalMeanOrderStatisticExpansion2(r, n));
            //Console.WriteLine(1.5 / Math.Sqrt(Math.PI));
            

        }

        /* INVERSE KS CDF */

        public static double InverseQ (double Q) {

            if (Q == 0.0) return (Double.PositiveInfinity);

            // This is an inversion technique suggested in Numerical Recipies, 3rd Edition, Section 6.14.12
            // Write the Q expansion as
            //   Q = 2 (y - y^4 + y^9 - y^{16} + y^{25} + \cdots)
            // where y = e^{-2 x^2}.  Rewrite this as
            //   y = \frac{Q}{2} + y^4 - y^9 + y^{16} - y^{25} + \cdots
            // and use this iteratively to generate y starting from y = \frac{Q}{2}. When you have y,
            //   x = \sqrt{-\log(x) / 2}
            // The next term, y^36, is below floating point accuracy for all y <~ 1/3, corresponding to Q <~ 2/3.
            // For Q <~ 1/2, 12 or fewer iterations are required for y to converge.

            double halfQ = Q / 2.0;
            double y = halfQ;
            for (int i = 0; i < 256; i++) {
                double y_old = y;
                double y4, y9, y16, y25;
                ComputePowerSet(y, out y4, out y9, out y16, out y25);
                y = halfQ + y4 - y9 + y16 - y25;
                if (y == y_old) {
                    return (Math.Sqrt(-Math.Log(y) / 2.0));
                }
            }

            throw new NonconvergenceException();

        }

        // This gives all the powers of y required in the Q series with the minimum number of operations.

        private static void ComputePowerSet (double y, out double y4, out double y9, out double y16, out double y25) {
            double y2, y8;
            y2 = y * y;
            y4 = y2 * y2;
            y8 = y4 * y4;
            y9 = y8 * y;
            y16 = y8 * y8;
            y25 = y16 * y9;
        }

        [TestMethod]
        public void TestInverseKS () {

            KolmogorovDistribution d = new KolmogorovDistribution();
            
            double x;
            
            Stopwatch s = Stopwatch.StartNew();
            for (double Q = 0.00001; Q <= 0.5; Q += 0.00001) {
                x = InverseQ(Q);
            }
            s.Stop();
            Console.WriteLine(s.ElapsedMilliseconds);

            s.Restart();
            for (double Q = 0.00001; Q <= 0.5; Q += 0.00001) {
                x = d.InverseRightProbability(Q);
            }
            s.Stop();
            Console.WriteLine(s.ElapsedMilliseconds);
            

            Console.WriteLine(InverseQ(0.5));
            Console.WriteLine(d.InverseRightProbability(0.5));

        }

        /* MINIMIZATION */

        private static void ParabolicFit (
            double x1, double f1, double x2, double f2, double x3, double f3,
            out double x0, out double fpp
        ) {

            // We want to find the parabola
            //   f = f0 + (f'' / 2) (x - x0)^2
            // that passes throught the points (x1,f1), (x2,f2), (x3,f3)

            // The solution, after much algebra is:
            //   x0 = x1 - \frac{1}{2} \frac{(f3-f1)(x1-x2)^2 - (f2-f1)(x3-x1)^2}{(f3-f1)(x1-x2) + (f2-f1)(x3-x1)}
            //   f0 = ?
            //   f'' = 2 \frac{(x1-x2)(f3-f1) + (x3-x1)(f2-f1)}{(x1-x2)(x3-x1)(x3-x2)}

            // compute the differences that appear in our expressions
            double x12 = x1 - x2; double f21 = f2 - f1;
            double x31 = x3 - x1; double f31 = f3 - f1;
            double x32 = x3 - x2;

            // compute the numerator and denominator in the expression for x0
            double t1 = f31 * x12;
            double t2 = f21 * x31;
            double p = t1 * x12 - t2 * x31;
            double q = 2.0 * (t1 + t2);

            if (q == 0.0) {
                // the denominator vanishes when the points are colinear; there is no parabolic fit
                x0 = Double.NaN;
                fpp = Double.NaN;
            } else {
                x0 = x1 - p / q;
                fpp = q / (x12 * x31 * x32);
            }

        }

        private static void CubicHermiteMinimum (double x0, double f0, double m0, double x1, double f1, double m1, out double x, out double fpp) {
            double dx = x1 - x0;
            double t;  CubicHermiteMinimum(f0, m0 * dx, f1, m1 * dx, out t, out fpp);
            x = x0 + dx * t;
            fpp = fpp / (dx * dx);
        }

        private static void CubicHermiteMinimum (double f0, double m0, double f1, double m1, out double t, out double fpp) {

            // Given f,f' = f0,m0 at t = 0 and f,f' = f1,m1 at t = 1, the cubic Hermite interpolating polynomial
            // (http://en.wikipedia.org/wiki/Cubic_Hermite_spline) is
            //   f = (2 t^3 - 3 t^2 + 1) f0 + (t^3 - 2 t^2 + t) m0 + (-2 t^3 + 3 t^2) f1 + (t^3 - t^2) m1
            //     = (1 + 2 t) (1 - t)^2 f0 + t (1 - t)^2 m0 + t^2 (3 - 2 t) f1 + t^2 (t - 1) m1
            //     = [2 (f0 - f1) + (m0 + m1)] t^3 - [3 (f0 - f1) + 2 m0 + m1] t^2 + m0 t + p0

            // The derivative of this function is
            //   f' = 3 [2 (f0 - f1) + (m0 + m1)] t^2 - 2 [3 (f0 - f1) + 2 m0 + m1] t + m0 = a t^2 - b t + c
            // Set this equal to zero to get a quadratic equation a t^2 - b t + c = 0 for minimum (and maximum).
            // The descriminant q^2 = b^2 - 4 a c = ?
            // If a > 0 the cubic is upward sloping, so minimum is rightmost solution. If a < 0, minimum is leftmost solution.

            // The second derivative is f'' = 2 a t - b.

            double df = f0 - f1;
            double sm = m0 + m1;
            double q2 = 4.0 * MoreMath.Pow(3.0 * df + sm, 2) - 4.0 * m0 * m1;

            // If q^2 < 0, there is no minimum

            if (q2 < 0.0) {
                t = Double.NaN;
                fpp = Double.NaN;
                return;
            }

            double q = Math.Sqrt(q2);

            double a = 3.0 * (2.0 * df + sm);
            double b = 2.0 * (3.0 * df + 2.0 * m0 + m1);
            double c = m0;

            Console.WriteLine("a={0} b={1} c={2}", a, b, c);

            // If a is very small or zero, our cubic becomes a parabola.
            // This happens, for example, given two points with equal function values and opposite slopes.
            if (Math.Abs(a) <= 1.0E-7 * Math.Abs(b) && Math.Abs(a * c) <= 1.0E-14 * b * b) {
                Console.WriteLine("parabola");
                if (b < 0.0) {
                    t = c / b;
                    fpp = -b;
                } else {
                    t = Double.NaN;
                    fpp = Double.NaN;
                }
                return;
            }

            double t1, t2;
            if (b > 0) {
                t1 = (b + q) / (2.0 * a);
                t2 = (2.0 * c) / (b + q);
            } else {
                t1 = (b - q) / (2.0 * a);
                t2 = (2.0 * c) / (b - q);
            }
            if (a > 0) {
                t = Math.Max(t1, t2);
            } else {
                t = Math.Min(t1, t2);
            }

            fpp = 2.0 * a * t - b;


            /*
            if (a > 0) {
                // rightmost solution is minimum
                // use expression with no cancelation
                if (b > 0) {
                    return ((b + q) / (2.0 * a));
                } else {
                    return ((2.0 * c) / (b - q));
                }
            } else {
                // leftmost solution is minimum
                // again, use expression with no cancelation
                if (b > 0) {
                    return ((2.0 * c) / (b + q));
                } else {
                    return ((b - q) / (2.0 * a));
                }
            }
            */
        }

        [TestMethod]
        public void TestCubicHermite () {

            double x, fpp; CubicHermiteMinimum(3.1428858, -0.99999916, 0.0012931460, 3.2067118, -0.997880, 0.0650731, out x, out fpp);

            //double x = CubicHermiteMinimum(1.0, -1.0, 1.0, +1.0);
            Console.WriteLine(x);

        }

        [TestMethod]
        public void TestParabolicFit () {

            double x0, fpp;
            ParabolicFit(2, 2, 1, 5, 0, 10, out x0, out fpp);

            Console.WriteLine(x0);
            Console.WriteLine(fpp);

        }

        private static double FindMinimum (
            Func<double, double> f,
            double a, double b
        ) {

            // evaluate three points within the bracket
            double u = (3.0 * a + b) / 4.0;
            double v = (a + b) / 2.0;
            double w = (a + 3.0 * b) / 4.0;

            double fu = f(u); double fv = f(v); double fw = f(w);

            Console.WriteLine("f({0})={1}  f({2})={3}  f({4})={5}", u, fu, v, fv, w, fw);

            // move in the bracket boundaries, if possible
            if (fv < fu) { a = u; if (fw < fv) a = v; }
            if (fv < fw) { b = w; if (fu < fv) b = v; }

            Console.WriteLine("a={0} b={1}", a, b);

            // sort u, v, w by fu, fv, fw values
            // these three comparisons are the most efficient three-item sort
            if (fv < fu) { double t = v; v = u; u = t; t = fv; fv = fu; fu = t; }
            if (fw < fu) { double t = w; w = u; u = t; t = fw; fw = fu; fu = t; }
            if (fw < fv) { double t = w; w = v; v = t; t = fw; fw = fv; fv = t; }

            // An evaluation budget of 32 is sufficient for all our test cases except for |x|, which requires 82 (!) evaluations to converge. Parabolic fitting just does a very poor job
            // for this function (at all scales, since it is scale invariant). We should look into cubic fitting.

            EvaluationSettings settings = new EvaluationSettings() { EvaluationBudget = 128, AbsolutePrecision = 0.0, RelativePrecision = 0.0 };
            return (FindMinimum(f, a, b, u, fu, v, fv, w, fw, settings, 3));

        }

        private static double FindMinimum (
            Func<double,double> f,
            double a, double b,
            double u, double fu, double v, double fv, double w, double fw,
            EvaluationSettings settings, int count
        ) {

            double tol = 0.0;

            while (count < settings.EvaluationBudget) {

                Console.WriteLine("n={0} tol={1}", count, tol);
                Console.WriteLine("[{0}  f({1})={2}  f({3})={4}  f({5})={6}  {7}]", a, u, fu, v, fv, w, fw, b);

                Debug.Assert(a < b);
                Debug.Assert((a <= u) && (u <= b));
                Debug.Assert((fu <= fv) && (fv <= fw));

                // Expected final situation is a<tol><tol>u<tol><tol>b, leaving no point left to evaluate that is not within tol of an existing point.

                if ((b - a) <= 4.0 * tol) return (u);

                // While a < u < b is guaranteed, a < v, w < b is not guaranteed, since the bracket can sometimes be made tight enough to exclude v or w.
                // For example, if u < v < w, then we can set b = v, placing w outside the bracket.

                double x, fpp;
                ParabolicFit(u, fu, v, fv, w, fw, out x, out fpp);
                Console.WriteLine("parabolic x={0} f''={1}", x, fpp);
 
                if (Double.IsNaN(fpp) || (fpp <= 0.0) || (x < a) || (x > b)) {

                    // the parabolic fit didn't work out, so do a golden section reduction instead

                    // to get the most reduction of the bracket, pick the larger of au and ub
                    // for self-similarity, pick a point inside it that divides it into two segments in the golden section ratio,
                    // i.e. 0.3820 = \frac{1}{\phi + 1} and 0.6180 = \frac{\phi}{\phi+1}
                    // put the smaller segment closer to u so that x is closer to u, the best minimum so far

                    double au = u - a;
                    double ub = b - u;

                    if (au > ub) {
                        x = u - au / (AdvancedMath.GoldenRatio + 1.0);
                    } else {
                        x = u + ub / (AdvancedMath.GoldenRatio + 1.0);
                    }

                    Console.WriteLine("golden section x={0}", x);

                }

                // ensure we don't evaluate within tolerance of an existing point
                if (Math.Abs(x - u) < tol) { Console.WriteLine("shift from u (x={0})", x); x = (x > u) ? u + tol : u - tol; }
                if ((x - a) < tol) { Console.WriteLine("shift from a (x={0})", x); x = a + tol; }
                if ((b - x) < tol) { Console.WriteLine("shift from b (x={0})", x); x = b - tol; }

                count++;
                double fx = f(x);
                Console.WriteLine("f({0}) = {1}", x, fx);

                Console.WriteLine("delta={0}", fu - fx);

                if (fx < fu) {

                    // the new point is lower than all the others; this is success 

                    // u now becomes a bracket point
                    if (u < x) {
                        a = u;
                    } else {
                        b = u;
                    }

                    // x -> u -> v -> w
                    w = v; fw = fv;
                    v = u; fv = fu;
                    u = x; fu = fx;

                } else {

                    // x now becomes a bracket point
                    if (x < u) {
                        a = x;
                    } else {
                        b = x;
                    }

                    if (fx < fv) {

                        // the new point is higher than u, but still lower than v and w
                        // this isn't what we expected, but we have lower points that before

                        // x -> v -> w
                        w = v; fw = fv;
                        v = x; fv = fx;

                    } else if (fx < fw) {

                        // x -> w
                        w = x; fw = fx;

                    } else {

                        // the new point is higher than all our other points; this is the worst case
                        // we might still want to replace w with x because
                        // (i) otherwise a parabolic fit will reproduce the same x and
                        // (ii) w is quite likely far outside the new bracket and not telling us much about the behavior near u

                        Console.WriteLine("bad point");
                        //throw new NotImplementedException();
                    }

                }

                // if the user has specified a tollerance, use it
                if ((settings.RelativePrecision > 0.0 || settings.AbsolutePrecision > 0.0)) {
                    tol = Math.Max(Math.Abs(u) * settings.RelativePrecision, settings.AbsolutePrecision);
                } else {
                    // otherwise, try to get the tollerance from the curvature
                    if (fpp > 0.0) {
                        tol = Math.Sqrt(2.0 * (Math.Abs(fu) * 1.0E-14 + 1.0E-28) / fpp);
                    } else {
                        // but if we don't have a useable curvature either, wing it
                    }
                }

            }

            throw new NonconvergenceException();

        }

        public delegate void FuncWithDerivative (double x, out double f, out double fp);

        private double FindMinimumWithDerivative (
            FuncWithDerivative f,
            double a, double b
        ) {

            // pick two points in the interval
            double u = 2.0 / 3.0 * a + 1.0 / 3.0 * b;
            double v = 1.0 / 3.0 * a + 2.0 / 3.0 * b;

            // evalue the function there
            double fu, fpu, fv, fpv;
            f(u, out fu, out fpu);
            f(v, out fv, out fpv);

            // move in the bound at the higher side
            if (fu > fv) {
                a = u;
            } else {
                b = v;
            }

            // if f(v) < f(u), swap the points to ensure that u and v are ordered as required
            if (fv < fu) {
                double t;
                t = u; u = v; v = t;
                t = fu; fu = fv; fv = t;
                t = fpu; fpu = fpv; fpv = t;
            }

            // An evaluation budget of 32 is sufficient for all our test cases except for |x|, which requires 82 (!) evaluations to converge. Parabolic fitting just does a very poor job
            // for this function (at all scales, since it is scale invariant). We should look into cubic fitting.

            return (FindMinimumWithDerivative(f, a, b, u, fu, fpu, v, fv, fpv, new EvaluationSettings() { EvaluationBudget = 64, AbsolutePrecision = 0.0, RelativePrecision = 0.0 }));
        }

        private double FindMinimumWithDerivative (
            FuncWithDerivative f,
            double a, double b,
            double u, double fu, double fpu,
            double v, double fv, double fpv,
            EvaluationSettings settings
        ) {

            double tol = 0.0;

            int count = 0;
            while (count < settings.EvaluationBudget) {

                Console.WriteLine("n = {0}, tol = {1}", count, tol);
                Console.WriteLine("[{0} f({1})={2}({3}) f({4})={5}({6}) {7}]", a, u, fu, fpu, v, fv, fpv, b);

                // a, b bracket minimum a < u, v < b and f(u), f(v) <= f(a), f(b)
                Debug.Assert(a < b);
                Debug.Assert((a <= u) && (u <= b));
                //Debug.Assert((a <= v) && (v <= b));
                Debug.Assert(fu <= fv);

                if ((b - a) <= 4.0 * tol) return (u);

                // compute the minimum of the interpolating Hermite cubic
                double x, fpp; CubicHermiteMinimum(u, fu, fpu, v, fv, fpv, out x, out fpp);

                Console.WriteLine("cubic x = {0}, fpp = {1}", x, fpp);

                // if the cubic had no minimum, or the minimum lies outside our bounds, fall back to bisection
                if (Double.IsNaN(x) || (x <= a) || (x >= b)) {

                    // the derivative tells us which side to choose
                    if (fpu > 0.0) {
                        x = (a + u) / 2.0;
                    } else {
                        x = (u + b) / 2.0;
                    }

                    Console.WriteLine("bisection x = {0}", x);

                }

                // ensure we don't evaluate within tolerance of an existing point
                if (Math.Abs(x - u) < tol) { Console.WriteLine("shift from u (x={0})", x); x = (x > u) ? u + tol : u - tol; }
                if ((x - a) < tol) { Console.WriteLine("shift from a (x={0})", x); x = a + tol; }
                if ((b - x) < tol) { Console.WriteLine("shift from b (x={0})", x); x = b - tol; }

                // evaluate the function plus derivative at the predicted minimum
                double fx, fpx;
                f(x, out fx, out fpx);
                count++;

                Console.WriteLine("f({0}) = {1}({2})", x, fx, fpx);

                // check if we have converged
                double df = fu - fx;
                Console.WriteLine("df={0}", df);
                if ((Math.Abs(df) < settings.AbsolutePrecision) || (2.0 * Math.Abs(df) < settings.RelativePrecision * (Math.Abs(fu) + Math.Abs(fx)))) {
                    Console.WriteLine("count = {0}", count);
                    return (x);
                }

                if (fx < fu) {

                    // x is the new lowest point: f(x) < f(u) < f(v)
                    // this is the expected outcome

                    // move the bracket
                    if (x < u) {
                        b = u;
                    } else {
                        a = u;
                    }

                    // x -> u -> v
                    v = u; fv = fu; fpv = fpu;
                    u = x; fu = fx; fpu = fpx;

                } else {

                    // move the bracket
                    if (x < u) {
                        a = x;
                    } else {
                        b = x;
                    }

                    if (fx < fv) {

                        // x lies between other two known points: f(u) < f(x) < f(v)

                        // x -> v
                        v = x; fv = fx; fpv = fpx;

                    } else {

                        // x is higher than both other points: f(u) < f(v) < f(x)
                        // this is a really poor outcome; we expected to get a point lower than our other two and we got a point higher than both
                        // next time we should bisect
                        Console.WriteLine("bad point");
                        //throw new NotImplementedException();

                        //v = x; fv = fx; fpv = fpx;

                    }

                }

                // if the user has specified a tollerance, use it
                if ((settings.RelativePrecision > 0.0 || settings.AbsolutePrecision > 0.0)) {
                    tol = Math.Max(Math.Abs(u) * settings.RelativePrecision, settings.AbsolutePrecision);
                } else {
                    // otherwise, try to get the tollerance from the curvature
                    if (fpp > 0.0) {
                        tol = Math.Sqrt(2.0 * 1.0E-14 * (Math.Abs(fu) + 1.0E-14) / fpp);
                    } else {
                        // but if we don't have a useable curvature either, wing it
                        if (tol == 0.0) tol = 1.0E-7;
                    }
                }


            }

            throw new NonconvergenceException();

        }

        public double FindMinimum (Func<double, double> f, double x, double d, EvaluationSettings settings) {

            // evaluate at x and x + d
            double fx = f(x);
            double y = x + d;
            double fy = f(y);
            int count = 2;

            // if we stepped uphill, reverse direction of steps and exchange x & y
            if (fy > fx) {
                double t = x; x = y; y = t;
                t = fx; fx = fy; fy = t;
                d = -d;
            }

            // we now know f(x) >= f(y) and we are stepping downhill
            // continue stepping until we step uphill
            double z, fz;
            while (true) {

                if (count >= settings.EvaluationBudget) throw new NonconvergenceException();

                z = y + d;
                fz = f(z);
                count++;

                Console.WriteLine("f({0})={1} f({2})={3} f({4})={5} d={6}", x, fx, y, fy, z, fz, d);

                if (fz > fy) break;

                // increase the step size each time
                d = AdvancedMath.GoldenRatio * d;

                // x <- y <- z
                x = y; fx = fy; y = z; fy = fz;


            }

            // we x and z now bracket a local minimum, with y the lowest point evaluated so far
            double a = Math.Min(x, z); double b = Math.Max(x, z);
            if (fz < fx) { double t = x; x = z; z = t; t = fx; fx = fz; fz = t; }

            return (FindMinimum(f, a, b, y, fy, x, fx, z, fz, settings, count));

        }

        [TestMethod]
        public void TestMinimizationWithoutDerivative () {

            //double t = FindMinimum(x => x * Math.Log(x), 0.1, 10.0); // 13
            //double t = FindMinimum(x => Math.Cos(x), 0.0, 5.0); // 7
            //double t = FindMinimum(x => AdvancedMath.Gamma(x), 0.0, 2.0); // 8
            //double t = FindMinimum(x => Math.Exp(-x), 0.0, 1.0); // 29 (not a minimum)
            //double t = FindMinimum(x => Math.Abs(x), -2.0, 3.0); // 82!
            //double t = FindMinimum(x => Math.Cosh(x), -4.0, 3.0); // 7
            //double t = FindMinimum(x => -1.0 / Math.Cosh(x), -1.5, 2.5); // 3
            //double t = FindMinimum(x => Math.Pow(x, 4.0), -2.0, 1.0); // 19

            //double t = FindMinimum(x => x * Math.Log(x), 1.0, -0.03125, new EvaluationSettings() { EvaluationBudget = 128, AbsolutePrecision = 0.0, RelativePrecision = 0.0 });
            double t = FindMinimum(AdvancedMath.Gamma, 1.0, 0.03125, new EvaluationSettings() { EvaluationBudget = 128, AbsolutePrecision = 0.0, RelativePrecision = 0.0 });
            Console.WriteLine(t);
        }

        [TestMethod]
        public void TestOldMinimization () {
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (x * Math.Log(x)); }, Interval.FromEndpoints(0.1, 10.0));
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (Math.Cos(x)); }, Interval.FromEndpoints(0.0, 5.0));
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (AdvancedMath.Gamma(x)); }, Interval.FromEndpoints(0.0, 2.0));
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (Math.Exp(-x)); }, Interval.FromEndpoints(0.0, 1.0));
            //FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (Math.Abs(x)); }, Interval.FromEndpoints(-2.0, 3.0));
            FunctionMath.FindMinimum(x => { Console.WriteLine(x); return (Math.Pow(x,4)); }, Interval.FromEndpoints(-2.0, 1.0));

        }

        [TestMethod]
        public void TestMinimizationWithDerivative () {

            //double x = FindMinimumWithDerivative(XLogXWithDerivative, 0.1, 10.0); // 9
            //double x = FindMinimumWithDerivative(CosineWithDerivative, 0.0, 5.0); // 5
            //double x = FindMinimumWithDerivative(GammaWithDerivative, 0.0, 2.0); // 9
            //double x = FindMinimumWithDerivative(ExponentialDecayWithDerivative, 0.0, 1.0); // 29 (not a minimum)
            //double x = FindMinimumWithDerivative(AbsoluteValueWithDerivative,-2.0, 3.0); // 19
            double x = FindMinimumWithDerivative(CoshWithDerivative, -4.0, 3.0); // 7
            //double x = FindMinimumWithDerivative(SechWithDerivative, -1.5, 2.5); // 5
            //double x = FindMinimumWithDerivative(FourthPowerWithDerivative, -2.0, 1.0); // 17
            Console.WriteLine(x);

        }

        private void CosineWithDerivative (double x, out double f, out double fp) {
            f = Math.Cos(x);
            fp = -Math.Sin(x);
        }

        private void GammaWithDerivative (double x, out double f, out double fp) {
            f = AdvancedMath.Gamma(x);
            fp = f * AdvancedMath.Psi(x);
        }

        private void XLogXWithDerivative (double x, out double f, out double fp) {
            f = x * Math.Log(x);
            fp = 1.0 + Math.Log(x);
        }

        private void ExponentialDecayWithDerivative (double x, out double f, out double fp) {
            f = Math.Exp(-x);
            fp = -f;
        }

        private void AbsoluteValueWithDerivative (double x, out double f, out double fp) {
            f = Math.Abs(x);
            fp = Math.Sign(x);
        }

        private void CoshWithDerivative (double x, out double f, out double fp) {
            f = Math.Cosh(x);
            fp = Math.Sinh(x);
        }

        private void SechWithDerivative (double x, out double f, out double fp) {
            f = -1.0 / Math.Cosh(x);
            fp = Math.Tanh(x) / Math.Cosh(x);
        }

        private void FourthPowerWithDerivative (double x, out double f, out double fp) {
            f = Math.Pow(x, 4);
            fp = 4.0 * Math.Pow(x, 3);
        }

        /* PRIME FACTORIZATION */

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

        /* INVERSE BETA */

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
        
        private static void WriteMatrix (AnyRectangularMatrix A) {
            for (int r = 0; r < A.RowCount; r++) {
                for (int c = 0; c < A.ColumnCount; c++) {
                    Console.Write("{0} ", A[r, c]);
                }
                Console.WriteLine();
            }
        }

    }

#endif

    }
