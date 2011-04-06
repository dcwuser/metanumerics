using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.ObjectModel;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.SignalProcessing;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace FutureTest {

    [TestClass]
    public class FutureTest {

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

        private Complex[] BluesteinCoefficients (int N) {

            Complex[] x = new Complex[N];
            double t = 2.0 * Math.PI / N;
            for (int i = 0; i < N; i++) {
                double ti = (i * i / 2.0) * t;
                x[i] = new Complex(Math.Cos(ti), -Math.Sin(ti));
            }
            return (x);
        }

        private void WriteSeries (Complex[] z) {
            for (int i = 0; i < z.Length; i++) {
                Console.WriteLine("{0} {1}", i, z[i]);
            }
        }

        [TestMethod]
        public void Bluestein () {

            Complex[] x = new Complex[5];
            x[1] = 1.0;
            Console.WriteLine("x=");
            WriteSeries(x);

            Complex[] b = BluesteinCoefficients(x.Length);

            FourierTransformer ft = new FourierTransformer(10);

            Complex[] x1 = new Complex[ft.Length];
            for (int i = 0; i < x.Length; i++) {
                x1[i] = b[i] * x[i];
            }
            Console.WriteLine("x1=");
            WriteSeries(x1);

            Complex[] b1 = new Complex[ft.Length];
            for (int i = 0; i < b.Length; i++) {
                b1[i] = b[i].Conjugate;
            }
            for (int i = 1; i < b.Length; i++) {
                b1[b1.Length - i] = b[i].Conjugate;
            }
            Console.WriteLine("b1=");
            WriteSeries(b1);

            Complex[] x1t = ft.Transform(x1);
            Complex[] b1t = ft.Transform(b1);

            Complex[] xb1t = new Complex[ft.Length];
            Console.WriteLine("xb1t=");
            for (int i = 0; i < ft.Length; i++) {
                xb1t[i] = x1t[i] * b1t[i];
                Console.WriteLine("{0} {1} * {2} = {3}", i, x1t[i], b1t[i], xb1t[i]);
            }

            Complex[] xb1 = ft.InverseTransform(xb1t);
            Console.WriteLine("xb1=");
            WriteSeries(xb1);

            Complex[] xt = new Complex[x.Length];
            for (int i = 0; i < xt.Length; i++) {
                xt[i] = b[i] * xb1[i];
            }
            Console.WriteLine("xt=");
            WriteSeries(xt);

        }


    }

}
