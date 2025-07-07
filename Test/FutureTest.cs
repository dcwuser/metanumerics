﻿using FluentAssertions;
using Meta.Numerics;
using Meta.Numerics.Data;
using Meta.Numerics.Extended;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Meta.Numerics.SignalProcessing;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Security.Policy;
using System.Xml.Linq;
using Assert = Microsoft.VisualStudio.TestTools.UnitTesting.Assert;
using ExpectedExceptionAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.ExpectedExceptionAttribute;
using TestClassAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestClassAttribute;
using TestMethodAttribute = Microsoft.VisualStudio.TestTools.UnitTesting.TestMethodAttribute;

namespace Test {


    [TestClass]
    public class FutureTest {

        public static void BadFirstPass (IEnumerable<double> xs, out int n, ref double m, out double c2, out double c3) {

            n = 0;
            //m = 0.5;

            c2 = 0.0;
            c3 = 0.0;
            foreach (double x in xs) {
                n++;
                double z = x - m;
                double z2 = z * z;
                double z3 = z * z2;
                c2 += z2;
                c3 += z3;
            }
            c2 /= n;
            c3 /= n;

        }

        public static void BadFirstPassPlusCorrection(IEnumerable<double> xs, out int n, out double m, out double c2, out double c3) {

            m = 0.5;
            BadFirstPass(xs, out n, ref m, out c2, out c3);

            double delta_m = 0.0;
            foreach (double x in xs) {
                double z = x - m;
                delta_m += z;
            }
            delta_m /= n;
            c3 += 3.0 * c2 * delta_m;
            c2 -= 1.0 * MoreMath.Sqr(delta_m);
            m += delta_m;

        }


        public static void NaiveTwoPass (IEnumerable<double> xs, out int n, out double m, out double c2, out double c3, out double c4) {

            n = 0;
            m = 0.0;

            foreach (double x in xs) {
                n++;
                m += x;
            }
            m /= n;

            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;
            foreach (double x in xs) {
                double z = x - m;
                double z2 = z * z;
                c2 += z2;
                c3 += z * z2;
                c4 += z2 * z2;
            }
            c2 /= n;
            c3 /= n;
            c4 /= n;

        }

        public static double[] CorrectionPass (IEnumerable<double> values, ref double mean, int maxMoment) {

            int n = 0;
            double delta = 0.0;
            double[] centralMoments = new double[maxMoment + 1];
            foreach (double value in values) {
                n++;
                double z = value - mean;
                delta += z;
                double zPower = z;
                for (int r = 2; r <= maxMoment; r++) {
                    zPower *= z;
                    centralMoments[r] += zPower;
                }
            }
            delta /= n;
            for (int r = 2; r <= maxMoment; r++) centralMoments[r] /= n;

            for (int r = maxMoment; r >= 3; r--) {
                centralMoments[r] += r * delta * (0.5 * (r - 1) * centralMoments[r - 2] * delta - centralMoments[r - 1]);
            }
            centralMoments[3] += delta * (2.0 * MoreMath.Sqr(delta) - 3.0 * centralMoments[2]);
            centralMoments[2] -= MoreMath.Sqr(delta);
            return centralMoments;
        }


        // Suppose we have computed central moments based on a value of the mean m (which will we suppose to be slightly wrong).
        //   c_r = n^{-1} \sum_i (x_i - m)^r
        // We will then find
        //   c_1 = ( n^{-1} \sum_i x_i - m ) \ne 0 = \delta m
        // is a correction to the mean. The true mean is m + \delta m.

        // We can then use Taylor expansion to get corrections to all the other moments.
        //   \frac{\partial c_r}{\partial m} = -r n^{-1} \sum_i (x_i - m)^{r-1} = -r c_{r-1}
        //   \frac{\partial^2 c_r}{\partial m^2} = r(r-1) n^{-1} \sum_i (x_i - m)^{r-2} = r(r-1) c_{r-2}
        // So for any r > 3, to 2nd order in (\delta m),
        //   \delta c_r = -r c_{r-1}  (\delta m) + 1/2 r(r-1) c_{r-2} (\delta m)^2 = (\delta m) [ 1/2 r (r-1) c_{r-2} (\delta m) - r c_{r-1} ]

        // For r = 2 or 3, c_1 = (\delta m) appears, which messes with the relative order of terms
        //   \delta c_2 = -2 c_1 (\delta m) + c_0 (\delta m)^2 = -2 (\delta m)^2 + (\delta m)^2 = - (\delta m)^2
        //   \delta c_3 = -3 c_2 (\delta m) + 3 c_1 (\delta m)^2 - 1 c_0 (\delta m)^3 = -3 c_2 (\delta m) + 3 (\delta m)^3 - (\delta m)^3 = (\delta m) [ 2 (\delta m)^2 - 3 c_2 ]

        public static void CorrectionPass (IEnumerable<double> xs, out int n, ref double m, out double c2, out double c3, out double c4) {

            n = 0;

            // Compute moments assuming the original m
            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;
            double delta_m = 0.0;
            foreach (double x in xs) {
                n++;
                double z = x - m;
                double z2 = z * z;
                delta_m += z;
                c2 += z2;
                c3 += z * z2;
                c4 += z2 * z2;
            }
            delta_m /= n;
            c2 /= n;
            c3 /= n;
            c4 /= n;

            // Add corrections based on computed delta_m
            // They must be ordered higher to lower moments, since original lower moments appear in corrections to higher
            c4 += delta_m * (6.0 * c2 * delta_m - 4.0 * c3);
            c3 += delta_m * (2.0 * MoreMath.Sqr(delta_m) - 3.0 * c2);
            c2 -= MoreMath.Sqr(delta_m);
            m += delta_m;

        }

        public static void KahanCorrectionPass(IEnumerable<double> xs, out int n, ref double m, out double c2, out double c3, out double c4) {

            n = 0;

            // Compute moments assuming original m
            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;
            double delta_m = 0.0;
            double delta_m_correction = 0.0;
            foreach (double x in xs) {
                n++;
                double z = x - m;
                delta_m = KahanSum(delta_m, z, ref delta_m_correction);
                double z2 = z * z;
                c2 += z2;
                c3 += z * z2;
                c4 += z2 * z2;
            }
            delta_m /= n;
            c2 /= n;
            c3 /= n;
            c4 /= n;

            // Add corrections
            c4 += delta_m * (6.0 * c2 * delta_m - 4.0 * c3);
            c3 += delta_m * (2.0 * MoreMath.Sqr(delta_m) - 3.0 * c2);
            c2 -= MoreMath.Sqr(delta_m);
            m += delta_m;

        }

        public static void DoubleDoubleMoments(IEnumerable<double>xs, out int n, out double m, out double c2, out double c3, out double c4) {

            n = 0;
            DoubleDouble m_dd = DoubleDouble.Zero;

            foreach (double x in xs) {
                n++;
                m_dd += x;
            }
            m_dd /= n;

            DoubleDouble c2_dd = DoubleDouble.Zero;
            DoubleDouble c3_dd = DoubleDouble.Zero;
            DoubleDouble c4_dd = DoubleDouble.Zero;

            foreach (double x in xs) {
                DoubleDouble z = x - m_dd;
                DoubleDouble z2 = z * z;
                c2_dd += z2;
                c3_dd += z * z2;
                c4_dd += z2 * z2;
            }

            m = (double) m_dd;
            c2 = (double) c2_dd / n;
            c3 = (double) c3_dd / n;
            c4 = (double) c4_dd / n;

        }

        public static void ThreePass(IEnumerable<double> xs, out int n, out double m, out double c2, out double c3, out double c4) {

            n = 0;
            m = 0.0;

            foreach (double x in xs) {
                n++;
                m += x;
            }
            m /= n;

            double delta_m = 0.0;
            foreach (double x in xs) delta_m += (x - m);
            m += delta_m / n;

            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;
            foreach (double x in xs) {
                double z = x - m;
                double z2 = z * z;
                c2 += z2;
                c3 += z * z2;
                c4 += z2 * z2;
            }
            c2 /= n;
            c3 /= n;
            c4 /= n;
        }


        public static void NaiveMeng (IEnumerable<double> xs, out int n, out double m, out double c2, out double c3, out double c4) {

            n = 0;
            m = 0.0;
            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;
            foreach (double x in xs) {
                n++;
                double d = x - m;
                double e = d / n;
                m += e;
                double ddme = d * (d - e);
                c2 += ddme;
                c3 += -3.0 * e * c2 + ddme * (d + e);
                c4 += -e * (4.0 * c3 + 6.0 * e * c2) + ddme * (e * e + d * e + d * d);
            }
            c2 /= n;
            c3 /= n;
            c4 /= n;
        }

        // In experiments it looks to me like the Kahan summation does a lot to improve m,
        // does less to improve c2, and does very little to improve higher moments.
        // Probably because more noisy operations are performed to compute what is being added to higher moments.
        // Because of this, we only do Kahan for m and c2.

        public static void KahanMeng(IEnumerable<double> xs, out int n, out double m, out double c2, out double c3, out double c4) {

            n = 0;
            m = 0.0;
            c2 = 0.0;
            c3 = 0.0;
            c4 = 0.0;

            double m_correction = 0.0;
            double c2_correction = 0.0;

            foreach (double x in xs) {
                n++;
                double d = x - m;
                double e = d / n;
                m = KahanSum(m, e, ref m_correction);
                double ddme = d * (d - e);
                c2 = KahanSum(c2, ddme, ref c2_correction);

                c3 += -3.0 * e * c2 + ddme * (d + e);
                c4 += -e * (4.0 * c3 + 6.0 * e * c2) + ddme * (e * e + d * e + d * d);
            }
            c2 /= n;
            c3 /= n;
            c4 /= n;
        }

        public static void RunningMean(IEnumerable<double> xs, out int n, out double m) {
            n = 0;
            m = 0.0;
            foreach (double x in xs) {
                n++;
                m += (x - m) / n;
            }
        }

        public static void KahanRunningMean (IEnumerable<double> xs, out int n, out double m) {

            n = 0;
            m = 0.0;
            double c = 0.0;
            
            foreach (double x in xs) {
                n++;
                double z = (x - m) / n;
                m = KahanSum(m, z, ref c);
            }

            // Naive mean needs 1 flop per cycle
            // Naive running mean needs 3 flops per cycle
            // Kahan running mean needs 6 flops per cycle

        }

        // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
        // Kahan suggested a way to do a running sum more accurately by keeping track of a correction term
        // This is essentially FastTwoSum with a re-used error term

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double KahanSum (double augend, double addend, ref double correction) {
            double y = addend - correction;
            double sum = augend + y;
            correction = (sum - augend) - y;
            return sum;
        }

        public static void NaiveWelford (IEnumerable<double> xs, out int n, out double m, out double c2) {

            n = 0;
            m = 0.0;
            c2 = 0.0;          
            foreach (double x in xs) {
                n++;
                double z = x - m;
                m += z / n;
                c2 += z * (x - m);
            }
            c2 /= n;
        }

        [TestMethod]
        public void OnlineMomentsTest () {

            Random rng = new Random(1);
            SummaryStatistics summary = new SummaryStatistics();
            List<double> values = new List<double>();
            ContinuousDistribution distribution = new LogisticDistribution(1.0, 4.0);

            for (int i = 0; i < 10000; i++) {
                double value = distribution.GetRandomValue(rng); //.NextDouble();
                values.Add(value);
                summary.Add(value);
            }

            int n;
            double m, c2, c3, c4;

            values.Count.Should().Be(summary.Count);
            values.Mean().Should().BeNearly(summary.Mean);
            values.Variance().Should().BeNearly(summary.Variance);
            //values.CentralMoment(3).Should().BeNearly(summary.ThirdCentralMoment);

            /*
            double mean = 0.0;
            foreach (double value in values) mean += value;
            mean = mean / values.Count;

            double meanDelta = 0.0;
            foreach (double value in values) meanDelta += (value - mean);
            meanDelta = meanDelta / values.Count;
            double meanPrime = mean + meanDelta;

            double meanPrimeDelta = 0.0;
            foreach (double value in values) meanPrimeDelta += (value - meanPrime);
            meanPrimeDelta = meanPrimeDelta / values.Count;
            double meanPrimePrime = meanPrime + meanPrimeDelta;

            double c3 = 0.0;
            foreach (double value in values) c3 += MoreMath.Pow(value - meanPrime, 3);
            c3 = c3 / summary.Count;

            UncertainValue m = values.PopulationMean();
            m.Value.Should().BeNearly(summary.PopulationMean.Value);
            m.Uncertainty.Should().BeNearly(summary.PopulationMean.Uncertainty);
            */

        }

        [TestMethod]
        public void CountGroupsThree () {

            int n = 5;

            int allEqual = 0;
            int twoEqual = 0;
            int allDifferent = 0;

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < n; k++) {

                        if ((i == j) && (j ==k)) {
                            allEqual++;
                            continue;
                        }

                        if ((i == j) || (i == k) || (j == k)) {
                            twoEqual++;
                            continue;
                        }

                        allDifferent++;
                    }
                }
            }

            Console.WriteLine($"All equal {allEqual} predicted {n}");
            Console.WriteLine($"Two equal {twoEqual} predicted {3 * n * (n - 1)}");
            Console.WriteLine($"All different {allDifferent} predicted {n * (n - 1) * (n - 2)}");
        }


        [TestMethod]
        public void CountGroupsFour() {

            int n = 6;

            int allEqual = 0;
            int threeEqual = 0;
            int twoPairs = 0;
            int twoEqual = 0;
            int allDifferent = 0;

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < n; k++) {
                        for (int l = 0; l < n; l++) {

                            if ((i == j) && (j == k) && (k == l)) {
                                allEqual++;
                                continue;
                            }

                            if (((i == j) && (j == k)) || ((i == j) && (j == l)) || ((i == k) && (k == l)) || ((j ==k) && (k ==l))) {
                                threeEqual++;
                                continue;
                            }

                            if (((i == j) && (k == l)) || ((i == k) && (j == l)) || ((i == l) && (j == k))) {
                                twoPairs++;
                                continue;
                            }

                            if ((i ==j) || (i == k) || (i == l) || (j == k) || (j == l) || (k == l)) {
                                twoEqual++;
                                continue;
                            }

                            allDifferent++;

                        }

                    }
                }
            }

            Console.WriteLine($"All equal {allEqual} predicted {n}");
            Console.WriteLine($"Two pairs {twoPairs} predicted {3 * n * (n - 1)}");
            Console.WriteLine($"Three equal {threeEqual} predicted {4 * n * (n - 1)}");
            Console.WriteLine($"Two equal {twoEqual} predicted {6 * n * (n - 1) * (n - 2)}");
            Console.WriteLine($"All different {allDifferent} predicted {n * (n - 1) * (n - 2) * (n - 3)}");
        }

        public void TestSample (ContinuousDistribution distribution, IReadOnlyList<double> sample) {

            Meta.Numerics.Statistics.TestResult result = sample.KolmogorovSmirnovTest(distribution);
            if (result.Probability < 0.01) throw new InvalidOperationException();

        }

        public void TestDistribution (ContinuousDistribution distribution, Random rng, int n) {

            List<double> sample = distribution.GetRandomValues(rng, n).ToList();
            TestSample(distribution, sample);

        }

        [TestMethod]
        public void MeanZeroMomentsTest () {

            ContinuousDistribution[] distributions = new ContinuousDistribution[] {
                //new NormalDistribution(),
                //new LogisticDistribution(),
                //new UniformDistribution(Interval.FromEndpoints(-1.0, 1.0)),
                new TriangularDistribution(-3.0, 1.0, 2.0)
            };

            int n = 10;
            Random rng = new Random(271828);

            foreach (ContinuousDistribution xDistribution in distributions) {

                xDistribution.Mean.Should().Be(0.0);

                List<double> m_32 = new List<double>();
                List<double> m_31_21_11 = new List<double>();
                List<double> m_22_12 = new List<double>();
                List<double> m_31_13 = new List<double>();
                List<double> m_21_14 = new List<double>();
                List<double> m_16 = new List<double>();

                for (int i = 0; i < 100000; i++) {
                    List<double> xSample = xDistribution.GetRandomValues(rng, n).ToList();

                    double m1 = xSample.RawMoment(1);
                    double m2 = xSample.RawMoment(2);
                    double m3 = xSample.RawMoment(3);
                    double m4 = xSample.RawMoment(4);
                    double m6 = xSample.RawMoment(6);

                    m_32.Add(m3 * m3);
                    m_31_21_11.Add(m3 * m2 * m1);
                    m_22_12.Add(m2 * m2 * m1 * m1);
                    m_31_13.Add(m3 * m1 * m1 * m1);
                    m_21_14.Add(m2 * m1 * m1 * m1 * m1);
                    m_16.Add(m1 * m1 * m1 * m1 * m1 * m1);
                }

                double C_61 = xDistribution.CentralMoment(6);
                double C_41_21 = xDistribution.CentralMoment(4) * xDistribution.CentralMoment(2);
                double C_32 = MoreMath.Sqr(xDistribution.CentralMoment(3));
                double C_23 = MoreMath.Pow(xDistribution.CentralMoment(2), 3);

                m_32.PopulationMean().ConfidenceInterval(0.99).Contains((C_61 + (n - 1) * C_32) / n).Should().BeTrue();
                m_31_21_11.PopulationMean().ConfidenceInterval(0.99).Contains((C_61 + (n-1) * (C_41_21 + C_32)) / MoreMath.Sqr(n)).Should().BeTrue();
                m_22_12.PopulationMean().ConfidenceInterval(0.99).Contains((C_61 + (n - 1) * (3.0 * C_41_21 + 2.0 * C_32) + (n - 1) * (n - 2) * C_23) / MoreMath.Pow(n, 3)).Should().BeTrue();
                m_31_13.PopulationMean().ConfidenceInterval(0.99).Contains((C_61 + (n - 1) * (3.0 * C_41_21 + 1.0 * C_32)) / MoreMath.Pow(n, 3)).Should().BeTrue();
                m_21_14.PopulationMean().ConfidenceInterval(0.99).Contains((C_61 + (n - 1) * (7.0 * C_41_21 + 4.0 * C_32) + 3 * (n - 1) * (n - 2) * C_23) / MoreMath.Pow(n, 4)).Should().BeTrue();
                m_16.PopulationMean().ConfidenceInterval(0.99).Contains((C_61 + (n-1) * (15.0 * C_41_21 + 10.0 * C_32) + 15 * (n-1) * (n-2) * C_23) / MoreMath.Pow(n, 5)).Should().BeTrue();
            }

        }

        [TestMethod]
        public void KStatisticTest() {

            ContinuousDistribution[] distributions = new ContinuousDistribution[] {
                new NormalDistribution(1.0, 2.0),
                new ExponentialDistribution(3.0),
                new UniformDistribution(Interval.FromEndpoints(-4.0, 5.0)),
                new BetaDistribution(2.0 / 3.0, 4.0),
                new GammaDistribution(5.0, 6.0)
            };

            foreach (ContinuousDistribution distribution in distributions) {
                KStatisticTest(distribution);
            }

        }

        // If I changed any factor n to n+1 or n-1, I saw this fail
        // Also, estimators don't agree for LogNormal cumulants; check generator / cumulant formulas
        // Finally, k3 estimator variance failed for Exponential and Gamma

        public void KStatisticTest(ContinuousDistribution xDistribution) {

            int n = 6;
            //int n = 20;

            Random rng = new Random(271828);

            List<double> m1s = new List<double>();
            List<double> m2s = new List<double>();
            List<double> m3s = new List<double>();
            List<double> m4s = new List<double>();

            List<double> m1vs = new List<double>();
            List<double> m2vs = new List<double>();
            List<double> m3vs = new List<double>();
            List<double> m4vs = new List<double>();

            List<double> c2s = new List<double>();
            List<double> c3s = new List<double>();
            List<double> c4s = new List<double>();

            List<double> c22s = new List<double>();
            List<double> c2c4s = new List<double>();
            List<double> c23s = new List<double>();
            List<double> c32s = new List<double>();

            List<double> h2s = new List<double>();
            List<double> h3s = new List<double>();
            List<double> h4s = new List<double>();

            List<double> h2vs = new List<double>();
            List<double> h3vs = new List<double>();
               
            List<double> k1s = new List<double>();
            List<double> k2s = new List<double>();
            List<double> k3s = new List<double>();
            List<double> k4s = new List<double>();

            List<double> k1vs = new List<double>();
            List<double> k2vs = new List<double>();

            List<double> s0s = new List<double>();
            List<double> s1s = new List<double>();
            List<double> s2s = new List<double>();

            List<double> g0s = new List<double>();
            List<double> g1s = new List<double>();


            for (int i = 0; i < 10000; i++) {

                List<double> xSample = xDistribution.GetRandomValues(rng, n).ToList();

                double m1 = xSample.RawMoment(1);
                double m2 = xSample.RawMoment(2);
                double m3 = xSample.RawMoment(3);
                double m4 = xSample.RawMoment(4);

                m1s.Add(m1);
                m2s.Add(m2);
                m3s.Add(m3);
                m4s.Add(m4);

                double m1v = (m2 - m1 * m1) / (n - 1);
                if (m1v < 0.0) throw new InvalidOperationException();
                double m2v = (m4 - m2 * m2) / (n - 1);
                if (m2v < 0.0) throw new InvalidOperationException();
                double m3v = (xSample.RawMoment(6) - m3 * m3) / (n - 1);
                if (m3v < 0.0) throw new InvalidOperationException();
                double m4v = (xSample.RawMoment(8) - m4 * m4) / (n - 1);
                if (m4v < 0.0) throw new InvalidOperationException();

                m1vs.Add(m1v);
                m2vs.Add(m2v);
                m3vs.Add(m3v);
                m4vs.Add(m4v);

                double c2 = xSample.CentralMoment(2);
                double c3 = xSample.CentralMoment(3);
                double c4 = xSample.CentralMoment(4);
                double c6 = xSample.CentralMoment(6);

                c2s.Add(c2);
                c3s.Add(c3);
                c4s.Add(c4);

                c22s.Add(c2 * c2);
                c23s.Add(c2 * c2 * c2);
                c2c4s.Add(c2 * c4);
                c32s.Add(c3 * c3);

                double h2 = c2 * n / (n - 1);
                double h3 = c3 * n * n / (n - 1) / (n - 2);
                double h4 = ((n * n - 2 * n + 3) * c4 - 3 * (2 * n - 3) * c2 * c2) * n / (n - 1) / (n - 2) / (n - 3);

                h2s.Add(h2);
                h3s.Add(h3);
                h4s.Add(h4);

                double h2v = (c4 - c2 * c2 * (n * n - 3) / MoreMath.Sqr(n - 1)) * n / (n - 2) / (n - 3);
                if (h2v < 0.0) Console.WriteLine($"{xDistribution.GetType().Name} h2v={h2v}");
                double h3v = (c6 * (n * n - n + 4) / (n * n)
                    - 3.0 * c4 * c2 * (2 * n * n * n - 5 * n * n - 5 * n + 20) / (n * n * (n - 1))
                    - c3 * c3 * (MoreMath.Pow(n, 5) + 4 * MoreMath.Pow(n, 4) - 41 * n * n * n + 40 * n * n + 100 * n - 80) / (n * n * (n - 1) * (n-1) * (n - 2))
                    + 3.0 * c2 * c2 * c2 * (3 * n * n - 15 * n + 20) / (n * (n - 1))) * n * n * n / (n-2) / (n-3) / (n-4) / (n-5);
                if (h3v < 0.0) Console.WriteLine($"{xDistribution.GetType().Name} h3v={h3v}");

                h2vs.Add(h2v);
                h3vs.Add(h3v);

                double k1 = xSample.Mean();
                double k2 = xSample.CentralMoment(2) * n / (n - 1);
                double k3 = xSample.CentralMoment(3) * n * n / (n - 1) / (n - 2);
                double k4 = ((n + 1) * xSample.CentralMoment(4) - 3.0 * (n - 1) * MoreMath.Sqr(xSample.CentralMoment(2))) * n * n / (n - 1) / (n - 2) / (n - 3);

                k1s.Add(k1);
                k2s.Add(k2);
                k3s.Add(k3);
                k4s.Add(k4);

                double k1v = k2 / n;
                if (k1v < 0.0) throw new InvalidOperationException();
                double k2v = ((n - 1) * k4 + 2 * n * k2 * k2) / n / (n + 1);
                if (k2v < 0.0) Console.WriteLine($"{xDistribution.GetType().Name} k2v={k2v}");

                k1vs.Add(k1v);
                k2vs.Add(k2v);

                double s0 = Math.Sqrt(h2);
                double s1 = s0 * (1.0 + (h4 - h2 * h2) / (8.0 * n * h2 * h2));

                s0s.Add(s0);
                s1s.Add(s1);

                double g0 = c3 / Math.Pow(c2, 3.0 / 2.0);

                g0s.Add(g0);

            }

            // Sample raw moments should on average agree with underlying raw mmoments
            m1s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.RawMoment(1)).Should().BeTrue();
            m2s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.RawMoment(2)).Should().BeTrue();
            m3s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.RawMoment(3)).Should().BeTrue();
            m4s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.RawMoment(4)).Should().BeTrue();

            // The variances of the sample raw moments (which are also predictors of the population raw moments) should agree with the predictions from the population raw moments
            m1s.PopulationVariance().ConfidenceInterval(0.99).Contains(xDistribution.CentralMoment(2) / n).Should().BeTrue();
            m2s.PopulationVariance().ConfidenceInterval(0.99).Contains((xDistribution.RawMoment(4) - MoreMath.Sqr(xDistribution.RawMoment(2))) / n).Should().BeTrue();
            // For an exponential distribution an n = 4 these failed
            //if (!(xDistribution is ExponentialDistribution)) {
                m3s.PopulationVariance().ConfidenceInterval(0.99).Contains((xDistribution.RawMoment(6) - MoreMath.Sqr(xDistribution.RawMoment(3))) / n).Should().BeTrue();
                m4s.PopulationVariance().ConfidenceInterval(0.99).Contains((xDistribution.RawMoment(8) - MoreMath.Sqr(xDistribution.RawMoment(4))) / n).Should().BeTrue();
            //} // Are our raw moments for exponential distribution wrong? Doesn't look like it. Maybe random generator? Maybe has to do with exponentiation/log expanding/compressing scale.

            // The variances of the sample raw moments (which are also predictors of the underlying raw moments) should agree with the predictions from the sample raw moments
            m1s.PopulationVariance().ConfidenceInterval(0.99).Contains(m1vs.Mean()).Should().BeTrue();
            m2s.PopulationVariance().ConfidenceInterval(0.99).Contains(m2vs.Mean()).Should().BeTrue();
            m3s.PopulationVariance().ConfidenceInterval(0.99).Contains(m3vs.Mean()).Should().BeTrue();
            m4s.PopulationVariance().ConfidenceInterval(0.99).Contains(m4vs.Mean()).Should().BeTrue();

            // Sample central moments should on average have values predicted from underlying central moments
            c2s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.CentralMoment(2) * (n-1) / n).Should().BeTrue();
            c3s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.CentralMoment(3) * (n-1) * (n-2) / MoreMath.Sqr(n)).Should().BeTrue();
            c4s.PopulationMean().ConfidenceInterval(0.99).Contains(((n * n - 3 * n + 3) * xDistribution.CentralMoment(4) + 3 * (2 * n - 3) * MoreMath.Sqr(xDistribution.CentralMoment(2))) * (n-1) / MoreMath.Pow(n, 3)).Should().BeTrue();

            // Sample central moments should exhibit the variance predicted from the underlying central moments
            c2s.PopulationVariance().ConfidenceInterval(0.99).Contains((
                (n - 1) * xDistribution.CentralMoment(4)
                - (n - 3) * MoreMath.Sqr(xDistribution.CentralMoment(2)))
                * (n - 1) / MoreMath.Pow(n, 3)).Should().BeTrue();
            c3s.PopulationVariance().ConfidenceInterval(0.99).Contains((
                xDistribution.CentralMoment(6)
                - 3 * xDistribution.CentralMoment(4) * xDistribution.CentralMoment(2) * (2 * n - 5) / (n - 1)
                - MoreMath.Sqr(xDistribution.CentralMoment(3)) * (n - 10) / (n - 1)
                + 3 * MoreMath.Pow(xDistribution.CentralMoment(2), 3) * (3 * n * n - 12 * n + 20) / (n - 1) / (n - 2))
                * (n - 1) * (n - 1) * (n - 2) * (n - 2) / MoreMath.Pow(n, 5)).Should().BeTrue();

            // Products of central moments should agree with predictions from underlying moments
            c22s.PopulationMean().ConfidenceInterval(0.99).Contains(((n - 1) * xDistribution.CentralMoment(4) + (n * n - 2 * n + 3) * MoreMath.Sqr(xDistribution.CentralMoment(2))) * (n - 1) / MoreMath.Pow(n, 3)).Should().BeTrue();
            c32s.PopulationMean().ConfidenceInterval(0.99).Contains((
                xDistribution.CentralMoment(6)
                - 3 * xDistribution.CentralMoment(4) * xDistribution.CentralMoment(2) * (2 * n - 5) / (n - 1)
                + MoreMath.Sqr(xDistribution.CentralMoment(3)) * (n * n - 2 * n + 10) / (n - 1)
                + 3 * MoreMath.Pow(xDistribution.CentralMoment(2), 3) * (3 * n * n - 12 * n + 20) / (n - 1) / (n - 2))
                * (n - 1) * (n - 1) * (n - 2) * (n - 2) / MoreMath.Pow(n, 5)).Should().BeTrue();

            // Estimated central moments shouild on average agree with underlying central moments
            h2s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.CentralMoment(2)).Should().BeTrue();
            h3s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.CentralMoment(3)).Should().BeTrue();
            h4s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.CentralMoment(4)).Should().BeTrue();

            // Variance of the estimated central moments should exhibit the variance predicted from the underlyling central moments
            h2s.PopulationVariance().ConfidenceInterval(0.99).Contains((
                xDistribution.CentralMoment(4)
                - MoreMath.Sqr(xDistribution.CentralMoment(2)) * (n - 3) / (n - 1)) / n).Should().BeTrue();
            h3s.PopulationVariance().ConfidenceInterval(0.99).Contains((
                xDistribution.CentralMoment(6)
                - 3 * xDistribution.CentralMoment(4) * xDistribution.CentralMoment(2) * (2 * n - 5) / (n-1)
                - MoreMath.Sqr(xDistribution.CentralMoment(3)) * (n - 10) / (n -1 )
                + 3 * MoreMath.Pow(xDistribution.CentralMoment(2), 3) * (3 * n * n - 12 * n + 20) / (n-1) / (n-2)) / n).Should().BeTrue();

            // Variance of the estimated central moments should on average agree with the estimated variance
            h2s.PopulationVariance().ConfidenceInterval(0.99).Contains(h2vs.Mean()).Should().BeTrue();
            h3s.PopulationVariance().ConfidenceInterval(0.99).Contains(h3vs.Mean()).Should().BeTrue();

            // Estimated cumulants should on average agree with the underlying cumulants 
            k1s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.Cumulant(1)).Should().BeTrue();
            k2s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.Cumulant(2)).Should().BeTrue();
            k3s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.Cumulant(3)).Should().BeTrue();
            k4s.PopulationMean().ConfidenceInterval(0.99).Contains(xDistribution.Cumulant(4)).Should().BeTrue();

            // Variance of estimated cumulants should agree with predictions from cumulants of underlying distribution
            k1s.PopulationVariance().ConfidenceInterval(0.99).Contains(xDistribution.Cumulant(2) / n).Should().BeTrue();
            k2s.PopulationVariance().ConfidenceInterval(0.99).Contains(xDistribution.Cumulant(4) / n + 2.0 * MoreMath.Sqr(xDistribution.Cumulant(2)) / (n - 1)).Should().BeTrue();
            if (!(xDistribution is ExponentialDistribution)) {
                k3s.PopulationVariance().ConfidenceInterval(0.999).Contains(xDistribution.Cumulant(6) / n + 9.0 * (xDistribution.Cumulant(4) * xDistribution.Cumulant(2) + MoreMath.Sqr(xDistribution.Cumulant(3))) / (n - 1) + 6.0 * MoreMath.Pow(xDistribution.Cumulant(2), 3) * n / (n - 1) / (n - 2)).Should().BeTrue();
            }

            // Variance of estimated cumulants should agree, on average, with predictions from sample-estimated cumulants
            k1s.PopulationVariance().ConfidenceInterval(0.99).Contains(k1vs.Mean()).Should().BeTrue();
            k2s.PopulationVariance().ConfidenceInterval(0.99).Contains(k2vs.Mean()).Should().BeTrue();

            // Square root approximations should get increasingly close
            // Observed: s1 does move us closer, but still consistently too small
            Math.Abs(s1s.Mean() - xDistribution.StandardDeviation).Should().BeLessThan(Math.Abs(s0s.Mean() - xDistribution.StandardDeviation));

            if (xDistribution.Skewness != 0.0) {
                double dg0 = Math.Abs(g0s.Mean() - xDistribution.Skewness);
                double dg = 33.0 / 8.0
                    + 15.0 / 8.0 * xDistribution.CentralMoment(4) / MoreMath.Sqr(xDistribution.CentralMoment(2))
                    - 3.0 / 2.0 * xDistribution.CentralMoment(5) / (xDistribution.CentralMoment(3) * xDistribution.CentralMoment(2));
                double dg1 = Math.Abs(g0s.Mean() - xDistribution.Skewness * (1.0 + dg / n));
            }

        }

        [TestMethod]
        public void SamplingTest () {

            int n = 4;

            ContinuousDistribution xDistribution = new NormalDistribution(0.0, 1.0);
            ContinuousDistribution sDistribution = new ChiDistribution(n - 1);

            Random rng = new Random(314159);

            List<double> sValues = new List<double>();
            for (int i = 0; i < 10000; i++) {

                List<double> xValues = xDistribution.GetRandomValues(rng, n).ToList();

                double m = xValues.Mean();

                double v = 0.0;
                foreach (double xValue in xValues) v += MoreMath.Sqr(xValue - m);
                sValues.Add(Math.Sqrt(v));
            }

            var result = sValues.KolmogorovSmirnovTest(sDistribution);

        }

        [TestMethod]
        public void TestGenerator () {

            //ContinuousDistribution distribution = new UniformDistribution(Interval.FromEndpoints(0.0, 1.0E-4));
            //ContinuousDistribution distribution = new ExponentialDistribution(1.0E-4);
            ContinuousDistribution distribution = new BetaDistribution(2.0, 3.0 / 2.0);
            Random rng = new Random(1);
            TestDistribution(distribution, rng, 1000000);

        }

        // Temme and Gil https://arxiv.org/pdf/1311.0681.pdf
        // See also https://arxiv.org/pdf/1306.1754.pdf for D(a, x) = x^a e^{-x} / \Gamma(a + 1)

        public double MarcumQ(double m, double x, double y) {

            double xPower = 1.0;
            double q = AdvancedMath.RightRegularizedGamma(m, y);
            double s = q;
            for (int k = 1; k < 100; k++) {
                double s_old = s;
                xPower *= x / k;
                //q = AdvancedMath.RightRegularizedGamma(m + k, y);
                q += Math.Pow(y, m) * Math.Exp(-y) / AdvancedMath.Gamma(m + 1);
                s += xPower * q;
                if (s == s_old) {
                    return Math.Exp(-x) * s;
                }
                m += 1;
            }

            throw new NonconvergenceException();
        }

        [TestMethod]
        public void TestMarcum () {

            MarcumQ(1, 10.0, 10.0);

        }


        [TestMethod]
        public void NonlinearRegressionLinearRegressionAgreementTest() {

            double[] x = new double[] { 1.0, 3.0, 3.0, 5.0 };
            double[] y = new double[] { 2.0, 2.0, 4.0, 4.0 };

            LinearRegressionResult linearFit = y.LinearRegression(x);
            NonlinearRegressionResult nonlinearFit = y.NonlinearRegression(x,
                (IReadOnlyList<double> p, double t) => p[0] + p[1] * t,
                new double[] { 1.0, 1.0 }
            );

            Assert.IsTrue(TestUtilities.IsNearlyEqual(linearFit.Parameters.ValuesVector, nonlinearFit.Parameters.ValuesVector, Math.Sqrt(TestUtilities.TargetPrecision)));

        }

        [TestMethod]
        public void TwoWayAnova() {

            List<double>[,] samples = new List<double>[,] {
                { new List<double> {54, 49, 59, 39, 55 }, new List<double> {25, 29, 47, 26, 28 } },
                { new List<double> {53, 72, 43, 56, 52 }, new List<double> {46, 51, 33, 47, 41 } },
                { new List<double> {33, 30, 26, 25, 29 }, new List<double> {18, 21, 34, 40, 24 } }
            };

            Univariate.TwoWayAnovaTest(samples);

        }


        [TestMethod]
        public void BinomialInverseCdf() {

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
        public void NegativeBinomialInverseCdf() {

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
        public void PoissonInverseCdf() {

            foreach (double mu in new double[] { 0.1, 2.0, 30.0 }) {
                PoissonDistribution d = new PoissonDistribution(mu);
                foreach (double P in new double[] { 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99 }) {

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

    }

    internal static class BesselAsymptotic {

        public static double ScaledBesselISeries (double nu, double z) {

            double t = MoreMath.Hypot(1.0, z);
            double y = 1.0 / (z + t) + Math.Log(z / (1.0 + t));

            double p = 1.0 / t;
            double pSquared = p * p;
            double s = 1.0;
            double u1 = p / 24.0 * (3.0 - 5.0 * pSquared);
            double s1 = u1 / nu;
            s += s1;
            double u2 = pSquared / 1152.0 * (81.0 - 462 * pSquared + 385.0 * pSquared * pSquared);
            double s2 = u2 / (nu * nu);
            s += s2;
            double u3 = p * pSquared / 414720.0 * (30375.0 - 369603.0 * pSquared + 765765 * pSquared * pSquared - 425425.0 * pSquared * pSquared * pSquared);
            double s3 = u3 / (nu * nu * nu);
            s += s3;
            double u4 = pSquared * pSquared / 39813120.0 * (4465125.0 - 94121676.0 * pSquared + 349922430.0 * pSquared * pSquared - 446185740.0 * pSquared * pSquared * pSquared + 185910725.0 * pSquared * pSquared * pSquared * pSquared);
            double s4 = u4 / (nu * nu * nu * nu);
            s += s4;

            return Math.Exp(nu * y) / Math.Sqrt(2.0 * Math.PI * nu * t) * s;

        }

    }


    // Uniform asymptotic exansion for Bessel functions of large order is summarized in A&S 9.3 and DLMF 10.20
    
    // The basic expression 
    //   J_{\nu}(\nu z) = \left(\frac{4 \zeta}{1 - z^2}\right)^{1/4} \left[
    //     \frac{Ai(\nu^{2/3} \zeta)}{\nu^{1/3}} \sum_{k=0}{\infty} \frac{a_k(\zeta)}{\nu^{2k}} +
    //     \frac{Ai'(\nu^{2/3} \zeta)}{\nu^{5/3}} \sum_{k=0}{\infty} \frac{b_k(\zeta)}{\nu^{2k}} \right]

    // Here \zeta is determined from z via a complicated expression that maps 0->\infinity, 1->0, and \infinity->-\infinity.

    // There are several complications that make this development difficult. The first complication is determing \zeta. \zeta
    // is related to z by
    //    \zeta z^2 = (1 - z^2) \left( \frac{dz}{d\zeta} \right)^2
    // which has the solution:
    //    0 < z < 1: 2/3 \zeta^{3/2} = \ln(\frac{1 + \sqrt{1-z^2}}{z}) - \sqrt{1-z^2}
    //    1 < z: 2/3 (-\zeta^{3/2}) = \sqrt{z^2 - 1} - \arccos(1/z)
    // which maps z < 1 to increasing positive \zeta, z > 1 to increasingly negative \zeta, and z = 1 to \zeta = 0.
    // These fromula are ugly, but also suffer from significant cancellation near z ~ 1, so to obtain accurate \zeta
    // in this region, it is required to develop the RHS in a series in e = z - 1.

    // The next complication is that the Airy functions must be evaluated with great precision, since they reproduce
    // the oscilations in J and Y. 

    // The next complication is that the a_k(\zeta) and b_k(\zeta) are even more complicated functions of \zeta. 
    
    
    // What makes this difficult is that the expressions for \zeta and a_k(\zeta) and b_k(\zeta) have canceling singularities
    // at z = 1 (\zeta = 0). While the functions are ultimimately analytic, the singularities that must cancel are quite
    // strong. (E.g. for b_0(\zeta) singularities of order \zeta^{-3} and \zeta^{-3/2} cancel near \zeta = 0, leaving a constant
    // value. The degree of cancelation gets higher with increasing k.) This makes the practical evaluation of these functions
    // very much more difficult than a cursory inspection leads one to expect.

    // We follow N. M. Temme, "Numerical algorithms for Airy-type asymptotic expansions", Numerical Algorithms 15 (1997) 207
    // (available at http://oai.cwi.nl/oai/asset/2178/2178A.pdf)

    internal static class BesselUae {

        private static readonly double OneOverSqrtTwoPi = 1.0 / Math.Sqrt(2.0 * Math.PI);

        public static SolutionPair Airy_Negative_Asymptotic (double x) {
            return (Airy_Negative_Asymptotic(x, 2.0 / 3.0 * Math.Pow(x, 3.0 / 2.0)));
        }

        public static SolutionPair Airy_Negative_Asymptotic (double x, double y) {

            double s = MoreMath.Sin(y);
            double c = MoreMath.Cos(y);

            int k = 1;
            double yInverse = 1.0 / y;
            double u = 5.0 / 72.0  * yInverse;
            double v = -7.0 / 72.0 * yInverse;
            double u_even = 1.0 + u;
            double u_odd = 1.0 - u;
            double v_even = 1.0 + v;
            double v_odd = 1.0 - v;
            while (k < 128) {
                double u_even_old = u_even;
                double u_odd_old = u_odd;
                double v_even_old = v_even;

                k++;
                u *= -(0.5 * (k - 1) + 5.0 / 72.0 / k) * yInverse;
                v = -(1.0 + 2.0 / (6 * k - 1)) * u;
                //u *= -(6 * k - 5) * (6 * k - 3) * (6 * k - 1) * z / ((2 * k - 1) * k);
                u_even += u;
                u_odd += u;
                v_even += v;
                v_odd += v;

                k++;
                u *= (0.5 * (k - 1) + 5.0 / 72.0 / k) * yInverse;
                v = -(1.0 + 2.0 / (6 * k - 1)) * u;
                //u *= (6 * k - 5) * (6 * k - 3) * (6 * k - 1) * z / ((2 * k - 1) * k);
                u_even += u;
                u_odd -= u;
                v_even += v;
                v_odd -= v;

                if ((u_even == u_even_old) && (u_odd == u_odd_old) && (v_even == v_even_old)) {
                    double xFourthRoot = Math.Pow(x, 1.0 / 4.0);
                    double n = OneOverSqrtTwoPi / xFourthRoot;
                    double m = OneOverSqrtTwoPi * xFourthRoot;
                    return new SolutionPair(
                        n * (s * u_even + c * u_odd), m * (s * v_odd - c * v_even),
                        n * (c * u_even - s * u_odd), m * (s * v_even + c * v_odd)
                    );
                }

            }
            throw new NonconvergenceException();

        }

        private static readonly double CbrtTwo = Math.Pow(2.0, 1.0 / 3.0);

        // The coefficeints of the series
        //   \zeta(z) = 2^{1/3} \sum_{k=1}^{\infty} c_k (z - 1)^k
        // which is a series devlopment of
        //   \zeta(z) = \sqrt{z^2 - 1} - \acos(1/z)
        // around z = 1.

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
            -19355768761912897160895124453273481957.0 / 720842364087730428996179297149658203125.0,
            20193853600897473400972883428162015393561.0 / 800135024137380776185759019836120605468750.0,
            -8998783481498741573746874301320903998981648.0 / 378063798904912416747771136872566986083984375.0,
            2355083733074596245284002631702400017597437.0 / 104592548954800331219019005477565765380859375.0,
            -8194181457446397986704636409997877371597271963852981.0 / 383635977269427222007301930021869110935955047607421875.0,
            46031971106083293475063974876599177738247302588282431.0 / 2266266583869589380969775880051968107073402404785156250.0,
            -445064482301791066659237430863379294049292971924644077632.0 / 22989385937870426278787568156560506472837106227874755859375.0,
            490476198188675169925160213749179549084836875610427400208.0 / 26526214543696645706293347872954430545581276416778564453125.0,
            -401690154573285433897214071223109161316172771519168044767966792.0 / 22702593348295492710959693243807414154588463327682018280029296875.0,
            2210868073858524888496964133636530234096101297055954328806560759492.0 / 130350723474796620648760238708194236270928760273107588291168212890625.0,
            -497213753307719277844728583245863479552772114438371686749150807854784.0 / 30532150229288900759651917451265496111152159617816354334354400634765625.0,
            193635801727377962771023777678202533839652872384146639842765327126560536.0 / 12365520842862004807659026567762525925016624645215623505413532257080078125.0
        };

        // If we could get analytic expression for terms of arccos(1/z), we wouldn't have to store these values.

        // We can evaluate
        //   \phi = \left( \frac{4 \zeta}{1-z^2} \right)^{1/4}
        // at the same time. 

        // Zeta is given by A&S 9.3.38 and 9.3.39:
        
        //

            /*
        public static double ZetaFromZ (double z) {
            if (z < 0.75) {
                double t = Math.Sqrt(1.0 - z * z);
                double zeta = Math.Pow(1.5 * Math.Log((1.0 + t) / z) - t, 2.0 / 3.0);
            } else if (z < 1.25) {
                double e = z - 1.0;
                double zeta = Math.Pow(2.0, 1.0 / 3.0) * e * (-1.0 + e * (3.0 / 10.0 + e * (-32.0 / 175.0 + e * (1037.0 / 7875.0 + e * (-103727.0 / 1010625.0 + e * 33060241.0 / 394143750.0)))));
            } else {
                double t = Math.Sqrt(z * z - 1.0);
                double zeta = -Math.Pow(1.5 * (t - Math.Acos(1.0 / t)), 2.0 / 3.0);
            }
        }
        */

        public static void ZetaFromZ (double z, out double zeta, out double phi, out double y, out double zeta_t) {

            if (z < 0.75) {
                double v = (1.0 - z) * (1.0 + z);
                double s = Math.Sqrt(v);
                y = Math.Log((1.0 + s) / z) - s;
                zeta = Math.Pow(3.0 / 2.0 * y, 2.0 / 3.0);
                phi = Math.Pow(4.0 * zeta / v, 1.0 / 4.0);
                zeta_t = zeta / v;
            } else if (z <= 1.25) {
                ZetaFromZSeries(z, out zeta, out phi, out zeta_t);
                y = 2.0 / 3.0 * Math.Pow(Math.Abs(zeta), 3.0 / 2.0);
            } else {
                double v = (z - 1.0) * (z + 1.0);
                double s = Math.Sqrt(v);
                y = s - Math.Acos(1.0 / z);
                zeta = Math.Pow(3.0 / 2.0 * y, 2.0 / 3.0);
                phi = Math.Pow(4.0 * zeta / v, 1.0 / 4.0);
                zeta_t = zeta / v;
                zeta = -zeta;
            }
        }

        public static void ZetaFromZSeries(double z, out double zeta, out double phi, out double zeta_t) {
            double e = z - 1.0;
            zeta = c[1];
            double ek = 1.0;
            for (int k = 2; k < c.Length; k++) {
                double zeta_old = zeta;
                ek *= e;
                zeta += c[k] * ek;
                if (zeta == zeta_old) {
                    phi = CbrtTwo * Math.Pow(-zeta / (1.0 + 0.5 * e), 0.25);
                    zeta_t = CbrtTwo * (-zeta) / (1.0 + z);
                    zeta = CbrtTwo * zeta * e;
                    return;
                }
            }
            throw new NonconvergenceException();
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
        }; // Last 3 are wrong!

        private static readonly double[] b0Coefficients = {
            1.0 / 70.0,
            2.0 / 225.0,
            138.0 / 67375.0,
            -296.0 / 511875.0,
            -38464.0 / 63669375.0,
            -117472.0 / 797671875.0,
            4835327648.0 / 90819588234375.0,
            196235264.0 / 3768999609375.0,
            1660313255231104.0 / 134643823706970703125.0,
            -66399236608.0 / 13534477597265625.0,
            -3601874963621281792.0 / 777568081907755810546875.0,
            -31198440677844992.0 / 28963218121581884765625.0,
            339118989445333019205632.0 / 757644590174265019544677734375.0,
            1778526038024190427136.0 / 4282935879728921209716796875.0,
            23549536309011999311770771456.0 / 246441833876147155212368990478515625.0,
            -93310845265220285037936508928.0 / 2297957529817804948349940032958984375.0,
            -2682064503630061540247431282688.0 / 71878868213876253603607622222900390625.0,
            -48526252331737774888068779474944.0 / 5687111848975905246448949416351318359375.0,
            4909761534958526661894859265996325658427392.0 / 1336710722193126209084675714361913278522491455078125.0,
            288827066125442828982554677743714304.0 / 86063906145641111265659573584224700927734375.0,
            395038580211589042747271815127586640142270464.0 / 516856005412127584113031674763288032342052459716796875.0,
            -3020453127690291145437890663054459443085312.0 / 9105227772572515212600628454363334487438201904296875.0,
            -415693811852901987488470406337119856435358859264.0 / 1376620400423164523223509095780671807415387630462646484375.0,
            -135785596536652911125363778092812910678351066219675648.0 / 1979960863899090463260579307483773619973095773875713348388671875.0,
            6752931067282482895206993683275844503156765416885321728.0 / 225597098936558707288627441254897451144546897822630405426025390625.0,
            38752618496816645289812974945044692001524843332763648.0 / 1425897741491114531154738184401894582285151277482509613037109375.0,
            1050643748528891762264542365939787384758785621457951663175397015552.0 / 170564061091910602507134320119758253348268386453085248322148621082305908203125.0,
            -198443088909204369513525850160250059674156051864171949512982528.0 / 73507181622110718234434843090222511683194777688575632981956005096435546875.0,
            -1594423240529392256942693018572814078726814644678067494840977850368.0 / 651728824538595214689523658919652300357741934836351014538295567035675048828125.0,
            -287452277097585176472798598230521879323692918524142193979974221824.0 / 519207610990370811865589760536068859381573041144862004953138530254364013671875.0,
            46158689539926545997946293538262538562766939122741043690093447245120111241429254144.0 / 189647829723697244011171944096524794232524042212212634869091488969100187900476157665252685546875.0,
            81263889690714659724831233098245316337847009260409370512461788086272.0 / 368967882186524096878157583155933589107213264279007752461438067257404327392578125.0
        };

        private static readonly double[][] aCoefficients;

        private static readonly double[][] bCoefficients;

        static BesselUae () {

            const int count = 6;
            aCoefficients = new double[count][];
            bCoefficients = new double[count][];

            aCoefficients[0] = new double[] { 1.0 };
            bCoefficients[0] = b0Coefficients;

            for (int n = 1; n < count; n++) {

                double[] bCoefficientsPrevious = bCoefficients[n - 1];

                double[] aCoefficientsNext = new double[bCoefficientsPrevious.Length - 1];

                // Determine zeroth a coefficient.
                double a0 = 0.0;
                for (int m = 1; m < n; m++) {
                    a0 += aCoefficients[m][0] * aCoefficients[n - m][0] - aCoefficients[m][1] * bCoefficients[n - 1 - m][0];
                }
                for (int m = 0; m < n; m++) {
                    a0 += aCoefficients[m][0] * bCoefficients[n - 1 - m][1];
                }
                aCoefficientsNext[0] = -0.5 * a0;
                
                // Determine a coefficients.
                for (int i = 1; i < aCoefficientsNext.Length; i++) {
                    double s = 0.0;
                    int jMax = i - 1;
                    for (int j = 0; j <= jMax; j++) {
                        s += (2 * j + 1) * b0Coefficients[j] * bCoefficientsPrevious[jMax - j];
                    }
                    aCoefficientsNext[i] = (2.0 * s - i * (i + 1) * bCoefficientsPrevious[i + 1]) / (2 * i);
                }

                // Determine b coefficients.
                double[] bCoefficientsNext = new double[aCoefficientsNext.Length - 2];
                for (int i = 0; i < bCoefficientsNext.Length; i++) {
                    double s = 0.0;
                    for (int j = 0; j <= i; j++) {
                        s += (2 * j + 1) * b0Coefficients[j] * aCoefficientsNext[i - j];
                    }
                    bCoefficientsNext[i] = (2.0 * s - (i + 1) * (i + 2) * aCoefficientsNext[i + 2]) / (2 * (2 * i + 1));
                }

                aCoefficients[n] = aCoefficientsNext;
                bCoefficients[n] = bCoefficientsNext;
            }

        }

        public static double EvaluateSeries(double[] coefficients, double x) {
            double s = coefficients[0] + coefficients[1] * x;
            double xk = x;
            for (int k = 2; k < coefficients.Length; k++) {
                double s_old = s;
                xk *= x;
                s += coefficients[k] * xk;
                if (s == s_old) return (s);
            }
            throw new NotImplementedException();
        }

        public static Tuple<double, double> Bessel_Uae_Series (double nu, double x) {

            double z = x / nu;
            ZetaFromZSeries(z, out double zeta, out double phi, out double zeta_t);
            SolutionPair p = AdvancedMath.Airy(Math.Pow(nu, 2.0 / 3.0) * zeta);

            double nu_oneThirdPower = Math.Pow(nu, 1.0 / 3.0);
            double nu_fiveThirdPower = Math.Pow(nu, 5.0 / 3.0);
            double ja = p.FirstSolutionValue / nu_oneThirdPower;
            double jb = p.FirstSolutionDerivative / nu_fiveThirdPower;
            double ya = p.SecondSolutionValue / nu_oneThirdPower;
            double yb = p.SecondSolutionDerivative / nu_fiveThirdPower;

            double a = 0.0;
            double b = 0.0;
            double J = 0.0;
            double Y = 0.0;
            double nu_power = 1.0;
            double nu_minusTwoPower = 1.0 / (nu * nu);
            double eta = zeta / CbrtTwo;
            for (int k = 0; k < 4; k++) {

                double j_old = J;
                double y_old = Y;

                double ak = k == 0 ? 1.0 : EvaluateSeries(aCoefficients[k], eta);
                a += ak * nu_power;
                double bk = CbrtTwo * EvaluateSeries(bCoefficients[k], eta);
                b += bk * nu_power;

                J = phi * (ja * a + jb * b);
                Y = -phi * (ya * a + yb * b);

                if ((J == j_old) && (Y == y_old)) {
                    return Tuple.Create(J, Y);
                }

                nu_power *= nu_minusTwoPower;
            }
            throw new NonconvergenceException();
        }

        public static double BesselJ (double nu, double x) {

            double z = x / nu;
            ZetaFromZ(z, out double zeta, out double phi, out double y, out double zeta_t);
            SolutionPair p = AdvancedMath.Airy(Math.Pow(nu, 2.0 / 3.0) * zeta);

            double an = p.FirstSolutionValue / Math.Pow(nu, 1.0 / 3.0);
            double bn = p.FirstSolutionDerivative / Math.Pow(nu, 5.0 / 3.0);

            double nui = 1.0 / (nu * nu);
            double cbrt2 = Math.Pow(2.0, 1.0 / 3.0);
            double eta = zeta / cbrt2;

            double a0 = 1.0;
            double v0a = phi * an;

            double b0 = cbrt2 * EvaluateSeries(bCoefficients[0], eta);
            double v0b = phi * (an * a0 + bn * b0);

            double a1 = EvaluateSeries(aCoefficients[1], eta);
            double v1a = phi * (an * (a0 + nui * a1) + bn * b0);

            double b1 = cbrt2 * EvaluateSeries(bCoefficients[1], eta);
            double v1b = phi * (an * (a0 + nui * a1) + bn * (b0 + nui * b1));

            double a2 = EvaluateSeries(aCoefficients[2], eta);
            double v2a = phi * (an * (a0 + nui * a1 + nui * nui * a2) + bn * (b0 + nui * b1));

            double b2 = cbrt2 * EvaluateSeries(bCoefficients[2], eta);
            double v2b = phi * (an * (a0 + nui * a1 + nui * nui * a2) + bn * (b0 + nui * b1 + nui * nui * b2));

            return (v2b);

        }

        private static Func<double, double>[] uFunctions = new Func<double, double>[] {
            t => 1.0,
            t => (3.0 - t * 5.0) / 24.0,
            t => (81.0 - t * (462.0 - t * 385.0)) / 1152.0,
            t => (30375.0 - t * (369603.0 - t * (765765.0 - t * 425425.0))) / 414720.0,
            t => (4465125.0 - t * (94121676.0 - t * (349922430.0 - t * (446185740.0 - t *  185910725.0)))) / 39813120.0,
            t => (1519035525.0 - t * (49286948607.0 - t * ( 284499769554.0 - t * (614135872350.0 - t * (566098157625.0 - t * 188699385875.0))))) / 6688604160.0,
            t => (2757049477875.0 - t * (127577298354750.0 - t * (1050760774457901.0 - t * (3369032068261860.0 - t * (5104696716244125.0 - t * (3685299006138750.0 - t * 1023694168371875.0)))))) / 4815794995200.0,
            t => (199689155040375.0 - t * (12493049053044375.0 - t * (138799253740521843.0 - t * (613221795981706275.0 - t * (1347119637570231525.0 - t * (1570320948552481125.0 - t * (931766432052080625.0 - t * 221849150488590625.0))))))) / 115579079884800.0
        };

        private static double[] mu;

        private static double[] lambda;

        private static void BesselUaePolynomialInit () {

            int nMax = 4;

            mu = new double[2 * nMax + 1];
            lambda = new double[mu.Length];

            mu[0] = 1.0;
            lambda[0] = 1.0;

            double q = 1.0;
            for (int i = 1; i < lambda.Length; i++) {
                double p = 1.0;
                int jMax = 6 * i - 1;
                for (int j = 2 * i + 1; j <= jMax; j += 2) {
                    p *= j;
                }
                q *= 144.0 * i;
                lambda[i] = p / q;
                mu[i] = -(jMax + 2.0) / jMax * lambda[i];
            }

        }

        public static double BesselJ_Uae_Polynomials (double nu, double x) {

            BesselUaePolynomialInit();

            double z = x / nu;
            ZetaFromZ(z, out double zeta, out double phi, out double y, out double zeta_t);

            SolutionPair p;
            double arg = Math.Pow(nu, 2.0 / 3.0) * zeta;
            if (arg < -32.0) {
                p = Airy_Negative_Asymptotic(arg, nu * y);
            } else {
                p = AdvancedMath.Airy(arg);
            }
            //SolutionPair p = AdvancedMath.Airy(Math.Pow(nu, 2.0 / 3.0) * zeta);

            double an = p.FirstSolutionValue / Math.Pow(nu, 1.0 / 3.0);
            double bn = p.FirstSolutionDerivative / Math.Pow(nu, 5.0 / 3.0);

            double t = 1.0 / ((1.0 - z) * (1.0 + z));
            double r = 0.5 * phi * phi; // = (\zeta t)^(1/2)

            double a0 = mu[0] * uFunctions[0](t);
            double b0 = - lambda[0] * r * MoreMath.Pow(zeta, -1) * uFunctions[1](t) - lambda[1] * MoreMath.Pow(zeta, -2) * uFunctions[0](t);

            double b0Book = -5.0 / 48.0 / (zeta * zeta) + 1.0 / Math.Sqrt(-zeta) * (5.0 / 24.0 / Math.Pow(z * z - 1.0, 1.5) + 1.0 / 8.0 / Math.Sqrt(z * z - 1.0));

            double J0 = phi * (an * a0 + bn * b0);

            double nui = 1.0 / (nu * nu);

            double a1 = mu[0] * r * r * MoreMath.Pow(zeta, -1) * uFunctions[2](t) + mu[1] * r * MoreMath.Pow(zeta, -2) * uFunctions[1](t) + mu[2] * MoreMath.Pow(zeta, -3) * uFunctions[0](t);
            double b1 = -lambda[0] * r * r * r * MoreMath.Pow(zeta, -2) * uFunctions[3](t) -
                lambda[1] * r * r * MoreMath.Pow(zeta, -3) * uFunctions[2](t) -
                lambda[2] * r * MoreMath.Pow(zeta, -4) * uFunctions[1](t) -
                lambda[3] * MoreMath.Pow(zeta, -5) * uFunctions[0](t);

            double J1 = phi * (an * (a0 + nui * a1) + bn * (b0 + nui * b1));

            return (J1);
        }

        public static Tuple<double, double> Bessel_Uae_Polynomials (double nu, double x) {

            BesselUaePolynomialInit();

            double z = x / nu;
            ZetaFromZ(z, out double zeta, out double phi, out double yy, out double zeta_t);

            SolutionPair p;
            double arg = Math.Pow(nu, 2.0 / 3.0) * zeta;
            if (arg < -32.0) {
                p = Airy_Negative_Asymptotic(-arg, nu * yy);
            } else {
                p = AdvancedMath.Airy(arg);
            }
            //SolutionPair p = AdvancedMath.Airy(Math.Pow(nu, 2.0 / 3.0) * zeta);

            double nu_oneThirdPower = Math.Pow(nu, 1.0 / 3.0);
            double nu_fiveThirdPower = Math.Pow(nu, 5.0 / 3.0);
            double ja = p.FirstSolutionValue / nu_oneThirdPower;
            double jb = p.FirstSolutionDerivative / nu_fiveThirdPower;
            double ya = p.SecondSolutionValue / nu_oneThirdPower;
            double yb = p.SecondSolutionDerivative / nu_fiveThirdPower;

            double t = 1.0 / ((1.0 - z) * (1.0 + z));
            double r = 0.5 * phi * phi;
            double nu_minusTwoPower = 1.0 / (nu * nu);

            double a = 0.0;
            double b = 0.0;
            double j = 0.0;
            double y = 0.0;
            double nu_power = 1.0;
            for (int k = 0; k < 4; k++) {

                double j_old = j;
                double y_old = y;

                double ak = 0.0;
                int iMax = 2 * k;
                for (int i = 0; i <= iMax; i++) {
                    ak += mu[i] * uFunctions[iMax - i](t) * MoreMath.Pow(r, iMax - i) * MoreMath.Pow(zeta, -(k + i));
                }
                a += ak * nu_power;

                double bk = 0.0;
                iMax += 1;
                for (int i = 0; i <= iMax; i++) {
                    bk -= lambda[i] * uFunctions[iMax - i](t) * MoreMath.Pow(r, iMax - i) * MoreMath.Pow(zeta, -(k + 1 + i));
                }
                b += bk * nu_power;

                j = phi * (ja * a + jb * b);
                y = -phi * (ya * a + yb * b);

                if ((j == j_old) && (y == y_old)) {
                    return Tuple.Create(j, y);
                }

                nu_power *= nu_minusTwoPower;
            }
            throw new NonconvergenceException();

        }

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
            ZetaFromZ(z, out zeta, out phi, out double y, out double zeta_t);
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


#if FUTURE    

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

    }
#endif

}