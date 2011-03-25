using System;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.SignalProcessing;

namespace Test {
    [TestClass]
    public class FourierTest {

        // Fibinacci numbers have a good mix of prime factorizations

        public static int[] sizes = new int[] { 2, 3, 5, 8, 13, 21, 34, 55, 89, 144 };

        [TestMethod]
        public void FourierTiming () {

            int n = 10020;

            FourierTransformer ft = new FourierTransformer(n);

            Complex[] x = new Complex[n];
            x[1] = 1.0;

            Stopwatch t = Stopwatch.StartNew();
            Complex[] xt = ft.Transform(x);
            t.Stop();

            Console.WriteLine("{0} {1}", n, t.ElapsedMilliseconds);

        }

        [TestMethod]
        public void FourierSpecialCases () {

            for (int n = 2; n <= 10; n++) {

                FourierTransformer ft = new FourierTransformer(n);
                Assert.IsTrue(ft.Length == n);

                // transform of uniform is zero frequency only
                Complex[] x1 = new Complex[n];
                for (int i = 0; i < n; i++) {
                    x1[i] = 1.0;
                }
                Complex[] y1 = ft.Transform(x1);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(y1[0], n));
                for (int i = 1; i < n; i++) {
                    Assert.IsTrue(ComplexMath.Abs(y1[i]) < TestUtilities.TargetPrecision);
                }
                Complex[] z1 = ft.InverseTransform(y1);
                for (int i = 0; i < n; i++) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(z1[i], 1.0));
                }

                // transform of pulse at 1 are nth roots of unity, read clockwise
                Complex[] x2 = new Complex[n];
                x2[1] = 1.0;
                Complex[] y2 = ft.Transform(x2);
                for (int i = 0; i < n; i++) {
                    double t = -2.0 * Math.PI / n * i;
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(y2[i], new Complex(Math.Cos(t), Math.Sin(t))));
                }

            }

        }

        [TestMethod]
        public void FourierInverse () {

            foreach (int n in sizes) {
                Console.WriteLine(n);
                
                FourierTransformer ft = new FourierTransformer(n);

                Complex[] x = TestUtilities.GenerateComplexValues(0.1, 10.0, n);

                Complex[] y = ft.Transform(x);
                Complex[] z = ft.InverseTransform(y);

                for (int i = 0; i < n; i++) {
                    Assert.IsTrue(TestUtilities.IsNearlyEqual(x[i], z[i]));
                }

            }

        }

        [TestMethod]
        public void FourierLinearity () {

            foreach (int n in sizes) {
                
                Console.WriteLine(n);
                FourierTransformer ft = new FourierTransformer(n);
                Random rng = new Random(1);

                Complex[] x = TestUtilities.GenerateComplexValues(0.1, 10.0, n, rng);
                Complex[] y = TestUtilities.GenerateComplexValues(0.1, 10.0, n, rng);
                Complex[] z = new Complex[n];
                for (int i = 0; i < n; i++) {
                    z[i] = x[i] + y[i];
                }

                Complex[] xt = ft.Transform(x);
                Complex[] yt = ft.Transform(y);
                Complex[] zt = ft.Transform(z);

                for (int i = 0; i < n; i++) {
                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(xt[i], yt[i], zt[i]));
                }

            }

        }

        [TestMethod]
        public void FourierParseval () {

            foreach (int n in sizes) {
                //int n = 15;
                Console.WriteLine(n);
                FourierTransformer ft = new FourierTransformer(n);
                Random rng = new Random(1);

                Complex[] x = TestUtilities.GenerateComplexValues(0.1, 10.0, n, rng);
                Complex[] y = TestUtilities.GenerateComplexValues(0.1, 10.0, n, rng);

                Complex s = 0.0;
                for (int i = 0; i < n; i++) {
                    s += x[i] * y[i].Conjugate;
                }

                Complex[] xt = ft.Transform(x);
                Complex[] yt = ft.Transform(y);

                Complex st = 0.0;
                for (int i = 0; i < n; i++) {
                    st += xt[i] * yt[i].Conjugate;
                }
                st /= n;

                Console.WriteLine("{0} {1}", s, st);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s, st));

            }

        }

        [TestMethod]
        public void TestMethod1 () {

            int N = 49 * 2;
            FourierTransformer ft = new FourierTransformer(N);

            Complex[] x = new Complex[N];
            x[1] = 1.0;

            Complex[] y = ft.Transform(x);

            Complex[] z = ft.InverseTransform(y);
            //Complex[] z = new Complex[N];

            for (int i = 0; i < x.Length; i++) {
                double t = - 2.0 * Math.PI / x.Length * i;
                Console.WriteLine("{0} -> {1} {2} -> {3}", x[i], y[i], new Complex(Math.Cos(t), Math.Sin(t)), z[i]);
            }

        }

    }
}
