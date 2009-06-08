using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;


namespace Test {

#if FUTURE

    /// <summary>
    /// Summary description for SpinTest
    /// </summary>
    [TestClass]
    public class SpinTest {
        public SpinTest () {
            //
            // TODO: Add constructor logic here
            //
        }

        private TestContext testContextInstance;

        /// <summary>
        ///Gets or sets the test context which provides
        ///information about and functionality for the current test run.
        ///</summary>
        public TestContext TestContext {
            get {
                return testContextInstance;
            }
            set {
                testContextInstance = value;
            }
        }

    #region Additional test attributes
        //
        // You can use the following additional attributes as you write your tests:
        //
        // Use ClassInitialize to run code before running the first test in the class
        // [ClassInitialize()]
        // public static void MyClassInitialize(TestContext testContext) { }
        //
        // Use ClassCleanup to run code after all tests in a class have run
        // [ClassCleanup()]
        // public static void MyClassCleanup() { }
        //
        // Use TestInitialize to run code before running each test 
        // [TestInitialize()]
        // public void MyTestInitialize() { }
        //
        // Use TestCleanup to run code after each test has run
        // [TestCleanup()]
        // public void MyTestCleanup() { }
        //
        #endregion

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SpinInvalid1 () {
            Spin s = new Spin(-1.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SpinInvalid2 () {
            Spin s = new Spin(1.25);
        }


        [TestMethod]
        public void SpinCast () {
            Spin s = 0.5;
        }

        [TestMethod]
        public void ClebschGordonSepcialCase () {

            // 0 x 0 => 0

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0, 0), new SpinState(0, 0), new SpinState(0, 0)), 1.0));

            // 1/2 x 1/2

            // => 1

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, 0.5), new SpinState(0.5, 0.5), new SpinState(1.0, 1.0)), 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, 0.5), new SpinState(0.5, -0.5), new SpinState(1.0, 0.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, -0.5), new SpinState(0.5, 0.5), new SpinState(1.0, 0.0)), Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, -0.5), new SpinState(0.5, -0.5), new SpinState(1.0, -1.0)), 1.0));

            // => 0

            // 1 x 1/2

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, 0.5), new SpinState(0.5, -0.5), new SpinState(0.0, 0.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, -0.5), new SpinState(0.5, 0.5), new SpinState(0.0, 0.0)), -Math.Sqrt(1.0 / 2.0)));

            // => 3/2

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(0.5, 0.5), new SpinState(1.5, 1.5)), 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(0.5, -0.5), new SpinState(1.5, 0.5)), Math.Sqrt(1.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(0.5, 0.5), new SpinState(1.5, 0.5)), Math.Sqrt(2.0 / 3.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(0.5, -0.5), new SpinState(1.5, -0.5)), Math.Sqrt(2.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(0.5, 0.5), new SpinState(1.5, -0.5)), Math.Sqrt(1.0 / 3.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(0.5, -0.5), new SpinState(1.5, -1.5)), 1.0));

            // => 1/2

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(0.5, -0.5), new SpinState(0.5, 0.5)), Math.Sqrt(2.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(0.5, 0.5), new SpinState(0.5, 0.5)), -Math.Sqrt(1.0 / 3.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(0.5, -0.5), new SpinState(0.5, -0.5)), Math.Sqrt(1.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(0.5, 0.5), new SpinState(0.5, -0.5)), -Math.Sqrt(2.0 / 3.0)));

            // 1 x 1

            // => 2

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, 1.0), new SpinState(2.0, 2.0)), 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, 0.0), new SpinState(2.0, 1.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 1.0), new SpinState(2.0, 1.0)), Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, -1.0), new SpinState(2.0, 0.0)), Math.Sqrt(1.0 / 6.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 0.0), new SpinState(2.0, 0.0)), Math.Sqrt(2.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 1.0), new SpinState(2.0, 0.0)), Math.Sqrt(1.0 / 6.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, -1.0), new SpinState(2.0, -1.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 0.0), new SpinState(2.0, -1.0)), Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0,-1.0), new SpinState(1.0,-1.0), new SpinState(2.0, -2.0)), 1.0));


            // => 1

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, 0.0), new SpinState(1.0, 1.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 1.0), new SpinState(1.0, 1.0)), -Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, -1.0), new SpinState(1.0, 0.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 0.0), new SpinState(1.0, 0.0)), 0.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 1.0), new SpinState(1.0, 0.0)), -Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, -1.0), new SpinState(1.0, -1.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 0.0), new SpinState(1.0, -1.0)), -Math.Sqrt(1.0 / 2.0)));

            // => 0

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, -1.0), new SpinState(0.0, 0.0)), Math.Sqrt(1.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 0.0), new SpinState(0.0, 0.0)), -Math.Sqrt(1.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 1.0), new SpinState(0.0, 0.0)), Math.Sqrt(1.0 / 3.0)));

        }

        [TestMethod]
        public void DebugTest () {
            Console.WriteLine(SpinMath.ThreeJ(new SpinState(12, 0), new SpinState(8, 0), new SpinState(14, 0)));
            Console.WriteLine(SpinMath.ThreeJ(new SpinState(11, 1), new SpinState(13, -5), new SpinState(10, 4)));
            //Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, 1.0), new SpinState(2.0, 2.0)), 1.0));
        }

        public void ThreeJExchangeSymmetry () {

        }

        [TestMethod]
        public void ThreeJOrthonormalityJM () {

            // pick two random spins to combine
            int tj1 = 4;
            int tj2 = 3;

            // loop over all possible values of j and j' and m and m'
            for (int tja = Math.Abs(tj1-tj2); tja <= (tj1+tj2); tja=tja+2) {
                for (int tjb = Math.Abs(tj1-tj2); tjb <= tja; tjb=tjb+2) {
                    for (int tma = -tja; tma <= tja; tma += 2) {
                        for (int tmb = -tjb; tmb <= tjb; tmb += 2) {

                            SpinState sa = new SpinState(tja / 2.0, tma / 2.0);
                            SpinState sb = new SpinState(tjb / 2.0, tmb / 2.0);

                            Console.WriteLine("tja = {0} tjb = {1} tma = {2} tmb = {3}:", tja, tjb, tma, tmb); 

                            // in each case, sum over all possible values of m1 and m2
                            List<double> t = new List<double>();
                            for (int tm1 = -tj1; tm1 <= tj1; tm1 = tm1 + 2) {
                                for (int tm2 = -tj2; tm2 <= tj2; tm2 = tm2 + 2) {
                                    SpinState s1 = new SpinState(tj1 / 2.0, tm1 / 2.0);
                                    SpinState s2 = new SpinState(tj2 / 2.0, tm2 / 2.0);

                                    double tt = SpinMath.ThreeJ(s1, s2, sa) * SpinMath.ThreeJ(s1, s2, sb);
                                    t.Add(tt);

                                    Console.WriteLine("tm1 = {0} tm2 = {1} tt = {2}", tm1, tm2, tt);

                                }
                            }

                            if ((tja == tjb) && (tma == tmb)) {
                                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(t, 1.0/(tja+1)));
                            } else {
                                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(t, 0.0));
                            }


                        }
                    }


                }
            }

        }

        public void ThreeJOrthonormalityMM () {
        }

        [TestMethod]
        public void ThreeJRacahSymmetry () {
            
            int tj1, tj2, tj3, tm1, tm2, tm3;

            tj1 = 8*2;
            tj2 = 10*2;
            tj3 = 12*2;
            tm1 = -3*2;
            tm2 = 2*2;
            tm3 = 1*2;

            //for (int i = 0; i < 5; i++) {

                SpinState s1 = new SpinState(tj1 / 2.0, tm1 / 2.0);
                SpinState s2 = new SpinState(tj2 / 2.0, tm2 / 2.0);
                SpinState s3 = new SpinState(tj3 / 2.0, tm3 / 2.0);
                double C = SpinMath.ThreeJ(s1, s2, s3);

                int tk1, tk2, tk3, tn1, tn2, tn3;
                SpinState t1, t2, t3;

                // rows 2<->3

                // rows 1->2->3

                tk1 = (tj2 + tm2 + tj3 + tm3) / 2;
                tk2 = (tj1 + tm1 + tj3 + tm3) / 2;
                tk3 = (tj1 + tm1 + tj2 + tm2) / 2;
                tn1 = tj1 - tm1 - tk1;
                tn2 = tj2 - tm2 - tk2;
                tn3 = tj3 - tm3 - tk3;

                t1 = new SpinState(tj1 / 2.0, tm1 / 2.0);
                t2 = new SpinState(tj2 / 2.0, tm2 / 2.0);
                t3 = new SpinState(tj3 / 2.0, tm3 / 2.0);
                double CC = SpinMath.ThreeJ(t1, t2, t3);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(C, CC));

                // transpose rows and columns

                tk1 = tj1;
                tk2 = (tj2 + tj3 + tm1) / 2;
                tk3 = (tj2 + tj3 - tm1) / 2;
                tn1 = tj2 - tj3;
                tn2 = tj3 - tm3 - tk2;
                tn3 = tj3 + tm3 - tk2;

                t1 = new SpinState(tj1 / 2.0, tm1 / 2.0);
                t2 = new SpinState(tj2 / 2.0, tm2 / 2.0);
                t3 = new SpinState(tj3 / 2.0, tm3 / 2.0);
                double CT = SpinMath.ThreeJ(t1, t2, t3);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(C, CT));

            //}

        }

        public void ThreeJRecursion () {

            /*
            SpinState jm1;
            SpinState jm2;
            SpinState jm3;

            SpinState jm1m1;
            SpinState jm2m1;
            SpinState jm3p1;

            SpinMath.ThreeJ(jm1, jm2, jm3) + SpinMath.ThreeJ(jm1m1, jm2, jm3p1) + SpinMath.ThreeJ(jm1, jm2m1, jm3p1) == 0;
            */

        }

        [TestMethod]
        public void ThreeJDebug () {

            Console.WriteLine(Stopwatch.IsHighResolution);
            Stopwatch watch = Stopwatch.StartNew();
            double c = SpinMath.ThreeJ(new SpinState(400.0, 0.0), new SpinState(600.0, 0.0), new SpinState(800.0, 0.0));
            watch.Stop();
            Console.WriteLine(c);
            Console.WriteLine(watch.ElapsedMilliseconds);

        }


    }

#endif

}
