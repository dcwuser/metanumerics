using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;


namespace Test {

    public struct ThreeJSymbol {

        public SpinState Column1;

        public SpinState Column2;

        public SpinState Column3;

    }

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

        private static ThreeJSymbol[] GenerateRandomThreeJSymbols (double j_max, int n) {

            int tj_max = (int) Math.Truncate(2.0 * j_max);

            ThreeJSymbol[] symbols = new ThreeJSymbol[n];
            Random rng = new Random(1);
            for (int i = 0; i < n; i++) {

                int tj1 = rng.Next(tj_max + 1);
                int tj2 = rng.Next(tj_max + 1);

                int tm1 = -tj1 + 2 * rng.Next(tj1);
                int tm2 = -tj2 + 2 * rng.Next(tj2);

                int tm3 = -(tm1 + tm2);

                int tj3_min = Math.Abs(tj1 - tj2);
                int tj3_max = tj1 + tj2;
                if (Math.Abs(tm3) > tj3_min) tj3_min = Math.Abs(tm3);
                int tj3 = tj3_min + 2 * rng.Next((tj3_max - tj3_min) / 2);

                ThreeJSymbol symbol = new ThreeJSymbol();
                symbol.Column1 = new SpinState(tj1 / 2.0, tm1 / 2.0);
                symbol.Column2 = new SpinState(tj2 / 2.0, tm2 / 2.0);
                symbol.Column3 = new SpinState(tj3 / 2.0, tm3 / 2.0);
                symbols[i] = symbol;

            }

            return (symbols);

        }

        private static Spin[] GenerateRandomSpins (double j_min, double j_max, int n) {

            int tj_min = (int) Math.Truncate(2.0 * j_min);
            int tj_max = (int) Math.Truncate(2.0 * j_max);

            Spin[] spins = new Spin[n];
            Random rng = new Random(1);
            for (int i = 0; i < n; i++) {
                spins[i] = new Spin(rng.Next(tj_min, tj_max) / 2.0);
            }

            return (spins);
        }

        private static Spin[] GenerateRandomCombinedSpins (Spin j1, Spin j2, int n) {

            int tj_min = (int) (2.0 * Math.Abs(j1.J - j2.J));
            int tj_max = (int) (2.0 * (j1.J + j2.J));

            Spin[] spins = new Spin[n];
            Random rng = new Random(1);
            for (int i = 0; i < n; i++) {
                int tj = tj_min + 2 * rng.Next((tj_max - tj_min) / 2);
                spins[i] = new Spin(tj / 2.0);
            }

            return (spins);

        }

        private static SpinState[] GenerateRandomSpinStates (Spin j, int n) {

            int tj = (int) (2.0 * j.J);

            SpinState[] result = new SpinState[n];
            Random rng = new Random(1);
            for (int i = 0; i < n; i++) {

                int tm = -tj + 2 * rng.Next(tj + 1);
                result[i] = new SpinState(j, tm / 2.0);

            }

            return (result);

        }

        // generates only spinors of the same parity as j_min

        private static SpinState[] GenerateRandomSpinStates (double j_min, double j_max, int n) {

            int tj_min = (int) Math.Truncate(2.0 * j_min);
            int tj_max = (int) Math.Truncate(2.0 * j_max);

            SpinState[] states = new SpinState[n];
            Random rng = new Random(1);
            for (int i = 0; i < n; i++) {
                int tj = tj_min + 2 * rng.Next((tj_max - tj_min) / 2);
                int tm = -tj + 2 * rng.Next(tj);
                states[i] = new SpinState(tj / 2.0, tm / 2.0);
            }

            return (states);

        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SpinInvalid1 () {
            // spin indices are non-negative
            Spin s = new Spin(-1.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SpinInvalid2 () {
            // spin indicates are integer or half-integer
            Spin s = new Spin(1.25);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SpinStateInvalid1 () {
            // spin state indices are smaller than spin indices
            SpinState s = new SpinState(1.0, -2.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SpinStateInvalid2 () {
            // spin states indices match the parity of spin indices
            SpinState s = new SpinState(1.0, 0.5);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SpinStateInvalid3 () {
            // spin state values are integer or half-integer
            SpinState s = new SpinState(3.0, 0.75);
        }

        [TestMethod]
        public void ClebschGordonSepcialCase () {

            // 0 x 0 => 0

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0, 0), new SpinState(0, 0), new SpinState(0, 0)), 1.0));

            // 1/2 x 1/2 => 1

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, 0.5), new SpinState(0.5, 0.5), new SpinState(1.0, 1.0)), 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, 0.5), new SpinState(0.5, -0.5), new SpinState(1.0, 0.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, -0.5), new SpinState(0.5, 0.5), new SpinState(1.0, 0.0)), Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, -0.5), new SpinState(0.5, -0.5), new SpinState(1.0, -1.0)), 1.0));

            // 1/2 x 1/2 => 0

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, 0.5), new SpinState(0.5, -0.5), new SpinState(0.0, 0.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(0.5, -0.5), new SpinState(0.5, 0.5), new SpinState(0.0, 0.0)), -Math.Sqrt(1.0 / 2.0)));

            // 1 x 1/2 => 3/2

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(0.5, 0.5), new SpinState(1.5, 1.5)), 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(0.5, -0.5), new SpinState(1.5, 0.5)), Math.Sqrt(1.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(0.5, 0.5), new SpinState(1.5, 0.5)), Math.Sqrt(2.0 / 3.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(0.5, -0.5), new SpinState(1.5, -0.5)), Math.Sqrt(2.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(0.5, 0.5), new SpinState(1.5, -0.5)), Math.Sqrt(1.0 / 3.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(0.5, -0.5), new SpinState(1.5, -1.5)), 1.0));

            // 1 x 1/2 => 1/2

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(0.5, -0.5), new SpinState(0.5, 0.5)), Math.Sqrt(2.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(0.5, 0.5), new SpinState(0.5, 0.5)), -Math.Sqrt(1.0 / 3.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(0.5, -0.5), new SpinState(0.5, -0.5)), Math.Sqrt(1.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(0.5, 0.5), new SpinState(0.5, -0.5)), -Math.Sqrt(2.0 / 3.0)));

            // 1 x 1 => 2

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, 1.0), new SpinState(2.0, 2.0)), 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, 0.0), new SpinState(2.0, 1.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 1.0), new SpinState(2.0, 1.0)), Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, -1.0), new SpinState(2.0, 0.0)), Math.Sqrt(1.0 / 6.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 0.0), new SpinState(2.0, 0.0)), Math.Sqrt(2.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 1.0), new SpinState(2.0, 0.0)), Math.Sqrt(1.0 / 6.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, -1.0), new SpinState(2.0, -1.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 0.0), new SpinState(2.0, -1.0)), Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0,-1.0), new SpinState(1.0,-1.0), new SpinState(2.0, -2.0)), 1.0));


            // 1 x 1=> 1

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, 0.0), new SpinState(1.0, 1.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 1.0), new SpinState(1.0, 1.0)), -Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, -1.0), new SpinState(1.0, 0.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 0.0), new SpinState(1.0, 0.0)), 0.0));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 1.0), new SpinState(1.0, 0.0)), -Math.Sqrt(1.0 / 2.0)));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, -1.0), new SpinState(1.0, -1.0)), Math.Sqrt(1.0 / 2.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 0.0), new SpinState(1.0, -1.0)), -Math.Sqrt(1.0 / 2.0)));

            // 1 x 1=> 0

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 1.0), new SpinState(1.0, -1.0), new SpinState(0.0, 0.0)), Math.Sqrt(1.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, 0.0), new SpinState(1.0, 0.0), new SpinState(0.0, 0.0)), -Math.Sqrt(1.0 / 3.0)));
            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.ClebschGodron(new SpinState(1.0, -1.0), new SpinState(1.0, 1.0), new SpinState(0.0, 0.0)), Math.Sqrt(1.0 / 3.0)));

        }

        [TestMethod]
        public void ClebschGordonOrthonormalityJM () {

            int t00 = 0;
            int t01 = 0;
            int t1 = 0;

            foreach (Spin j1 in GenerateRandomSpins(0.0, 50.0, 4)) {
                foreach (Spin j2 in GenerateRandomSpins(0.0, 50.0, 4)) {

                    foreach (Spin j3a in GenerateRandomCombinedSpins(j1, j2, 4)) {
                        foreach (Spin j3b in GenerateRandomCombinedSpins(j1, j2, 4)) {
                            foreach (SpinState s3a in GenerateRandomSpinStates(j3a, 4)) {
                                foreach (SpinState s3b in GenerateRandomSpinStates(j3b, 4)) {

                                    // skip the trivial zeros; if the sums of the Ms for each coefficient
                                    // are not the same, no term can be non-zero
                                    if (s3a.M != s3b.M) continue;


                                    // sum over m1 and m2
                                    double s = 0.0;
                                    bool nonZero = false;

                                    foreach (SpinState s1 in j1.States()) {
                                        foreach (SpinState s2 in j2.States()) {

                                            double ds = SpinMath.ClebschGodron(s1, s2, s3a) * SpinMath.ClebschGodron(s1, s2, s3b);
                                            if (ds != 0.0) nonZero = true;
                                            s += ds;

                                        }
                                    }

                                    // check orthonormality
                                    if ((s3a.J == s3b.J) && (s3a.M == s3b.M)) {
                                        Assert.IsTrue(TestUtilities.IsNearlyEqual(s, 1.0));
                                        t1++;
                                    } else {
                                        Assert.IsTrue(Math.Abs(s) <= TestUtilities.TargetPrecision);
                                        if (nonZero) {
                                            t01++;
                                        } else {
                                            t00++;
                                        }
                                    }



                                }

                            }

                        }

                    }

                }
            }

            Console.WriteLine("Trivial zeros: {0}", t00);
            Console.WriteLine("Non-trivial zerios: {0}", t01);
            Console.WriteLine("Ones: {0}", t1);


        }

        [TestMethod]
        public void ClebschGordonOrthonormalityMM () {

            int t00 = 0;
            int t01 = 0;
            int t1 = 1;

            foreach (Spin j1 in GenerateRandomSpins(0.0, 50.0, 4)) {
                foreach (Spin j2 in GenerateRandomSpins(0.0, 50.0, 4)) {
                    foreach (SpinState s1a in GenerateRandomSpinStates(j1, 4)) {
                        foreach (SpinState s1b in GenerateRandomSpinStates(j1, 4)) {
                            foreach (SpinState s2a in GenerateRandomSpinStates(j2, 4)) {
                                foreach (SpinState s2b in GenerateRandomSpinStates(j2, 4)) {

                                    // if the sum of M's are equal, all terms will be trivially zero
                                    if ((s1a.M + s2a.M) != (s1b.M + s2b.M)) continue;

                                    // sum over j and m
                                    double s = 0.0;
                                    bool nonZero = false;

                                    double j_min = Math.Abs(j1.J - j2.J);
                                    double j_max = j1.J + j2.J;
                                    for (double j = j_min; j <= j_max; j = j + 1.0) {
                                        Spin j3 = new Spin(j);
                                        foreach (SpinState s3 in j3.States()) {
                                            double ds = SpinMath.ClebschGodron(s1a, s2a, s3) * SpinMath.ClebschGodron(s1b, s2b, s3);
                                            if (ds != 0.0) nonZero = true;
                                            s += ds;
                                        }
                                    }

                                    if ((s1a.M == s1b.M) && (s2a.M == s2b.M)) {
                                        Assert.IsTrue(TestUtilities.IsNearlyEqual(s, 1.0));
                                        t1++;
                                    } else {
                                        Assert.IsTrue(Math.Abs(s) < TestUtilities.TargetPrecision);
                                        if (nonZero) {
                                            t01++;
                                        } else {
                                            t00++;
                                        }
                                    }

                                }
                            }
                        }
                    }
                }
            }

            Console.WriteLine("Trivial zeros: {0}", t00);
            Console.WriteLine("Non-trivial zerios: {0}", t01);
            Console.WriteLine("Ones: {0}", t1);

        }

        [TestMethod]
        public void ThreeJExchangeSymmetry () {

            foreach (ThreeJSymbol symbol in GenerateRandomThreeJSymbols(50.0, 10)) {

                SpinState s1 = symbol.Column1;
                SpinState s2 = symbol.Column2;
                SpinState s3 = symbol.Column3;

                Console.WriteLine("( {0} {1} {2} )", s1.J, s2.J, s3.J);
                Console.WriteLine("( {0} {1} {2} )", s1.M, s2.M, s3.M);

                double s123 = SpinMath.ThreeJ(s1, s2, s3);
                double s231 = SpinMath.ThreeJ(s2, s3, s1);
                double s312 = SpinMath.ThreeJ(s3, s1, s2);

                Console.WriteLine("{0},{1},{2}", s123, s231, s312);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(s123, s231));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s123, s312));

                int P;
                if (((int) (s1.J + s2.J + s3.J)) % 2 == 0) {
                    P = 1;
                } else {
                    P = -1;
                }

                double s132 = SpinMath.ThreeJ(s1, s3, s2);
                double s213 = SpinMath.ThreeJ(s2, s1, s3);
                double s321 = SpinMath.ThreeJ(s3, s2, s1);

                Console.WriteLine("{0},{1},{2}", s132, s213, s321);


                Assert.IsTrue(TestUtilities.IsNearlyEqual(s123, P * s132));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s123, P * s213));
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s123, P * s321));

            }

        }

        [TestMethod]
        public void ThreeJRacahSymmetry () {

            foreach (ThreeJSymbol symbol in GenerateRandomThreeJSymbols(50.0, 10)) {

                // compute the 3j symbol
                SpinState s1 = symbol.Column1;
                SpinState s2 = symbol.Column2;
                SpinState s3 = symbol.Column3;
                double C = SpinMath.ThreeJ(s1, s2, s3);

                // form other 3j symbols related by Racah symmetry

                double k1, k2, k3, n1, n2, n3;
                SpinState t1, t2, t3;

                // rows 1->2->3

                k1 = (s2.J + s2.M + s3.J + s3.M) / 2.0;
                k2 = (s1.J + s1.M + s3.J + s3.M) / 2.0;
                k3 = (s1.J + s1.M + s2.J + s2.M) / 2.0;
                n1 = s1.J - s1.M - k1;
                n2 = s2.J - s2.M - k2;
                n3 = s3.J - s3.M - k3;

                t1 = new SpinState(k1, n1);
                t2 = new SpinState(k2, n2);
                t3 = new SpinState(k3, n3);
                double CC = SpinMath.ThreeJ(t1, t2, t3);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(C, CC));

                // transpose rows and columns

                /*
                k1 = s1.J;
                k2 = (s2.J + s3.J + s1.M) / 2.0;
                k3 = (s2.J + s3.J - s1.M) / 2.0;
                n1 = s2.J - s3.J;
                n2 = s3.J - s3.M - k2;
                n3 = s3.J + s3.M - k2;

                t1 = new SpinState(k1, n1);
                t2 = new SpinState(k2, n2);
                t3 = new SpinState(k3, n3);
                double CT = SpinMath.ThreeJ(t1, t2, t3);

                Assert.IsTrue(TestUtilities.IsNearlyEqual(C, CT));
                */



            }

        }

        [TestMethod]
        public void ThreeJRecursion () {

            foreach (ThreeJSymbol symbol in GenerateRandomThreeJSymbols(50.0, 10)) {

                double j1 = symbol.Column1.J;
                double m1 = symbol.Column1.M;
                double j2 = symbol.Column2.J;
                double m2 = symbol.Column2.M;
                double j3 = symbol.Column3.J;
                double m3 = symbol.Column3.M;

                double s1 = SpinMath.ThreeJ(new SpinState(j1, m1), new SpinState(j2, m2), new SpinState(j3, m3));

                double s2 = 0.0;
                if ((Math.Abs(m2-1.0) <= j2) && (Math.Abs(m3+1.0) <= j3))
                    s2 = SpinMath.ThreeJ(new SpinState(j1, m1), new SpinState(j2, m2 - 1.0), new SpinState(j3, m3 + 1.0));
                double s3 = 0.0;
                if ((Math.Abs(m1-1.0) <= j1) && (Math.Abs(m3+1.0) <= j3))
                    s3 = SpinMath.ThreeJ(new SpinState(j1, m1 - 1.0), new SpinState(j2, m2), new SpinState(j3, m3 + 1.0));

                double c1 = Math.Sqrt((j3 + m3 + 1.0) * (j3 - m3));
                double c2 = Math.Sqrt((j2 - m2 + 1.0) * (j2 + m2));
                double c3 = Math.Sqrt((j1 - m1 + 1.0) * (j1 + m1));

                Assert.IsTrue(TestUtilities.IsSumNearlyEqual(c1 * s1, c2 * s2, -c3 * s3));

            }

        }

        [TestMethod]
        public void ThreeJLegendreIntegral () {

            // pick three legendre polynomials

            foreach (int l1 in TestUtilities.GenerateUniformIntegerValues(0, 16, 8)) {
                foreach (int l2 in TestUtilities.GenerateUniformIntegerValues(0, 16, 4)) {
                    foreach (int l3 in TestUtilities.GenerateUniformIntegerValues(0, 16, 2)) {

                        Console.WriteLine("{0} {1} {2}", l1, l2, l3);

                        // do the integral over their product

                        Function<double, double> f = delegate(double x) {
                            return (
                                OrthogonalPolynomials.LegendreP(l1, x) *
                                OrthogonalPolynomials.LegendreP(l2, x) *
                                OrthogonalPolynomials.LegendreP(l3, x)
                            );
                        };
                        double I = FunctionMath.Integrate(f, Interval.FromEndpoints(-1.0, 1.0));

                        // it should be the same as 2 ( l1 l2 l3 )^2
                        //                            ( 0  0  0  )

                        double S = SpinMath.ThreeJ(new SpinState(l1, 0), new SpinState(l2, 0), new SpinState(l3, 0));

                        Console.WriteLine("  {0} {1}", I, 2.0 * S * S);

                        if (Math.Abs(S) < TestUtilities.TargetPrecision) {
                            Assert.IsTrue(I < TestUtilities.TargetPrecision);
                        } else {
                            Assert.IsTrue(TestUtilities.IsNearlyEqual(I, 2.0 * S * S));
                        }


                    }
                }
            }

        }

    }

}
