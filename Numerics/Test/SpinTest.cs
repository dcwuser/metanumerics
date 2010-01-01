using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Functions;


namespace Test {

    public struct SpinRange {

        public Spin Minimum;

        public Spin Maximum;

    }

    public struct ThreeJSymbol {

        public SpinState Column1;

        public SpinState Column2;

        public SpinState Column3;

    }

    public struct SixJSymbol {

        public Spin J1;

        public Spin J2;

        public Spin J3;

        public Spin J4;

        public Spin J5;

        public Spin J6;

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

        private static SixJSymbol[] GenerateRandomSixJSymbols (double j_max, int n) {

            int tj_max = (int) Math.Truncate(2.0 * j_max);

            SixJSymbol[] symbols = new SixJSymbol[n];
            Random rng = new Random(1);

            int i = 0;
            while (i < n) {
            //for (int i = 0; i < n; i++) {

                // pick 1, 4, and 2 randomly
                int tj1 = rng.Next(tj_max + 1);
                int tj4 = rng.Next(tj_max + 1);
                int tj2 = rng.Next(tj_max + 1);

                // make sure 5 has appropriate wholess so that 1+2 has same wholeness as 4+5
                int tj5 = rng.Next(tj_max + 1);
                if (((tj1 + tj2) % 2) != ((tj4 + tj5) % 2)) tj5++;

                // make sure 3 can be formed from 1+2 or 4+5

                int tj12_min = Math.Abs(tj1 - tj2);
                int tj45_min = Math.Abs(tj4 - tj5);
                int tj3_min;
                if (tj12_min > tj45_min) {
                    tj3_min = tj12_min;
                } else {
                    tj3_min = tj45_min;
                }

                int tj12_max = tj1 + tj2;
                int tj45_max = tj4 + tj5;
                int tj3_max;
                if (tj12_max < tj45_max) {
                    tj3_max = tj12_max;
                } else {
                    tj3_max = tj45_max;
                }

                if (tj3_min > tj3_max) continue;

                int tj3 = tj3_min + 2 * rng.Next((tj3_max - tj3_min) / 2 + 1);

                // make sure 6 can be formed from 1+5 or 4+2

                int tj15_min = Math.Abs(tj1 - tj5);
                int tj42_min = Math.Abs(tj4 - tj2);
                int tj6_min;
                if (tj15_min > tj42_min) {
                    tj6_min = tj15_min;
                } else {
                    tj6_min = tj42_min;
                }

                int tj15_max = tj1 + tj5;
                int tj42_max = tj4 + tj2;
                int tj6_max;
                if (tj15_max < tj42_max) {
                    tj6_max = tj15_max;
                } else {
                    tj6_max = tj42_max;
                }

                if (tj6_min > tj6_max) continue;

                int tj6 = tj6_min + 2 * rng.Next((tj6_max - tj6_min) / 2 + 1);

                // make the symbol
                SixJSymbol symbol = new SixJSymbol();
                symbol.J1 = new Spin(tj1 / 2.0);
                symbol.J2 = new Spin(tj2 / 2.0);
                symbol.J3 = new Spin(tj3 / 2.0);
                symbol.J4 = new Spin(tj4 / 2.0);
                symbol.J5 = new Spin(tj5 / 2.0);
                symbol.J6 = new Spin(tj6 / 2.0);

                symbols[i] = symbol;
                i++;

            }

            return (symbols);
        }


        private static SpinRange CombinedSpinRange (Spin j1, Spin j2) {

            int tj1 = (int) Math.Round(2 * j1.J);
            int tj2 = (int) Math.Round(2 * j2.J);

            int tj_min = Math.Abs(tj1 - tj2);
            int tj_max = tj1 + tj2;

            SpinRange range = new SpinRange();
            range.Minimum = new Spin(tj_min / 2.0);
            range.Maximum = new Spin(tj_max / 2.0);

            return (range);

        }

        private static SpinRange CombinedSpinRange (SpinRange r1, SpinRange r2) {

            SpinRange range = new SpinRange();

            if (r1.Minimum.J > r2.Minimum.J) {
                range.Minimum = r1.Minimum;
            } else {
                range.Minimum = r2.Minimum;
            }

            if (r1.Maximum.J < r2.Maximum.J) {
                range.Maximum = r1.Maximum;
            } else {
                range.Maximum = r2.Maximum;
            }

            return (range);

        }

        private static Spin[] AllSpinsInRange (SpinRange range) {

            int tj_min = (int) Math.Round(2 * range.Minimum.J);
            int tj_max = (int) Math.Round(2 * range.Maximum.J);

            List<Spin> spins = new List<Spin>();
            for (int tj = tj_min; tj <= tj_max; tj+=2) {
                spins.Add(new Spin(tj / 2.0));
            }

            return (spins.ToArray());
        }

        private static Spin[] RandomSpinsInRange (SpinRange range, int n) {

            int tj_min = (int) Math.Truncate(2.0 * range.Minimum.J);
            int tj_max = (int) Math.Truncate(2.0 * range.Maximum.J);

            Spin[] spins = new Spin[n];
            Random rng = new Random(1);
            for (int i = 0; i < n; i++) {
                spins[i] = new Spin((tj_min + 2 * rng.Next((tj_max - tj_min) / 2 + 1))/2.0);
            }

            return (spins);

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

        private static Spin[] GenerateCombinedSpins (Spin j1, Spin j2) {

            int tj_min = (int) (2.0 * Math.Abs(j1.J - j2.J));
            int tj_max = (int) (2.0 * (j1.J + j2.J));

            List<Spin> spins = new List<Spin>();
            for (int tj = tj_min; tj <= tj_max; tj += 2) {
                spins.Add(new Spin(tj / 2.0));
            }

            return (spins.ToArray());

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

        [TestMethod]
        public void SixJSpecialCase () {

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.SixJ(Spin.SpinZero, Spin.SpinZero, Spin.SpinZero, Spin.SpinZero, Spin.SpinZero, Spin.SpinZero), 1.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.SixJ(Spin.SpinZero, Spin.SpinZero, Spin.SpinZero, Spin.SpinOneHalf, Spin.SpinOneHalf, Spin.SpinOneHalf), - 1.0 / Math.Sqrt(2.0) ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.SixJ(Spin.SpinOneHalf, Spin.SpinOneHalf, Spin.SpinZero, Spin.SpinOneHalf, Spin.SpinOneHalf, Spin.SpinZero), - 1.0 / 2.0)); 

        }

        [TestMethod]
        public void SixJReggeSymmetry () {

            SixJSymbol[] symbols = GenerateRandomSixJSymbols(15.0, 15);
            foreach (SixJSymbol symbol in symbols) {

                Spin sa1 = symbol.J1;
                Spin sa2 = symbol.J2;
                Spin sa3 = symbol.J3;
                Spin sa4 = symbol.J4;
                Spin sa5 = symbol.J5;
                Spin sa6 = symbol.J6;
                
                Console.WriteLine("{0} {1} {2}", sa1.J, sa2.J, sa3.J);
                Console.WriteLine("{0} {1} {2}", sa4.J, sa5.J, sa6.J);

                double sa = SpinMath.SixJ(sa1, sa2, sa3, sa4, sa5, sa6);

                Spin sb1 = new Spin(sa1.J);
                Spin sb2 = new Spin((sa3.J + sa5.J + sa6.J - sa2.J) / 2.0);
                Spin sb3 = new Spin((sa2.J + sa5.J + sa6.J - sa3.J) / 2.0);
                Spin sb4 = new Spin(sa4.J);
                Spin sb5 = new Spin((sa2.J + sa3.J + sa6.J - sa5.J) / 2.0);
                Spin sb6 = new Spin((sa2.J + sa3.J + sa5.J - sa6.J) / 2.0);

                Console.WriteLine("{0} {1} {2}", sb1.J, sb2.J, sb3.J);
                Console.WriteLine("{0} {1} {2}", sb4.J, sb5.J, sb6.J);

                double sb = SpinMath.SixJ(sb1, sb2, sb3, sb4, sb5, sb6);

                Console.WriteLine("{0} vs. {1}", sa, sb);
                Console.WriteLine("---");

                Assert.IsTrue(TestUtilities.IsNearlyEqual(sa, sb));

            }


        }


        [TestMethod]
        public void SixJOrthonormality () {

            SixJSymbol[] sas = GenerateRandomSixJSymbols(50.0, 5);
            foreach (SixJSymbol sa in sas) {

                Spin j1 = sa.J1;
                Spin j2 = sa.J2;
                Spin j4 = sa.J4;
                Spin j5 = sa.J5;

                Spin j6a = sa.J6;

                SpinRange r15 = CombinedSpinRange(j1, j5);
                SpinRange r42 = CombinedSpinRange(j4, j2);
                SpinRange r6b = CombinedSpinRange(r15, r42);
                Spin[] j6bs = RandomSpinsInRange(r6b, 5);
                j6bs[0] = j6a;

                SpinRange r12 = CombinedSpinRange(j1, j2);
                SpinRange r45 = CombinedSpinRange(j4, j5);
                SpinRange r3 = CombinedSpinRange(r12, r45);
                Spin[] j3s = AllSpinsInRange(r3);

                foreach (Spin j6b in j6bs) {

                    double s = 0.0;
                    foreach (Spin j3 in j3s) {

                        double t = (2.0 * j3.J + 1.0) *
                            SpinMath.SixJ(j1, j2, j3, j4, j5, j6a) *
                            SpinMath.SixJ(j1, j2, j3, j4, j5, j6b);
                        s += t;

                    }

                    Console.WriteLine("j1={0} j2={1} j4={2} j5={3} j6a={4} j6b={5} s={6}", j1.J, j2.J, j4.J, j5.J, j6a.J, j6b.J, s);

                    if (j6a == j6b) {
                        Assert.IsTrue(TestUtilities.IsNearlyEqual(s, 1.0 / (2.0 * j6a.J + 1.0)));
                    } else {
                        Assert.IsTrue(Math.Abs(s) < TestUtilities.TargetPrecision);
                    }


                }

            }

        }

        //                         { a b x }
        // \sum_x (-1)^{2x} (2x+1) { a b f } = 1

        [TestMethod]
        public void SixJSum () {

            Spin[] spins = GenerateRandomSpins(0.0, 50.0, 5);
            foreach (Spin a in spins) {
                foreach (Spin b in spins) {

                    Spin[] combined_spins = GenerateRandomCombinedSpins(a, b, 5);
                    foreach (Spin f in combined_spins) {

                        double s = 0.0;
                        Spin[] xs = GenerateCombinedSpins(a, b);
                        foreach (Spin x in xs) {

                            double t = (2 * x.J + 1) * SpinMath.SixJ(a, b, x, a, b, f);
                            if ((((int) Math.Round(2 * x.J)) % 2) != 0) t = -t;
                            s += t;
                        }

                        Assert.IsTrue(TestUtilities.IsNearlyEqual(s, 1.0));

                    }


                }
            }
        }

        // \sum_{x} (-1)^{f+g+x} (2x+1) { a b x } { c d x } = { a d f }
        //                              { c d f } { b a g }   { b c g }

        [TestMethod]
        public void SixJProductSum () {

            SixJSymbol[] symbols = GenerateRandomSixJSymbols(50.0, 5);

            foreach (SixJSymbol symbol in symbols) {

                Spin a = symbol.J1;
                Spin b = symbol.J2;
                Spin c = symbol.J4;
                Spin d = symbol.J5;

                Spin f = symbol.J6;

                SpinRange ac = CombinedSpinRange(a, c);
                Console.WriteLine("ac: [{0},{1}]", ac.Minimum.J, ac.Maximum.J);
                SpinRange bd = CombinedSpinRange(b, d);
                Console.WriteLine("bd: [{0},{1}]", bd.Minimum.J, bd.Maximum.J);
                SpinRange ac_bd = CombinedSpinRange(ac, bd);
                Console.WriteLine("g: [{0},{1}]", ac_bd.Minimum.J, ac_bd.Maximum.J);
                Spin[] gs = RandomSpinsInRange(ac_bd, 3);

                SpinRange ab = CombinedSpinRange(a, b);
                SpinRange cd = CombinedSpinRange(c, d);
                SpinRange ab_cd = CombinedSpinRange(ab, cd);
                Spin[] xs = AllSpinsInRange(ab_cd);

                foreach (Spin g in gs) {

                    //double s = 0.0;
                    List<double> s = new List<double>();
                    foreach (Spin x in xs) {

                        double t = (2.0 * x.J + 1.0) *
                            SpinMath.SixJ(a, b, x, c, d, f) *
                            SpinMath.SixJ(c, d, x, b, a, g);
                        int p = (int) Math.Round(f.J + g.J + x.J);
                        if (p % 2 != 0) t = -t;
                        Console.WriteLine(" {0} {1}", x.J, t);
                        //s += t;
                        s.Add(t);

                    }

                    double r = SpinMath.SixJ(a, d, f, b, c, g);

                    Console.WriteLine("a={0} b={1} c={2} d={3} f={4} g={5}: {6} {7}", a.J, b.J, c.J, d.J, f.J, g.J, s, r);

                    Assert.IsTrue(TestUtilities.IsSumNearlyEqual(s, r));
                    //Assert.IsTrue(TestUtilities.IsNearlyEqual(s, r));

                }

            }

        }

        // ( a b c ) { a b c }
        // ( A B C ) { d e f } = \sum_{D, E, F} (-1)^{d + e + f + D + E + F}
        //
        // ( d  e c ) ( e  f a ) ( f  d b )
        // ( D -E C ) ( E -F A ) ( F -D B )

        [TestMethod]
        public void SixJThreeJRelation () {

            SixJSymbol[] symbols = GenerateRandomSixJSymbols(50.0, 5);
            foreach (SixJSymbol symbol in symbols) {

                Spin a = symbol.J1;
                Spin b = symbol.J2;
                Spin c = symbol.J3;
                Spin d = symbol.J4;
                Spin e = symbol.J5;
                Spin f = symbol.J6;

                SpinState[] ams = GenerateRandomSpinStates(a, 2);
                SpinState[] bms = GenerateRandomSpinStates(b, 2);
                foreach (SpinState am in ams) {
                    foreach (SpinState bm in bms) {

                        if (Math.Abs(am.M + bm.M) > c.J) continue;
                        SpinState cm = new SpinState(c, -(am.M + bm.M));

                        double g1 = SpinMath.ThreeJ(am, bm, cm);
                        double g2 = SpinMath.SixJ(a, b, c, d, e, f);
                        double p = g1 * g2;

                        double q = 0.0;
                        List<double> ts = new List<double>();
                        SpinState[] dms = d.States();
                        foreach (SpinState dm in dms) {

                            if (Math.Abs(dm.M + cm.M) > e.J) continue;
                            SpinState em = new SpinState(e, dm.M + cm.M);
                            SpinState mem = new SpinState(e, -em.M);
                            double f1 = SpinMath.ThreeJ(dm, mem, cm);

                            if (Math.Abs(em.M + am.M) > f.J) continue;
                            SpinState fm = new SpinState(f, em.M + am.M);
                            SpinState mfm = new SpinState(f, -fm.M);
                            double f2 = SpinMath.ThreeJ(em, mfm, am);

                            SpinState mdm = new SpinState(d, -dm.M);
                            double f3 = SpinMath.ThreeJ(fm, mdm, bm);

                            double t = f1 * f2 * f3;

                            int s = (int) Math.Round(dm.J + dm.M + em.J + em.M + fm.J + fm.M);
                            if (s % 2 != 0) t = -t;

                            q += t;
                            ts.Add(t);

                        }

                        Console.WriteLine("{0} v. {1}", p, q);
                        //Assert.IsTrue(TestUtilities.IsNearlyEqual(p, q));
                        Assert.IsTrue(TestUtilities.IsSumNearlyEqual(ts, p));
                    }
                }

            }

        }

        [TestMethod]
        public void SixJExchangeSymmetry () {

            SixJSymbol[] symbols = GenerateRandomSixJSymbols(50.0, 5);
            foreach (SixJSymbol symbol in symbols) {

                Spin j1 = symbol.J1;
                Spin j2 = symbol.J2;
                Spin j3 = symbol.J3;
                Spin j4 = symbol.J4;
                Spin j5 = symbol.J5;
                Spin j6 = symbol.J6;

                double s1 = SpinMath.SixJ(j1, j2, j3, j4, j5, j6);

                // odd column permutation
                double s2 = SpinMath.SixJ(j2, j1, j3, j5, j4, j6);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s1, s2));

                // even column permutation
                double s3 = SpinMath.SixJ(j2, j3, j1, j5, j6, j4);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s1, s3));

                // flip
                double s4 = SpinMath.SixJ(j1, j5, j6, j4, j2, j3);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(s1, s4));

            }

        }

        [TestMethod]
        public void SixJZeroMinimum () {

            // we include these special cases in order to exercise code in the recursion routine that is active
            // only when the minimum spin of the recursion is zero

            Spin s1 = new Spin(1.0);
            Spin s2 = new Spin(2.0);
            Spin s3 = new Spin(3.0);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.SixJ(s2, s2, s2, s1, s2, s2), -1.0 / 10.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.SixJ(s2, s2, s2, s2, s2, s2), -3.0 / 70.0));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(SpinMath.SixJ(s2, s2, s2, s3, s2, s2), 4.0 / 35.0));

        }

        /*
        [TestMethod]
        public void SixJZeroMinimumTest () {

            Spin sa1 = new Spin(1.5);
            Spin sa2 = new Spin(1.5);
            Spin sa3 = new Spin(3.0);
            Spin sa4 = new Spin(3.0);
            Spin sa5 = new Spin(3.0);
            Spin sa6 = new Spin(1.5);

            double c = SpinMath.SixJ_ShultenGorton_Recurse(4, 6, 6, 6, 4, 4);
            //double c = SpinMath.SixJ(sa1, sa2, sa3, sa4, sa5, sa6);
            Console.WriteLine(c);


        }
        */

    }

}
