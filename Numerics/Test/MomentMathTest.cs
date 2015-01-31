using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class MomentMathTest {

        [TestMethod]
        public void MomentMathConsistency () {

            // We can't be too demanding here, since there can be strong cancelations.
            // We take a low number of simple integer cumulants and try to verify consistency.

            double[] K0 = new double[] { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
            double mu = K0[1];

            double[] C1fromK0 = MomentMath.CumulantToCentral(K0);
            Assert.IsTrue(C1fromK0[0] == 1.0);
            Assert.IsTrue(C1fromK0[1] == 0.0);
            Assert.IsTrue(C1fromK0[2] == K0[2]);

            double[] M1fromK0 = MomentMath.CumulantToRaw(K0);
            Assert.IsTrue(M1fromK0[0] == 1.0);
            Assert.IsTrue(M1fromK0[1] == K0[1]);

            double[] K2fromC1 = MomentMath.CentralToCumulant(mu, C1fromK0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(K2fromC1, K0));

            double[] M2fromC1 = MomentMath.CentralToRaw(mu, C1fromK0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(M2fromC1, M1fromK0));

            double[] K2fromM1 = MomentMath.RawToCumulant(M1fromK0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(K2fromM1, K0));

            double[] C2fromM1 = MomentMath.RawToCentral(M1fromK0);
            Assert.IsTrue(TestUtilities.IsNearlyEqual(C2fromM1, C1fromK0));

        }

    }
}
