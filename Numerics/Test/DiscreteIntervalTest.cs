using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics.Distributions;

namespace Test {
    [TestClass]
    public class DiscreteIntervalTest {

        [TestMethod]
        public void DiscreteIntervalEquality () {

            DiscreteInterval i = DiscreteInterval.FromEndpoints(1, 2);
            DiscreteInterval j = DiscreteInterval.FromEndpoints(1, 3);

            Assert.IsFalse(i == j); Assert.IsTrue(i != j);

            DiscreteInterval ic = i;

            Assert.IsTrue(i == ic); Assert.IsFalse(i != ic);
            Assert.IsTrue(i.GetHashCode() == ic.GetHashCode());
        }

    }
}
