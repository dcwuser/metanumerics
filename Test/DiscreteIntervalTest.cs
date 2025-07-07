using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using FluentAssertions;

namespace Test {


    [TestClass]
    public class DiscreteIntervalTest {

        [TestMethod]
        public void DiscreteIntervalEquality () {

            DiscreteInterval a = DiscreteInterval.FromEndpoints(1, 2);
            DiscreteInterval aPrime = DiscreteInterval.FromEndpoints(2, 1);
            DiscreteInterval b = DiscreteInterval.FromEndpoints(1, 3);

            (a == aPrime).Should().BeTrue();
            (a != aPrime).Should().BeFalse();
            a.Equals(aPrime).Should().BeTrue();
            a.Equals((object)aPrime).Should().BeTrue();
            DiscreteInterval.Equals(a, aPrime).Should().BeTrue();
            a.GetHashCode().Should().Be(aPrime.GetHashCode());

            (a == b).Should().BeFalse();
            (a != b).Should().BeTrue();
            a.Equals(b).Should().BeFalse();
            a.Equals((object)b).Should().BeFalse();
            DiscreteInterval.Equals(a, b).Should().BeFalse();
        }

    }

}
