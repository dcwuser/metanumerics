using System;
using FluentAssertions;
using Meta.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {
    
    
    /// <summary>
    ///This is a test class for IntervalTest and is intended
    ///to contain all IntervalTest Unit Tests
    ///</summary>
    [TestClass]
    public class IntervalTest {

        private double a = -3.1;
        private double b = -2.0;
        private double c = 2.7;
        private double d = 5.1;

        [TestMethod]
        public void IntervalZeroWidth () {
            Interval aa = Interval.FromEndpoints(a, a);
            Assert.IsTrue(aa.Width == 0.0);
            Assert.IsTrue(aa.LeftEndpoint == aa.RightEndpoint);
            Assert.IsTrue(aa.Midpoint == a);
        }

        [TestMethod]
        public void IntervalWidth () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.IsTrue(ab.Width == Math.Abs(a - b));
        }

        [TestMethod]
        public void IntervalEndpoints () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.IsTrue(ab.LeftEndpoint == a);
            Assert.IsTrue(ab.RightEndpoint == b);
        }

        [TestMethod]
        public void IntervalMidpoint () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.IsTrue(ab.Midpoint == (a + b) / 2.0);
            Assert.IsTrue(ab.Contains(ab.Midpoint));
        }

        [TestMethod]
        public void IntervalOpenContains () {
            Interval ac = Interval.FromEndpoints(a, c);
            Assert.IsFalse(ac.OpenContains(a));
            Assert.IsFalse(ac.OpenContains(c));
            Assert.IsTrue(ac.OpenContains(b));
            Assert.IsFalse(ac.OpenContains(d));
        }

        [TestMethod]
        public void IntervalContainsEndpoints () {
            Interval ab = Interval.FromEndpoints(a, b);
            Assert.IsTrue(ab.Contains(a, IntervalType.Closed));
            Assert.IsFalse(ab.Contains(a, IntervalType.Open));
            Assert.IsTrue(ab.Contains(a, IntervalType.Closed, IntervalType.Open));
            Assert.IsTrue(ab.Contains(b, IntervalType.Closed));
            Assert.IsFalse(ab.Contains(b, IntervalType.Open));
            Assert.IsFalse(ab.Contains(b, IntervalType.Closed, IntervalType.Open));
        }

        [TestMethod]
        public void IntervalFromMidpointAndWidth () {
            double m = 1.0;
            double w = 4.0;
            Interval ab = Interval.FromMidpointAndWidth(m, w);
            Assert.IsTrue(ab.Width == w);
            Assert.IsTrue(ab.Midpoint == m);
            Assert.IsTrue(ab.LeftEndpoint == m - w / 2.0);
            Assert.IsTrue(ab.RightEndpoint == m + w / 2.0);
        }

        [TestMethod]
        public void IntervalFromEndpoints () {
            Interval ab = Interval.FromEndpoints(a, b);
            Interval ba = Interval.FromEndpoints(b, a);
            Assert.IsTrue(ab.LeftEndpoint == ba.LeftEndpoint);
            Assert.IsTrue(ab.RightEndpoint == ba.RightEndpoint);
            Assert.IsTrue(ab.Width == ba.Width);
            Assert.IsTrue(ab == ba);
            Assert.IsFalse(ab != ba);
            Assert.IsTrue(ab.GetHashCode() == ba.GetHashCode());
        }

        [TestMethod]
        public void IntervalClosedContains () {
            Interval ac = Interval.FromEndpoints(a, c);
            Assert.IsTrue(ac.ClosedContains(a));
            Assert.IsTrue(ac.ClosedContains(c));
            Assert.IsTrue(ac.ClosedContains(b));
            Assert.IsFalse(ac.ClosedContains(d));
        }

        [TestMethod]
        public void IntervalEquality () {
            Interval ab = Interval.FromEndpoints(a, b);
            Interval ac = Interval.FromEndpoints(a, c);
            Assert.IsTrue(ab.Equals(ab));
            Assert.IsTrue(ab.Equals((object) ab));
            Assert.IsFalse(ab.Equals(ac));
            Assert.IsFalse(ab.Equals((object) ac));
            Assert.IsFalse(ab.Equals(new Object()));
        }

        [TestMethod]
        public void IntervalToString() {
            Interval ab = Interval.FromEndpoints(a, b);
            ab.ToString().Should().Be($"[{a},{b}]");
        }

    }
}
