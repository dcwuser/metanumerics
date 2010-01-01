using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Statistics;

namespace Test {

    [TestClass]
    public class BugTests {

        [TestMethod]
        public void Bug2811 () {

            ChiSquaredDistribution d = new ChiSquaredDistribution(1798);
            double x = d.InverseLeftProbability(0.975);

            //ChiSquaredDistribution d = new ChiSquaredDistribution(4);
            //double x = d.InverseLeftProbability(1.0E-10);
            Console.WriteLine(x);
        }

    }

}