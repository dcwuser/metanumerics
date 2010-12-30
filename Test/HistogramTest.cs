using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;

namespace Test {
    [TestClass]
    public class HistogramTest {

        [TestMethod]
        public void TestMethod1 () {

            Histogram h = new Histogram(6);
            for (int j = 0; j < 3; j++) {
                for (int i = 0; i < 6; i++) {
                    h.Add(i);
                }
            }
            for (int i = 0; i < 6; i++) {
                Console.WriteLine(h[i]);
            }

        }

    }
}
