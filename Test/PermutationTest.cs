using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Test {

    [TestClass]
    public class PermutationTest {

        [TestMethod]
        public void PermutationCount () {

            foreach (int n in TestUtilities.GenerateIntegerValues(2, 8, 4)) {

                int totalCount = 0;
                int evenCount = 0;
                int derangementCount = 0;

                foreach (Permutation p in Permutation.Permutations(n)) {
                    totalCount++;
                    if (p.IsEven) evenCount++;
                    if (p.IsDerangement) derangementCount++;
                }

                Assert.IsTrue(totalCount == (int) AdvancedIntegerMath.Factorial(n));
                Assert.IsTrue(evenCount == totalCount / 2);
                //Assert.IsTrue(derangementCount == AdvancedIntegerMath.Subfactorial(n));

            }

        }

        [TestMethod]
        public void PermutationEquality () {

            Permutation a1 = Permutation.Parse("(0 1)(2 3)");
            Permutation a2 = Permutation.Parse("[1 0 3 2]");

            Permutation b = Permutation.Parse("(0 1 2)(3)");

            Assert.IsTrue(a1 == a2);
            Assert.IsTrue(a2 == a1);
            Assert.IsFalse(a1 == b);
            Assert.IsFalse(b == a1);
            Assert.IsFalse(b == null);
            Assert.IsFalse(null == b);
            Assert.IsTrue((Permutation) null == (Permutation) null);

            Assert.IsFalse(a1 != a2);
            Assert.IsFalse(a2 != a1);
            Assert.IsTrue(a1 != b);
            Assert.IsTrue(b != a1);
            Assert.IsTrue(b != null);
            Assert.IsTrue(null != b);
            Assert.IsFalse((Permutation) null != (Permutation) null);

            Assert.IsTrue(b.Equals(b));
            Assert.IsFalse(b.Equals(a1));
            Assert.IsFalse(b.Equals(new object()));
            Assert.IsFalse(b.Equals(null));

            Assert.IsTrue(a1.GetHashCode() == a2.GetHashCode());

        }

        [TestMethod]
        public void PermutationProduct () {

            Random rng = new Random(4);
            string s = "abcdefghijklmnopqrstuvwxyz";

            for (int i = 0; i < 16; i++) {

                int n = rng.Next(2, s.Length);

                char[] c0 = s.Substring(0, n).ToCharArray();
                char[] c1 = (char[]) c0.Clone();
                char[] c2 = (char[]) c0.Clone();

                Permutation p = Permutation.GetRandomPermutation(n, rng);
                Permutation q = Permutation.GetRandomPermutation(n, rng);

                q.Apply(c1);
                p.Apply(c1);
                string s1 = new String(c1);

                Permutation pq = p * q;
                pq.Apply(c2);
                string s2 = new String(c2);

                Console.WriteLine("{0} {1}", s1, s2);

                Assert.IsTrue(s1 == s2);

            }

        }


        [TestMethod]
        public void PermutationInverse () {

            Random rng = new Random(1);
            for (int i = 0; i < 16; i++) {

                int n = rng.Next(2, 32);
                Permutation p = Permutation.GetRandomPermutation(n, rng);

                Permutation pi = p.Inverse();

                Assert.IsTrue((p * pi).IsIdentity);
                Assert.IsTrue((pi * p).IsIdentity);

                if (p.IsInvolution) {
                    Assert.IsTrue(p == pi);
                    if (p.IsIdentity) {
                        Assert.IsTrue(p.Order == 1);
                    } else {
                        Assert.IsTrue(p.Order == 2);
                    }
                } else {
                    Assert.IsTrue(p != pi);
                    Assert.IsTrue(p.Order > 2);
                }

            }


        }


        [TestMethod]
        public void PermutationOrder () {

            Random rng = new Random(37498327);

            for (int i = 0; i < 16; i++) {
                int n = rng.Next(1, 16);
                Permutation p = Permutation.GetRandomPermutation(n, rng);

                int order = (int) p.Order;
                Assert.IsTrue(order > 0);

                Permutation q = p;
                for (int j = 1; j < order; j++) {
                    q *= p;
                }
                Assert.IsTrue(q.IsIdentity);

                Console.WriteLine("{0} {1} {2}", p.ToString("M"), p.ToString("C"), order);
            }

        }

        [TestMethod]
        public void PermutationDistribution () {

            int n = 4;
            Dictionary<Permutation, int> index = new Dictionary<Permutation, int>();
            int count = 0;
            foreach (Permutation p in Permutation.Permutations(n)) {
                index.Add(p, count);
                count++;
            }

            Histogram bins = new Histogram(count);

            Random rng = new Random(2);
            for (int i = 0; i < 8 * count; i++) {
                Permutation p = Permutation.GetRandomPermutation(n, rng);
                bins[index[p]].Increment();
            }

            for (int i = 0; i < count; i++) {
                Console.WriteLine("{0} {1}", i, bins[i].Counts);
            }

            TestResult result = bins.ChiSquaredTest(new DiscreteUniformDistribution(0, count - 1));
            Console.WriteLine(result.RightProbability);

        }

        [TestMethod]
        public void PermutationFormat () {

            Random rng = new Random(5);

            for (int i = 0; i < 16; i++) {

                int n = rng.Next(1, 16);
                Permutation p = Permutation.GetRandomPermutation(n, rng);

                string m = p.ToString("M");
                Permutation pm = Permutation.Parse(m);
                Assert.IsTrue(pm == p);

                string c = p.ToString("C");
                Permutation pc = Permutation.Parse(c);
                Assert.IsTrue(pc == p);


            }

        }

    }
}
