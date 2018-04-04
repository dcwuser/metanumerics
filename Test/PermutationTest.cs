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

                foreach (Permutation p in Permutation.GetPermutations(n)) {
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

            // These two perumations are the same. One is just in cycle notation and the other in list notation.
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

            // We want to test that GetRandomPermutation actually samples all permutations equally.

            // Don't let n get too big or we will have a ridiculously large number of bins.
            for (int n = 2; n < 8; n++) {

                // Build a mapping that assigns each permutation a unique integer index from 0 to (n! - 1). 
                Dictionary<Permutation, int> index = new Dictionary<Permutation, int>();
                int count = 0;
                foreach (Permutation p in Permutation.GetPermutations(n)) {
                    index.Add(p, count);
                    count++;
                }

                // Create a histogram of randomly generated permutation indexes.
                Histogram histogram = new Histogram(count);
                Random rng = new Random(2);
                for (int i = 0; i < 8 * count; i++) {
                    Permutation p = Permutation.GetRandomPermutation(n, rng);
                    histogram.Bins[index[p]].Increment();
                }

                // Test the uniformity of the distribution.
                TestResult result = histogram.ChiSquaredTest(new DiscreteUniformDistribution(0, count - 1));
                Assert.IsTrue(result.Probability > 0.01);

            }

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

        [TestMethod]
        public void GenerateGroups () {

            Permutation a = Permutation.Parse("(0 1 2)(3)(4)");
            Permutation b = Permutation.Parse("(0 1 2 3 4)");
            HashSet<Permutation> A5 = GenerateGroup(new Permutation[] { a, b });
            Console.WriteLine(A5.Count);

            Permutation c = Permutation.Parse("(0 1 2 3 4 5 6 7 8 9 10)");
            Permutation d = Permutation.Parse("(0)(1)(2 6 10 7)(3 9 4 5)(8)");
            HashSet<Permutation> M11 = GenerateGroup(new Permutation[] { c, d });
            Console.WriteLine(M11.Count);

        }

        private HashSet<Permutation> GenerateGroup (IList<Permutation> generators) {

            // The basic idea is to keep multiplying by generators until we stop generating new elements

            // Validate the generators. There should be at least one, all should have the same dimension, and none should be the identity.
            int d = generators[0].Dimension;

            // Start with the identity, add it to the group
            Permutation e = Permutation.Identity(d);
            HashSet<Permutation> group = new HashSet<Permutation>();
            group.Add(e);

            // We will generate elements in "generations". The first generation is just the identity.
            List<Permutation> parents = new List<Permutation>();
            parents.Add(e);

            // Keep going as long as we have some elements in the parent generation.
            while (parents.Count > 0) {

                // Start with an empty list of children.
                List<Permutation> childern = new List<Permutation>();

                // Multiply every parent by a generator.
                foreach (Permutation parent in parents) {
                    foreach (Permutation generator in generators) {
                        Permutation child = parent * generator;

                        // If it creates a new element, add that element to the group and to the children.
                        // The children list will thus contain only elements that were first created in this generation.
                        if (!group.Contains(child)) {
                            group.Add(child);
                            childern.Add(child);
                        }
                        // The first generation will thus be the generators themselves, and will end with the group set containing
                        // the identity and the generators.

                    }
                }

                // The children become the parents of the next generation.
                // If no new elements are generated, there will be no children and the algorithm will terminate.
                parents = childern;

            }

            return (group);

        }

    }
}
