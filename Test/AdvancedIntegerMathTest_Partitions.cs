using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;


using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Functions;
using FluentAssertions;

namespace Test {

    [TestClass]
    public class AdvancedIntegerMathTest_Partitions {

        [TestMethod]
        public void IntegerPartitionChecks () {
            foreach (int n in TestUtilities.GenerateIntegerValues(2, 32, 4)) {
                CheckPartitions(n);
            }
        }

        private void CheckPartitions (int n) {

            // Some sets to make sure partition and conjugate only appear once.
            // Also tests equality, hash set.
            HashSet<IntegerPartition> set = new HashSet<IntegerPartition>();
            HashSet<IntegerPartition> conjugateSet = new HashSet<IntegerPartition>();

            foreach(IntegerPartition partition in IntegerPartition.GetPartitions(n)) {

                Assert.IsTrue(partition != null);

                // Values should add to number.
                int vCount = 0;
                foreach (int v in partition.Values) {
                    Assert.IsTrue(v > 0);
                    vCount += v;
                }
                Assert.IsTrue(vCount == n);

                // Elements should add to number.
                int eCount = 0;
                foreach (Element e in partition.Elements) {
                    Assert.IsTrue(e.Value > 0);
                    Assert.IsTrue(e.Multiplicity > 0);
                    eCount += e.Value * e.Multiplicity;
                }
                Assert.IsTrue(eCount == n);

                // Partition should be generated only once
                Assert.IsTrue(!set.Contains(partition));
                set.Add(partition);

                IntegerPartition conjugate = partition.Conjugate();

                // Conjugate values should add to same number
                int cCount = 0;
                    foreach (int c in conjugate.Values) {
                    Assert.IsTrue(c > 0);
                    cCount += c;
                }
                Assert.IsTrue(cCount == n);

                // Conjugate should be generated only once
                Assert.IsTrue(!conjugateSet.Contains(conjugate));
                conjugateSet.Add(conjugate);

                // Rank should fulfill inequality and be related to conjugate rank
                Assert.IsTrue((-n < partition.Rank) && (partition.Rank < n));
                Assert.IsTrue(conjugate.Rank == -partition.Rank);

                // Double-conjugating should return us to the original
                Assert.IsTrue(conjugate.Conjugate() == partition);

                // Equality methods should work
                (partition == null).Should().BeFalse();
                (partition != null).Should().BeTrue();
                partition.Equals(partition).Should().BeTrue();
                partition.Equals((object)partition).Should().BeTrue();

            }

        }

        [TestMethod]
        public void IntegerPartitionSmallCounts () {

            // These counts are from Table 21.5 of Abromowitz & Stegun
            Assert.IsTrue(PartitionFunction(1) == 1);
            Assert.IsTrue(PartitionFunction(2) == 2);
            Assert.IsTrue(PartitionFunction(3) == 3);
            Assert.IsTrue(PartitionFunction(4) == 5);
            Assert.IsTrue(PartitionFunction(5) == 7);
            Assert.IsTrue(PartitionFunction(6) == 11);
            Assert.IsTrue(PartitionFunction(7) == 15);
            Assert.IsTrue(PartitionFunction(8) == 22);

        }

        [TestMethod]
        public void IntegerPartitionRamanujanCongruences () {

            // https://en.wikipedia.org/wiki/Ramanujan%27s_congruences

            for (int k = 0; k < 3; k++) {
                Assert.IsTrue(PartitionFunction(5 * k + 4) % 5L == 0);
                Assert.IsTrue(PartitionFunction(7 * k + 5) % 7L == 0);
                Assert.IsTrue(PartitionFunction(11 * k + 6) % 11L == 0);
            }

        }

        private long PartitionFunction (int n) {
            return (IntegerPartition.GetPartitions(n).LongCount());
        }

    }
}
