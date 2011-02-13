using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;

namespace Test {

    [TestClass]
    public class ContingencyTableTest {

        [TestMethod]
        public void ContingencyTableOperations () {

            ContingencyTable t = new ContingencyTable(4, 3);
            Assert.IsTrue(t.RowCount == 4);
            Assert.IsTrue(t.ColumnCount == 3);

            Assert.IsTrue(t.RowTotal(2) == 0);
            Assert.IsTrue(t.ColumnTotal(1) == 0);
            Assert.IsTrue(t.Total == 0);

            t[1, 1] = 2;
            Assert.IsTrue(t.RowTotal(2) == 0);
            Assert.IsTrue(t.ColumnTotal(1) == 2);
            Assert.IsTrue(t.Total == 2);

        }

        [TestMethod]
        public void ContingencyTableNamedOperations () {

            ContingencyTable t = new ContingencyTable(2, 3);
            t.RowNames[0] = "Male";
            t.RowNames[1] = "Female";
            t.ColumnNames[0] = "Party 1";
            t.ColumnNames[1] = "Party 2";
            t.ColumnNames[2] = "Party 3";
            t["Male", "Party 1"] = 10;
            t["Male", "Party 2"] = 20;
            t["Male", "Party 3"] = 30;
            t["Female", "Party 1"] = 30;
            t["Female", "Party 2"] = 20;
            t["Female", "Party 3"] = 10;

        }

    }
}
