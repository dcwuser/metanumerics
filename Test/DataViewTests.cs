using System;
using System.Collections.Generic;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Data;
using Meta.Numerics.Statistics;

namespace DataTest
{
    [TestClass]
    public class DataViewTests
    {

        private DataFrame GetTestFrame ()
        {
            DataFrame frame = new DataFrame(
                new DataHeader<string>("name"),
                new DataHeader<double>("height"),
                new DataHeader<bool?>("male")
            );
            frame.AddRow("a", 7.0, false);
            frame.AddRow(null, 6.5, true);
            frame.AddRow("c", 6.0, false);
            frame.AddRow("d", 5.5, true);
            frame.AddRow("e", 5.0, null);
            frame.AddRow("f", 4.5, true);

            return (frame);
        }

        [TestMethod]
        public void DataViewSelectClause ()
        {
            DataView original = GetTestFrame();

            // We should get the columns we specify, in the order we specify
            DataView selected = original.Select(original.Columns[1].Name, original.Columns[0].Name);
            Assert.IsTrue(selected.Columns.Count == 2);
            Assert.IsTrue(selected.Rows.Count == original.Rows.Count);
            Assert.IsTrue(selected.Columns[0].Name == original.Columns[1].Name);
            Assert.IsTrue(selected.Columns[1].Name == original.Columns[0].Name);

        }

        [TestMethod]
        public void DataViewOrderByColumn ()
        {
            DataView original = GetTestFrame();

            string columnName = "height";

            DataView reordered = original.OrderBy(columnName);
            Assert.IsTrue(reordered.Columns.Count == original.Columns.Count);
            Assert.IsTrue(reordered.Rows.Count == original.Rows.Count);

            for (int i = 1; i < reordered.Rows.Count; i++)
            {
                IComparable previous = (IComparable) reordered.Rows[i-1][columnName];
                IComparable current = (IComparable) reordered.Rows[i][columnName];
                Assert.IsTrue(previous.CompareTo(current) <= 0);
            }

        }

        [TestMethod]
        public void DataViewOrderByColumnFunc () {

            DataView original = GetTestFrame();

            // This should support nulls
            DataView reordered = original.OrderBy("name");
            //DataView reordered = original.OrderBy<string>("name", String.Compare);
            Assert.IsTrue(reordered.Columns.Count == original.Columns.Count);
            Assert.IsTrue(reordered.Rows.Count == original.Rows.Count);

            List<string> reorderedCopy = original.Column<string>("name").ToList();
            reorderedCopy.Sort();

            DataColumn<string> reorderedColumn = reordered.Column<string>("name");

            for (int i = 0; i < reorderedColumn.Count; i++) {
                Assert.IsTrue(reorderedColumn[i] == reorderedCopy[i]);
            }

        }

        [TestMethod]
        public void DataViewWhereClause ()
        {
            DataView original = GetTestFrame();

            DataView filtered = original.Where<double>("height", h => (h < 5.5));
            Assert.IsTrue(filtered.Columns.Count == original.Columns.Count);
            Assert.IsTrue(filtered.Rows.Count <= original.Rows.Count);
        }

        [TestMethod]
        public void DataViewGroupByClause ()
        {
            DataView original = GetTestFrame();//.WhereNotNull("male");

            DataFrame grouped = original.GroupBy("male", v => v.Column<double>("height").Mean(), "meanHeight");

            foreach (DataRow row in grouped.Rows) {
                object value = row["male"];
                // r["male"] == value doesn't work, because this is an object comparison,
                // and two equal values boxed to objects are not equal. This is a problem.
                DataView selected = original.Where<bool?>("male", m => m.Equals(value));
                //DataView selected = original.Where(r => (r["male"] == value));
                double height = selected.Column<double>("height").Mean();
                Assert.IsTrue((double) row["meanHeight"] == height);
            }

        }

    }
}
