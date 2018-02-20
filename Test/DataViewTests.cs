using System;
using System.Collections.Generic;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Data;
using Meta.Numerics.Statistics;

namespace Test
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
        public void DataViewDiscardClause () {
            DataView original = GetTestFrame();

            DataView filtered = original.Discard(original.Columns[0].Name);
            Assert.IsTrue(filtered.Columns.Count == original.Columns.Count - 1);
            Assert.IsTrue(filtered.Rows.Count == original.Rows.Count);

            for (int i = 0; i < filtered.Columns.Count; i++) {
                Assert.IsTrue(filtered.Columns[i].Name == original.Columns[i + 1].Name);
            }
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
        public void DataViewComputedColumn () {

            DataView original = GetTestFrame();

            original.AddComputedColumn("sex", r => {
                bool? isMale = (bool?) r["male"];
                if (isMale.HasValue) {
                    return (isMale.Value ? "male" : "female");
                } else {
                    return (null);
                }
            });

            Assert.IsTrue(original.Columns["sex"].Name == "sex");
            Assert.IsTrue(original.Columns["sex"].StorageType == typeof(string));
            Assert.IsTrue(original.Columns["sex"].Count == original.Rows.Count);
            for (int i = 0; i < original.Rows.Count; i++) {
                bool? isMale = (bool?) original.Rows[i]["male"];
                string sex = (string) original.Rows[i]["sex"];
                Assert.IsTrue(
                    ((isMale == null) && (sex == null)) ||
                    ((isMale == true) && (sex == "male")) ||
                    ((isMale == false) && (sex == "female"))
                );
            }

        }

        public void DataViewOrderBy () {

            DataView original = GetTestFrame();

            DataView reordered = original.OrderBy((ra, rb) => {
                int heightComparison = ((double) ra["height"]).CompareTo((double) rb["height"]);
                return (heightComparison);
            });

            Assert.IsTrue(reordered.Rows.Count == original.Rows.Count);
            Assert.IsTrue(reordered.Columns.Count == original.Columns.Count);


        }

        [TestMethod]
        public void DataViewWhereClause ()
        {
            DataView original = GetTestFrame();
            Func<DataRow, bool> filter = row => (double) row["height"] < 5.5;
            DataView filtered = original.Where(filter);
            Assert.IsTrue(filtered.Columns.Count == original.Columns.Count);
            Assert.IsTrue(filtered.Rows.Count <= original.Rows.Count);
            foreach (DataRow row in filtered.Rows) {
                Assert.IsTrue(filter(row));
            }
        }

        [TestMethod]
        public void DataViewWhereNotNull () {

            DataView original = GetTestFrame();

            DataView filteredAll = original.WhereNotNull();
            Assert.IsTrue(filteredAll.Columns.Count == original.Columns.Count);
            Assert.IsTrue(filteredAll.Rows.Count <= original.Rows.Count);
            foreach (DataRow row in filteredAll.Rows) {
                foreach (KeyValuePair<string, object> entry in (IReadOnlyDictionary<string,object>) row) {
                    Assert.IsTrue(entry.Value != null);
                }
            }

            foreach(DataColumn<object> column in (IReadOnlyList<object>) original.Columns) {
                DataView filtered = original.WhereNotNull(column.Name);
                Assert.IsTrue(filtered.Columns.Count == original.Columns.Count);
                Assert.IsTrue(filtered.Rows.Count <= original.Rows.Count);
                foreach(DataRow row in filtered.Rows) {
                    Assert.IsTrue(row[column.Name] != null);
                }
            }

        }

        [TestMethod]
        public void DataViewGroupBy () {

            DataView original = GetTestFrame();

            HashSet<bool?> values = new HashSet<bool?>(original.Column<bool?>("male").Distinct());

            DataFrame grouped = original.GroupBy(
                "male",
                new AggregateColumn<double>("heightMean", v => v.Column<double>("height").Mean()),
                new AggregateColumn<double>("heightStandardDeviation", v => v.Rows.Count < 2 ? Double.NaN : v.Column<double>("height").StandardDeviation())
            );
            Assert.IsTrue(grouped.Rows.Count == values.Count);
            Assert.IsTrue(grouped.Columns.Count == 3);

            for (int i = 0; i < grouped.Rows.Count; i++) {

                bool? value = grouped.Column<bool?>("male")[i];
                Assert.IsTrue(values.Contains(value));

                DataView selected = original.Where(r => (bool?) r["male"] == value);
                Assert.IsTrue(selected.Rows.Count > 0);

                double mean = selected.Column<double>("height").Mean();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(grouped.Column<double>("heightMean")[i], mean));

                double standardDeviation = selected.Rows.Count < 2 ? Double.NaN : selected.Column<double>("height").StandardDeviation();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(grouped.Rows.Count < 2 ? Double.NaN : grouped.Column<double>("heightStandardDeviation")[i], standardDeviation));
            }

        }

        [TestMethod]
        public void DataViewGroupByClause ()
        {
            DataView original = GetTestFrame();//.WhereNotNull("male");

            HashSet<bool?> values = new HashSet<bool?>(original.Column<bool?>("male").Distinct());

            DataFrame grouped = original.GroupBy("male", v => v.Column<double>("height").Mean(), "meanHeight");
            Assert.IsTrue(grouped.Rows.Count == values.Count);

            foreach (DataRow row in grouped.Rows) {
                bool? value = (bool?) row["male"];
                Assert.IsTrue(values.Contains(value));
                // r["male"] == value doesn't work, because this is an object comparison,
                // and two equal values boxed to objects are not equal. This is a problem.
                //DataView selected = original.Where<bool?>("male", m => m.Equals(value));
                DataView selected = original.Where(r => (bool?) r["male"] == value);
                double height = selected.Column<double>("height").Mean();
                Assert.IsTrue((double) row["meanHeight"] == height);
            }

        }

    }
}
