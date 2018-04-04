using System;
using System.Collections.Generic;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Data;
using Meta.Numerics.Statistics;

namespace Test
{
    [TestClass]
    public class FrameViewTests
    {

        private FrameTable GetTestFrame ()
        {
            FrameTable frame = new FrameTable();
            frame.AddColumn<string>("name");
            frame.AddColumn<double>("height");
            frame.AddColumn<double>("weight");
            frame.AddColumn<bool?>("male");
            
            frame.AddRow("a", 7.0, 10.0, false);
            frame.AddRow(null, 6.5, 11.0, true);
            frame.AddRow("c", 6.0, 12.0, false);
            frame.AddRow("d", 5.5, 11.0, true);
            frame.AddRow("e", 5.0, 12.0, null);
            frame.AddRow("f", 4.5, 13.0, true);
            frame.AddRow(null, 4.0, 12.0, false);

            frame.AddComputedColumn("bmi", r => ((double) r["weight"]) / MoreMath.Sqr((double) r["height"]));

            return (frame);
        }

        [TestMethod]
        public void FrameViewRowInterfaces () {

            FrameView original = GetTestFrame();
            FrameRow row = original.Rows[0];

            IReadOnlyDictionary<string, object> rowAsDictionary = (IReadOnlyDictionary<string, object>) row;
            Assert.IsTrue(rowAsDictionary.Count == original.Columns.Count);
            Assert.IsTrue(rowAsDictionary.ContainsKey(original.Columns[0].Name));
            Assert.IsFalse(rowAsDictionary.ContainsKey("NoName"));

            object value = null;
            Assert.IsTrue(rowAsDictionary.TryGetValue(original.Columns[0].Name, out value));
            Assert.IsTrue(original.Columns[0][0].Equals(value));
            Assert.IsFalse(rowAsDictionary.TryGetValue("NoName", out value));

            Assert.IsTrue(rowAsDictionary.Keys.Count() == original.Columns.Count);
            Assert.IsTrue(rowAsDictionary.Values.Count() == original.Columns.Count);

            IReadOnlyList<object> rowAsList = (IReadOnlyList<object>) row;
            Assert.IsTrue(rowAsList.Count == original.Columns.Count);
            Assert.IsTrue(rowAsList[0].Equals(original.Columns[0][0]));
        }

        [TestMethod]
        public void FrameViewColumnInterfaces () {

            FrameView original = GetTestFrame();

            FrameColumnCollection columns = original.Columns;

            IReadOnlyList<FrameColumn> columnsAsList = (IReadOnlyList<FrameColumn>) columns;
            Assert.IsTrue(columnsAsList.Count == original.Columns.Count);

            /*
            IReadOnlyDictionary<string, FrameColumn> columnsAsDictionary = (IReadOnlyDictionary<string, FrameColumn>) columns;
            Assert.IsTrue(columnsAsDictionary.ContainsKey(original.Columns[1].Name));
            Assert.IsFalse(columnsAsDictionary.ContainsKey("NoName"));

            foreach(KeyValuePair<string, FrameColumn> entry in columnsAsDictionary) {
                Assert.IsTrue(entry.Key == entry.Value.Name);
                FrameColumn column;
                bool found = columnsAsDictionary.TryGetValue(entry.Key, out column);
                Assert.IsTrue(found);
                Assert.IsTrue(column != null);
            }

            foreach (string name in columnsAsDictionary.Keys) {
                Assert.IsTrue(columns.GetIndexOf(name) >= 0);
            }

            foreach (FrameColumn column in columnsAsDictionary.Values) {
                Assert.IsTrue(columns.GetIndexOf(column.Name) >= 0);
            }
            */

        }

        [TestMethod]
        public void FrameViewSelectClause ()
        {
            FrameView original = GetTestFrame();

            // We should get the columns we specify, in the order we specify
            FrameView selected = original.Select(original.Columns[1].Name, original.Columns[0].Name);
            Assert.IsTrue(selected.Columns.Count == 2);
            Assert.IsTrue(selected.Rows.Count == original.Rows.Count);
            Assert.IsTrue(selected.Columns[0].Name == original.Columns[1].Name);
            Assert.IsTrue(selected.Columns[1].Name == original.Columns[0].Name);

        }

        [TestMethod]
        public void FrameViewDiscardClause () {
            FrameView original = GetTestFrame();

            FrameView filtered = original.Discard(original.Columns[0].Name);
            Assert.IsTrue(filtered.Columns.Count == original.Columns.Count - 1);
            Assert.IsTrue(filtered.Rows.Count == original.Rows.Count);

            for (int i = 0; i < filtered.Columns.Count; i++) {
                Assert.IsTrue(filtered.Columns[i].Name == original.Columns[i + 1].Name);
            }
        }

        [TestMethod]
        public void FrameViewOrderByColumn ()
        {
            FrameView original = GetTestFrame();

            string columnName = "bmi";

            FrameView reordered = original.OrderBy(columnName, SortOrder.Descending);
            Assert.IsTrue(reordered.Columns.Count == original.Columns.Count);
            Assert.IsTrue(reordered.Rows.Count == original.Rows.Count);

            int columnIndex = reordered.GetColumnIndex(columnName);
            for (int i = 1; i < reordered.Rows.Count; i++)
            {
                IComparable previous = (IComparable) reordered[i-1, columnIndex];
                IComparable current = (IComparable) reordered[i, columnIndex];
                Assert.IsTrue(previous.CompareTo(current) >= 0);
            }

        }

        [TestMethod]
        public void FrameViewOrderByColumnWithNulls () {

            FrameView original = GetTestFrame();

            // This should support nulls
            FrameView reordered = original.OrderBy("name");
            Assert.IsTrue(reordered.Columns.Count == original.Columns.Count);
            Assert.IsTrue(reordered.Rows.Count == original.Rows.Count);

            List<string> reorderedCopy = original["name"].As<string>().ToList();
            reorderedCopy.Sort();

            IReadOnlyList<string> reorderedColumn = reordered["name"].As<string>();

            for (int i = 0; i < reorderedColumn.Count; i++) {
                Assert.IsTrue(reorderedColumn[i] == reorderedCopy[i]);
            }

        }

        [TestMethod]
        public void FrameViewOrderByRowComparer () {

            FrameView original = GetTestFrame();

            // Define a comparison which sorts first by sex and then by height
            Comparison<FrameRow> comparer = (a, b) => {
                bool? aMale = (bool?) a["male"];
                bool? bMale = (bool?) b["male"];
                int aMaleValue = aMale.HasValue ? aMale == true ? 0 : 1 : 2;
                int bMaleValue = bMale.HasValue ? bMale == true ? 0 : 1 : 2;
                int maleComparison = aMaleValue.CompareTo(bMaleValue);
                if (maleComparison != 0) {
                    return (maleComparison);
                } else {
                    double aHeight = (double) a["height"];
                    double bHeight = (double) b["height"];
                    return (aHeight.CompareTo(bHeight));
                }
            };

            FrameView reordered = original.OrderBy(comparer);
            Assert.IsTrue(reordered.Columns.Count == original.Columns.Count);
            Assert.IsTrue(reordered.Rows.Count == original.Rows.Count);

            for (int i = 1; i < reordered.Rows.Count; i++) {
                Assert.IsTrue(comparer(reordered.Rows[i - 1], reordered.Rows[i]) <= 0);
            }

        }

        [TestMethod]
        public void FrameViewOrderByAgreement () {

            FrameView view = GetTestFrame();

            FrameView sorted1 = view.OrderBy("bmi");
            FrameView sorted2 = view.OrderBy("bmi", SortOrder.Descending);
            FrameView sorted3 = view.OrderBy<double>("bmi", (v1, v2) => v1.CompareTo(v2));
            FrameView sorted4 = view.OrderBy((r1, r2) => ((double) r1["bmi"]).CompareTo((double) r2["bmi"]));

            List<double> sorted = view["bmi"].As<double>().ToList();
            sorted.Sort();

            for (int i = 0; i < sorted.Count; i++) {
                Assert.IsTrue(sorted1["bmi"].As<double>()[i] == sorted[i]);
                Assert.IsTrue(sorted2["bmi"].As<double>()[sorted.Count - 1 - i] == sorted[i]);
                Assert.IsTrue(sorted3["bmi"].As<double>()[i] == sorted[i]);
                Assert.IsTrue(sorted4["bmi"].As<double>()[i] == sorted[i]);
            }

        }

        [TestMethod]
        public void FrameViewComputedColumn () {

            FrameView original = GetTestFrame();

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

        [TestMethod]
        public void FrameViewWhereClause ()
        {
            FrameView original = GetTestFrame();
            Func<FrameRow, bool> filter = row => (double) row["height"] < 5.5;
            FrameView filtered = original.Where(filter);
            Assert.IsTrue(filtered.Columns.Count == original.Columns.Count);
            Assert.IsTrue(filtered.Rows.Count <= original.Rows.Count);
            foreach (FrameRow row in filtered.Rows) {
                Assert.IsTrue(filter(row));
            }
        }

        [TestMethod]
        public void FrameViewWhereNotNull () {

            FrameView original = GetTestFrame();

            FrameView filteredAll = original.WhereNotNull();
            Assert.IsTrue(filteredAll.Columns.Count == original.Columns.Count);
            Assert.IsTrue(filteredAll.Rows.Count <= original.Rows.Count);
            foreach (FrameRow row in filteredAll.Rows) {
                foreach (KeyValuePair<string, object> entry in (IReadOnlyDictionary<string,object>) row) {
                    Assert.IsTrue(entry.Value != null);
                }
            }

            foreach(FrameColumn column in original.Columns) {
                FrameView filtered = original.WhereNotNull(column.Name);
                Assert.IsTrue(filtered.Columns.Count == original.Columns.Count);
                Assert.IsTrue(filtered.Rows.Count <= original.Rows.Count);
                foreach(FrameRow row in filtered.Rows) {
                    Assert.IsTrue(row[column.Name] != null);
                }
            }

        }

        [TestMethod]
        public void FrameViewGroupBy () {

            FrameView original = GetTestFrame();

            HashSet<bool?> values = new HashSet<bool?>(original["male"].As<bool?>().Distinct());

            FrameTable grouped = original.GroupBy("male", v => {
                SummaryStatistics summary = new SummaryStatistics(v["height"].As<double>());
                return (new Dictionary<string, object>() {
                    {"count", summary.Count },
                    {"heightMean", summary.Mean },
                    {"heightStandardDeviation", summary.StandardDeviation }
                });
            });
            Assert.IsTrue(grouped.Rows.Count == values.Count);
            Assert.IsTrue(grouped.Columns.Count == 4);

            for (int i = 0; i < grouped.Rows.Count; i++) {

                bool? value = grouped["male"].As<bool?>()[i];
                Assert.IsTrue(values.Contains(value));

                FrameView selected = original.Where(r => (bool?) r["male"] == value);
                Assert.IsTrue(selected.Rows.Count > 0);

                double mean = selected["height"].As<double>().Mean();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(grouped["heightMean"].As<double>()[i], mean));

                double standardDeviation = selected["height"].As<double>().StandardDeviation();
                Assert.IsTrue(TestUtilities.IsNearlyEqual(grouped["heightStandardDeviation"].As<double>()[i], standardDeviation));
            }

        }

        [TestMethod]
        public void FrameViewGroupByClause ()
        {
            FrameView original = GetTestFrame();

            HashSet<bool?> values = new HashSet<bool?>(original["male"].As<bool?>().Distinct());

            FrameTable grouped = original.GroupBy("male", "meanHeight", v => v["height"].As<double>().Mean());
            Assert.IsTrue(grouped.Rows.Count == values.Count);

            foreach (FrameRow row in grouped.Rows) {
                bool? value = (bool?) row["male"];
                Assert.IsTrue(values.Contains(value));
                // r["male"] == value doesn't work, because this is an object comparison,
                // and two equal values boxed to objects are not equal. This is a problem.
                //DataView selected = original.Where<bool?>("male", m => m.Equals(value));
                FrameView selected = original.Where(r => (bool?) r["male"] == value);
                double height = selected["height"].As<double>().Mean();
                Assert.IsTrue((double) row["meanHeight"] == height);
            }

        }

    }
}
