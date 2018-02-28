using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Security.Cryptography;
using System.Text;
using System.Net;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Data;
using Meta.Numerics.Statistics;

namespace DataTest
{
    [TestClass]
    public class DataFrameTests
    {

        [TestMethod]
        public void DataFrameColumnManipulations ()
        {
            ColumnDefinition column0 = new ColumnDefinition<int>("Integer");
            ColumnDefinition column1 = new ColumnDefinition<double>("Double");
            ColumnDefinition column2 = new ColumnDefinition<DateTime>("Timestamp");

            DataFrame frame = new DataFrame(column0, column1);
            Assert.IsTrue(frame.Columns.Count == 2);

            Assert.IsTrue(frame.Columns[0].Name == column0.Name);
            Assert.IsTrue(frame.Columns[0].StorageType == column0.StorageType);
            Assert.IsTrue(frame.Columns[column0.Name].Name == column0.Name);
            Assert.IsTrue(frame.Column<object>(column0.Name).Name == column0.Name);

            frame.AddColumn(column2);
            Assert.IsTrue(frame.Columns.Count == 3);

            Assert.IsTrue(frame.Columns[2].Name == column2.Name);
            Assert.IsTrue(frame.Columns[2].StorageType == column2.StorageType);
            Assert.IsTrue(frame.Column<object>(column2.Name).Name == column2.Name);

            frame.RemoveColumn(column0.Name);
            Assert.IsTrue(frame.Columns.Count == 2);

            Assert.IsTrue(frame.Columns[0].Name == column1.Name);
            Assert.IsTrue(frame.Columns[0].StorageType == column1.StorageType);
            Assert.IsTrue(frame.Column<object>(column1.Name).Name == column1.Name);

        }

        [TestMethod]
        public void DataFrameRowManipulations () {

            // In future, test with computed columns
            DataFrame frame = new DataFrame(new ColumnDefinition<double>("Height"), new ColumnDefinition<string>("Name"));
            //DataFrame frame = new DataFrame(new ColumnDefinition<double>("Height"), new ColumnDefinition<string>("Name"), new ComputedColumnDefinition<int>("NameLength", r => ((string) r["Name"]).Length));
            Assert.IsTrue(frame.Columns.Count == 2);
            Assert.IsTrue(frame.Rows.Count == 0);

            // Insert a row
            Dictionary<string, object> row = new Dictionary<string, object>() { { "Name", "John" }, { "Height", 1.1 } };
            frame.AddRow(row);
            Assert.IsTrue(frame.Rows.Count == 1);

            // Try to insert a row with missing values
            Dictionary<string, object> smallRow = new Dictionary<string, object>() { { "Name", "Mark" } };
            try {
                frame.AddRow(smallRow);
                Assert.Fail();
            } catch (Exception) { }

            // Try to insert a row with too many values
            Dictionary<string, object> bigRow = new Dictionary<string, object>() { { "Name", "Luke" }, { "Height", 1.2 }, { "Weight", 60.0 } };
            try {
                frame.AddRow(bigRow);
                Assert.Fail();
            } catch (Exception) { }

        }

        private const string csvFileName = "train.csv";

        [TestMethod]
        [DeploymentItem(csvFileName)]
        public void DataFrameCsvRoundtrip()
        {
            DataFrame frame;
            using (TextReader reader = File.OpenText(csvFileName))
            {
                frame = DataFrame.ReadCsv(reader);
            }

            Assert.IsTrue(frame != null);
            Assert.IsTrue(frame.Columns.Count > 0);
            Assert.IsTrue(frame.Rows.Count > 0);

            string outputPath = Path.GetTempFileName();
            try
            {
                using (FileStream stream = File.OpenWrite(outputPath))
                {
                    using (TextWriter writer = new StreamWriter(stream))
                    {
                        frame.WriteCsv(writer);
                    }
                }

                Guid inputHash = ComputeMD5Hash(csvFileName);
                Guid outputHash = ComputeMD5Hash(outputPath);
                Assert.IsTrue(inputHash == outputHash);
            }
            finally
            {
                File.Delete(outputPath);
            }


        }

        private static Guid ComputeMD5Hash (string path)
        {
            using (MD5 hasher = new MD5CryptoServiceProvider())
            {
                using (FileStream stream = File.OpenRead(path))
                {
                    byte[] hash = hasher.ComputeHash(stream);
                    return (new Guid(hash));
                }
            }
        }

        [TestMethod]
        public void DataFrameDictionariesRoundtrip()
        {
            DataFrame frame = new DataFrame(
                new ColumnDefinition<string>("name"),
                new ColumnDefinition<double>("height"),
                new ColumnDefinition<bool?>("male")
            );
            frame.AddRow("a", 5.0, false);
            frame.AddRow("b", 6.0, true);
            frame.AddRow("c", 5.5, null);

            List<Dictionary<string, object>> dictionaries = frame.ToDictionaries().ToList();
      
            Assert.IsTrue(dictionaries.Count == frame.Rows.Count);
            Assert.IsTrue(dictionaries[0].Count == frame.Columns.Count);

            DataFrame frame2 = DataFrame.FromDictionaries(dictionaries);

            Assert.IsTrue(frame2.Rows.Count == frame.Rows.Count);
            Assert.IsTrue(frame2.Columns.Count == frame.Columns.Count);
            Assert.IsTrue(frame2.Columns[0].Name == frame.Columns[0].Name);
            Assert.IsTrue(frame2.Columns[1].StorageType == frame.Columns[1].StorageType);
            Assert.IsTrue(frame2.Rows[2]["male"] == frame2.Rows[2]["male"]);

        }

        [TestMethod]
        public void DataFrameCsvRoundtrip2 () {

            // Let's exercise all our data adaptors
            DataFrame original = new DataFrame(
                new ColumnDefinition<string>("String"),
                new ColumnDefinition<double?>("Double?"),
                new ColumnDefinition<int>("Int"),
                new ColumnDefinition<DateTime?>("DateTime?"),
                new ColumnDefinition<TimeSpan>("TimeSpan"),
                new ColumnDefinition<Boolean?>("Boolean?")
            );
            original.AddRow("z", null, 1, DateTime.Today, TimeSpan.FromMinutes(5.0), true);
            original.AddRow("y", 4.3, 2, null, TimeSpan.FromHours(4.0), null);
            original.AddRow("x", 2.0, 3, DateTime.UtcNow.Date, TimeSpan.FromDays(3.0), false);

            TextWriter storage = new StringWriter();
            original.WriteCsv(storage);

            DataFrame copy = DataFrame.ReadCsv(new StringReader(storage.ToString()));
            for(int i = 0; i < original.Columns.Count; i++) {
                Assert.IsTrue(original.Columns[i].Name == copy.Columns[i].Name);
                Assert.IsTrue(original.Columns[i].StorageType == copy.Columns[i].StorageType);
            }

            for (int i = 0; i < original.Rows.Count; i++) {
                for (int j = 0; j < original.Columns.Count; j++) {
                    // This awkwardness is necessary because == resolves to a static method,
                    // so object == object does a reference check which will fail even if
                    // both sides are equal structures. Equals, on the other hand, is a
                    // virtual method, so it will do the appropriate comparison, but will
                    // fail if the instance is null.
                    if (original.Rows[i][j] == null) {
                        Assert.IsTrue(original.Rows[i][j] == null);
                    } else {
                        Assert.IsTrue(original.Rows[i][j].Equals(copy.Rows[i][j]));
                    }
                }
            }
        }

        [TestMethod]
        public void Smoketest () {

            DataFrame frame;
            string url = "https://raw.github.com/pandas-dev/pandas/master/pandas/tests/data/tips.csv";
            WebRequest request = WebRequest.Create(url);
            using (WebResponse response = request.GetResponse()) {
                using (Stream responseStream = response.GetResponseStream()) {
                    using (TextReader reader = new StreamReader(responseStream)) {
                        frame = DataFrame.ReadCsv(reader);
                    }
                }
            }
            frame.AddComputedColumn("tip_fraction", r => ((double) r["tip"]) / ((double) r["total_bill"]));

            DataView counts = frame.GroupBy("day", v => v.Rows.Count, "total").OrderBy("day");
            DataView means = frame.GroupBy("sex", v => v.Column<double>("tip_fraction").Mean(), "mean_tip_fraction");

        }

        [TestMethod]
        public void SmokeTest2 () {

            DataFrame frame;
            string path = @"C:\Users\dcw-b\Desktop\DataSets\551184489_52017_210_airline_delay_causes\551184489_52017_210_airline_delay_causes.csv";
            using (StreamReader stream = File.OpenText(path)) {
                frame = DataFrame.ReadCsv(stream);
            }

            DataView view = frame.GroupBy("carrier", (DataView q) => q.Rows.Count, "count");

        }
    }
}
