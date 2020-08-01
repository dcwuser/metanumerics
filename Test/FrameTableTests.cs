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

namespace Test
{
    [TestClass]
    public class FrameTableTests
    {

        [TestMethod]
        public void FrameTableColumnManipulations ()
        {

            FrameTable frame = new FrameTable();
            frame.AddColumn<int>("Integer");
            frame.AddColumn<double>("Double");
            Assert.IsTrue(frame.Columns.Count == 2);

            Assert.IsTrue(frame.Columns[0].Name == frame[frame.Columns[0].Name].Name);
            Assert.IsTrue(frame.Columns[0].StorageType == frame[frame.Columns[0].Name].StorageType);

            frame.AddColumn<DateTime>("Timestamp");
            Assert.IsTrue(frame.Columns.Count == 3);

            frame.RemoveColumn(frame.Columns[0].Name);
            Assert.IsTrue(frame.Columns.Count == 2);
        }

        [TestMethod]
        public void FrameTableRowManipulations () {

            // In future, test with computed columns
            FrameTable frame = new FrameTable();
            frame.AddColumn<double>("Height");
            frame.AddColumn<string>("Name");

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
        public void FrameTableCsvRoundtrip()
        {
            FrameTable frame;
            using (TextReader reader = File.OpenText(csvFileName))
            {
                frame = FrameTable.FromCsv(reader);
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
                        frame.ToCsv(writer);
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
        public void FrameTableDictionariesRoundtrip()
        {
            FrameTable frame = new FrameTable();
            frame.AddColumn<string>("name");
            frame.AddColumn<double>("height");
            frame.AddColumn<bool?>("male");

            frame.AddRow("a", 5.0, false);
            frame.AddRow("b", 6.0, true);
            frame.AddRow("c", 5.5, null);

            List<Dictionary<string, object>> dictionaries = frame.ToDictionaries().ToList();
      
            Assert.IsTrue(dictionaries.Count == frame.Rows.Count);
            Assert.IsTrue(dictionaries[0].Count == frame.Columns.Count);

            FrameTable frame2 = FrameTable.FromDictionaries(dictionaries);

            Assert.IsTrue(frame2.Rows.Count == frame.Rows.Count);
            Assert.IsTrue(frame2.Columns.Count == frame.Columns.Count);
            Assert.IsTrue(frame2.Columns[0].Name == frame.Columns[0].Name);
            Assert.IsTrue(frame2.Columns[1].StorageType == frame.Columns[1].StorageType);
            Assert.IsTrue(frame2.Rows[2]["male"] == frame2.Rows[2]["male"]);

        }

        [TestMethod]
        public void FrameTableCsvRoundtrip2 () {

            // Let's exercise all our data adaptors
            FrameTable original = new FrameTable();
            original.AddColumn<string>("String");
            original.AddColumn<double?>("Double?");
            original.AddColumn<int>("Int");
            original.AddColumn<DateTime?>("DateTime?");
            original.AddColumn<TimeSpan>("TimeSpan");
            original.AddColumn<Boolean?>("Boolean?");

            original.AddRow("z", null, 1, DateTime.Today, TimeSpan.FromMinutes(5.0), true);
            original.AddRow("y", 4.3, 2, null, TimeSpan.FromHours(4.0), null);
            original.AddRow("x", 2.0, 3, DateTime.UtcNow.Date, TimeSpan.FromDays(3.0), false);

            TextWriter storage = new StringWriter();
            original.ToCsv(storage);

            FrameTable copy = FrameTable.FromCsv(new StringReader(storage.ToString()));
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

            FrameTable frame;
            string url = "https://raw.githubusercontent.com/pandas-dev/pandas/master/doc/data/tips.csv";
            WebRequest request = WebRequest.Create(url);
            using (WebResponse response = request.GetResponse()) {
                using (Stream responseStream = response.GetResponseStream()) {
                    using (TextReader reader = new StreamReader(responseStream)) {
                        frame = FrameTable.FromCsv(reader);
                    }
                }
            }
            frame.AddComputedColumn("tip_fraction", r => ((double) r["tip"]) / ((double) r["total_bill"]));

            FrameView counts = frame.GroupBy("day", v => v.Rows.Count, "total").OrderBy("day");
            FrameView means = frame.GroupBy("sex", v => v["tip_fraction"].As<double>().Mean(), "mean_tip_fraction");

        }

        [TestMethod]
        public void SmokeTest2 () {

            FrameTable frame;
            string path = @"C:\Users\dcw-b\Desktop\DataSets\551184489_52017_210_airline_delay_causes\551184489_52017_210_airline_delay_causes.csv";
            using (StreamReader stream = File.OpenText(path)) {
                frame = FrameTable.FromCsv(stream);
            }

            FrameView view = frame.GroupBy("carrier", (FrameView q) => q.Rows.Count, "count");

        }


        [TestMethod]
        public void CsvWhitespaceParsing () {
            StringBuilder text = new StringBuilder();
            text.AppendLine(" c0 ,\tc1  ");
            text.AppendLine(",");
            text.AppendLine(" , ");
            text.AppendLine("  ,\t");
            text.AppendLine("\t\" \" ,  \"\t\" ");
            text.AppendLine(" \" a\t\"\t,\t\"\t a \"");
            text.AppendLine("\t \" \"\"\",  \" , \" ");

            FrameTable table = FrameTable.FromCsv(new StringReader(text.ToString()));

            Assert.IsTrue(table.Columns[0].Name == "c0");
            Assert.IsTrue(table.Columns[1].Name == "c1");

            Assert.IsTrue((string) table.Rows[0][0] == null);
            Assert.IsTrue((string) table.Rows[0][1] == null);
            Assert.IsTrue((string) table.Rows[1][0] == null);
            Assert.IsTrue((string) table.Rows[1][1] == null);
            Assert.IsTrue((string) table.Rows[2][0] == null);
            Assert.IsTrue((string) table.Rows[2][1] == null);
            Assert.IsTrue((string) table.Rows[3][0] == " ");
            Assert.IsTrue((string) table.Rows[3][1] == "\t");
            Assert.IsTrue((string) table.Rows[4][0] == " a\t");
            Assert.IsTrue((string) table.Rows[4][1] == "\t a ");
            Assert.IsTrue((string) table.Rows[5][0] == " \"");
            Assert.IsTrue((string) table.Rows[5][1] == " , ");
        }
    }
}
