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
            DataList column1 = new DataList<int>("Integer");
            DataList column2 = new DataList<double>("Double");
            DataList column3 = new DataList<DateTime>("Timestamp");

            DataFrame frame = new DataFrame(column1, column2);
            Assert.IsTrue(frame.Columns.Count == 2);

            Assert.IsTrue(frame.Columns[0].Name == column1.Name);
            Assert.IsTrue(frame.Columns[0].StorageType == column1.StorageType);
            Assert.IsTrue(frame.Column<object>(column1.Name).Name == column1.Name);

            frame.AddColumn(column3);
            Assert.IsTrue(frame.Columns.Count == 3);

            Assert.IsTrue(frame.Columns[2].Name == column3.Name);
            Assert.IsTrue(frame.Columns[2].StorageType == column3.StorageType);
            Assert.IsTrue(frame.Column<object>(column3.Name).Name == column3.Name);

            frame.RemoveColumn(column1.Name);
            Assert.IsTrue(frame.Columns.Count == 2);

            Assert.IsTrue(frame.Columns[0].Name == column2.Name);
            Assert.IsTrue(frame.Columns[0].StorageType == column2.StorageType);
            Assert.IsTrue(frame.Column<object>(column2.Name).Name == column2.Name);

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
                new DataHeader<string>("name"),
                new DataHeader<double>("height"),
                new DataHeader<bool?>("male")
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
