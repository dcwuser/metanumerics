using System;
using System.Collections.Generic;
using System.IO;
using System.Net;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Data;
using Meta.Numerics.Statistics;

namespace Test {

    [TestClass]
    public class FrameTableTest {

        [TestMethod]
        public void FrameTableManipulation () {

            FrameTable table = new FrameTable(
                new ColumnDefinition<int>("Id"),
                new ColumnDefinition<DateTime>("Birthdate"),
                new ColumnDefinition<List<int>>("Results")
            );

            Assert.IsTrue(table.GetColumnIndex("Birthdate") >= 0);
            Assert.IsTrue(table.GetColumnIndex("None") < 0);

            // Add rows
            Assert.IsTrue(table.Rows.Count == 0);
            table.AddRow(1, DateTime.Parse("1990-01-01"), new List<int>() { 1, 2, 4 });
            table.AddRow(2, DateTime.Parse("2000-02-02"), new List<int>() { 2, 3 });
            Assert.IsTrue(table.Rows.Count == 2);

            // Adding a new column with the wrong length should fail
            try {
                table.AddColumn(new ColumnDefinition<string>("Name"));
                Assert.Fail();
            } catch (Exception) { }

            // Adding a new computed column should work
            table.AddColumn(
                new ComputedColumnDefinition<TimeSpan?>("Age", r => {
                    DateTime? b = (DateTime?) r["Birthdate"];
                    if (b.HasValue) {
                        return (DateTime.Now - b.Value);
                    } else {
                        return (null);
                    }
                })
            );
            Assert.IsTrue(table.GetColumnIndex("Age") >= 0);
            table.AddRow(new Dictionary<string, object>() { { "Id", 3 }, { "Birthdate", DateTime.Now.AddYears(-3) }, { "Results", null } });

            Assert.IsTrue((int) table.Columns["Id"][0] == 1);

            table.Clear();
            Assert.IsTrue(table.Columns.Count > 0);
            Assert.IsTrue(table.Rows.Count == 0);

        }
        
        [TestMethod]
        public void InternetCsvDownload () {
            FrameTable table = DownloadFrameTable(new Uri("https://raw.githubusercontent.com/Dataweekends/zero_to_deep_learning_udemy/master/data/weight-height.csv"));
            FrameView view = table.WhereNotNull().Where("Gender", (string s) => (s == "Male"));
            view.AddComputedColumn("Bmi", (FrameRow r) => {
                double h = (double) r["Height"];
                double w = (double) r["Weight"];
                return (w / (h * h));
            });
            SampleSummary summary = new SampleSummary(view["Bmi"].As<double>());
            PolynomialRegressionResult result1 = view["Height"].As<double>().PolynomialRegression(view["Weight"].As<double>(), 1);
            PolynomialRegressionResult result2 = view["Height"].As<double>().PolynomialRegression(view["Weight"].As<double>(), 2);
            PolynomialRegressionResult result3 = view["Height"].As<double>().PolynomialRegression(view["Weight"].As<double>(), 3);
        }

        public FrameTable DownloadFrameTable (Uri url) {
            FrameTable frame;
            WebRequest request = WebRequest.Create(url);
            using (WebResponse response = request.GetResponse()) {
                using (Stream responseStream = response.GetResponseStream()) {
                    using (TextReader reader = new StreamReader(responseStream)) {
                        frame = FrameTable.FromCsv(reader);
                    }
                }
            }
            return (frame);
        }

    }
}
