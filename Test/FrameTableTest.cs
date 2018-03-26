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
        public void FrameViewColumnCoercion () {

            // Create nullable double and integer columns.
            FrameTable table = new FrameTable();
            table.AddColumn<double?>("one");
            table.AddColumn<int>("two");

            table.AddRow(1.1, 2);
            table.AddRow(null, 3);

            // Coerce the nullable double into a non-nullable double
            // Should work when value is non-null, and fail with value is null
            IReadOnlyList<double> one = table["one"].As<double>();
            Assert.IsTrue(one[0] == 1.1);
            try {
                double v = one[1];
                Assert.Fail();
            } catch (Exception) { }

            // Coerce the integer to a double.
            IReadOnlyList<double> two = table.Columns[1].As<double>();
            Assert.IsTrue(two[0] == 2.0);

        }

        [TestMethod]
        public void FrameTableManipulation () {

            FrameTable table = new FrameTable();
            table.AddColumn<int>("Id");
            table.AddColumn<DateTime?>("Birthdate");
            table.AddColumns<string>("FirstName", "LastName");
            Assert.IsTrue(table.Columns.Count == 4);

            // Index lookup should work
            Assert.IsTrue(table.GetColumnIndex("Birthdate") >= 0);
            Assert.IsTrue(table.GetColumnIndex("None") < 0);

            // Add rows
            Assert.IsTrue(table.Rows.Count == 0);
            table.AddRow(1, DateTime.Parse("1990-01-01"), "a", "p");
            table.AddRow(2, DateTime.Parse("2000-02-02"), null, null);
            table.AddRow(new Dictionary<string, object>() {
                {"Id", 3}, {"Birthdate", null }, { "FirstName", "c" }, { "LastName", "r" }
            });
            Assert.IsTrue(table.Rows.Count == 3);

            // Adding rows with the wrong types and/or entries should fail
            // Careful, some of these will leave the table in a bad state
            //try {
            //    table.AddRow(4, DateTime.Parse("2010-04-04"), 1.0, "s");
            //    Assert.Fail();
            //} catch (Exception) { }
            try {
                table.AddRow(4, DateTime.Parse("2010-04-04"));
                Assert.Fail();
            } catch (Exception) { }
            //try {
            //    table.AddRow(new Dictionary<string, object>() {
            //        {"Id", 4}, { "FirstName", "d" }, { "LastName", "r" }
            //    });
            //    Assert.Fail();
            //} catch (Exception) { }
            //try {
            //    table.AddRow(new Dictionary<string, object>() {
            //        {"Id", 4}, { "Birthdate", null }, { "FirstName", "d" }, { "LastName", "r" }, { "MiddleName", "u" }
            //    });
            //    Assert.Fail();
            //} catch (Exception) { }

            // Adding a new column with the wrong length should fail
            try {
                table.AddColumn<double>("Score");
                Assert.Fail();
            } catch (Exception) { }
            Assert.IsTrue(table.GetColumnIndex("Score") < 0);

            // Adding a new column with the right length should work
            List<double> scores = new List<double>() { 1.1, 1.2, 1.3 };
            table.AddColumn("Score", scores);
            Assert.IsTrue(table.GetColumnIndex("Score") >= 0);

            // Adding a new computed column should work
            table.AddComputedColumn<TimeSpan?>("Age", r=> {
                DateTime? b = (DateTime?) r["Birthdate"];
                if (b.HasValue) {
                    return (DateTime.Now - b.Value);
                } else {
                    return (null);
                }
            });
            Assert.IsTrue(table.GetColumnIndex("Age") >= 0);

            // Changing a value should change the result of the computed column that depends on it
            int birthdateIndex = table.GetColumnIndex("Birthdate");
            int ageIndex = table.GetColumnIndex("Age");
            TimeSpan age1 = (TimeSpan) table[0, ageIndex];
            table[0, birthdateIndex] = DateTime.Parse("2010-01-01");
            TimeSpan age2 = (TimeSpan) table[0, ageIndex];
            Assert.IsTrue(age2 != age1);

            // Clearing a table should work
            table.Clear();
            Assert.IsTrue(table.Columns.Count > 0);
            Assert.IsTrue(table.Rows.Count == 0);

        }
        
        [TestMethod]
        public void InternetSampleDownload () {
            FrameTable table = DownloadFrameTable(new Uri("https://raw.githubusercontent.com/Dataweekends/zero_to_deep_learning_udemy/master/data/weight-height.csv"));
            FrameView view = table.WhereNotNull();
            view.AddComputedColumn("Bmi", (FrameRow r) => {
                double h = (double) r["Height"];
                double w = (double) r["Weight"];
                return (w / (h * h));
            });

            FrameView males = view.Where("Gender", (string s) => (s == "Male"));
            FrameView females = view.Where("Gender", (string s) => (s == "Female"));

            SampleSummary maleSummary = new SampleSummary(males["Height"].As<double>());
            SampleSummary femaleSummary = new SampleSummary(females["Height"].As<double>());

            TestResult allNormal = view["Height"].As<double>().ShapiroFranciaTest();
            TestResult maleNormal = males["Height"].As<double>().ShapiroFranciaTest();
            TestResult femaleNormal = females["Height"].As<double>().ShapiroFranciaTest();

            TestResult tTest = Univariate.StudentTTest(males["Height"].As<double>(), females["Height"].As<double>());
            TestResult mwTest = Univariate.MannWhitneyTest(males["Height"].As<double>(), females["Height"].As<double>());

            LinearRegressionResult result0 = males["Weight"].As<double>().LinearRegression(males["Height"].As<double>());
            PolynomialRegressionResult result1 = males["Height"].As<double>().PolynomialRegression(males["Weight"].As<double>(), 1);
            PolynomialRegressionResult result2 = males["Height"].As<double>().PolynomialRegression(males["Weight"].As<double>(), 2);
            PolynomialRegressionResult result3 = males["Height"].As<double>().PolynomialRegression(males["Weight"].As<double>(), 3);

            //MultiLinearRegressionResult multi = view["Weight"].As<double>().MultiLinearRegression(view["Height"].As<double>(), view["Gender"].As<string>().Select(s => (s == "Male") ? 1.0 : 0.0).ToList());

        }

        [TestMethod]
        public void InternetTimeSeriesDownload () {

            FrameTable table = DownloadFrameTable(new Uri("https://timeseries.weebly.com/uploads/2/1/0/8/21086414/sea_ice.csv"));

            double[] powerSpectrum = table["Arctic"].As<double>().PowerSpectrum();
            double v12 = table["Arctic"].As<double>().Autocovariance(12);
            TestResult lbTest = table["Arctic"].As<double>().LjungBoxTest();
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
