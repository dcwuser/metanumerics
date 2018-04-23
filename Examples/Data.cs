using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;

using Meta.Numerics;
using Meta.Numerics.Data;
using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;

namespace Examples {
    
    public static class Data {

        [ExampleMethod]
        public static void ManipulatingData () {

            FrameTable table;
            Uri url = new Uri("https://raw.githubusercontent.com/dcwuser/metanumerics/master/Examples/Data/example.csv");
            WebRequest request = WebRequest.Create(url);
            using (WebResponse response = request.GetResponse()) {
                using (StreamReader reader = new StreamReader(response.GetResponseStream())) {
                    table = FrameTable.FromCsv(reader);
                }
            }

            FrameView selected = table.Select("Height", "Weight", "Sex");

            FrameView discarded = table.Discard("Name");

            table.AddComputedColumn("Bmi", r => ((double) r["Weight"])/MoreMath.Sqr((double) r["Height"] / 100.0));
            Console.WriteLine($"Bmi of first subject is {table["Bmi"][0]}.");

            FrameView noNulls = table.WhereNotNull();
            FrameView noNullWeights = table.WhereNotNull("Weight");
            FrameView noNullWeightsOrHeights = table.WhereNotNull("Weight", "Height");

            double meanWeight = table.WhereNotNull("Weight").Columns["Weight"].As<double>().Mean();

            FrameView men = table.Where<string>("Sex", s => s == "M");

            FrameView shortMen = table.Where(
                r => ((string) r["Sex"]) == "M" && ((double) r["Height"] < 175.0)
            );

            FrameView ordered = table.OrderBy("Height");

            FrameView reversed = table.OrderBy("Height", SortOrder.Descending);

            FrameView alsoOrdered = table.OrderBy<double>("Height", (h1, h2) => h1.CompareTo(h2));

            FrameView sorted = table.OrderBy((r1, r2) => {
                int first = ((string) r1["Sex"]).CompareTo((string ) r2["Sex"]);
                int second = ((double) r1["Height"]).CompareTo((double) r2["Height"]);
                return first != 0 ? first : second;
            });

            List<string> sexes = table["Sex"].As<string>().Distinct().ToList();

            FrameTable counts = table.GroupBy("Sex", v => v.Rows.Count, "Count");

            FrameTable summarize = table.GroupBy("Sex", v => { 
                SummaryStatistics summary = new SummaryStatistics(v["Height"].As<double>());
                return (new Dictionary<string, object>() {
                    {"Count", summary.Count},
                    {"Mean", summary.Mean},
                    {"StdDev", summary.StandardDeviation}
                });
            });

        }

        [ExampleMethod]
        public static void AnalyzingData () {

            FrameTable table;
            Uri url = new Uri("https://raw.githubusercontent.com/dcwuser/metanumerics/master/Examples/Data/example.csv");
            WebRequest request = WebRequest.Create(url);
            using (WebResponse response = request.GetResponse()) {
                using (StreamReader reader = new StreamReader(response.GetResponseStream())) {
                    table = FrameTable.FromCsv(reader);
                }
            }

            SummaryStatistics summary = new SummaryStatistics(table["Height"].As<double>());
            Console.WriteLine($"Count = {summary.Count}");
            Console.WriteLine($"Mean = {summary.Mean}");
            Console.WriteLine($"Standard Deviation = {summary.StandardDeviation}");
            Console.WriteLine($"Skewness = {summary.Skewness}");
            Console.WriteLine($"Estimated population mean = {summary.PopulationMean}");
            Console.WriteLine($"Estimated population standard deviation = {summary.PopulationStandardDeviation}");

            IReadOnlyList<double> maleHeights =
                table.Where<string>("Sex", s => s == "M").Columns["Height"].As<double>();
            IReadOnlyList<double> femaleHeights =
                table.Where<string>("Sex", s => s == "F").Columns["Height"].As<double>();
            TestResult test = Univariate.StudentTTest(maleHeights, femaleHeights);
            Console.WriteLine($"{test.Statistic.Name} = {test.Statistic.Value}, P = {test.Probability}");

            TestResult maleHeightNormality = maleHeights.ShapiroFranciaTest();
            TestResult totalHeightNormality = table["Height"].As<double>().ShapiroFranciaTest();
            TestResult heightCompatibility = Univariate.KolmogorovSmirnovTest(maleHeights, femaleHeights);

            LinearRegressionResult fit =
                table["Weight"].As<double>().LinearRegression(table["Height"].As<double>());
            Console.WriteLine($"Model weight = ({fit.Slope}) * height + ({fit.Intercept}).");
            Console.WriteLine($"Model explains {fit.RSquared * 100.0}% of variation.");

            ContingencyTable<string, bool> contingency =
                Bivariate.Crosstabs(table["Sex"].As<string>(), table["Result"].As<bool>());
            Console.WriteLine($"Male incidence: {contingency.ProbabilityOfColumnConditionalOnRow(true, "M")}");
            Console.WriteLine($"Female incidence: {contingency.ProbabilityOfColumnConditionalOnRow(false, "F")}");
            Console.WriteLine($"Log odds ratio = {contingency.Binary.LogOddsRatio}");

            table.AddComputedColumn("Bmi", r => ((double) r["Weight"])/MoreMath.Sqr((double) r["Height"] / 100.0));
            table.AddComputedColumn("Age", r=> (DateTime.Now - (DateTime) r["Birthdate"]).TotalDays / 365.24);

            MultiLinearLogisticRegressionResult result = 
                table["Result"].As<bool>().MultiLinearLogisticRegression(
                    table["Bmi"].As<double>(),
                    table["Sex"].As<string, double>(s => s == "M" ? 1.0 : 0.0)
                );
            foreach (Parameter parameter in result.Parameters) {
                Console.WriteLine($"{parameter.Name} = {parameter.Estimate}");
            }

            //TestResult ageResultPearson = Bivariate.PearsonRTest(table["Age"].As<double>(), table["Result"].As<double>());
            TestResult spearman = Bivariate.SpearmanRhoTest(table["Age"].As<double>(), table["Result"].As<double>());
            Console.WriteLine($"{spearman.Statistic.Name} = {spearman.Statistic.Value} P = {spearman.Probability}");

        }        

        public static void ConstructData () {

            FrameTable table = new FrameTable();
            table.AddColumn<int>("Id");
            table.AddColumn<string>("Name");
            table.AddColumn<string>("Sex");
            table.AddColumn<DateTime>("Birthdate");
            table.AddColumns<double>("Height", "Weight");
            table.AddColumn<bool>("Result");

            //Random rng = new Random(3);
            //Random rng = new Random(314159);
            // Random rng = new Random(271828);
            Random rng = new Random(1000001);

            //string[] maleNames = new string[1024];
            string[] maleNames = new string[] {"Alex", "Chris", "David", "Eric", "Frederic", "George", "Hans", "Igor", "John", "Kevin", "Luke", "Mark", "Oscar", "Peter", "Richard", "Stephan", "Thomas", "Vincent" };
            AddRows(table, maleNames, "M", 175.0, 12.0, 24.0, 3.0, 1, rng);

            //string[] femaleNames = new string[1024];
            string[] femaleNames = new string[] {"Anne", "Belle", "Dorothy", "Elizabeth", "Fiona", "Helen", "Julia", "Kate", "Louise", "Mary", "Natalie", "Olivia", "Ruth", "Sarah", "Theresa", "Viola" };
            AddRows(table, femaleNames, "F", 160.0, 10.0, 24.0, 3.0, 0, rng);

            string path = @"C:\Users\dawright\Documents\example.csv";
            using (StreamWriter writer = new StreamWriter(File.OpenWrite(path))) {
                table.ToCsv(writer);
            }
            Console.WriteLine(File.Exists(path));

        }

       private static void AddRows(FrameTable table, IReadOnlyList<string> names, string sex, double meanHeight, double stddevHeight, double meanBmi, double stddevBmi, int flag, Random rng) {

            NormalDistribution gauss = new NormalDistribution();
            UniformDistribution ages = new UniformDistribution(Interval.FromEndpoints(15.0, 75.0));

            foreach (string name in names) {

                double zHeight = gauss.GetRandomValue(rng);
                double height = meanHeight + stddevHeight * zHeight;

                double zBmi = gauss.GetRandomValue(rng);
                double bmi = meanBmi + stddevBmi * zBmi;

                double weight = MoreMath.Sqr(height / 100.0) * bmi;
                
                double t = -0.4 + 0.6 * zBmi + 0.8 * flag;
                double p = 1.0 / (1.0 + Math.Exp(-t));
                bool r = rng.NextDouble() < p;

                int id = table.Rows.Count;

                TimeSpan age = TimeSpan.FromDays(365.24 * ages.GetRandomValue(rng));
                DateTime birthdate = (DateTime.Now - age).Date;

                table.AddRow(id, name, sex, birthdate, height, weight, r);

            }

        }

    }

}