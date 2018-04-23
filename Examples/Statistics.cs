using System;
using System.Collections.Generic;

using Meta.Numerics;
using Meta.Numerics.Matrices;
using Meta.Numerics.Statistics;

namespace Examples {
    
    public static class Statistics {

        [ExampleMethod]
        private static void CompareSamples () {

            List<double> a = new List<double>() { 130.0, 140.0, 150.0, 150.0, 160.0, 190.0 };
            List<double> b = new List<double>() { 120.0, 150.0, 180.0, 170.0, 185.0, 175.0, 190.0, 200.0 };

            TestResult student = Univariate.StudentTTest(a, b);
            Console.WriteLine($"{student.Statistic.Name} = {student.Statistic.Value}");
            Console.WriteLine($"{student.Type} P = {student.Probability}");

            student.Type = TestType.LeftTailed;
            Console.WriteLine($"{student.Type} P = {student.Probability}");

            TestResult mannWhitney = Univariate.MannWhitneyTest(a, b);
            Console.WriteLine($"{mannWhitney.Statistic.Name} = {mannWhitney.Statistic.Value}");
            Console.WriteLine($"{mannWhitney.Type} P = {mannWhitney.Probability}");

            TestResult kolmogorov = Univariate.KolmogorovSmirnovTest(a, b);
            Console.WriteLine($"{kolmogorov.Statistic.Name} = {kolmogorov.Statistic.Value}");
            Console.WriteLine($"{kolmogorov.Type} P = {kolmogorov.Probability}");
        }

        [ExampleMethod]
        public static void Association () {

            double[] x = new double[] {-0.58, 0.92, 1.41, 1.62, 2.72, 3.14 };
            double[] y = new double[] {1.00, 0.00, 2.00, 16.00, 18.0, 20.0 };

            TestResult pearson = Bivariate.PearsonRTest(x, y);
            Console.WriteLine($"Pearson {pearson.Statistic.Name} = {pearson.Statistic.Value}");
            Console.WriteLine($"{pearson.Type} P = {pearson.Probability}");

            TestResult spearman = Bivariate.SpearmanRhoTest(x, y);
            Console.WriteLine($"Spearman {spearman.Statistic.Name} = {spearman.Statistic.Value}");
            Console.WriteLine($"{spearman.Type} P = {spearman.Probability}");

            TestResult kendall = Bivariate.KendallTauTest(x, y);
            Console.WriteLine($"Kendall {kendall.Statistic.Name} = {kendall.Statistic.Value}");
            Console.WriteLine($"{kendall.Type} P = {kendall.Probability}");

        }

        [ExampleMethod]
        public static void LinearRegression () {

            List<double> x = new List<double>() {-1.1, 2.2, 1.4, 0.5, 3.7, 2.8};
            List<double> y = new List<double>() {-2.9, 3.4, 0.9, 0.1, 6.8, 5.7}; 

            LinearRegressionResult result = y.LinearRegression(x);
            Console.WriteLine($"y = ({result.Intercept}) + ({result.Slope}) x");

            Console.WriteLine($"Fit explains {result.RSquared * 100.0}% of the variance");

            Console.WriteLine($"Probability of no dependence {result.R.Probability}.");

            OneWayAnovaResult anova = result.Anova;
            Console.WriteLine("Fit        dof = {0} SS = {1}", anova.Factor.DegreesOfFreedom, anova.Factor.SumOfSquares);
            Console.WriteLine("Residual   dof = {0} SS = {1}", anova.Residual.DegreesOfFreedom, anova.Residual.SumOfSquares);
            Console.WriteLine("Total      dof = {0} SS = {1}", anova.Total.DegreesOfFreedom, anova.Total.SumOfSquares);
            Console.WriteLine($"Probability of no dependence {anova.Result.Probability}.");

            // Print a 95% confidence interval on the slope
            Console.WriteLine($"slope is in {result.Slope.ConfidenceInterval(0.95)} with 95% confidence");

            IReadOnlyList<double> residuals = result.Residuals;

            ColumnVector parameters = result.Parameters.ValuesVector;
            SymmetricMatrix covariance = result.Parameters.CovarianceMatrix;

            result.Parameters.CovarianceOf("Intercept", "Slope");

            double x1 = 3.0;
            UncertainValue y1 = result.Predict(x1);
            Console.WriteLine($"Predicted y({x1}) = {y1}.");

        }

        [ExampleMethod]
        public static void ContingencyTable () {
            ContingencyTable<string, bool> contingency = new ContingencyTable<string, bool>(
                new string[] { "P", "N" }, new bool[] { true, false }
            );
            contingency["P", true] = 35;
            contingency["P", false] = 65;
            contingency["N", true] = 4;
            contingency["N", false] = 896;

            IReadOnlyList<string> x = new string[] { "N", "P", "N", "N", "P", "N", "N", "N", "P" };
            IReadOnlyList<bool> y = new bool[] { false, false, false, true, true, false, false, false, true };
            ContingencyTable<string, bool> contingencyFromLists = Bivariate.Crosstabs(x, y);

            foreach (string row in contingency.Rows) {
                Console.WriteLine($"Total count of {row}: {contingency.RowTotal(row)}");
            }
            foreach (bool column in contingency.Columns) {
                Console.WriteLine($"Total count of {column}: {contingency.ColumnTotal(column)}");
            }
            Console.WriteLine($"Total counts: {contingency.Total}");

            foreach (string row in contingency.Rows) {
                UncertainValue probability = contingency.ProbabilityOfRow(row);
                Console.WriteLine($"Estimated probability of {row}: {probability}");
            }
            foreach (bool column in contingency.Columns) {
                UncertainValue probability = contingency.ProbabilityOfColumn(column);
                Console.WriteLine($"Estimated probablity of {column}: {probability}");
            }

            UncertainValue sensitivity = contingency.ProbabilityOfRowConditionalOnColumn("P", true);
            Console.WriteLine($"Chance of P result given true condition: {sensitivity}");
            UncertainValue precision = contingency.ProbabilityOfColumnConditionalOnRow(true, "P");
            Console.WriteLine($"Chance of true condition given P result: {precision}");

            UncertainValue logOddsRatio = contingency.Binary.LogOddsRatio;
            Console.WriteLine($"log(r) = {logOddsRatio}");

            TestResult pearson = contingency.PearsonChiSquaredTest();
            Console.WriteLine($"Pearson χ² = {pearson.Statistic.Value} has P = {pearson.Probability}.");

            TestResult fisher = contingency.Binary.FisherExactTest();
            Console.WriteLine($"Fisher exact test has P = {fisher.Probability}.");

        }

    }

}