using System;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;
using Meta.Numerics.Statistics.Distributions;
using System.Linq;

namespace Test {

    [TestClass]
    public class AnovaTest {


        [TestMethod]
        public void TwoWayAnova () {

            // We will construct a 3 X 2 two factor model, with row and column effects
            // but no interaction effect. We should detect this with a two-way ANOVA.

            Random rng = new Random(1);

            double[,][] samples = new double[3, 2][]; 
            for (int r = 0; r < 3; r++) {
                for (int c = 0; c < 2; c++) {
                    double mu = 1.0;

                    if (c == 0) {
                        mu -= 2.0;                        
                    } else if (c == 1) {
                        mu += 2.0;                        
                    }

                    if (r == 1) {
                        mu -= 3.0;
                    } else if (r == 2) {
                        mu += 3.0;
                    }

                    NormalDistribution sDistribution = new NormalDistribution(mu, 4.0);
                    samples[r, c] = sDistribution.GetRandomValues(rng, 24).ToArray();
                }
            }

            TwoWayAnovaResult result = Univariate.TwoWayAnovaTest(samples);
            Assert.IsTrue(result.RowFactor.Result.Probability < 0.05);
            Assert.IsTrue(result.ColumnFactor.Result.Probability < 0.05);
            Assert.IsTrue(result.Interaction.Result.Probability > 0.05);

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                result.RowFactor.SumOfSquares + result.ColumnFactor.SumOfSquares + result.Interaction.SumOfSquares + result.Residual.SumOfSquares,
                result.Total.SumOfSquares
            ));

            Assert.IsTrue(TestUtilities.IsNearlyEqual(
                result.RowFactor.DegreesOfFreedom + result.ColumnFactor.DegreesOfFreedom + result.Interaction.DegreesOfFreedom + result.Residual.DegreesOfFreedom,
                result.Total.DegreesOfFreedom
            ));

        }

    }
}
