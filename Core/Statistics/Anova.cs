using System;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// The result of a one-way ANOVA test.
    /// </summary>
    public class OneWayAnovaResult {

        internal OneWayAnovaResult (AnovaRow factor, AnovaRow residual, AnovaRow total) {

            this.Factor = new AnovaFactor(factor.SumOfSquares, factor.DegreesOfFreedom, this);
            this.Residual = residual;
            this.Total = total;

        }

        /// <summary>
        /// Gets design factor variance data.
        /// </summary>
        public AnovaFactor Factor { get; private set; }

        /// <summary>
        /// Gets residual variance data.
        /// </summary>
        public AnovaRow Residual { get; private set; }

        /// <summary>
        /// Gets total variance data.
        /// </summary>
        public AnovaRow Total { get; private set; }

        /// <summary>
        /// Gets the result of the F test for the influence of the factor.
        /// </summary>
        public TestResult FTest {
            get {
                return (Factor.FTest);
            }
        }

    }

    /// <summary>
    /// A row in an ANOVA table.
    /// </summary>
    public class AnovaRow {

        internal AnovaRow (double sumOfSquares, int degreesOfFreedom) {
            this.SumOfSquares = sumOfSquares;
            this.DegreesOfFreedom = degreesOfFreedom;
        }

        /// <summary>
        /// Gets the sum of squares contributed by the row.
        /// </summary>
        public double SumOfSquares { get; private set; }

        /// <summary>
        /// Gets the degrees of freedom associated with the row.
        /// </summary>
        public int DegreesOfFreedom { get; private set; }

    }

    /// <summary>
    /// A 
    /// </summary>
    public class AnovaFactor : AnovaRow {

        internal AnovaFactor (double sumOfSquares, int degreesOfFreedom, OneWayAnovaResult result) : base(sumOfSquares, degreesOfFreedom) {
            this.result = result;
        }

        private OneWayAnovaResult result;

        /// <summary>
        /// Gets the result of an F-test measuring the significance of the factor.
        /// </summary>
        public TestResult FTest {
            get {
                double F = (this.SumOfSquares / this.DegreesOfFreedom) / (result.Residual.SumOfSquares / result.Residual.DegreesOfFreedom);
                Distribution D = new FisherDistribution(this.DegreesOfFreedom, result.Residual.DegreesOfFreedom);
                return (new TestResult(F, D));
            }
        }

    }

}