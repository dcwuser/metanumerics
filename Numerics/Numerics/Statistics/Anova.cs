using System;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    /// <summary>
    /// The result of a one-way ANOVA test.
    /// </summary>
    /// <remarks>
    /// <para>A one way ANOVA test detects the influence of a categorical factor on the mean of a measured variable, which is assumed
    /// to be normally distributed.</para>
    /// <para>A one way ANOVA result is returned by the static <see cref="Sample.OneWayAnovaTest(System.Collections.Generic.IList{Sample})"/>
    /// method.</para>
    /// <para>While, fundamentally, a one-way ANOVA is a simple statistical test like any other, with a single test statistic (F) and
    /// a single associated distribution (the F distribution), some ANOVA users like to examine and report some intermediate quantities
    /// used in the computation of the test. In particular, the sum of square deviations and corresponding degrees of freedom associated
    /// with the design factor and the residual, and their sum may be of interest. Each of these appear as rows in the common tabular
    /// representation of an ANOVA test. To enable this, the class makes this information available as <see cref="AnovaRow"/> objects
    /// returned by the <see cref="Factor"/>, <see cref="Residual"/>, and <see cref="Total"/> properties. This has the unfortunate
    /// side-effect of making the AVOVA look more complicated than it really is. If you just want the test result, you can get it
    /// from the <see cref="Result"/> property.</para>
    /// </remarks>
    /// <example>
    /// <para>Suppose you have sampled the heights of aliens from three planets. Heights are approximately normally distributed
    /// on each planet. You want to know whether planet-of-origin affects average height. You can do a one-way ANOVA to determine
    /// if the planet factor affects mean height.</para>
    /// <code lang="C#">
    /// Sample group1 = new Sample(4, 5, 6);
    /// Sample group2 = new Sample(3, 4, 5);
    /// Sample group3 = new Sample(5, 6, 8, 9);
    /// OneWayAnovaResult result = Sample.OneWayAnovaTest(group1, group2, group3);
    /// return(result.Result.RightProbability);
    /// </code>
    /// </example>
    /// <seealso cref="Sample.OneWayAnovaTest(System.Collections.Generic.IList{Sample})"/>
    public class OneWayAnovaResult {

        internal OneWayAnovaResult (AnovaRow factor, AnovaRow residual, AnovaRow total) {

            this.Factor = new AnovaTestRow(factor.SumOfSquares, factor.DegreesOfFreedom, this);
            this.Residual = residual;
            this.Total = total;

        }

        /// <summary>
        /// Gets design factor variance data.
        /// </summary>
        public AnovaTestRow Factor { get; private set; }

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
        public TestResult Result {
            get {
                return (Factor.Result);
            }
        }

    }

    /// <summary>
    /// A row in an ANOVA table.
    /// </summary>
    /// <remarks>
    /// <para>An ANOVA seperates the variance associated with one or more sources into "rows", each of which has an associated
    /// sum of square deviations and number of degrees of freedom.</para>
    /// <para>It is used in properties of the <see cref="OneWayAnovaResult"/> class.</para>
    /// </remarks>
    /// <seealso cref="OneWayAnovaResult"/>
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
    /// A row in an ANOVA table for which an F-test is available.
    /// </summary>
    /// <seealso cref="OneWayAnovaResult"/>
    public class AnovaTestRow : AnovaRow {

        internal AnovaTestRow (double sumOfSquares, int degreesOfFreedom, OneWayAnovaResult result) : base(sumOfSquares, degreesOfFreedom) {
            this.result = result;
        }

        private OneWayAnovaResult result;

        /// <summary>
        /// Gets the result of an F-test measuring the significance of the row.
        /// </summary>
        public TestResult Result {
            get {
                double F = (this.SumOfSquares / this.DegreesOfFreedom) / (result.Residual.SumOfSquares / result.Residual.DegreesOfFreedom);
                Distribution D = new FisherDistribution(this.DegreesOfFreedom, result.Residual.DegreesOfFreedom);
                return (new TestResult("F", F, TestType.TwoTailed, D));
            }
        }

    }

}