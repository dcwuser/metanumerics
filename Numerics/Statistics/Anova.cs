using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Statistics.Distributions;

namespace Meta.Numerics.Statistics {

    // This interface returns that data that an AnovaTestRow needs from other rows in order
    // to do an F-test. 
    internal interface IAnovaResult {

        AnovaRow Residual { get; }

        AnovaRow Total { get; }

    }

    /// <summary>
    /// The result of a one-way ANOVA test.
    /// </summary>
    /// <remarks>
    /// <para>A one way ANOVA test detects the influence of a single factor on the mean of a measured variable, which is assumed
    /// to be normally distributed.</para>
    /// <para>A one way ANOVA result is returned by the static <see cref="Univariate.OneWayAnovaTest(System.Collections.Generic.IReadOnlyCollection{double}[])"/>
    /// method.</para>
    /// <para>Fundamentally, a one-way ANOVA is a simple statistical test like any other, with a single test statistic (F) and
    /// a single associated null distribution (the <see cref="FisherDistribution"/>), but some ANOVA users like to examine and
    /// report intermediate quantities used in the computation of the test.
    /// In particular, the sum of square deviations and degrees of freedom associated
    /// with the design factor and the residual, and their sum may be of interest. Each of these appear as rows in the common tabular
    /// representation of an ANOVA. To enable this, the class makes this information available as <see cref="AnovaRow"/> objects
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
    /// OneWayAnovaResult anova = Sample.OneWayAnovaTest(group1, group2, group3);
    /// Console.WriteLine("P = {0}", anova.Result.Probability);
    /// </code>
    /// </example>
    /// <seealso cref="Univariate.OneWayAnovaTest(System.Collections.Generic.IReadOnlyCollection{double}[])"/>
    public sealed class OneWayAnovaResult : IAnovaResult {

        internal OneWayAnovaResult (AnovaRow factor, AnovaRow residual, AnovaRow total) {
            Debug.Assert(factor != null);
            Debug.Assert(residual != null);
            Debug.Assert(total != null);
            Debug.Assert(factor.DegreesOfFreedom + residual.DegreesOfFreedom == total.DegreesOfFreedom);
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
                return Factor.Result;
            }
        }

    }

    /// <summary>
    /// A row in an analysis of variance (ANOVA) table.
    /// </summary>
    /// <remarks>
    /// <para>An ANOVA separates the variance associated with one or more sources into "rows", each of which has an associated
    /// sum of square deviations and number of degrees of freedom.</para>
    /// <para>It is used in properties of the <see cref="OneWayAnovaResult"/> class.</para>
    /// </remarks>
    /// <seealso cref="OneWayAnovaResult"/>
    public class AnovaRow {

        internal AnovaRow (double sumOfSquares, int degreesOfFreedom) {
            Debug.Assert(sumOfSquares >= 0.0);
            Debug.Assert(degreesOfFreedom >= 0);
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

        internal AnovaTestRow (double sumOfSquares, int degreesOfFreedom, IAnovaResult result) : base(sumOfSquares, degreesOfFreedom) {
            Debug.Assert(result != null);
            this.result = result;
        }

        private IAnovaResult result;

        /// <summary>
        /// Gets the result of an F-test measuring the significance of the row.
        /// </summary>
        public TestResult Result {
            get {
                double F = (this.SumOfSquares / this.DegreesOfFreedom) / (result.Residual.SumOfSquares / result.Residual.DegreesOfFreedom);
                ContinuousDistribution D = new FisherDistribution(this.DegreesOfFreedom, result.Residual.DegreesOfFreedom);
                return new TestResult("F", F, D, TestType.RightTailed);
            }
        }

    }

    /// <summary>
    /// Represents the result of a two-factor analysis of variance.
    /// </summary>
    /// <seealso cref="Univariate.TwoWayAnovaTest(IReadOnlyCollection{double}[,])"/>
    public sealed class TwoWayAnovaResult : IAnovaResult {

        internal TwoWayAnovaResult (AnovaRow row, AnovaRow column, AnovaRow interaction, AnovaRow residual) {
            Debug.Assert(row != null);
            Debug.Assert(column != null);
            Debug.Assert(interaction != null);
            Debug.Assert(residual != null);
            this.row = row;
            this.column = column;
            this.interaction = interaction;
            this.residual = residual;
            this.total = new AnovaRow(
                row.SumOfSquares + column.SumOfSquares + interaction.SumOfSquares + residual.SumOfSquares,
                row.DegreesOfFreedom + column.DegreesOfFreedom + interaction.DegreesOfFreedom + residual.DegreesOfFreedom
            );
        }

        private readonly AnovaRow row;
        private readonly AnovaRow column;
        private readonly AnovaRow interaction;
        private readonly AnovaRow residual;
        private readonly AnovaRow total;

        /// <summary>
        /// Gets the variance associated with the effect of the row factor.
        /// </summary>
        public AnovaTestRow RowFactor {
            get {
                return(new AnovaTestRow(row.SumOfSquares, row.DegreesOfFreedom, this));
            }
        }

        /// <summary>
        /// Gets the variance associated with the effect of the column factor. 
        /// </summary>
        public AnovaTestRow ColumnFactor {
            get {
                return (new AnovaTestRow(column.SumOfSquares, column.DegreesOfFreedom, this));
            }
        }

        /// <summary>
        /// Gets the variance associated with the effect of the interaction of the row and column factors.
        /// </summary>
        public AnovaTestRow Interaction {
            get {
                return (new AnovaTestRow(interaction.SumOfSquares, interaction.DegreesOfFreedom, this));
            }
        }

        /// <summary>
        /// Gets the variance associated with the unpredicted residuals.
        /// </summary>
        public AnovaRow Residual {
            get {
                return (residual);
            }
        }

        /// <summary>
        /// Gets the total variance.
        /// </summary>
        public AnovaRow Total {
            get {
                return (total);
            }
        }

    } 

}