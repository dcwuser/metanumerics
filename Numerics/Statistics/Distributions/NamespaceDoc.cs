using System.Runtime.CompilerServices;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Contains types that describe probability distributions.
    /// </summary>
    /// <remarks>
    /// <para>Distributions are assignments of a probability-weight to each of the elements in a set. Most commonly,
    /// those sets are subsets of the integers or real numbers.</para>
    /// <para>Distribution on the integers inherit from the abstract <see cref="DiscreteDistribution"/> class. For any discrete
    /// distribution, you can determine its range (called <see cref="DiscreteDistribution.Support"/>), the
    /// probability weight of each value (using <see cref="DiscreteDistribution.ProbabilityMass(int)"/>),
    /// and many other properties. You can generate pseduo-random integers distributed according to a
    /// discrete distribution (using <see cref="DiscreteDistribution.GetRandomValue(System.Random)"/>).
    /// Many discrete distributions are defined, including <see cref="PoissonDistribution"/> and
    /// <see cref="BinomialDistribution"/>.</para>
    /// <para>Distributions on the real numbers inherit from the abstract <see cref="ContinuousDistribution"/> class. For
    /// any continuous distribution, you can determine its range (called <see cref="ContinuousDistribution.Support"/>),
    /// the probability density at each value (using <see cref="ContinuousDistribution.ProbabilityDensity(double)"/>),
    /// the cumulative distribution function (using <see cref="ContinuousDistribution.LeftProbability(double)"/>),
    /// and many other properties. You can generate pseudo-random floating-point values distributed according to
    /// a continuous distribution (using <see cref="ContinuousDistribution.GetRandomValue(System.Random)"/>).
    /// Many continuous distributions are defined, including <see cref="NormalDistribution"/>,
    /// <see cref="BetaDistribution"/>, <see cref="GammaDistribution"/>, and <see cref="WeibullDistribution"/>.</para>
    /// <para>All one-dimensional distibutions, continuous and discrete, inherit from the abstract <see cref="UnivariateDistribution"/> class.
    /// Using the properties and methods of this class, you can determine raw moments (<see cref="UnivariateDistribution.RawMoment(int)"/>)
    /// such as the <see cref="UnivariateDistribution.Mean"/>, central moments (<see cref="UnivariateDistribution.CentralMoment(int)"/>)
    /// such as the <see cref="UnivariateDistribution.Variance"/>, or cumulants (<see cref="UnivariateDistribution.Cumulant(int)"/>).</para>
    /// <para>Many distributions also offer methods that allow you to find the parameters that best fit a given set of data points
    /// and measure the quality of the fit.</para>
    /// <para>You can add your own continous and discrete distributions by inheriting from <see cref="ContinuousDistribution"/> or
    /// <see cref="DiscreteDistribution"/> and implementing only a few abstract methods. All the remaining properties and methods
    /// are then automatically determined for your distribution.</para>
    /// </remarks>
    /// <seealso href="https://en.wikipedia.org/wiki/Probability_distribution"/>
    [CompilerGenerated]
    internal class NamespaceDoc {
    }
}
