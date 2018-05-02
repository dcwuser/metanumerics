using System;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics.Distributions;

namespace Examples {

    public static class Distributions {

        [ExampleMethod]
        public static void DistributionFunctions () {

            ContinuousDistribution gumbel = new GumbelDistribution();

            // Use PDF to compute absolute deviation
            IntegrationResult r = FunctionMath.Integrate(
                z => gumbel.ProbabilityDensity(z) * Math.Abs(z - gumbel.Mean),
                gumbel.Support
            );
            Console.WriteLine($"mean absolute deviation = {r.Value}");

            // Shorter form
            double gumbelMad = gumbel.ExpectationValue(z => Math.Abs(z - gumbel.Mean));
            Console.WriteLine($"mean absolute deviation = {gumbelMad}");

            double x = 1.5;

            // PDF
            Console.WriteLine($"p({x}) = {gumbel.ProbabilityDensity(x)}");

            // CDF, aka percentile
            double P = gumbel.LeftProbability(x);
            Console.WriteLine($"P({x}) = {P}");

            // Right CDF 
            double Q = gumbel.RightProbability(x);
            Console.WriteLine($"Q({x}) = {Q}");

            Console.WriteLine($"P + Q = {P + Q}");

            // Far tail
            double xt = 100.0;
            double qt = gumbel.RightProbability(xt);
            Console.WriteLine($"Q({xt}) = {qt}");

             // Inverse CDF, aka quantile
            Console.WriteLine($"PI({P}) = {gumbel.InverseLeftProbability(P)}");
            Console.WriteLine($"QI({qt} = {gumbel.InverseRightProbability(qt)}");


            DiscreteDistribution binomial = new BinomialDistribution(0.4, 8);

            Console.WriteLine($"support {binomial.Support}");

            int k = 4;
            Console.WriteLine($"P({k}) = {binomial.ProbabilityMass(k)}");

            double binomialMad = binomial.ExpectationValue(i => Math.Abs(i - binomial.Mean));
            Console.WriteLine($"mean absolute deviation = {binomialMad}");

            Console.WriteLine($"P(k < {k}) = {binomial.LeftExclusiveProbability(k)}");
            Console.WriteLine($"P(k <= {k}) = {binomial.LeftInclusiveProbability(k)}");
            Console.WriteLine($"P(k > {k}) = {binomial.RightExclusiveProbability(k)}");

            int k0 = binomial.InverseLeftProbability(0.5);
            Console.WriteLine($"min k0 to achieve P(k <= k0) > 0.5: {k0}");
            Console.WriteLine($"P(k < {k0}) = {binomial.LeftExclusiveProbability(k0)}");
            Console.WriteLine($"P(k <= {k0}) = {binomial.LeftInclusiveProbability(k0)}");

         }


        [ExampleMethod]
        public static void DistributionMoments () {

            //ContinuousDistribution d = new GumbelDistribution();
            DiscreteDistribution d = new PoissonDistribution(5);
            Console.WriteLine($"support = {d.Support}");

            Console.WriteLine($"mean = {d.Mean}");
            Console.WriteLine($"mean as expectation = {d.ExpectationValue(x => x)}");

            Console.WriteLine($"variance = {d.Variance}");
            Console.WriteLine($"variance as expectation = {d.ExpectationValue(x => MoreMath.Sqr(x - d.Mean))}");

            Console.WriteLine($"standard deviation = {d.StandardDeviation}");
            Console.WriteLine($"skewness = {d.Skewness}");
            Console.WriteLine($"excess kuritosis = {d.ExcessKurtosis}");

            for (int r = 0; r <= 4; r++) {
                Console.WriteLine($"M_{r} = {d.RawMoment(r)}");
                Console.WriteLine($"C_{r} = {d.CentralMoment(r)}");
                Console.WriteLine($"K_{r} = {d.Cumulant(r)}");
            }

        }

    }

}