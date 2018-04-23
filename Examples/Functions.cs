using System;

using Meta.Numerics;
using Meta.Numerics.Functions;

namespace Examples {
    
    public static class Functions {

        [ExampleMethod]
        public static void ComputeSpecialFunctions () {

            // Compute the value x at which erf(x) is just 10^{-15} from 1.
            double x = AdvancedMath.InverseErfc(1.0E-15);

            // The Gamma function at 1/2 is sqrt(pi)
            double y = AdvancedMath.Gamma(0.5);

            // Compute a Coulomb Wave Function in the quantum tunneling region
            SolutionPair s = AdvancedMath.Coulomb(2, 4.5, 3.0);

            // Compute the Reimann Zeta function at a complex value
            Complex z = AdvancedComplexMath.RiemannZeta(new Complex(0.75, 6.0));

            // Compute the 100th central binomial coefficient
            double c = AdvancedIntegerMath.BinomialCoefficient(100, 50);            
        }

    }

}