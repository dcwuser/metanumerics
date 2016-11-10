using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;

namespace Test {

    [TestClass]
    public class AdvancedMathTest_Bessel {

        [TestMethod]
        public void BesselAtZero () {

            // Normal Bessel \nu = 0

            Assert.IsTrue(AdvancedMath.BesselJ(0, 0.0) == 1.0);
            Assert.IsTrue(AdvancedMath.BesselJ(0.0, 0.0) == 1.0);

            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(0, 0.0)));
            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(0.0, 0.0)));

            SolutionPair jy0 = AdvancedMath.Bessel(0.0, 0.0);
            Assert.IsTrue(jy0.FirstSolutionValue == 1.0);
            Assert.IsTrue(jy0.FirstSolutionDerivative == 0.0);
            Assert.IsTrue(Double.IsNegativeInfinity(jy0.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(jy0.SecondSolutionDerivative));

            // Normal Bessel 0 < \nu < 1

            Assert.IsTrue(AdvancedMath.BesselJ(0.1, 0.0) == 0.0);

            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(0.1, 0.0)));

            SolutionPair jyf = AdvancedMath.Bessel(0.9, 0.0);
            Assert.IsTrue(jyf.FirstSolutionValue == 0.0);
            Assert.IsTrue(jyf.FirstSolutionDerivative == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNegativeInfinity(jyf.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(jyf.SecondSolutionDerivative));

            // Normal Bessel \nu = 1

            Assert.IsTrue(AdvancedMath.BesselJ(1, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.BesselJ(1.0, 0.0) == 0.0);

            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(1, 0.0)));
            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(1.0, 0.0)));

            SolutionPair jy1 = AdvancedMath.Bessel(1.0, 0.0);
            Assert.IsTrue(jy1.FirstSolutionValue == 0.0);
            Assert.IsTrue(jy1.FirstSolutionDerivative == 0.5);
            Assert.IsTrue(Double.IsNegativeInfinity(jy1.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(jy1.SecondSolutionDerivative));

            // Normal Bessel \nu > 1

            Assert.IsTrue(AdvancedMath.BesselJ(2, 0.0) == 0.0);
            Assert.IsTrue(AdvancedMath.BesselJ(1.2, 0.0) == 0.0);

            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(2, 0.0)));
            Assert.IsTrue(Double.IsNegativeInfinity(AdvancedMath.BesselY(1.2, 0.0)));

            SolutionPair jy2 = AdvancedMath.Bessel(1.7, 0.0);
            Assert.IsTrue(jy2.FirstSolutionValue == 0.0);
            Assert.IsTrue(jy2.FirstSolutionDerivative == 0.0);
            Assert.IsTrue(Double.IsNegativeInfinity(jy2.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(jy2.SecondSolutionDerivative));

        }

        [TestMethod]
        public void ModifiedBesselAtZero () {

            // Modified Bessel \nu = 0

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(0.0, 0.0) == 1.0);

            Assert.IsTrue(AdvancedMath.ModifiedBesselK(0.0, 0.0) == Double.NegativeInfinity); 

            SolutionPair ik0 = AdvancedMath.ModifiedBessel(0.0, 0.0);
            Assert.IsTrue(ik0.FirstSolutionValue == 1.0);
            Assert.IsTrue(ik0.FirstSolutionDerivative == 0.0);
            Assert.IsTrue(Double.IsNegativeInfinity(ik0.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(ik0.SecondSolutionDerivative));

            // Modified Bessel 0 < \nu < 1

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(0.1, 0.0) == 0.0);

            Assert.IsTrue(AdvancedMath.ModifiedBesselK(0.1, 0.0) == Double.NegativeInfinity);

            SolutionPair ikf = AdvancedMath.ModifiedBessel(0.9, 0.0);
            Assert.IsTrue(ikf.FirstSolutionValue == 0.0);
            Assert.IsTrue(ikf.FirstSolutionDerivative == Double.PositiveInfinity);
            Assert.IsTrue(Double.IsNegativeInfinity(ikf.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(ikf.SecondSolutionDerivative));

            // Modified Bessel \nu = 1

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(1.0, 0.0) == 0.0);

            Assert.IsTrue(AdvancedMath.ModifiedBesselK(1.0, 0.0) == Double.NegativeInfinity);

            SolutionPair ik1 = AdvancedMath.ModifiedBessel(1.0, 0.0);
            Assert.IsTrue(ik1.FirstSolutionValue == 0.0);
            Assert.IsTrue(ik1.FirstSolutionDerivative == 0.5);
            Assert.IsTrue(Double.IsNegativeInfinity(ik1.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(ik1.SecondSolutionDerivative));

            // Modified Bessel \nu > 1

            Assert.IsTrue(AdvancedMath.ModifiedBesselI(1.2, 0.0) == 0.0);

            Assert.IsTrue(AdvancedMath.ModifiedBesselK(1.2, 0.0) == Double.NegativeInfinity);

            SolutionPair ik2 = AdvancedMath.ModifiedBessel(1.7, 0.0);
            Assert.IsTrue(ik2.FirstSolutionValue == 0.0);
            Assert.IsTrue(ik2.FirstSolutionDerivative == 0.0);
            Assert.IsTrue(Double.IsNegativeInfinity(ik2.SecondSolutionValue));
            Assert.IsTrue(Double.IsPositiveInfinity(ik2.SecondSolutionDerivative));

        }

    }

}