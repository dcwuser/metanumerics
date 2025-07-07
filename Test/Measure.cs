using System;
using System.Collections.Generic;
using System.Linq;

using Microsoft.VisualStudio.TestTools.UnitTesting;
using FluentAssertions;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Extended;
using Meta.Numerics.Statistics;
using System.Data;


namespace Test {
    
    [TestClass]
    public class Measure {

        [TestMethod]
        public void MeasureGammaSeries () {

            //GammaSeries.LogGammaTwoPlus(0.25);
            //Stirling.LogGamma(24);

            double e_abs_g = 0.0;
            double e_abs_z = 0.0;
            double e_rel_g = 0.0;
            double e_rel_z = 0.0;

            List<double> e_rel_s = new List<double>();
            List<double> e_abs_s = new List<double>();

            Random rng = new Random(1);
            foreach (double z in TestUtilities.GenerateRealValues(1.5, 2.5, rng).Take(4000)) {

                //double z = (rng.Next() < 0.5 ? -zz : zz);

                //double g1 = GammaSeries.LogGammaTwoPlus(z);
                //double g1 = Stirling.Gamma(z);
                //double g1 = Lanczos.Gamma(z);
                double g1 = AdvancedMath.LogGamma(z);

                DoubleDouble g2 = AdvancedDoubleDoubleMath.LogGamma(z);
                //if (Math.Abs(z) < 0.5) {
                //    g2 = AdvancedDoubleDoubleMath.LogGammaTwoPlus(z);
                //} else {
                //    g2 = AdvancedDoubleDoubleMath.LogGamma(2.0 + z);
                //}

                DoubleDouble e_abs = DoubleDouble.Abs(g1 - g2);
                DoubleDouble e_rel = DoubleDouble.Abs(e_abs / g2);

                e_abs_s.Add((double) e_abs);
                e_rel_s.Add((double) e_rel);

                if (e_abs > e_abs_g) {
                    e_abs_g = (double) e_abs;
                    e_abs_z = z;
                }

                if (e_rel > e_rel_g) {
                    e_rel_g = (double)e_rel;
                    e_rel_z = z;
                }

            }

            double e_rel_m = e_rel_s.Mean();
            double e_abs_m = e_abs_s.Mean();


        }

    }
}
