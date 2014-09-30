using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics;
using Meta.Numerics.Analysis;
using Meta.Numerics.Functions;
using Meta.Numerics.Statistics;

namespace Test {


    internal class TestIntegral {

        public TestIntegral (Func<double, double> f, Interval r, double I) {
            this.Integrand = f;
            this.Range = r;
            this.Result = I;
        }

        public Func<double, double> Integrand { get; private set; }

        public Interval Range { get; private set; }

        public double Result { get; private set; }


    }

    /// <summary>
    /// Summary description for IntegrateTest
    /// </summary>
    [TestClass]
    public class IntegrateTest {
        
        [TestMethod]
        public void IntegrateTestIntegrals () {

            Console.WriteLine(integrals.Length);

            for (int i = 0; i < integrals.Length; i++) {
                TestIntegral integral = integrals[i];
                double result = FunctionMath.Integrate(integral.Integrand, integral.Range);
                Assert.IsTrue(TestUtilities.IsNearlyEqual(result, integral.Result));
            }

        }

        private TestIntegral[] integrals = new TestIntegral[] {

            // smooth functions

            // to add:
            // low-order polynomial
            // high-order polynomial
            // rational function
            // rational function involving roots
            // multiple scales

            new TestIntegral(
                delegate(double x) {
                    return (Math.Sqrt(x));
                },
                Interval.FromEndpoints(0.0, 1.0),
                2.0 / 3.0
            ),

            new TestIntegral(
                delegate (double x) {
			        return (1.0 / (1.0 + x));
                },
		        Interval.FromEndpoints(0,1),
                Math.Log(2.0)
            ),

            new TestIntegral(
                delegate (double x) {
			        return(1.0 / (1.0 + x * x));
		        },
		        Interval.FromEndpoints(-1, 1),
                Math.PI / 2.0
            ),

            new TestIntegral(
                delegate (double x) {
			        return( Math.Exp(x) * Math.Cos(x) );
		        },
                Interval.FromEndpoints(0.0, Math.PI),
                -(Math.Exp(Math.PI) + 1.0) / 2.0
            ),

            new TestIntegral(
                delegate (double x) {
			        return( Math.Sqrt( 1.0 - x * x) );
		        },
                Interval.FromEndpoints(0.0, 1.0),
                Math.PI / 4.0
            ),
            
            new TestIntegral(
                delegate (double x) {
                    return( Math.Cos(Math.PI * x / 2.0) );
                },
                Interval.FromEndpoints(-1.0,1.0),
                4.0 / Math.PI
            ),

            new TestIntegral(
                delegate (double t) {
                    return( Math.Log ( 1.0 + Math.Tan(t) ) );
                },
                Interval.FromEndpoints(0,Math.PI/4.0),
                Math.PI * Math.Log(2.0) / 8.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return(1.0/x);
                },
                Interval.FromEndpoints(1.0, Math.E),
                1.0
            ),

            new TestIntegral(
                delegate (double x) {
			        return( Math.Pow(x, -x) );
		        },
                Interval.FromEndpoints(0.0, 1.0),
                1.2912859970626635404 // 1 + 1/2^2 + 1/3^3 + 1/4^4 + 1/5^5 + ...
            ),

            new TestIntegral(
                delegate (double x) {
                    return( Math.Pow(x, x) ); // has minimum at 
                },
                Interval.FromEndpoints(0.0, 1.0),
                0.78343051071213440706 // 1 - 1/2^2 + 1/3^3 - 1/4^4 + 1/5^5 + ...
            ),

            // nondifferentiable points

            new TestIntegral(
                delegate (double x) {
			        double s = Math.Sin(x);
			    return(s*s);
		        },
		        Interval.FromEndpoints(0.0, 2.0 * Math.PI),
                Math.PI
            ),

            new TestIntegral(
                delegate (double x) {
			        return( Math.Sqrt( Math.Abs( x + 0.5 ) ) );
		        },
		        Interval.FromEndpoints(-1.0,1.0),
                (Math.Sqrt(2.0) + 3.0 * Math.Sqrt(6.0)) / 6.0
            ),

            // discontinuities

            new TestIntegral(
                delegate (double x) {
			        if (x < 0) {
				        return(0.0);
			        } else {
				        return(1.0);
		    	    }
		        },
		        Interval.FromEndpoints(-2.0,1.0),
                1.0
            ),

            // log singularities

            new TestIntegral(
		        delegate (double x) {
			        return( Math.Log(x) );
		        },
		        Interval.FromEndpoints(0,1),
                -1.0
            ),

            new TestIntegral(
		        delegate (double x) {
			        return( Math.Sqrt(x) * Math.Log(x) );
		        },
		        Interval.FromEndpoints(0,1),
                -4.0 / 9.0
            ),

            new TestIntegral(
		        delegate (double x) {
			        return( Math.Sqrt( - Math.Log(x) ) );
		        },
		        Interval.FromEndpoints(0,1),
                Math.Sqrt(Math.PI) / 2.0
            ),

            new TestIntegral(
                delegate (double x) {
			        return( Math.Log(x) / (1.0 + x * x) );
		        },
		        Interval.FromEndpoints(0,1),
                -0.915965594177219015054  // Catalan constant = DirichletBeta(2)
            ),

            new TestIntegral(
                delegate (double t) {
			        return( t * Math.Log ( Math.Sin(t) ) );
		        },
		        Interval.FromEndpoints(0,Math.PI),
                -Math.PI * Math.PI * Math.Log(2.0) / 2.0
            ),

            new TestIntegral(
                delegate (double t) {
			        return( Math.Sin(t) * Math.Log ( Math.Sin(t) ) );
		        },
		        Interval.FromEndpoints(0,Math.PI / 2.0),
                Math.Log(2.0) - 1.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return(Math.Log(x) * Math.Log(1.0-x));
                },
                Interval.FromEndpoints(0.0,1.0),
                2.0 - Math.PI * Math.PI / 6
            ),

            new TestIntegral(
                delegate (double x) {
                    return( Math.Log(x) / (1.0 + x) );
                },
                Interval.FromEndpoints(0.0, 1.0),
                -Math.PI * Math.PI / 12.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( Math.Log(x) / (1.0 - x) );
                },
                Interval.FromEndpoints(0.0, 1.0),
                -Math.PI * Math.PI / 6.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( Math.Log(x) / (1.0 - x*x) );
                },
                Interval.FromEndpoints(0.0, 1.0),
                -Math.PI * Math.PI / 8.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( Math.Log( -Math.Log(x) ) );
                },
                Interval.FromEndpoints(0.0, 1.0),
                -AdvancedMath.EulerGamma
            ),

            // power law singularities

            new TestIntegral(
                delegate (double x) {
                    return( Math.Pow(x,-1.0/6.0));
                },
                Interval.FromEndpoints(0.0, 1.0),
                6.0 / 5.0
            ),

            // multiple scales

            new TestIntegral(
                delegate (double x) {
                    return( Math.Exp(-x) * Math.Sin(20.0 * x) );
                },
                Interval.FromEndpoints(0.0, 2.0 * Math.PI),
                20.0 / 401.0 * ( 1.0 - Math.Exp(-2.0 * Math.PI) )
            ),

            // semi-infinite

            // power law fall-off

            new TestIntegral(
                delegate (double x) {
                    return( 1.0 / (1.0 + x * x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                Math.PI / 2.0
            ),

            // exponential fall-off
            
            new TestIntegral(
                delegate (double x) {
                    return( Math.Exp(-x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                1.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( Math.Sqrt(x) * Math.Exp(-x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                Math.Sqrt(Math.PI) / 2.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( x / (Math.Exp(x) + 1.0) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                Math.PI * Math.PI / 12.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( x / (Math.Exp(x) - 1.0) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                Math.PI * Math.PI / 6.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( x * x * x / (Math.Exp(x) - 1.0) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                Math.Pow(Math.PI, 4.0) / 15.0
            ),

            // e^(-x^2) fall-off

            new TestIntegral(
                delegate (double x) {
                    return( Math.Exp(-x*x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                Math.Sqrt(Math.PI) / 2.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( x * Math.Exp(-x*x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                1.0 / 2.0
            ),

            // oscilatory

            /*
            // nonconvergent
            new TestIntegral(
                delegate (double x) {
                    return( Math.Cos(x) / (1.0 + x*x) );                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                Math.PI / Math.E / 2.0
            ),
            */

            new TestIntegral(
                delegate (double x) {
                    return( Math.Exp(-x) * Math.Sin(x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                1.0 / 2.0
            ), 

            new TestIntegral(
                delegate (double x) {
                    return( Math.Exp(-x) * Math.Cos(x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                1.0 / 2.0
            ),

            new TestIntegral(
                delegate (double x) {
                    return( Math.Exp(-x*x/2.0) * Math.Cos(x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                Math.Sqrt(Math.PI / Math.E / 2.0)
            ),

            // singular

            new TestIntegral(
                delegate (double x) {
                    return( Math.Exp(-x) * Math.Log(x) );
                },
                Interval.FromEndpoints(0.0, Double.PositiveInfinity),
                -AdvancedMath.EulerGamma
            ),

            // fully infinite

            new TestIntegral(
                delegate (double x) {
                    double z = x - 10.0;
                    return( 1.0 / (1.0 + z*z) );
                },
                Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity),
                Math.PI
            ),

            new TestIntegral(
                delegate (double x) {
                    return( 1.0 / (1.0 + Math.Pow(x, 4.0)) );
                },
                Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity),
                Math.PI / Math.Sqrt(2.0)
            ),

            new TestIntegral(
                delegate (double x) {
                    return( Math.Exp(-Math.Abs(x)) );
                },
                Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity),
                2.0
            ),

            new TestIntegral(
                delegate (double x) {
                    double z = (x + 100.0) / 10.0;
                    return( Math.Exp(-z*z/2.0) );
                },
                Interval.FromEndpoints(Double.NegativeInfinity, Double.PositiveInfinity),
                Math.Sqrt(2.0 * Math.PI) * 10.0
            ),

        };

    }
}
