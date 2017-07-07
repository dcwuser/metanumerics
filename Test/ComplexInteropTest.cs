using System;
using Meta.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using NetComplex = System.Numerics.Complex;


namespace Test
{
    [TestClass]
    public class ComplexInteropTest
    {

        [TestMethod]
        public void ComplexNetToMy ()
        {
            NetComplex netComplex = new NetComplex(1.0, -2.0);
            Complex myComplex = netComplex;
            Assert.IsTrue(netComplex.Real == myComplex.Re);
            Assert.IsTrue(netComplex.Imaginary == myComplex.Im);
            Assert.IsTrue(ComplexMath.Sqrt(netComplex) == ComplexMath.Sqrt(myComplex));
        }

    }
}
