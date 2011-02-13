using System;
using System.Text;
using System.IO;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {
    [TestClass]
    public class RandomTest {

        //[TestMethod]
        public void TestMethod1 () {

            FileStream s = File.OpenWrite(@"C:\Users\dawright\downloads\diehard\random.bin");

            Random rng = new Random(1);
            byte[] buffer = new byte[1024];
            for (int i = 0; i < 16384; i++) {
                rng.NextBytes(buffer);
                s.Write(buffer, 0, buffer.Length);
            }

            s.Close();
        }

    }
}
