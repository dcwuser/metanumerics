using System;
using System.Diagnostics;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Data;

namespace DataTest
{
    [TestClass]
    public class DataListTests
    {

        [TestMethod]
        public void Timings()
        {
            int n = 10000;
            int m = 10000;
            Random rng = new Random(1);
            double mean = 0.0;

            double[] array = new double[n];
            for (int j = 0; j < array.Length; j++)
            {
                array[j] = rng.NextDouble();
            }
            Stopwatch arrayTimer = Stopwatch.StartNew();
            for (int i = 0; i < m; i++)
            {
                mean += array.Mean();
            }
            arrayTimer.Stop();
            Console.WriteLine(arrayTimer.ElapsedMilliseconds);

            DataList<double> list = new DataList<double>("list");
            for (int j = 0; j < list.Count; j++)
            {
                list.Add(rng.NextDouble());
            }
            Stopwatch listTimer = Stopwatch.StartNew();
            for (int i = 0; i < m; i++)
            {
                mean += list.Mean();
            }
            listTimer.Stop();
            Console.WriteLine(listTimer.ElapsedMilliseconds);

            DataFrame frame = new DataFrame(list);
            Stopwatch frameTimer = Stopwatch.StartNew();
            for (int i = 0; i < m; i++)
            {
                mean += frame.Column<double>("list").Mean();
            }
            frameTimer.Stop();
            Console.WriteLine(frameTimer.ElapsedMilliseconds);



        }



        [TestMethod]
        public void DataListManipulations ()
        {
            string name = "test";
            DataList<double> list = new DataList<double>(name);
            Assert.IsTrue(list.Name == name);
            Assert.IsTrue(list.Count == 0);
            Assert.IsTrue(list.StorageType == typeof(double));

            list.Add(1.1);
            list.Add(2.2);

            Assert.IsTrue(list.Count == 2);
            Assert.IsTrue(list[0] == 1.1);
            Assert.IsTrue(list[1] == 2.2);

            Assert.IsTrue(list[list.IndexOf(2.2)] == 2.2);

            Assert.IsFalse(list.Contains(0.0));
            Assert.IsFalse(list.Remove(0.0));

            Assert.IsTrue(list.Contains(1.1));
            Assert.IsTrue(list.Remove(1.1));
            Assert.IsFalse(list.Contains(1.1));
            Assert.IsTrue(list.Count == 1);

            list.InsertAt(0, 3.3);
            Assert.IsTrue(list.Contains(3.3));
            Assert.IsTrue(list.Count == 2);

            list.Clear();
            Assert.IsFalse(list.Contains(2.2));
            Assert.IsTrue(list.Count == 0);

        }
    }
}
