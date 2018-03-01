using System;
using System.Diagnostics;

using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Data;
using Meta.Numerics.Statistics;

namespace DataTest
{
    [TestClass]
    public class DataListTests
    {

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
