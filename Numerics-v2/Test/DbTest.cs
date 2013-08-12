using System;
using System.Text;
using System.Collections.Generic;
using System.Data;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;

namespace Test {

    [TestClass]
    public class DbTest {

        private static DataTable CreateTable () {

            DataTable table = new DataTable();
            table.Columns.Add("Name", typeof(string));
            table.Columns.Add("Age", typeof(int));
            table.Columns.Add("Day", typeof(DayOfWeek));
            table.Columns.Add("Male", typeof(bool));
            table.Rows.Add("A", 1, DayOfWeek.Monday, false);
            table.Rows.Add("B", 2, DayOfWeek.Tuesday, false);
            table.Rows.Add("C", 3, DayOfWeek.Wednesday, true);
            table.Rows.Add(null, null, null, null);
            return (table);

        }

        [TestMethod]
        public void LoadSample () {

            DataTable table = CreateTable();

            Sample sample = new Sample();
            sample.Load(table.CreateDataReader(), 1);
            Console.Write("{0} {1} {2}", sample.Count, sample.Mean, sample.StandardDeviation);
            Assert.IsTrue(sample.Count == table.Rows.Count - 1);


        }

        [TestMethod]
        public void LoadBivariateSample () {

            DataTable table = CreateTable();

            BivariateSample sample = new BivariateSample();
            sample.Load(table.CreateDataReader(), 1, 3);
            Assert.IsTrue(sample.Count == table.Rows.Count - 1);

        }

        [TestMethod]
        public void LoadMultivariateSample () {

            DataTable table = CreateTable();

            MultivariateSample sample = new MultivariateSample(3);
            sample.Load(table.CreateDataReader(), 1, 2, 3);
            Assert.IsTrue(sample.Count == table.Rows.Count - 1);

        }

    }
}
