using System;
using System.Text;
using System.Collections.Generic;
using System.Data;
using Microsoft.VisualStudio.TestTools.UnitTesting;

using Meta.Numerics.Statistics;

namespace Test {
    [TestClass]
    public class DbTest {

        [TestMethod]
        public void CastTest () {
            int i = 32;
            object o = (object) i;
            double d = Convert.ToDouble(o);
        }

        [TestMethod]
        public void SampleLoad () {

            DataTable table = new DataTable();
            table.Columns.Add("Name", typeof(string));
            table.Columns.Add("Age", typeof(int));
            table.Columns.Add("Day", typeof(DayOfWeek));
            table.Columns.Add("Male", typeof(bool));
            table.Rows.Add("A", 1, DayOfWeek.Monday, false);
            table.Rows.Add("B", 2, DayOfWeek.Tuesday, false);
            table.Rows.Add("C", 3, DayOfWeek.Wednesday, true);
            table.Rows.Add(null, null, null, null);

            Sample sample = new Sample();
            sample.Load(table.CreateDataReader(), 1);
            Console.Write("{0} {1} {2}", sample.Count, sample.Mean, sample.StandardDeviation);


        }

    }
}
