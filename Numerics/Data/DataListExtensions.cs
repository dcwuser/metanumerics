using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Data
{
    /*
    public static class Univariate
    {

        public static double Mean (this IReadOnlyCollection<double> sample)
        {
            double sum = 0.0;
            foreach (double value in sample)
            {
                sum += value;
            }
            return (sum / sample.Count);
        }

        public static double Variance (this IReadOnlyCollection<double> sample)
        {
            throw new NotImplementedException();
        }

        public static double PopulationMean (this IReadOnlyCollection<double> sample)
        {
            throw new NotImplementedException();
        }

    }
    */

    public static class DataListExtensions
    {
        // frame.Column("foo").Distinct()
        // frame.GroupBy<T>("a", Func<IEnumerable<DataFrameRow>, T>)
        // frame.GroupBy<T>("a", "b", Func<IEnumerable<DataFrameRow>, T>)
        // frame.OrderBy("a", Comparer)

        public static HashSet<T> Distinct<T> (this DataList<T> list)
        {
            HashSet<T> values = new HashSet<T>();
            foreach (T value in list)
            {
                values.Add(value);
            }
            return (values);
        }

        public static DataList<U> Transform<T, U>(this DataList<T> list, Func<T, U> transformation)
        {
            DataList<U> result = new DataList<U>(list.Name);
            foreach(T t in list)
            {
                U u = transformation(t);
                result.Add(u);
            }
            return (result);
        }
        
    }
}
