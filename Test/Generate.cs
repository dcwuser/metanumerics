using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test {
    public abstract class Generate<T> : IEnumerable<T> {

        public abstract IEnumerator<T> GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() {
            return GetEnumerator();
        }

        public T Minimum { get; set; }

        public T Maximum { get; set; }

        public Random Rng { get; set; }

        public Generate<T> From (T min) {
            Minimum = min;
            return this;
        }

        public Generate<T> To(T max) {
            Maximum = max;
            return this;
        }

        public Generate<T> Using (Random rng) {
            Rng = rng;
            return this;
        }

        public static Generate<double> Doubles() {
            return new DoubleGenerator();
        }

    }

    internal class DoubleGenerator : Generate<Double> {

        public override IEnumerator<double> GetEnumerator() {
            while (true) {
                double u = Rng.NextDouble();
                yield return (1.0 - u) * Minimum + u * Maximum;
            }
        }
    }
}
