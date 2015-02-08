using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace Meta.Numerics.Functions {

#if FUTURE

    public sealed class Partition {

        internal Partition (int[] list) {
            Debug.Assert(list != null);
            this.list = list;
        }

        private int[] list;

        public ReadOnlyCollection<int> GetListRepresentation () {
            return (new ReadOnlyCollection<int>(list));
        }

        public ReadOnlyCollection<Element> GetElementRepresentation () {

            List<Element> elements = new List<Element>();
            int v = 0;
            int m = 0;
            for (int i = 0; i < list.Length; i++) {
                if (v == list[i]) {
                    m++;
                } else {
                    elements.Add(new Element(v, m));
                    v = list[i];
                    m = 1;
                }
            }
            return (new ReadOnlyCollection<Element>(elements));
        }


    }

#endif

}
