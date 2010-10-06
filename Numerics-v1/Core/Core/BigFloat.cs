using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics.Core {

#if FUTURE

    public struct BigFloat {

        private byte[] mantissa;
        private int exponent;
        private bool sign;

        private static byte[] Add (IList<byte> n1, IList<byte> n2) {
            if (n1.Count != n2.Count) throw new InvalidOperationException();
            byte[] s = new byte[n1.Count + 1];
            for (int k = n1.Count - 1; k >= 0; k--) {
                int t = ((int) n1[k]) + ((int) n2[k]);

            }
            return (s);
        }

    }

#endif

}
