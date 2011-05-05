using System;
using System.Collections.Generic;
using System.Text;

namespace Meta.Numerics {


    public struct BigFloat {

        public BigFloat (byte[] mantissa, int exponent) {
            this.mantissa = mantissa;
            this.exponent = exponent;
            this.sign = true;
        }

        private byte[] mantissa;
        private int exponent;
        private bool sign;

        public static BigFloat operator - (BigFloat a) {
            BigFloat ma = new BigFloat();
            ma.mantissa = new byte[a.mantissa.Length];
            Array.Copy(a.mantissa, ma.mantissa, a.mantissa.Length);
            ma.exponent = a.exponent;
            ma.sign = !a.sign;
            return (ma);
        }

        public static BigFloat operator + (BigFloat a, BigFloat b) {

            // compute the shift
            int s = a.exponent - b.exponent;

            // if the shift is negative, do the operation in the opposite way
            if (s < 0) {
                Global.Swap<BigFloat>(ref a, ref b);
                s = -s;
            }

            BigFloat ab = new BigFloat();
            ab.mantissa = new byte[Math.Min(a.mantissa.Length, b.mantissa.Length - s)+1];
            ab.exponent = a.exponent;

            int c = 0;
            for (int k = ab.mantissa.Length - 1; k >= 0; k--) {
                int q = c + (a.mantissa[k] + b.mantissa[k - s]);
                if (q < 8) {
                    ab.mantissa[k+1] = (byte) q;
                    c = 0;
                } else {
                    ab.mantissa[k+1] = (byte)(q - 8);
                    c = 1;
                }
            }
            if (c > 0) {
                ab.mantissa[0] = 1;
                ab.exponent++;
            }

            return (ab);

        }

        public override string ToString () {

            StringBuilder text = new StringBuilder();
            for (int k = 0; k < mantissa.Length; k++) {
                text.AppendFormat("{0} ", mantissa[k]);
            }
            text.AppendFormat("X 256^{0}", exponent);

            return (text.ToString());

        }

    }


}
