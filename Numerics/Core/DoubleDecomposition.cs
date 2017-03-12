using System;


namespace Meta.Numerics {

    /// <summary>
    /// Represents the decomposition of a double into its sign, mantissa, and exponent.
    /// </summary>
    public struct DoubleDecomposition {

        public DoubleDecomposition (double value) : this(BitConverter.ToInt64(BitConverter.GetBytes(value), 0)) {

        }

        internal DoubleDecomposition (long bits) {

            this.bits = bits;

            exponent = (int) ((bits >> 52) & 0x7ffL);
            mantissa = bits & 0xfffffffffffffL;

            //if (mantissa == 0L) return;

            exponent -= 1075;

            mantissa = mantissa | (1L << 52);

            while ((mantissa & 1) == 0) {
                mantissa >>= 1;
                exponent++;
            }
        }

        private long mantissa;

        private int exponent;

        private long bits;

        public long Mantissa {
            get {
                return (mantissa);
            }
        }

        public int Exponent {
            get {
                return (exponent);
            }
        }
        
        public bool IsNegative {
            get {
                return (bits < 0);
            }
        }

        public double Value {
            get {
                return (BitConverter.ToDouble(BitConverter.GetBytes(bits), 0));
            }
        }

        public DoubleDecomposition Next() {
            long nextBits = bits + 1L;
            return (new DoubleDecomposition(nextBits));
        }

        public string ToString () {
            return (String.Format("{0}{1} X 2^({2})", this.IsNegative ? "-" : String.Empty, this.Mantissa, this.Exponent));
        }

    }
}
