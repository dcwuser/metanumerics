using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Meta.Numerics.Extended {

    // There are basically two ways to structure DoubleInfo: 1. Just store the bit pattern. This keeps the
    // storage cost to 64 bits, the same as the orginal double. You will need to expend some computational effort
    // to get at the data you want when each propery is invoked, but most of that computational effort is very
    // fast bitmask operations. 2. Store the extracted data. This multiplies the storage size by about 4. At
    // very least you need to store a sign bit, the exponent, and the mantissa. To easily distinguish subnormals
    // and allow fast Next, Previous, and Value implementations you probably also want to store the original bit
    // pattern.


    // double format is 64b: (1b sign)(11b exponent)(52b fraction)
    // for normals, value is (-1)^s (1.ffffff) X 2^(e-1023) = (-1)^s (1ffffff) X 2^(e-(1023 + 52))
    // for subnormals, e is all zeros and value is (-1)^s (0.ffffff) X 2^(1-1023)
    // for infinities and NaNs, 

    /// <summary>
    /// Contains information on the stored represenation of a double value.
    /// </summary>
    public readonly struct DoubleInfo {

        /// <summary>
        /// Initializes a new double info object for the given double value.
        /// </summary>
        /// <param name="value">The value to analyze.</param>
        public DoubleInfo(double value) : this(unchecked((ulong) BitConverter.DoubleToInt64Bits(value))) {

        }

        internal DoubleInfo (ulong storage) {
            this.storage = storage;
            GetExponentAndMantissa(storage, out exponent, out mantissa);
        }

        private readonly ulong storage;
        private readonly int exponent;
        private readonly long mantissa;

        /// <summary>
        /// Creates a new double info object for the given double value.
        /// </summary>
        /// <param name="value">The value to analyze.</param>
        /// <returns>The corresponding <see cref="DoubleInfo"/> data.</returns>
        /// <remarks>
        /// <para>This is equivilent to the constructor <see cref="DoubleInfo.DoubleInfo(double)"/>.</para>
        /// </remarks>
        public static DoubleInfo For (double value) {
            return new DoubleInfo(value);
        }

        /// <summary>
        /// Gets a value indicating whether the number is negative.
        /// </summary>
        /// <remarks>
        /// <para>Note that -0.0 is negative for the purposes of this test,
        /// as are NaNs with the negative bit set.</para>
        /// </remarks>
        public bool IsNegative {
            get {
                return (storage & 0x8000000000000000UL) != 0;
            }
        }

        /// <summary>
        /// Gets a value indicating whether the floating point value is a finite number.
        /// </summary>
        /// <remarks>
        /// <para>Infinities and NaNs are not finite.</para>
        /// </remarks>
        public bool IsFinite {
            get {
                return (storage & 0x7FF0000000000000UL) < 0x7FF0000000000000UL;
            }
        }

        /// <summary>
        /// Gets a value indicating whether the floating point value is zero.
        /// </summary>
        /// <remarks>
        /// <para>Both +0.0 and -0.0 are zero for the purposes of this test.</para>
        /// </remarks>
        public bool IsZero {
            get {
                // All exponent and mantissa bits are 0.
                return (storage & 0x7FFFFFFFFFFFFFFFUL) == 0UL;
            }
        }

        /// <summary>
        /// Gets a value indicating whether the value is infinite.
        /// </summary>
        public bool IsInfinite {
            get {
                // All exponent bits are 1 and all mantissa bits are zero.
                return (storage & 0x7FFFFFFFFFFFFFFFUL) == 0x7FF0000000000000UL;
            }
        }

        /// <summary>
        /// Gets a value indicating whether the value is not-a-number.
        /// </summary>
        public bool IsNaN {
            get {
                // All exponent bits are one and there are non-zero mantissa bits
                return (storage & 0x7FFFFFFFFFFFFFFFUL) > 0x7FF0000000000000UL;
            }
        }

        /// <summary>
        /// Gets a value indicating whether the floating point value is sub-normal.
        /// </summary>
        public bool IsSubnormal {
            get {
                // All exponent bits are zero and there are non-zero mantissa bits.
                return ((storage & 0x7FF0000000000000UL) == 0) && ((storage & 0x000FFFFFFFFFFFFFUL) != 0);
            }
        }

        /// <summary>
        /// Gets the double value.
        /// </summary>
        public double Value {
            get {
                unchecked {
                    return BitConverter.Int64BitsToDouble((long) storage);
                }
            }
        }

        private static void GetExponentAndMantissa (ulong storage, out int exponent, out long mantissa) {
            exponent = (int) ((storage & 0x7FF0000000000000) >> 52);
            mantissa = (long) (storage & 0x000FFFFFFFFFFFFF);

            if (exponent == 0) {
                // subnormals and zeros
                if (mantissa == 0L) return;
                exponent++;
            } else {
                // normals, infinities, nans
                mantissa = (1L << 52) | mantissa;
            }

            exponent -= (1023 + 52);
            while ((mantissa & 1L) == 0) {
                exponent++;
                mantissa = mantissa >> 1;
            }

        }

        /// <summary>
        /// Gets the base-2 exponent of the floating point value.
        /// </summary>
        public int Exponent {
            get {
                return exponent;
            }
        }

        /// <summary>
        /// Gets the mantissa of the floating point value.
        /// </summary>
        public long Mantissa {
            get {
                return mantissa;
            }
        }

        /// <summary>
        /// Gets the next higher floating point value.
        /// </summary>
        public DoubleInfo Next {
            get {
                // For a generic positive number, you get the next number by adding 1.
                // The placement of the sign bit complicates things just a bit.

                if (this.IsNegative) {
                    if (storage == 0x8000000000000000UL) {
                        return new DoubleInfo(0UL);
                    } else {
                        return new DoubleInfo(storage - 1UL);
                    }
                } else {
                    if (storage == 0x8000000000000000UL) {
                        return new DoubleInfo(0UL);
                    } else {
                        return new DoubleInfo(storage + 1UL);
                    }
                }
            }
        }

        /// <summary>
        /// Gets the next lower floating point value.
        /// </summary>
        public DoubleInfo Previous {
            get {
                if (this.IsNegative) {
                    if (storage == 0xFFFFFFFFFFFFFFFFUL) {
                        return new DoubleInfo(0x7FFFFFFFFFFFFFFFUL);
                    } else {
                        return new DoubleInfo(storage + 1UL);
                    }
                } else {
                    if (storage == 0UL) {
                        return new DoubleInfo(0x8000000000000000UL);
                    } else {
                        return new DoubleInfo(storage - 1UL);
                    }
                }
            }
        }

        /// <summary>
        /// Gets the internal representation of the floating point value.
        /// </summary>
        [CLSCompliant(false)]
        public ulong Bits {
            get {
                return storage;
            }
        }

        /// <inheritdoc/>
        public override string ToString () {
            return ToString(CultureInfo.CurrentCulture);
        } 

        private string ToString (IFormatProvider provider) {
            Debug.Assert(provider != null);
            return String.Format(provider, "{0}{1} X 2^({2})", this.IsNegative ? "-" : String.Empty, this.Mantissa, this.Exponent);
        }
    }

}
