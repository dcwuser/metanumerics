using System;
using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    // Spearman's rho is just Pearson's r evaluated on the ranks of the values, i.e. the numbers 1 through n.
    // As such, it is completely invariant to any monotonic remapping of the values.
    //   \rho = \frac{\sum_i (x_i - \bar{x})(y_i - \bar{y}}{\sqrt{\sum_i (x_i - \bar{x})^2 \cdot \sum_i (y_i - \bar{y})^2}}
    // There is a straightforward relation between \rho as defined here and
    //   S = \sum_i x_i y_i
    //   D = \sum_i (x_i - y_i)^2
    // It is
    //   \rho = 1 - \frac{6 D}{n (n+1) (n-1)}
    //   D = n (n+1) (2n + 1) / 3 - 2 S
    // These can be derived by plugging the expressions for the sum of the first n integers and the squares of the first
    // n integers into the definitions of \rho and D above. The only remaining terms invovle S.

    // The limits of \rho, D, and S are easy to derive
    //   max   \rho = 1    D = 0                    S = n (n+1) (2n + 1) / 6
    //   min   \rho = -1   D = n (n+1) (n-1) / 3    S = n (n+1) (n+2) / 6
    //   mid   \rho = 0    D = n (n+1) (n-1) / 6    S = n (n+1)^2 / 4
    // Note the midpoint of S is not an integer for even n.

    // Under null hypothesis of independence, all relative orderings are equally likely. To compute null distribution, enumerate
    // all orderings, i.e. all n! permutations of n integers. E.g. n=3 has 3!=6 permutations
    //   1 2 3   S = 1*1 + 2*2 + 3*3 = 14
    //   1 3 2   S = 1*1 + 2*3 + 3*2 = 13
    //   2 1 3   S = 1*2 + 2*1 + 3*3 = 13
    //   2 3 1   S = 1*2 + 2*3 + 3*1 = 11
    //   3 1 2   S = 1*3 + 2*1 + 3*2 = 11
    //   3 2 1   S = 1*3 + 2*2 + 3*1 = 10
    // Thus the distribution of S is 10->1, 11->2, 13->2, 14->1

    // It is also true that the distribution of S is given by the coefficients of the polynomial perm(A) where A_{ij} = x^{ij}. E.g.
    //        | x    x^2  x^3 |
    //   perm | x^2  x^4  x^6 | = x^10 + 2 x^11 + 2 x^13 + x^14
    //        | x^3  x^6  x^9 |
    // For large n, it's fastest to evaluate the permanent via Ryser's formula.
    //   mask  sign  contribution
    //   001   +     x * x^2 * x^3 = x^6
    //   010   +     x^2 * x^4 * x^6 = x^12
    //   011   -     (x + x^2) * (x^2 + x^4) * (x^3 + x^6) = x^6 + x^7 + x^8 + 2 x^9 + x^10 + x^11 + x^12
    //   100   +     x^3 * x^6 * x^9 = x^18
    //   101   -     (x + x^3) * (x^2 + x^6) * (x^3 + x^9) = x^6 + x^8 + x^10 + 2 x^12 + x^14 + x^16 + x^18
    //   110   -     (x^2 + x^3) * (x^4 + x^6) * (x^6 + x^9) = x^12 + x^13 + x^14 + 2 x^15 + x^16 + x^17 + x^18
    //   111   +     (x + x^2 + x^3) * (x^2 + x^4 + x^6) * (x^3 + x^6 + x^9) =
    //                 x^6 + x^7 + 2 x^8 + 2 x^9 + 3 x^10 + 3 x^11 + 3 x^12 + 3 x^13 + 3 x^14 + 2 x^15 + 2 x^16 + x^17 + x^18
    // Adding all gives x^10 + 2 x^11 + 0 x^12 + 2 x^13 + x^14, same as above. (Yeah, for n=3 this is more work than enumeration, but
    // for higher n it's less. Enumeration work grows like n! while work for Ryser's permanent grows like 2^n.)
    // Note some things that help us:
    //    1. We know from the start that coefficients of powers below S_min will vanish, so we needn't compute them if we can avoid it.
    //    2. We know from the start that coefficients of powers above S_mid will mirror those below, so we needn't compute them either.
    //    3. Contributions of left-shifted bitmasks just add n(n+1)/2 powers for each shift. In our example, compare 001 to 010 to 100 and 011 to 110.
    //       So after we compute the contribution of each odd bitmask we can get the contributions of higher members trivially.
    // This method was summarized by Luke Gustafson at http://www.luke-g.com/math/spearman/index.html, referencing
    //   van de Wiel and Bucchianico, "Fast Computation of the Exact Null Distribution of Spearman's rho and Page's L Statistic
    //   for Samples with and without Ties", 1998

    internal sealed class SpearmanExactDistribution : DiscreteDistribution {

        public SpearmanExactDistribution (int n) {
            if (n < 2) throw new ArgumentOutOfRangeException("n");
            this.n = n;

            // n = 20 is maximum n for which n! fits in long
            this.totalCounts = 1;
            for (long i = 2; i <= n; i++) this.totalCounts *= i;

            ComputeCounts();
        }

        private readonly int n;

        int sMin, sMax, sMid;
        long[] counts;
        long totalCounts;

        private long GetCount (int s) {
            if (s <= sMid) {
                return (counts[s - sMin]);
            } else {
                return (counts[sMax - s]);
            }
        }

        public long GetLeftExclusiveCountSum (int s) {
            if (s <= sMin) {
                return (0);
            } else if (s <= sMid) {
                long sum = 0;
                for (int i = 0; i < (s - sMin); i++) sum += counts[i];
                return (sum);
            } else if (s <= sMax) {
                long sum = totalCounts;
                for (int i = 0; i <= (sMax - s); i++) sum -= counts[i];
                return (sum);
            } else {
                throw new InvalidOperationException();
            }
        }

        // We expect to get P = 0 at left boundary and P = 1 at right boundary
        // If we use exclusive summation, we get P = 0 at left boundary but P < 1 at right boundary
        // If we use inclusive summation, we get P = 1 at right bountary but P > 1 at left boundary

        public override double LeftExclusiveProbability (int k) {
            if (k < sMin) {
                return (0.0);
            } else if (k > sMax) {
                return (1.0);
            } else {
                return ((double) GetLeftExclusiveCountSum(k) / totalCounts);
            }
        }

        public override double RightExclusiveProbability (int k) {
            if (k < sMin) {
                return (1.0);
            } else if (k > sMax) {
                return (0.0);
            } else {
                return ((double) GetLeftExclusiveCountSum(sMax - (k - sMin)) / totalCounts);
            }
        }

        public override double ProbabilityMass (int k) {
            if ((k < sMin) || (k > sMax)) {
                return(0.0);
            } else {
                return((double) GetCount(k) / totalCounts);
            }
        }

        public override DiscreteInterval Support {
            get {
                return (new DiscreteInterval(sMin, sMax));
            }
        }

        public override double Mean {
            get {
                return (n * (n + 1) * (n + 1) / 4.0);
            }
        }

        public override double Variance {
            get {
                // V(S) = n^2 (n + 1)^2 (n - 1) / 144
                // Since V(S) / DS^2 = V(\rho) / 2^2 and DS = n (n + 1) (n - 1) / 6, V(\rho) = 1 / (n - 1)
                return ((double) n * n * (n + 1) * (n + 1) * (n - 1) / 144.0);
            }
        }


        public override double Skewness {
            get {
                return (0.0);
            }
        }

        private void ComputeCounts () {

            // our results must fit into a long, which takes us up to about n~20
            if (AdvancedIntegerMath.LogFactorial(n) > Math.Log(Int64.MaxValue)) throw new InvalidOperationException();

            // pre-compute some values we will use
            int shiftIncrement = n * (n + 1) / 2;
            sMin = n * (n + 1) * (n + 2) / 6;
            sMax = n * (n + 1) * (2 * n + 1) / 6;
            //Console.WriteLine("sMin = {0} sMax = {1}", sMin, sMax);

            // to improve performance, we calculate counts only up to the middle bin and use symmetry to obtain the higher ones
            // the midpoint is of course (sMin+sMax)/2, but for even n this is a half-integer
            // since we want to include the middle bin but integer arithmetic rounds down, we
            // add one before halving; this gives the right maximum bin to compute for both even and odd n
            sMid = (sMin + sMax + 1) / 2;

            bool baseSign = (n % 2 == 0);

            // array to accumulate permanent polynomial coefficients
            // only the entries between sMin and sMax matter, as per (1) and (2) above
            long[] totalP = new long[sMid + 1];

            // as per Ryser's formula, loop over all cominations of rows
            // use b as a length 2^n bitmask indicating which columns to include in a combination
            ulong bmax = (((ulong) 1) << n);

            // skip b = 0 since that means include no columns, resulting in a zero contribution

            // skip even b values as we will compute their contributions by using (3) above

            // do the b=1 case specially
            // when only one column enters, the product over all rows is a simply the monomial x^{(c+1) n n(+1) / 2}
            // knowing this, we can avoid all the multiplication to get us there
            // all the one column cases are produced by left-shifting the b=1 case
            // only those for which the power lies in [sMin,sMax] need be recorded
            int p = shiftIncrement;
            while (p < sMin) p += shiftIncrement;
            while (p <= sMid) {
                totalP[p] += baseSign ? -1 : 1;
                //totalP[p] += 1;
                p += shiftIncrement;
            }

            // okay, now we are ready for the remaining combinations of columns
            for (ulong b = 3; b < bmax; b += 2) {
                //Console.WriteLine("k={0}", k);

                // compute the row sum polynomial for the first row
                // while we are at it, we compute the sign of the combination
                // and note the degree of the polynomial as well
                long[] productP = new long[n + 1];
                //bool sign = true;
                bool sign = baseSign;
                int d = 0;
                for (int c = 0; c < n; c++) {
                    ulong bp = ((ulong) 1) << c;
                    if ((b & bp) != 0) {
                        /* productP[c + 1] = 1; */
                        productP[c] = 1;
                        /* d = c + 1; */
                        d = c;
                        sign = !sign;
                    }
                }

                // row sum polynomials for higher rows are obtained from those from
                // the first row by increasing the stride. just continue to multiply
                // the product polynomial by each row sum polynomial to get the total
                // contribution of the set
                for (int r = 1; r < n; r++) {
                    int da = productP.Length - 1;
                    productP = LongPolynomialMultiply(productP, da, b, d, r + 1, sMid);
                }

                //LongPolynomialMath.Write(productP, 0, productP.Length-1);

                // add the contribution that particular k,
                // then the contributions from left-shifting it,
                // which as per (3) are obtained by shifting the powers by n(n+1)/2
                // for each left-shift

                /* int shift = 0; */
                int shift = shiftIncrement;
                ulong k2 = b;
                while (k2 < bmax) {
                    LongPolynomialAdd(totalP, productP, shift, sMin, sMid, sign);
                    k2 = k2 << 1;
                    shift += shiftIncrement;
                }

            }

            // include the overall (-1)^n
            //if (n % 2 == 0) {
            //    LongPolynomialNegate(totalP, sMin, sMid);
            //}

            counts = new long[sMid - sMin + 1];
            Array.Copy(totalP, sMin, counts, 0, counts.Length);
            //for (int i = 0; i < counts.Length; i++) {
            //    counts[i] = totalP[sMin + i];
            //}
            //for (int i = sMin; i <= sMid; i++) {
            //    Console.WriteLine("{0} {1}", i, totalP[i]);
            //}

        }

        // a is a polynomial a[0] + a[1] x + a[2] x^2 + \cdots + a[d_a] x^{d_a} of degree d_a, so a.Length >= da + 1
        // b is a Z_2 polynomial of degree d_b * s, encoded as a binary number (b_db \cdots b_1 b_0)
        // representing b_0 + b_1 x^s + b_2 x^{2s} + \cdots + b_{db} x^{d_b s}
        // we require that b_0 = 1, i.e. b is an odd binary number

        // the product polynomial will have degree d_a + d_b s, but we will limit our computation up to terms of degree max

        private static long[] LongPolynomialMultiply (long[] a, int da, ulong b, int db, int s, int max) {
            int pbmax = db * s;
            int dab = Math.Min(da + pbmax, max);
            long[] ab = new long[dab + 1];
            /* for (int pb = s; pb <= pbmax; pb += s) { */
            for (int pb = 0; pb <= pbmax; pb += s) {
                if (b == 0) break;
                if ((b & 1) != 0) {
                    int pamax = Math.Min(dab - pb, da);
                    for (int pa = 0; pa <= pamax; pa++) {
                        ab[pa + pb] += a[pa];
                    }
                }
                b = b >> 1;
            }
            return (ab);
        }

        private static void LongPolynomialAdd (long[] c, long[] d, int shift, int min, int max, bool invert) {
            int diMin = Math.Max(0, min - shift);
            int diMax = Math.Min(d.Length - 1, max - shift);
            if (!invert) {
                for (int di = diMin; di <= diMax; di++) {
                    c[di + shift] -= d[di];
                }
            } else {
                for (int di = diMin; di <= diMax; di++) {
                    c[di + shift] += d[di];
                }
            }
        }


    }

}
