using System;
using System.Diagnostics;

namespace Meta.Numerics.Spin {

    /// <summary>
    /// Contains methods for computing functions of spin and spin states.
    /// </summary>
    public static class SpinMath {

        /// <summary>
        /// Enumerates all the spins that can be obtained by combining two spins.
        /// </summary>
        /// <param name="j1">The first spin.</param>
        /// <param name="j2">The second spin.</param>
        /// <returns>A list of spins which may be obtained.</returns>
        public static Spin[] Combine (Spin j1, Spin j2) {

            int tj_min = Math.Abs(j1.TwoJ - j2.TwoJ);
            int tj_max = j1.TwoJ + j2.TwoJ;

            Spin[] spins = new Spin[(tj_max - tj_min) / 2 + 1];
            for (int i = 0; i < spins.Length; i++) {
                spins[i] = new Spin((tj_min + 2 * i) / 2.0);
            }

            return (spins);
        }

        /// <summary>
        /// Enumerates all spin states that may be obtained by combining two spin states.
        /// </summary>
        /// <param name="s1">The first spin state.</param>
        /// <param name="s2">The second spin state.</param>
        /// <returns>A list of spin states which may  be obtained.</returns>
        public static SpinState[] Combine (SpinState s1, SpinState s2) {

            int tj_max = s1.TwoJ + s2.TwoJ;
            int tj_min = Math.Abs(s1.TwoJ - s2.TwoJ);
            int tm = s1.TwoM + s2.TwoM;
            if (Math.Abs(tm) > tj_min) tj_min = Math.Abs(tm);

            SpinState[] states = new SpinState[(tj_max - tj_min) / 2 + 1];
            for (int i = 0; i < states.Length; i++) {
                states[i] = new SpinState((tj_min + 2 * i) / 2.0, tm / 2.0);
            }

            return (states);

        }

        // ( J1 J2 M1 M2 | J M ) =
        // (-1)^{J1-J2+M} * Sqrt(2J+1) *
        // ( J1  J2   J )
        // ( M1  M1  -M )

        /// <summary>
        /// Computes a Clebsch-Gordon coefficient.
        /// </summary>
        /// <param name="s1">The first spin state.</param>
        /// <param name="s2">The second spin state.</param>
        /// <param name="s">The total spin state.</param>
        /// <returns>The Clebsch-Gordon coefficient measuring the contribution of the given first and
        /// second spin states to the given total spin state.</returns>
        public static double ClebschGodron (SpinState s1, SpinState s2, SpinState s) {

            double f = Math.Sqrt(s.TwoJ + 1);
            if (s1.JPlusM % 2 != 0) f = -f;
            if (s2.JMinusM % 2 != 0) f = -f;
            return (f * ThreeJ(s1, s2, s.Invert()));

        }

        /// <summary>
        /// Computes a 3j symbol.
        /// </summary>
        /// <param name="s1">The first column spin state.</param>
        /// <param name="s2">The second column spin state.</param>
        /// <param name="s3">The third column spin state.</param>
        /// <returns>The 3j symbol.</returns>
        public static double ThreeJ (SpinState s1, SpinState s2, SpinState s3) {

            // require that Js satisfy triangle inequality
            if (s3.TwoJ < Math.Abs(s1.TwoJ - s2.TwoJ)) return (0.0);
            if (s3.TwoJ > (s1.TwoJ + s2.TwoJ)) return (0.0);

            // require that the parity of the spins agree
            // for example, spin-1 and spin-2 cannot combine to spin 3/2, even though 1 < 3/2 < 3
            if ((s1.TwoJ + s2.TwoJ + s3.TwoJ) % 2 != 0) return (0.0);

            // require that Ms add to zero
            if ((s1.TwoM + s2.TwoM + s3.TwoM) != 0) return (0.0);

            // determine permutation factor for row switching or inverting
            int P;
            if (((s1.TwoJ + s2.TwoJ + s3.TwoJ) / 2) % 2 == 0) {
                P = 1;
            } else {
                P = -1;
            }

            // shuffle so J1 <= J2 <= J3
            // keeping track of the sign induced by the permutation
            int s = 1;
            if (s1.TwoJ > s3.TwoJ) {
                Exchange(ref s1, ref s3);
                s = s * P;
            }
            if (s1.TwoJ > s2.TwoJ) {
                Exchange(ref s1, ref s2);
                s = s * P;
            }
            if (s2.TwoJ > s3.TwoJ) {
                Exchange(ref s2, ref s3);
                s = s * P;
            }

            // now we can test for special cases,
            // knowing where the Js will fall in each case

            if (s1.TwoJ == 0) {
                // J1 = 0
                return (s * ThreeJ_ZeroJ(s2));
            } else if (s1.TwoJ == 1) {
                // J1 = 1/2
                return (s * ThreeJ_HalfJ(s1, s2));
            } else if (s1.TwoJ == 2) {
                // J1 = 1
                return (s * ThreeJ_OneJ(s1, s2, s3));
            } else {
                // arbitrary Js
                if ((s1.TwoM == 0) && (s2.TwoM == 0) && (s3.TwoM == 0)) {
                    // special case of all M's zero
                    return (s * ThreeJ_ZeroM_Int(s1, s2, s3));
                } else {
                    // general case
                    return (s * ThreeJ_ShultenGordon_RecurseJ(s1, s2, s3));
                }
            }

            // other possible special cases include:
            // J1 = 3/2, J1 = 2
            // Ms = 1/2 -1/2 0, 1/2 1/2 1, 1 -1 0 and permutations
            // maximum coupling J3 = J1 + J2, minimum coupling J3 = |J1 - J2|

        }

        // exchange two references

        private static void Exchange (ref SpinState s1, ref SpinState s2) {
            SpinState t = s1;
            s1 = s2;
            s2 = t;
        }

        private static void Swap<T> (ref T a, ref T b) where T : struct {
            T t = a;
            a = b;
            b = t;
        }

        // ( 0  J   J )
        // ( 0  M  -M )
        private static double ThreeJ_ZeroJ (SpinState s) {
            double r = 1.0 / Math.Sqrt(s.TwoJ + 1);
            if (s.JMinusM % 2 != 0) r = -r;
            return (r);
        }

        // ( 1/2  J2   J+1/2 )
        // ( M1   M2  -M1-M2 )
        private static double ThreeJ_HalfJ (SpinState s1, SpinState s2) {

            Debug.Assert(s1.TwoJ == 1);

            if (s1.TwoM < 0) {
                // m1 = -1/2
                // invert to turn into m1 = +1/2 case
                double r = ThreeJ_HalfJ(s1.Invert(), s2.Invert());
                if (s2.TwoJ % 2 == 0) r = -r;
                return (r);
            } else {
                // m1 = +1/2
                double r = -Math.Sqrt((s2.JPlusM + 1.0) / (s2.TwoJ + 1.0) / (s2.TwoJ + 2.0));
                if (s2.JMinusM % 2 != 0) r = -r;
                return (r);
            }

        }

        // ( 1   J2  J3 )
        // ( M1  M2  M3 )
        private static double ThreeJ_OneJ (SpinState s1, SpinState s2, SpinState s3) {

            Debug.Assert(s1.TwoJ == 2);

            if (s1.TwoM < 0) {
                // m1 = -1
                double r = ThreeJ_OneJ(s1.Invert(), s2.Invert(), s3.Invert());
                if ((s1.TwoJ + s2.TwoJ + s3.TwoJ) / 2 % 2 != 0) r = -r;
                return (r);
            } else if (s1.TwoM == 0) {
                // m1 = 0
                if (s3.TwoJ == s2.TwoJ) {
                    // j3 = j2
                    double r = s2.TwoM / Math.Sqrt(s2.TwoJ * (s2.TwoJ + 1) * (s2.TwoJ + 2));
                    if ((s2.JMinusM % 2) != 0) r = -r;
                    return (r);
                } else {
                    // j3 = j2 + 1
                    double r = -Math.Sqrt(2.0 * (s2.JPlusM + 1) * (s2.JMinusM + 1) / (s2.TwoJ + 1) / (s2.TwoJ + 2) / (s2.TwoJ + 3));
                    if ((s2.JMinusM % 2) != 0) r = -r;
                    return (r);
                }
            } else {
                // m1 = 1
                if (s3.TwoJ == s2.TwoJ) {
                    // j3 = j2
                    double r = Math.Sqrt(2.0 * s2.JMinusM * (s2.JPlusM + 1) / s2.TwoJ / (s2.TwoJ + 1) / (s2.TwoJ + 2));
                    if ((s2.JMinusM % 2) != 0) r = -r;
                    return (r);
                } else {
                    // j3 = j2 + 1
                    double r = Math.Sqrt(1.0 * (s2.JPlusM + 1) * (s2.JPlusM + 2) / (s2.TwoJ + 1) / (s2.TwoJ + 2) / (s2.TwoJ + 3));
                    if ((s2.JMinusM % 2) != 0) r = -r;
                    return (r);
                }
            }

        }

#if FUTURE

        // ( 3/2  J2  J3 )
        // ( M1   M2  M3 )
        private static double ThreeJ_ThreeHalvesJ (SpinState s1, SpinState s2, SpinState s3) {

            Debug.Assert(s1.TwoJ == 3);

            if (s1.TwoM < 0) {
                // m1 is negative; invert the m's
                double r = ThreeJ_ThreeHalvesJ(s1.Invert(), s2.Invert(), s3.Invert());
                // fix sign!
                return (r);
            } else if (s1.TwoM == 1) {
                // m1 = 1/2
                if (s3.TwoJ == (s2.TwoJ + 1)) {
                    // j3 = j2 + 1/2
                    double r = (s2.JMinusM-s2.TwoM) * Math.Sqrt(1.0 * (s2.JPlusM + 1) / s2.TwoJ / (s2.TwoJ + 1) / (s2.TwoJ + 2) / (s2.TwoJ + 3));
                    if (s2.JMinusM % 2 != 0) r = -r;
                    return (r);
                } else {
                    // j3 = j2 + 3/2
                    double r = Math.Sqrt(3.0 * (s2.JPlusM + 1) * (s2.JPlusM + 2) * (s2.JMinusM + 1) / (s2.TwoJ + 1) / (s2.TwoJ + 2) / (s2.TwoJ + 3) / (s2.TwoJ + 4));
                    if (s2.JMinusM % 2 != 0) r = -r;
                    return (r);
                }
            } else {
                // m1 = 3/2
                if (s3.TwoJ == (s2.TwoJ + 1)) {
                    // j3 = j2 + 1/2
                    double r = -Math.Sqrt(3.0 * (s2.JPlusM+1) * (s2.JPlusM+2) * (s2.JMinusM) / s2.TwoJ / (s2.TwoJ+1) / (s2.TwoJ+2) / (s2.TwoJ+3) );
                    if (s2.JMinusM % 2 != 0) r = -r;
                    return (r);
                } else {
                    // j3 = j2 + 3/2
                    double r = -Math.Sqrt(1.0 * (s2.JPlusM + 1) * (s2.JPlusM + 2) * (s2.JPlusM + 3) / (s2.TwoJ + 1) / (s2.TwoJ + 2) / (s2.TwoJ + 3) / (s2.TwoJ + 4));
                    if (s2.JMinusM % 2 != 0) r = -r;
                    return (r);
                }
            }

        }

        // ( J1  J2    J1+J2  )
        // ( M1  M2  -(M1+M2) )
        private static double ThreeJ_MaxJ (SpinState j1, SpinState j2) {
            double r = AdvancedIntegerMath.LogFactorial(j1.TwoJ) +
                AdvancedIntegerMath.LogFactorial(j2.TwoJ) +
                AdvancedIntegerMath.LogFactorial(j1.JPlusM + j2.JPlusM) +
                AdvancedIntegerMath.LogFactorial(j1.JMinusM + j2.JMinusM) -
                AdvancedIntegerMath.LogFactorial(j1.TwoJ + j2.TwoJ + 1) -
                AdvancedIntegerMath.LogFactorial(j1.JPlusM) -
                AdvancedIntegerMath.LogFactorial(j1.JMinusM) -
                AdvancedIntegerMath.LogFactorial(j2.JPlusM) -
                AdvancedIntegerMath.LogFactorial(j2.JMinusM);
            r = Math.Exp(r / 2.0);
            return (r);
        }

        private static double ThreeJ_Any (SpinState j1, SpinState j2, SpinState j3) {

            // compute prefactor

            double f1 = 0.0
                + AdvancedIntegerMath.LogFactorial((j1.TwoJ + j2.TwoJ - j3.TwoJ) / 2)
                + AdvancedIntegerMath.LogFactorial((j1.TwoJ - j2.TwoJ + j3.TwoJ) / 2)
                + AdvancedIntegerMath.LogFactorial((-j1.TwoJ + j2.TwoJ + j3.TwoJ) / 2)
                - AdvancedIntegerMath.LogFactorial((j1.TwoJ + j2.TwoJ + j3.TwoJ) / 2 + 1);
            double f2 = 0.0
                + AdvancedIntegerMath.LogFactorial(j1.JPlusM)
                + AdvancedIntegerMath.LogFactorial(j1.JMinusM)
                + AdvancedIntegerMath.LogFactorial(j2.JPlusM)
                + AdvancedIntegerMath.LogFactorial(j2.JMinusM)
                + AdvancedIntegerMath.LogFactorial(j3.JPlusM)
                + AdvancedIntegerMath.LogFactorial(j3.JMinusM);
            double f = Math.Exp((f1 + f2) / 2.0);
            if ((j1.JPlusM - j2.JMinusM) % 2 != 0) f = -f;
            Console.WriteLine("f={0}", f);

            // determine maximum and minimum values of k in sum

            int kmin = 0;
            int k23 = j2.JPlusM - j3.JMinusM;
            if (kmin < k23) kmin = k23;
            int k13 = j1.JMinusM - j3.JPlusM;
            if (kmin < k13) kmin = k13;

            int k123 = (j1.TwoJ + j2.TwoJ - j3.TwoJ) / 2;
            int kmax = k123;
            int k1 = j1.JMinusM;
            if (k1 < kmax) kmax = k1;
            int k2 = j2.JPlusM;
            if (k2 < kmax) kmax = k2;


            Console.WriteLine("{0} <= k <= {1}", kmin, kmax);
            Debug.Assert(kmin <= kmax);

            // compute the sum

            double g = 0.0;
            for (int k = kmin; k <= kmax; k++) {
                double gt = AdvancedIntegerMath.LogFactorial(k)
                    + AdvancedIntegerMath.LogFactorial(k123 - k)
                    + AdvancedIntegerMath.LogFactorial(k1 - k)
                    + AdvancedIntegerMath.LogFactorial(k2 - k)
                    + AdvancedIntegerMath.LogFactorial(k - k13)
                    + AdvancedIntegerMath.LogFactorial(k - k23);
                gt = Math.Exp(-gt);
                if (k % 2 != 0) gt = -gt;
                g += gt;
            }
            Console.WriteLine("g={0}", g);

            return (f * g);
        }

#endif

        // ( J1  J2  J3 )
        // ( 0   0   0  )
        private static double ThreeJ_ZeroM_Int (SpinState s1, SpinState s2, SpinState s3) {

            Debug.Assert(s1.TwoM == 0);
            Debug.Assert(s2.TwoM == 0);
            Debug.Assert(s3.TwoM == 0);

            // this method depends on the ordering j1 <= j2 <= j3
            Debug.Assert(s1.TwoJ <= s2.TwoJ);
            Debug.Assert(s2.TwoJ <= s3.TwoJ);

            int j = (s1.TwoJ + s2.TwoJ + s3.TwoJ) / 2;

            if (j % 2 != 0) {
                return (0.0);
            } else {

                // use the method formula derived in 

                int A = j - s1.TwoJ;
                int B = j - s2.TwoJ;
                int C = j - s3.TwoJ;

                //Console.WriteLine("A={0} B={1} C={2}", A, B, C);

                // note that A >= B >= C, because of the ordering of js; therefore we eliminate factors involving C, which are largest
                // note also that A, B, and C are even, because j is even and we subtract even numbers from j to get them

                int hA = A / 2;
                int hB = B / 2;
                int hC = C / 2;

                // calculate first factor
                double f = 1.0 / (j + 1);
                int fp1 = B / 2;
                int fp2 = A / 2;
                int fq1 = A;
                int fq2 = A + B / 2;
                for (int i = 1; i <= hB; i++) {
                    //fp1++; fp2++; fq1++; fq2++;
                    //f = f * fp1 * fp2 * fp2 / fq1 / fq2 / i;
                    f = f * (hB + i) * (hA + i) * (hA + i) / (A + i) / (A + hB + i) / i;
                }
                f = Math.Sqrt(f);

                // calculate second factor
                double g = 1.0;
                //int gp1 = C / 2;
                //int gp2 = (A + B) / 2;
                //int gq1 = A + B;
                //int gq2 = A + B + C / 2;
                for (int i = 1; i <= hC; i++) {
                    g = g * (hC + i) * (hA + hB + i) * (hA + hB + i) / (A + B + i) / (A + B + hC + i) / i;
                }
                //Console.WriteLine("g={0}", g);
                g = Math.Sqrt(g);

                if ((j / 2) % 2 != 0) f = -f;

                return (f * g);
            }


        }

#if FUTURE

        // ( J1  J2  J3 )
        // ( 0   0   0  )
        private static double ThreeJ_ZeroM (SpinState s1, SpinState s2, SpinState s3) {

            Debug.Assert(s1.TwoM == 0);
            Debug.Assert(s2.TwoM == 0);
            Debug.Assert(s3.TwoM == 0);

            int j = (s1.TwoJ + s2.TwoJ + s3.TwoJ) / 2;

            if (j % 2 != 0) {
                // zero by symmetry
                return (0.0);
            } else {

                double f = -AdvancedIntegerMath.LogFactorial(j + 1)
                    + AdvancedIntegerMath.LogFactorial(j - s1.TwoJ)
                    + AdvancedIntegerMath.LogFactorial(j - s2.TwoJ)
                    + AdvancedIntegerMath.LogFactorial(j - s3.TwoJ);
                f = Math.Exp(f / 2.0);

                int k = j / 2;
                double g = AdvancedIntegerMath.LogFactorial(k)
                    - AdvancedIntegerMath.LogFactorial(k - s1.TwoJ / 2)
                    - AdvancedIntegerMath.LogFactorial(k - s2.TwoJ / 2)
                    - AdvancedIntegerMath.LogFactorial(k - s3.TwoJ / 2);
                g = Math.Exp(g);

                if (k % 2 != 0) {
                    return (-f * g);
                } else {
                    return (f * g);
                }

            }

        }

        private static double ThreeJ_ZeroHalfM (SpinState s1, SpinState s2, SpinState s3) {

            Debug.Assert(s1.TwoM == 0);
            Debug.Assert(Math.Abs(s2.TwoM) == 0);

            int j = (s1.TwoJ + s2.TwoJ + s3.TwoJ) / 2;

            // if m2 < 0 invert and re-call

            if (j % 2 == 0) {

                double f = - Math.Sqrt(1.0 * (s1.TwoJ - s2.TwoJ + s3.TwoJ + 2) * (s1.TwoJ + s2.TwoJ - s3.TwoJ) / (s2.TwoJ + 1) / (s3.TwoJ + 1)) / 2.0;
                double g = ThreeJ_ZeroM(new SpinState(s1.TwoJ / 2.0, 0.0), new SpinState((s2.TwoJ - 1) / 2.0, 0.0), new SpinState((s3.TwoJ + 1) / 2.0, 0.0));
                return (f * g);

            } else {

                double f = Math.Sqrt(1.0 * (j + 2) / (s2.TwoJ + 1) / (s3.TwoJ + 1)) / 2.0;
                throw new NotImplementedException();
            }

        }
#endif

        // method described by Shulten & Gordon in J. Math. Phys. 16 (1975) 1961
        // uses the recurrsion
        //  j3 A_(j3+1) ( j1 j2 j3+1 ) + B_j3 ( j1 j2 j3 ) + (j3+1) A_j3 ( j1 j2 j3-1 )
        //              ( m1 m2  m3  )        ( m1 m2 m3 )               ( m1 m2  m3  )
        // plus the normalization condition
        // \sum_j3 (2 j3 + 1) ( j1 j2 j3 )^2 = 1
        //                    ( m1 m2 m3 )

        // we recurse on j3 because, for j1 <= j2 <= j3, the typical number of recursion steps
        // will be 2 * j1, the lowest of the three possible spins to recurse on

        private static double ThreeJ_ShultenGordon_RecurseJ (SpinState s1, SpinState s2, SpinState s3) {

            // find the range of j
            int tj_min = Math.Abs(s1.TwoJ - s2.TwoJ);
            if (tj_min < Math.Abs(s3.TwoM)) tj_min = Math.Abs(s3.TwoM);
            int tj_max = s1.TwoJ + s2.TwoJ;

            // find the midpoint of the range
            int tj_mid = (tj_min + tj_max) / 2;
            if ((tj_mid % 2) != (tj_max % 2)) tj_mid++;

            // variables to store recurrence coefficients and (unscaled) 3j values
            double ap, a0, b0;
            double cm, c0, cp;
            int tj;

            // a variable to store the sought value
            double c = 0.0;

            // iterate in from the right
            double NR = 0.0;
            cp = 0.0;
            ap = 0.0;
            c0 = 0.5 / tj_max;
            if (((s1.TwoJ - s2.TwoJ - s3.TwoM) / 2) % 2 != 0) c0 = -c0; // sign convention
            tj = tj_max;
            while (true) {

                // keep track of the normalization constant
                NR += (tj + 1) * c0 * c0;

                // remember the (unscaled) value of the desired 3j symbol
                if (tj == s3.TwoJ) {
                    c = c0;
                }

                // check whether we have reached the center
                if (tj <= tj_mid) {
                    if (Math.Abs(c0) <= Math.Pow(2.0, -40.0) * Math.Abs(cp)) {
                        // if the middle value is zero, the left and right values cannot be matached
                        // (acutally, the numbers will be two random tiny numbers, and their ratio will
                        // be a random number of order one); in this case, go down one more step
                        tj_mid -= 2;
                    } else {
                        break;
                    }
                }


                // compute required coefficients
                b0 = ShultenGodron_B(s1.TwoJ, s2.TwoJ, tj, s1.TwoM, s2.TwoM, s3.TwoM);
                a0 = ShultenGordon_A(s1.TwoJ, s2.TwoJ, tj, s3.TwoM);

                // use recursion to compute the tj-2 symbol
                cm = (-b0 * c0 - (tj / 2.0) * ap * cp) / ((tj + 2) / 2.0) / a0;

                //Console.WriteLine("{0}: {1} {2} + {3} {4} + {5} {6}", tj, ap, cp, b0, c0, a0, cm);

                // prepare for the next cycle
                cp = c0;
                c0 = cm;
                ap = a0;
                tj -= 2;

            }

            // remember the (unscaled) middle value
            double c1 = c0;

            // iterate in from the left
            double NL = 0.0;
            cm = 0.0;
            a0 = 0.0;
            c0 = 0.5 / tj_max;
            tj = tj_min;

            // j_min = 0 iff j1 = j2 and m3 = 0
            // in this case for j3=0, b0 = 0 and all three coefficients vanish, so the recursion is indeterminate
            // in this case we do the first recursion step "by hand" using relationships

            if (tj == 0) {
                NL += c0 * c0;
                cm = c0;
                c0 = s1.TwoM / Math.Sqrt(s1.TwoJ * (s1.TwoJ + 2)) * cm;
                a0 = ShultenGordon_A(s1.TwoJ, s2.TwoJ, 2, 0);
                tj = 2;
            }

            while (true) {

                // check whether we have reached the center
                if (tj >= tj_mid) {
                    break;
                }

                // notice that we terminate before checking whether we are at the desired value,
                // or adding to the normalization constant; this is different from the previous loop
                // we do this so that the middle element does not double contribute to the normalization
                // and because if the middle element is desired, it has already been stored from the previous loop

                // remember the desired (unscaled) value
                if (tj == s3.TwoJ) {
                    c = c0;
                }

                // add current value to normalization
                NL += (tj + 1) * c0 * c0;

                // check for end
                if (tj >= tj_max) break;

                // compute required coeficients
                b0 = ShultenGodron_B(s1.TwoJ, s2.TwoJ, tj, s1.TwoM, s2.TwoM, s3.TwoM);
                ap = ShultenGordon_A(s1.TwoJ, s2.TwoJ, tj + 2, s3.TwoM);

                // use recursion relation to compute tj+2
                cp = (-b0 * c0 - ((tj + 2) / 2.0) * a0 * cm) / (tj / 2.0) / ap;

                //Console.WriteLine("{0}: {1} {2} + {3} {4} + {5} {6}", tj, ap, cp, b0, c0, a0, cm);

                // prepare for next iteration
                cm = c0;
                c0 = cp;
                a0 = ap;
                tj += 2;

            }

            // match in the middle
            double r = c0 / c1;
            NL = NL / (r * r);
            if (s3.TwoJ < tj_mid) c = c / r;

            // normalize
            double N = NL + NR;
            c = c / Math.Sqrt(N);

            return (c);

        }

        private static double ShultenGordon_A (int tj1, int tj2, int tj3, int tm3) {
            double f1 = (tj3 - (tj1 - tj2)) * (tj3 + (tj1 - tj2)) / 4.0;
            double f2 = ((tj1 + tj2 + 2) - tj3) * ((tj1 + tj2 + 2) + tj3) / 4.0;
            double f3 = (tj3 + tm3) * (tj3 - tm3) / 4.0;
            return (Math.Sqrt(f1 * f2 * f3));
        }

        private static double ShultenGodron_B (int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {
            double f = (tm3 * (tj1 * (tj1 + 2) - tj2 * (tj2 + 2)) + (tm1 - tm2) * tj3 * (tj3 + 2)) / 8.0;
            return (-(tj3 + 1) * f);
        }

        // SixJ

        /// <summary>
        /// Computes the value of the 6j symbol for the six given spins.
        /// </summary>
        /// <param name="j1">Upper left spin.</param>
        /// <param name="j2">Upper middle spin.</param>
        /// <param name="j3">Upper right spin.</param>
        /// <param name="j4">Lower left spin.</param>
        /// <param name="j5">Lower middle spin.</param>
        /// <param name="j6">Lower right spin.</param>
        /// <returns>The value of {{j1,j2,j3},{j4,j5,j6}}.</returns>
        public static double SixJ (Spin j1, Spin j2, Spin j3, Spin j4, Spin j5, Spin j6) {

            // check for required triangle relations
            bool t = (Triangle(j1, j2, j3) && Triangle(j1, j5, j6) && Triangle(j4, j2, j6) && Triangle(j4, j5, j3));
            if (!t) return (0.0);


            // move the smallest entry to the lower-right corner
            // move lowest entry in first column to bottom
            if (j1.TwoJ < j4.TwoJ) { Swap(ref j1, ref j4); Swap(ref j2, ref j5); }
            // move lowest entry in second column to bottom
            if (j2.TwoJ < j5.TwoJ) { Swap(ref j2, ref j5); Swap(ref j3, ref j6); }
            // move lowest of first two entries in lower row to right
            if (j4.TwoJ < j5.TwoJ) { Swap(ref j1, ref j2); Swap(ref j4, ref j5); }
            // move lowest entry in third column to bottom
            if (j3.TwoJ < j6.TwoJ) { Swap(ref j3, ref j6); Swap(ref j1, ref j4); }
            // move lowest of last two entries in lower row to right
            if (j5.TwoJ < j6.TwoJ) { Swap(ref j2, ref j3); Swap(ref j5, ref j6); }


            if (j6.TwoJ == 0) {
                // special case for 0 entry
                return (SixJ_Zero(j1, j2, j3));
            } else if (j6.TwoJ == 1) {
                // special case for 1/2 entry
                return (SixJ_OneHalf(j1, j2, j3, j4, j5));
            } else {
                // general case
                return (SixJ_ShultenGorton_Recurse(j1.TwoJ, j2.TwoJ, j3.TwoJ, j4.TwoJ, j5.TwoJ, j6.TwoJ));
            }

        }

        private static bool Triangle (Spin a, Spin b, Spin c) {

            // check integer relationships
            int s = (a.TwoJ + b.TwoJ + c.TwoJ);
            if (s % 2 != 0) return (false);

            // check triangle inequality
            return ((Math.Abs(a.TwoJ - b.TwoJ) <= c.TwoJ) && (c.TwoJ <= (a.TwoJ + b.TwoJ)));

        }

        // j1 j2 j3
        // j2 j1 0

        private static double SixJ_Zero (Spin j1, Spin j2, Spin j3) {
            double f = 1.0 / Math.Sqrt((j1.TwoJ + 1) * (j2.TwoJ + 1));
            if ((j1.TwoJ + j2.TwoJ + j3.TwoJ) / 2 % 2 != 0) f = -f;
            return (f);
        }

        // j1 j2 j3
        // j4 j5 1/2

        private static double SixJ_OneHalf (Spin j1, Spin j2, Spin j3, Spin j4, Spin j5) {

            Debug.Assert(Math.Abs(j1.TwoJ - j5.TwoJ) == 1); Debug.Assert(Math.Abs(j2.TwoJ - j4.TwoJ) == 1);

            // make j5 = j1-1/2, reducing the number of formulae required
            if (j5.TwoJ > j1.TwoJ) { Swap(ref j1, ref j5); Swap(ref j2, ref j4); }

            double f;
            if (j4.TwoJ < j2.TwoJ) {
                //   a     b    c
                // b-1/2 a-1/2 1/2 
                f = (j1.TwoJ + j2.TwoJ - j3.TwoJ) * (j1.TwoJ + j2.TwoJ + j3.TwoJ + 2) / 4.0 / (j1.TwoJ * (j1.TwoJ + 1)) / (j2.TwoJ * (j2.TwoJ + 1));

            } else {
                //   a     b    c
                // b+1/2 a-1/2 1/2
                f = (j1.TwoJ - j2.TwoJ + j3.TwoJ) * (j2.TwoJ + j3.TwoJ - j1.TwoJ + 2) / 4.0 / (j1.TwoJ * (j1.TwoJ + 1)) / ((j2.TwoJ + 1) * (j2.TwoJ + 2));

            }
            f = Math.Sqrt(f);
            if ((j1.TwoJ + j2.TwoJ + j3.TwoJ) % 4 != 0) f = -f;
            return (f);

            /*
            // ensure that J1, J2 are integer-spins; J3, J4 are half-spins
            if (j1.TwoJ % 2 != 0) { Swap(ref j1, ref j4); Swap(ref j2, ref j5); }

            // apply relation to 3j symbols
            // re-write to call directly into specialized 3j routines
            if ((j1.TwoJ + j2.TwoJ + j3.TwoJ) % 4 == 0) {
                double p = ThreeJ(new SpinState(j4, 0.5), new SpinState(j5, -0.5), new SpinState(j3, 0.0));
                double q = ThreeJ(new SpinState(j1, 0.0), new SpinState(j2, 0.0), new SpinState(j3, 0.0));
                double r = Math.Sqrt((j1.TwoJ + 1) * (j2.TwoJ + 1));
                return (- p / q / r);
            } else {
                return (0.0);
            }
            */
        }

        private static double SixJ_ShultenGorton_Recurse (int tj1, int tj2, int tj3, int tj4, int tj5, int tj6) {

            // determine minimum
            int tk_min;
            int tk_min_23 = Math.Abs(tj2 - tj3);
            int tk_min_56 = Math.Abs(tj5 - tj6);
            if (tk_min_23 > tk_min_56) {
                tk_min = tk_min_23;
            } else {
                tk_min = tk_min_56;
            }

            // determine maximum
            int tk_max;
            int tk_max_23 = tj2 + tj3;
            int tk_max_56 = tj5 + tj6;
            if (tk_max_23 < tk_max_56) {
                tk_max = tk_max_23;
            } else {
                tk_max = tk_max_56;
            }

            // find the midpoint of the range
            int tk_mid = (tk_min + tk_max) / 2;
            if ((tk_mid % 2) != (tk_max % 2)) tk_mid++;

            //Console.WriteLine("{0} <= {1} <= {2}", tk_min, tk_mid, tk_max);

            // variables to store recurrence coefficients and (unscaled) 6j values
            double ap, a0, b0;
            double cm, c0, cp;
            int tj;

            // a variable to store the sought value
            double c = 0.0;

            // iterate in from the right
            double NR = 0.0;
            cp = 0.0;
            ap = 0.0;
            c0 = 0.5 / tk_max / tk_max;

            // fix sign
            if ((tj2 + tj3 + tj5 + tj6) % 4 != 0) { c0 = -c0; }

            tj = tk_max;
            while (true) {

                // keep track of the normalization constant
                NR += (tj + 1) * c0 * c0;

                // remember the (unscaled) value of the desired 3j symbol
                if (tj == tj1) {
                    c = c0;
                }

                // check whether we have reached the center
                if (tj <= tk_mid) {
                    if (Math.Abs(c0) <= Math.Pow(2.0, -40.0) * Math.Abs(cp)) {
                        // if the middle value is zero, the left and right values cannot be matached
                        // (acutally, the numbers will be two random tiny numbers, and their ratio will
                        // be a random number of order one); in this case, go down one more step
                        tk_mid -= 2;
                    } else {
                        break;
                    }
                }


                // compute required coefficients
                b0 = ShultenGordon_F(tj, tj2, tj3, tj4, tj5, tj6);
                a0 = ShultenGordon_E(tj, tj2, tj3, tj4, tj5, tj6);

                // use recursion to compute the tj-2 symbol
                cm = -(b0 * c0 + (tj / 2.0) * ap * cp) / ((tj + 2) / 2.0) / a0;

                //Console.WriteLine("{0}: {1} {2} + {3} {4} + {5} {6}", tj, a0, cm, b0, c0, ap, cp);

                // prepare for the next cycle
                cp = c0;
                c0 = cm;
                ap = a0;
                tj -= 2;

            }

            // remember the (unscaled) middle value
            double c1 = c0;

            // iterate in from the left
            double NL = 0.0;
            cm = 0.0;
            a0 = 0.0;
            c0 = 0.5 / tk_max / tk_max;
            tj = tk_min;

            // j_min = 0 iff j2 == j3 && j5 == j6
            // in this case all three coefficients vanish, so the recursion is indeterminate,
            // so we do the first recursion step "by hand" using explicit expressions
            if (tj == 0) {
                Debug.Assert(tj2 == tj3); Debug.Assert(tj5 == tj6);
                NL += c0 * c0;
                cm = c0;
                c0 = (tj4 * (tj4 + 2) - tj3 * (tj3 + 2) - tj6 * (tj6 + 2)) / Math.Sqrt(tj3 * (tj3 + 2)) / Math.Sqrt(tj6 * (tj6 + 2)) / 2.0 * cm;
                a0 = ShultenGordon_E(2, tj2, tj3, tj4, tj5, tj6);
                tj = 2;
            }

            while (true) {

                // check whether we have reached the center
                if (tj >= tk_mid) {
                    break;
                }

                // notice that we terminate before checking whether we are at the desired value,
                // or adding to the normalization constant; this is different from the previous loop
                // we do this so that the middle element does not double contribute to the normalization
                // and because if the middle element is desired, it has already been stored from the previous loop

                // remember the desired (unscaled) value
                if (tj == tj1) {
                    c = c0;
                    //Console.WriteLine("{0}: {1}", tj1, c);
                }

                // add current value to normalization
                NL += (tj + 1) * c0 * c0;

                // check for end
                // this should never happen, no? we just checked before.
                if (tj >= tk_max) break;

                // compute required coeficients
                b0 = ShultenGordon_F(tj, tj2, tj3, tj4, tj5, tj6);
                ap = ShultenGordon_E(tj + 2, tj2, tj3, tj4, tj5, tj6);
                //b0 = ShultenGodron_B(s1.TwoJ, s2.TwoJ, tj, s1.TwoM, s2.TwoM, s3.TwoM);
                //ap = ShultenGordon_A(s1.TwoJ, s2.TwoJ, tj + 2, s3.TwoM);

                // use recursion relation to compute tj+2
                //cp = (-b0 * c0 - ((tj + 2) / 2.0) * a0 * cm) / (tj / 2.0) / ap;
                cp = -(b0 * c0 + ((tj + 2) / 2.0) * a0 * cm) / (tj / 2.0) / ap;

                //Console.WriteLine("{0}: {1} {2} + {3} {4} + {5} {6}", tj, a0, cm, b0, c0, ap, cp);

                // prepare for next iteration
                cm = c0;
                c0 = cp;
                a0 = ap;
                tj += 2;

            }

            // match in the middle
            double r = c0 / c1;
            NL = NL / (r * r);
            if (tj1 < tk_mid) c = c / r;

            // normalize
            double N = NL + NR;
            c = c / Math.Sqrt(N * (tj4 + 1));

            return (c);

        }

        private static double ShultenGordon_E (int tj1, int tj2, int tj3, int tl1, int tl2, int tl3) {
            double f1 = (tj1 + tj2 - tj3) * (tj1 - tj2 + tj3) / 4.0;
            double f2 = (tj2 + tj3 + 2 + tj1) * (tj2 + tj3 + 2 - tj1) / 4.0;
            double f3 = (tj1 - tl2 + tl3) * (tj1 + tl2 - tl3) / 4.0;
            double f4 = (tl2 + tl3 + 2 + tj1) * (tl2 + tl3 + 2 - tj1) / 4.0;
            return (Math.Sqrt(f1 * f2 * f3 * f4));
        }

        private static double ShultenGordon_F (int tj1, int tj2, int tj3, int tj4, int tj5, int tj6) {
            double pj1 = tj1 * (tj1 + 2);
            double pj2 = tj2 * (tj2 + 2);
            double pj3 = tj3 * (tj3 + 2);
            double pj4 = tj4 * (tj4 + 2);
            double pj5 = tj5 * (tj5 + 2);
            double pj6 = tj6 * (tj6 + 2);
            double s = pj1 * (-pj1 + pj2 + pj3) + pj5 * (pj1 + pj2 - pj3) + pj6 * (pj1 - pj2 + pj3) - 2 * pj1 * pj4;
            return ((tj1 + 1) * s / 16.0);
            // don't to this all in integer math, or you will silently overflow!
        }

    }

}
