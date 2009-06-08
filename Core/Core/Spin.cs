using System;
using System.Diagnostics;

using Meta.Numerics.Functions;

namespace Meta.Numerics {

#if FUTURE

    public struct Spin {

        // construction

        public Spin (double j) {

            // no negative (or too large) spins
            double tj = 2.0 * j;
            if ((tj < 0.0) || (tj > Int32.MaxValue)) throw new ArgumentOutOfRangeException("j");

            // spin must be integer or half-integer
            double tjt = Math.Truncate(tj);
            if (tjt != tj) throw new ArgumentOutOfRangeException("j");

            // store 2*j
            twoJ = (int) tjt;
            
        }

        // stored data

        private int twoJ;

        // accessors

        internal int TwoJ {
            get {
                return (TwoJ);
            }
        }

        // casting

        public static implicit operator Spin (double j) {
            return (new Spin(j));
        }

        public static implicit operator double (Spin s) {
            return (s.TwoJ / 2.0);
        }

    }

    public struct SpinState {

        public SpinState (double j, double m) {

            // no negative (or too large) spins
            double tj = 2.0 * j;
            if ((tj < 0.0) || (tj > Int32.MaxValue)) throw new ArgumentOutOfRangeException("j");

            // spin must be integer or half-integer
            double tjt = Math.Truncate(tj);
            if (tjt != tj) throw new ArgumentOutOfRangeException("j");

            // store 2*j
            twoJ = (int) tjt;

            // M must be an integer
            double tm = 2.0 * m;
            double tmt = Math.Truncate(tm);
            if (tmt != tm) throw new ArgumentOutOfRangeException("m");

            twoM = (int) tmt;

            // -J <= M <= J
            if (Math.Abs(twoM) > Math.Abs(twoJ)) throw new ArgumentOutOfRangeException("m");

            // half-integer J requires half-integer M; integer J requires integer M
            if ((twoJ % 2) != Math.Abs(twoM % 2)) throw new ArgumentOutOfRangeException("m");
            
        }

        private int twoJ;

        private int twoM;

        internal int TwoJ {
            get {
                return (twoJ);
            }
        }

        internal int TwoM {
            get {
                return (twoM);
            }
        }

        internal int JPlusM {
            get {
                return ((twoJ + twoM) / 2);
            }
        }

        internal int JMinusM {
            get {
                return ((twoJ - twoM) / 2);
            }
        }

        internal SpinState Invert () {
            twoM = -twoM;
            return (this);
        }

    }


    /// <summary>
    /// Contains methods for computing functions of spin and spin states.
    /// </summary>
    public static class SpinMath {

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
            if (s3.TwoJ < Math.Abs(s1.TwoJ - s2.TwoJ)) throw new InvalidOperationException();
            if (s3.TwoJ > (s1.TwoJ + s2.TwoJ)) throw new InvalidOperationException();

            // require that Ms add to zero
            if ((s1.TwoM + s2.TwoM + s3.TwoM) != 0) return(0.0);

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
            } else if (s1.TwoJ == 3) {
                // J1 = 3/2
                return (s * ThreeJ_Any(s1, s2, s3));
                // replace this with special case!
            } else if ((s1.TwoJ + s2.TwoJ) == s3.TwoJ) {
                // J1+J2=J3, maximum coupling
                return (s * ThreeJ_Any(s1, s2, s3));
                //return (ThreeJ_MaxJ(s1, s2));
                // should this really be special cased? there is only one term in sum.
            } else {
                // arbitrary Js
                if ((s1.TwoM == 0) && (s2.TwoM == 0) && (s3.TwoM == 0)) {
                    // special case of all M's zero
                    return (s * ThreeJ_ZeroM_Int(s1, s2, s3));
                    // also for other low Ms, because there are many terms to compute these cases
                } else {
                    // general case
                    return (s * ThreeJ_Any(s1, s2, s3));
                }
            }

            throw new NotImplementedException();
        }

        // exchange two references

        private static void Exchange (ref SpinState s1, ref SpinState s2) {
            SpinState t = s1;
            s1 = s2;
            s2 = t;
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
                double r = ThreeJ_OneJ(s1.Invert(),s2.Invert(),s3.Invert());
                if ((s1.TwoJ + s2.TwoJ + s3.TwoJ)/2 % 2 != 0) r = -r;
                return(r);
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

                Console.WriteLine("A={0} B={1} C={2}", A, B, C);

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
                Console.WriteLine("g={0}", g);
                g = Math.Sqrt(g);

                if ((j/2) % 2 != 0) f = -f;

                return (f * g);
            }


        }

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

        // SixJ

        public static double SixJ (Spin j1, Spin j2, Spin j3, Spin j4, Spin J5, Spin J6) {
            throw new NotImplementedException();
        }

    }

#endif

}
