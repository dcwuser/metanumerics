using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Statistics.Distributions {

    /// <summary>
    /// Contains methods for converting between different kinds of moments.
    /// </summary>
    /// <remarks>
    /// Distributions can be described by by sets of various moments. Raw moments
    /// and central moments are encountered most commonly, but cumulants (also
    /// called semi-invariants) and factorial moments are also seen. Given any
    /// set of one kind of moment for a distribution up to a given order
    /// (e.g. the 1st, 2nd, and 3rd central moments), the methods of this
    /// class return the values of other kinds of moments, up to the same order
    /// (e.g. the 1st, 2nd, and 3rd raw moments). 
    /// </remarks>
    public static class MomentMath {

        /// <summary>
        /// Converts raw moments to central moments.
        /// </summary>
        /// <param name="M">A set of raw moments.</param>
        /// <returns>The corresponding set of central moments.</returns>
        /// <remarks>The computation of central moments from raw moments is often subject
        /// to significant cancellation errors, so you should be wary of the higher
        /// digits of central moments returned by this method.</remarks>
        /// <exception cref="ArgumentNullException"><paramref name="M"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The zeroth raw moment is not one.</exception>
        public static double[] RawToCentral (double[] M) {

            if (M == null) throw new ArgumentNullException(nameof(M));

            double[] C = new double[M.Length];

            if (C.Length == 0) return (C);

            if (M[0] != 1.0) throw new ArgumentOutOfRangeException(nameof(M));

            C[0] = 1.0;

            if (C.Length == 1) return (C);

            C[1] = 0.0;

            for (int r = 2; r < C.Length; r++) {
                C[r] = RawToCentral(M, r);
            }

            return (C);
        }

        // Start from definition C_r = < (x-M1)^r > and use the binomial formula to write (x - M1)^r in terms of powers of x and M1.
        //   C_r = \sum_{k=0}^{n} {n \choose k} (-\mu)^k M_{r-k}

        internal static double RawToCentral (double[] M, int r) {

            Debug.Assert(M.Length > 1);
            Debug.Assert(0 <= r && r < M.Length);
            Debug.Assert(M[0] == 1.0);

            double mmu = -M[1];
            IEnumerator<double> B = AdvancedIntegerMath.BinomialCoefficients(r).GetEnumerator();

            // k = 0 term is M_r
            B.MoveNext();
            int i = r; // tracks r - k
            double mmuk = 1.0; // tracks (-\mu)^k
            double C = M[r];

            while (B.MoveNext()) {
                i--;
                mmuk *= mmu;
                C += mmuk * B.Current * M[i];
            }

            // final term is (-\mu)^r {r \choose r} M_0 = (-\mu)^r; this will be wrong if M[0] != 1 

            return (C);

        }

        /// <summary>
        /// Converts central moments to raw moments.
        /// </summary>
        /// <param name="mu">The mean.</param>
        /// <param name="C">A set of central moments.</param>
        /// <returns>The corresponding set of central moments.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="C"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The zeroth central moment is not one, or the first central moment is not zero.</exception>
        public static double[] CentralToRaw (double mu, double[] C) {

            if (C == null) throw new ArgumentNullException(nameof(C));

            double[] M = new double[C.Length];

            if (M.Length == 0) return (M);

            if (C[0] != 1.0) throw new ArgumentOutOfRangeException(nameof(C));

            M[0] = 1.0;

            if (M.Length == 1) return (M);

            if (C[1] != 0.0) throw new ArgumentOutOfRangeException(nameof(C));

            M[1] = mu;

            for (int r = 2; r < M.Length; r++) {
                M[r] = CentralToRaw(mu, C, r);
            }

            return (M);

        }

        // M_r = \sum_{k=0}^{r} {r \choose k} \mu^k C_{r-k}

        internal static double CentralToRaw (double mu, double[] C, int r) {

            Debug.Assert(C.Length > 0);
            Debug.Assert(0 <= r && r < C.Length);
            Debug.Assert(C[0] == 1.0);
            //Debug.Assert(C[1] == 0.0);

            IEnumerator<double> B = AdvancedIntegerMath.BinomialCoefficients(r).GetEnumerator();

            // k = 0 term
            B.MoveNext();
            int kp = r; // tracks k' = k - r
            double muk = 1.0; // tracks \mu^k
            double M = C[r];

            while (B.MoveNext()) {
                kp--;
                muk *= mu;
                M += B.Current * muk * C[kp];
            }

            // final term is \mu^r {r \choose r} C_0 = \mu^r; this will be wrong if C[0] != 1
            // second to final term is 0; this will be wrong if C[1] != 0

            return (M);

        }

        /// <summary>
        /// Converts cumulants to raw moments.
        /// </summary>
        /// <param name="K">A set of cumulants.</param>
        /// <returns>The corresponding set of raw moments.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="K"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The zeroth cumulant is not zero.</exception>
        public static double[] CumulantToRaw (double[] K) {

            if (K == null) throw new ArgumentNullException(nameof(K));

            double[] M = new double[K.Length];

            if (M.Length == 0) return (M);

            if (K[0] != 0.0) throw new ArgumentOutOfRangeException(nameof(K));

            M[0] = 1.0;

            for (int r = 1; r < M.Length; r++) {
                M[r] = CumulantToRaw(K, r);
            }

            return (M);
        }


        internal static double CumulantToRaw (double[] K, int r) {
            return (CumulantToMoment(K, r, false));
        }

        /// <summary>
        /// Converts cumulants to central moments.
        /// </summary>
        /// <param name="K">A set of cumulants.</param>
        /// <returns>The corresponding set of central moments.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="K"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The zeroth cumulant is not zero.</exception>
        public static double[] CumulantToCentral (double[] K) {

            if (K == null) throw new ArgumentNullException(nameof(K));

            double[] C = new double[K.Length];

            if (C.Length == 0) return (C);

            if (K[0] != 0.0) throw new ArgumentOutOfRangeException(nameof(K));

            C[0] = 1.0;

            if (C.Length == 1) return (C);

            C[1] = 0.0;

            for (int r = 2; r < C.Length; r++) {
                C[r] = CumulantToCentral(K, r);
            }

            return (C);
        }

        internal static double CumulantToCentral (double[] K, int r) {
            return (CumulantToMoment(K, r, true));
        }

        // Faa di Bruno's formula expresses raw moments in terms of cumulants.

        // M_r = \sum_{m_1 + \cdots + m_k = r} \frac{r!}{m_1 \cdots m_k} K_1 \cdots K_k

        // That is: take all partitions of r. (E.g. for r = 4, thre are 5 partitions: 1 + 1 + 1 + 1 = 4, 1 + 1 + 2 = 4, 2 + 2 = 4, 1 + 3 = 4, and 4 = 4).
        // Each partition will contribute one term. Each appearance of an integer k in the partition will contribute one factor of the kth cumulant
        // to that term. (E.g. K_1^4, K_1^2 K_2, K_2^2, K_1 K_3, and K_4.)
        // Each term has a combinatoric factor of r! divided by the product of all the integers in the partition. (E.g. 1^4, 1^2 * 2, 2^2, 1 * 3, and 4.)

        internal static double CumulantToMoment (double[] K, int r, bool central) {
            Debug.Assert(K != null);
            Debug.Assert(r > 0);
            Debug.Assert(K.Length >= r);
            double M = 0.0;
            foreach (int[] partition in AdvancedIntegerMath.InternalPartitions(r)) {
                double dM = AdvancedIntegerMath.Factorial(r);
                int u = 0; // tracks the last observed partition member
                int m = 1; // tracks the multiplicity of the current partition member
                foreach (int v in partition) {

                    // if we are computing central moments, ignore all terms with a factor of K_1
                    if (central && (v == 1)) {
                        dM = 0.0;
                        break;
                    }

                    if (v == u) {
                        m++;
                    } else {
                        m = 1;
                    }
                    dM *= K[v] / AdvancedIntegerMath.Factorial(v) / m;
                    u = v;
                }
                M += dM;
            }
            return (M);

            // This method of reading out a partition is a bit klunky. We really want to know the multiplicity of each element,
            // but that isn't straightforwardly available. In fact, this method will break if identical elements are not adjacent
            // in the array (e.g. 4 = 1 + 2 + 1 + 1). Would be nice to address this, perhaps by actually returning a multiplicity
            // representation of parititions.

        }

        /// <summary>
        /// Converts raw moments to cumulants.
        /// </summary>
        /// <param name="M">A set of raw moments.</param>
        /// <returns>The corresponding set of cumulants.</returns>
        /// <remarks>The computation of cumulants from raw moments is often subject
        /// to significant cancellation errors, so you should be wary of the higher
        /// digits of cumulants returned by this method.</remarks>
        /// <exception cref="ArgumentNullException"><paramref name="M"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The zeroth raw moment is not one.</exception>
        public static double[] RawToCumulant (double[] M) {

            if (M == null) throw new ArgumentNullException(nameof(M));

            double[] K = new double[M.Length];

            if (K.Length == 0) return (K);

            if (M[0] != 1.0) throw new ArgumentOutOfRangeException(nameof(M));

            K[0] = 0.0;

            if (K.Length == 1) return (K);

            for (int r = 0; r < K.Length - 1; r++) {
                double t = M[r + 1];
                for (int s = 0; s < r; s++) {
                    t -= AdvancedIntegerMath.BinomialCoefficient(r, s) * K[s + 1] * M[r - s];
                }
                K[r + 1] = t;
            }

            return (K);
        }

        /// <summary>
        /// Converts central moments to cumulants.
        /// </summary>
        /// <param name="mu">The mean.</param>
        /// <param name="C">A set of central moments.</param>
        /// <returns>The corresponding set of cumulants.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="C"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The zeroth central moment is not one, or the first central moment is not zero.</exception>
        public static double[] CentralToCumulant (double mu, double[] C) {

            if (C == null) throw new ArgumentNullException(nameof(C));

            double[] K = new double[C.Length];
            if (K.Length == 0) return (K);

            // C0 = 1 and K0 = 0
            if (C[0] != 1.0) throw new ArgumentOutOfRangeException(nameof(C));
            K[0] = 0.0;
            if (K.Length == 1) return (K);

            // C1 = 0 and K1 = M1
            if (C[1] != 0.0) throw new ArgumentOutOfRangeException(nameof(C));
            K[1] = mu;
            if (K.Length == 2) return (K);

            // Determine higher K
            // s = 0 term involves K1 = M1, so ignore
            // s = r - 1 term involves C1 = 0, so ignore
            for (int r = 1; r < K.Length - 1; r++) {
                double t = C[r + 1];
                for (int s = 1; s < r - 1; s++) {
                    t -= AdvancedIntegerMath.BinomialCoefficient(r, s) * K[s + 1] * C[r - s];
                }
                K[r + 1] = t;
            }

            return (K);
        }

    }

}
