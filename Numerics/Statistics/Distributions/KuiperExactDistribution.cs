using System;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;

namespace Meta.Numerics.Statistics.Distributions {

    // w = n V, support 1 < w < n

    // n = 2
    // 1 < w < 2    P = 1 (w - 1) = 1 - 1 (2 - w)

    // n = 3
    // 1 < w < 2    P = 2/3 (w - 1)^2
    // 2 < w < 3    P = 2/3 (-3 + 3 w - 1/2 w^2) = 1 - 1/3 (3 - w)^2

    // n = 4
    // 1 < w < 2    P = 3/8 (w - 1)^3
    // 2 < w < 3    P = 3/8 (1 - 4 w + 3 w^2 - 1/2 w^3)
    // 3 < w < 4    P = 3/8 (-8 + 8 w - 2 w^2 + 1/6 w^3) = 1 - 1/16 (4 - w)^3

    // n = 5
    // 1 < w < 2    P = 24/125 (w - 1)^4
    // 2 < 2 < 3    P = 24/125 (5 - 7 w + 3/2 w^2 + w^3 - 1/4 w^4)
    // 3 < w < 4    P = 24/125 (1/2 - 17/2 w + 31/4 w^2 - 2 w^3 + 1/6 w^4)
    // 4 < w < 5    P = 24/125 (-125/6 + 125/6 w - 25/4 w^2 + 5/6 w^3 - 1/24 w^4) = 1 - 1/125 (5 - w)^4

    // n = 6
    // 1 < w < 2    P = 5/54 (w - 1)^5
    // 2 < w < 3    P = 5/54 (-7 + 22 w - 23 w^2 + 19/2 w^3 - 5/4 w^4 + 0 w^5)
    // 3 < w < 4    P = 5/54 (35/4 - 65/4 w + 11/2 w^2 + 5/3 w^3 - 5/6 w^4 + 1/12 w^5)
    // 4 < w < 5    P = 5/54 (-23/12 - 227/12 w + 39/2 w^2 - 37/6 w^3 + 5/6 w^4 - 1/24 w^5)
    // 5 < w < 6    P = 5/54 (-54 + 54 w - 18 w^2 + 3 w^3 - 1/4 w^4 + 1/120 w^5) = 1 - 1/1296 (6 - w)^5

    // n = 7
    // 1 < w < 2    P = 720/16807 (w - 1)^6
    // 2 < w < 3    P = 720/16807 (-3 - 5 w + 51/2 w^2 - 28 w^3 + 25/2 w^4 - 9/4 w^5 + 1/8 w^6)
    // 3 < w < 4    P = 720/16807 (105/4 - 29 w + 5/2 w^2 + 9/2 w^3 - 5/8 w^4 - 1/6 w^5 + 1/36 w^6)
    // 4 < w < 5    P = 720/16807 (187/12 - 115/3 w + 97/6 w^2 + 91/36 w^3 - 35/16 w^4 + 3/8 w^5 - 1/48 w^6)
    // 5 < w < 6    P = 720/16807 (-251/24 - 1045/24 w + 2351/48 w^2 - 629/36 w^3 + 143/48 w^4 - 1/4 w^5 + 1/120 w^6)
    // 6 < w < 7    P = 1 - 1/16807 (7 - w)^6

    // At the ends P_n(1 < w < 2) = n! / n^(n-1) (w-1)^(n-1) and Q_n(n-1 < w < n) = (n-w)^(n-1) / n^(n-2)

    // Taking a derivative wrt w gives p(w) and integrating p(w) * w^m gives the mth moment.

    // n = 2: <w> = 3/2, <w^2> = 7/3, <w^3> = 15/4 => C2 = 1/12, C3 = 0
    // n = 3: <w> = 17/9 , <w^2> = 67/18, <w^3> = 229/30 => C2 = 25/162, C3 = 71/3645
    // n = 4: <w> = 71/32, <w^2> = 103/20, <w^3> = 1997/160 => C2 = 1163/5120, C3 = 3827/81920
    // n = 5: <w> = 1569/625, <w^2> = 2476/375, <w^3> = 11352/625 => C2 = 352217/1171875
    // n = 6: <w> = 899/324, <w^2> = 7847/972, <w^3> = 223099/9072 => C2 = 39275/104976
    // n = 7: <w> = 355081/117649, <w^2> = 1124366/117649, <w^3> = 11189840/352947 => C2 = 6198018973/13841287201

    // Moments provide an additional check. We should have <1> = 1 and <w> = \frac{n!}{n^n} \sum_{k=0}^{n-1} \frac{n^k}{k!}.
    // So far so good.

    internal class KuiperExactDistribution : ContinuousDistribution {

        public KuiperExactDistribution(int n) {
            if (n < 2) throw new ArgumentOutOfRangeException(nameof(n));
            this.n = n;
        }

        private readonly int n;

        public override Interval Support {
            get {
                return Interval.FromEndpoints(1.0, n);
            }
        }

        public override double Mean {
            get {
                return ComputeMean();
            }
        }

        // As per Stephens, Biopmetrika (1965) 52, p. 309-, the mean
        //   <w> = \frac{n!}{n^n} \sum_{k=0}^{n-1} \frac{n^k}{k!}
        // This is derived by taking Birnbaum's expression of the mean of D_{\pm} and noting that <V> = <(D_+ + D_-)> = <D_+> + < D_->.
        // It's pretty cool that we can get this in simple analytic form. I'd love to do the same for <D> or <V^2>, but D=\max(D_+,D_-) does
        // not seperate into terms each involving only D_+ or D_- and V^2 involves the cross-term D_+ D_-, so both depend on the joint
        // distribution p(D_+,D_-).

        private double ComputeMean() {
            double t = 1.0;
            double s = t;
            for (int k = 1; k < n; k++) {
                t = t * n / k;
                s += t;
            }
            // Canceling factors of n, \frac{n!}{n^n} = \frac{(n-1)!}{n^{n-1}}, so the prefactor is actually just the inverse of the last term.
            return s / t;
        }

        public override double Variance {
            get {
                switch (n) {
                    case 2:
                        // 1/12 = 0.041667 * 2
                        return 1.0 / 12.0;
                    case 3:
                        // 25/162 = 0.051440 * 3
                        return 25.0 / 162.0;
                    case 4:
                        // 1163/5120 = 0.056787 * 4
                        return 1163.0 / 5120.0;
                    case 5:
                        // 352217/1171875 = 0.060111 * 5
                        return 352217.0 / 1171875.0;
                    case 6:
                        // 39275/104976 = 0.062355 * 6
                        return 39275.0 / 104976.0;
                    case 7:
                        // 6198018973/13841287201 = 0.063970 * 7
                        return 6198018973.0 / 13841287201.0;
                    default:
                        // \frac{\pi}{2} \left( \frac{\pi}{3} - 1 \right) = 0.074138
                        // this will perform numerical integration, which will be very slow
                        return base.CentralMoment(2);
                }
            }
        }

        public override double LeftProbability(double w) {
            if (w <= 1.0) {
                return 0.0;
            } else if (w < n) {
                return DurbinMatrixP(w);
            } else {
                return 1.0;
            }
        }

        public override double RightProbability(double w) {
            if (w <= 1.0) {
                return 1.0;
            } else if (w < n) {
                return 1.0 - DurbinMatrixP(w);
            } else {
                return 0.0;
            }
        }

        public override double ProbabilityDensity(double w) {
            if (w <= 1.0) {
                return 0.0;
            } else if (w < n) {
                return DurbinMatrixPPrime(w);
            } else {
                return 0.0;
            }
        }

        ///<inheritdoc />
        public override double RawMoment(int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else if (r == 1) {
                return Mean;
            } else {
                return base.RawMoment(r);
            }
        }

        ///<inheritdoc />
        public override double CentralMoment(int r) {
            if (r < 0) {
                throw new ArgumentOutOfRangeException(nameof(r));
            } else if (r == 0) {
                return 1.0;
            } else if (r == 1) {
                return 0.0;
            } else if (r == 2) {
                return this.Variance;
            } else {
                return base.CentralMoment(r);
            }
        }

        // Durbin derives a matrix result for the Kuiper statistic analogous to his result for the Kolmogorov-Smirnov that was used by Marsaglia.
        // In this case, if w = n V = k - h, H is k X k with entries of the form
        //       { 0  1           0           0           0        }
        //       { 0  1/1!        1           0           0        }
        //   H = { 0  1/2!        1/1!        1           0        }
        //       { 0  1/3!        1/2!        1/1!        1        }
        //       { 0  (1-h^4)/4!  (1-h^3)/3!  (1-h^2)/2!  (1-h)/1! }
        // and
        //   P = \frac{n!}{n^{n-1}} \left( H^n \right)_{1,2}
        // Note there is no discontinuity at half-integer values of w in the Kuiper case.
        // See Durbin, "Distribution Theory for Tests Based on the Sample Distribution Function", 1973, p. 12-13, 35.

        private double DurbinMatrixP(double w) {

            int k = (int) Math.Ceiling(w);
            double h = k - w;

            SquareMatrix H = GetDurbinMatrix(k, h);

            SquareMatrix Hn = H.Power(n);

            double f = AdvancedIntegerMath.Factorial(n - 1) / MoreMath.Pow(n, n - 2);
            return f * Hn[0, 1];

        }

        private double DurbinMatrixPPrime(double w) {

            int k = (int) Math.Ceiling(w);
            double h = k - w;

            // Compute derivative of H
            SquareMatrix DH = GetDurbinMatrixPrime(k, h);

            // Compute powers of H
            SquareMatrix[] PowerH = new SquareMatrix[n];
            PowerH[1] = GetDurbinMatrix(k, h);
            for (int i = 2; i < n; i++) {
                PowerH[i] = PowerH[1] * PowerH[i - 1];
            }

            // Use D(H^n) = (DH) H^(n-1) + H (DH) H^(n-2) + H^2 (DH) H^(n-3) + \cdots + H^(n-2) (DH) H + H^(n-1) (DH)
            SquareMatrix HnP = DH * PowerH[n - 1];
            for (int i = 1; i < (n - 1); i++) {
                HnP += PowerH[i] * DH * PowerH[n - 1 - i];
            }
            HnP += PowerH[n - 1] * DH;

            double f = AdvancedIntegerMath.Factorial(n - 1) / MoreMath.Pow(n, n - 2);
            return f * HnP[0, 1];


        }

        private static SquareMatrix GetDurbinMatrix(int k, double h) {

            // dimension of matrix
            int m = k;
            SquareMatrix H = new SquareMatrix(m);

            // populate the matrix along diagonals, since they share factorial factors

            // first superdiagonal is all 1s
            for (int j = 1; j < m; j++) {
                H[j - 1, j] = 1.0;
            }

            // bottom row and diagonals
            double Fi = 1.0; // variable for 1/i!
            double hi = h; // variable for h^i
            for (int i = 1; i < m; i++) {
                // bottom row
                H[m - 1, m - i] = Fi * (1.0 - hi);
                // diagonal
                for (int j = i + 1; j < m; j++) {
                    H[j - 1, j - i] = Fi;
                }
                // prepare for the next recurrsion
                hi = hi * h;
                Fi = Fi / (i + 1);
            }

            return H;

        }

        private static SquareMatrix GetDurbinMatrixPrime(int k, double h) {

            // dimension of matrix
            int m = k;
            SquareMatrix DH = new SquareMatrix(m);

            // only lower row of H is non-constant, hence only lower row of H' is non-zero
            double Fi = 1.0;
            double hi = 1.0;
            for (int i = 1; i < m; i++) {
                DH[m - 1, m - i] = Fi * hi;
                Fi /= i;
                hi *= h;
            }

            return DH;

        }

    }

}
