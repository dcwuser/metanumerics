using System;

using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;


namespace Meta.Numerics.Statistics.Distributions {

    // Useful to have some symbolic results for low orders. These can be derived in several ways. One is to form the Durbin matrix symbolicly, take powers, and expand
    // the results into polynomials. The results in the region 1 < t < n/2 are all formed this way, are those for n/2 < t < 5, which could be obtained with our largest
    // 7X7 symbolic Durbin matrix. The region 1/2 < t < 1 corresponds to a 1X1 Durbin matrix of which we can take trivial powers, or 

    // For example, in the region 5/2 < t < 3, k = 3 and h = 3 - t and the Durbin matrix is:
    //      { [1-(3-t)]       1              0              0              0         }
    //      { [1-(3-t)^2]/2   1              1              0              0         }
    //  H = { [1-(3-t)^3]/3!  1/2            1              1              0         }
    //      { [1-(3-t)^4]/4!  1/3!           1/2            1              1         }
    //      { [1-2(3-t)^5]/5! [1-(3-t)^4]/4! [1-(3-t)^3]/3! [1-(3-t)^2]/2  [1-(3-t)] }
    // and P(t) = (n!/n^n) (H^n)_{3,3}.

    // n = 1
    // 1/2 < t < 1      P = 1 (2 t - 1) =  1 - 2 (1 - t) = 2 t - 1

    // n = 2
    // 1/2 < t < 1      P = 1/2 (2 t - 1)^2
    //   1 < t < 2      P = 1 - 1/2 (2 - t)^2 = 1/2 (-2 + 4 t - t^2)

    // n = 3
    // 1/2 < t < 1      P = 2/9 (2 t - 1)^3
    //   1 < t < 3/2    P = 2/9 (0 - 4 t + 7 t^2 - 2 t^3)
    // 3/2 < t < 2      P = 1 - 2/27 (3 - t)^3 - 2/9 t (2 - t)^2 = 2/9 (-9/2 + 5 t + t^2 - 2/3 t^3)
    //   2 < t < 3      P = 1 - 2/27 (3 - t)^3 = 2/9 (-9/2 + 9 t - 3 t^2 + 1/3 t^3)

    // n = 4
    // 1/2 < t < 1      P = 3/32 (2 t - 1)^4
    //   1 < t < 3/2    P = 3/32 (4 - 12 t + 7 t^2 + 4 t^3 - 2 t^4)
    // 3/2 < t < 2      P = 3/32 (4 - 21 t + 25 t^2 - 8 t + 2/3 t^4)
    //   2 < t < 3      P = 1 - 1/128 (4 - t)^4 - 1/32 t (3 - t)^3 = 3/32 (-32/3 + 37/3 t + t^2 - 5/3 t^3 + 1/4 t^4)
    //   3 < t < 4      P = 1 - 1/128 (4 - t)^4 = 3/32 (-32/3 + 64/3 t - 8 t^2 + 4/3 t^3 - 1/12 t^4)

    // n = 5
    // 1/2 < t < 1      P = 24/625 (2 t - 1)^5
    //   1 < t < 3/2    P = 24/625 (-4 + 28 t - 61 t^2 + 50 t^3 - 12 t^4 + 0 t^5)
    // 3/2 < t < 2      P = 24/625 (14 - 35 t + 25/2 t^2 + 53/3 t^3 - 10 t^4 + 4/3 t^5)
    //   2 < t < 5/2    P = 24/625 (0 - 91/3 t + 140/3 t^2 - 19 t^3 + 37/12 t^4 - 1/6 t^5)
    // 5/2 < t < 3      P = 1 - 2/3125 (5 - t)^5 - 2/625 t (4 - t)^4 - 4/625 t (t + 2) (3 - t)^3 = 24/625 (-625/24 + 87/4 t + 5 t^2 - 7/3 t^3 - 1/4 t^4 + 1/10 t^5)
    //   3 < t < 4      P = 1 - 2/3125 (5 - t)^5 - 2/625 t (4 - t)^4 = 24/625 (-625/24 + 123/4 t + 1/2 t^2 - 23/6 t^3 + 11/12 t^4 - 1/15 t^5)
    //   4 < t < 5      P = 1 - 2/3125 (5 - t)^5

    // n = 6
    // 1/2 < t < 1      P = 5/324 (2 t - 1)^6
    //   1 < t < 3/2    P = 5/324 (-4 + 4 t + 47 t^2 - 128 t^3 + 118 t^4 - 40 t^5 + 4 t^6)
    // 3/2 < t < 2      P = 5/324 (-7/4 + 58 t - 157 t^2 + 424/3 t^3 - 130/3 t^4 + 8/3 t^5 + 4/9 t^6)
    //   2 < t < 5/2    P = 5/324 (81/4 - 226/3 t + 305/6 t^2 + 103/6 t^3 - 223/12 t^4 + 14/3 t^5 - 7/18 t^6)
    // 5/2 < t < 3      P = 5/324 (81/4 - 1529/12 t + 155 t^2 - 397/6 t^3 + 59/4 t^4 - 2 t^5 + 13/90 t^6)
    //   3 < t < 4      P = 1 - 1/23328 (6 - t)^6 - 1/3888 t (5 - t)^5 - 5/7776 t (t + 2) (4 - t)^4 = 5/324 (-324/5 + 3371/60 t + 35/4 t^2 - 37/6 t^3 + 4/15 t^5 - 1/36 t^6)
    //   4 < t < 5      P = 1 - 1/23328 (6 - t)^6 - 1/3888 t (5 - t)^5
    //   5 < t < 6      P = 1 - 1/23328 (6 - t)^6

    // n = 7
    // 1/2 < t < 1      P = 720/117649 (2 t - 1)^7
    //   1 < t < 3/2    P = 720/117649 (12 - 84 t + 203 t^2 - 178 t^3 - 28 t^4 + 128 t^5 - 60 t^6 + 8 t^7)
    // 3/2 < t < 2      P = 720/117649 (-177/4 + 195 t - 491/2 t^2 + 65/3 t^3 + 416/3 t^4 - 80 t^5 + 140/9 t^6 - 8/9 t^7)
    //   2 < t < 5/2    P = 720/117649 (303/4 - 167 t + 151/3 t^2 + 563/9 t^3 - 1595/72 t^4 - 59/12 t^5 + 53/18 t^6 - 1/3 t^7)
    // 5/2 < t < 3      P = 720/117649 (303/4 - 2629/12 t + 1229/12 t^2 + 1501/18 t^3 - 5195/72 t^4 + 87/4 t^5 - 287/90 t^6 + 1/5 t^7)
    //   3 < t < 7/2    P = 720/117649 (0 - 12821/60 t + 9191/30 t^2 - 2575/18 t^3 + 2639/72 t^4 - 73/12 t^5 + 13/20 t^6 - 1/30 t^7)
    // 7/2 < t < 4      P = 720/117649 (-117649/720 + 40723/360 t + 105/4 t^2 - 29/3 t^3 - 35/24 t^4 + 9/20 t^5 + 1/36 t^6 - 1/126 t^7)
    //   4 < t < 5      P = 1 - 2/823543 (7 - t)^7 - 2/117649 t (6 - t)^6 - 6/117649 t (t + 2) (5 - t)^5
    //   5 < t < 6      P = 1 - 2/823543 (7 - t)^7 - 2/117649 t (6 - t)^6
    //   6 < t < 7      P = 1 - 2/823543 (7 - t)^7

    // These results can be checked for consistency by verifying that P is continuous at the region boundaries, and by verifying that the matrix and Q-series results agree
    // in the regions where they overlap.

    // Differentiation wrt t gives piecewise expressions for the PDF p(t). Unlike the CDF P(t), p(t) is not continuous, although it appears to get quite close as n increases.
    // (The largest discontinuity appears to be consistently at t=1.) By doing piecewise integration of p(t) * t^m, we can get raw moments, and from these we can get central moments.
    // Checking that <1> = 1 provides a sanity check, as does verifying that the moments approach their asymptotic 

    // n = 1: <t> = 3/4, <t^2> = 7/12, <t^3> = 15/32 => C2 = 1/48, C3 = 0
    // n = 2: <t> = 13/12, <t^2> = 61/48, <t^3> = 257/160 => C2 = 7/72
    // n = 3: <t> = 293/216, <t^2> = 2167/1080, <t^3> = 
    // n = 4: <t> = 813/512, <t^2> = 4231/1536, <t^3> =
    // n = 5: <t> = 134377/75000, <t^2> = 614267/175000, <t^3> = 524179/70000 => C2 = 11809828097/39375000000, C3 = 183803886921931/1476562500000000
    // n = 6: <t> = 1290643/653184, <t^2> = 1594285/373248, <t^3> = 11686771/1161216 => C2 = 156623190071/426649337856, C3 = 24188301396000379/139340260549066752
    // n = 7: <t> = 7067335/3294172, <t^2> = 298688525/59295096, <t^3> = 2039651983/158120256, => C2 = 42440671868125/97664122490256, C3 = 8138432234802217189/35746935301330176448

    public class KolmogorovExactDistribution : Distribution {

        public KolmogorovExactDistribution (int size) {
            if (size < 1) throw new ArgumentOutOfRangeException(nameof(size));

            n = size;
        }

        private readonly int n;

        /// <inheritdoc />
        public override Interval Support {
            get {
                return (Interval.FromEndpoints(0.5, n));
            }
        }

        /// <inheritdoc />
        public override double ProbabilityDensity (double t) {
            if (t <= 0.5) {
                return (0.0);
            } else if (t < n / 2.0) {
                return (DurbinMatrixPPrime(t));
            } else if (t < n) {
                return (DurbinSeriesQPrime(t));
            } else {
                return (0.0);
            }
        }

        /// <inheritdoc />
        public override double LeftProbability (double t) {
            if (t <= 0.5) {
                return (0.0);
            } else if (t < n / 2.0) {
                return (DurbinMatrixP(t));
            } else if (t < n) {
                return (1.0 - DurbinSeriesQ(t));
            } else {
                return (1.0);
            }
        }

        /// <inheritdoc />
        public override double RightProbability (double t) {
            if (t <= 0.5) {
                return (1.0);
            } else if (t < n / 2.0) {
                return (1.0 - DurbinMatrixP(t));
            } else if (t < n) {
                return (DurbinSeriesQ(t));
            } else {
                return (0.0);
            }
        }

        // Durbin generalized results of Massey to obtain a formula for the right tail that reads, in the t > n / 2 region:
        //   Q = \frac{2t}{n} \sum_{k=0}^{\floor{n-t}} \frac{n!}{k! (n-k)!} \left( \frac{t+k}{n} \right)^{k-1} \left( \frac{n - k - t}{n} \right)^{n - k}
        // This formula consists of a single term for n-1 < t < n, two terms for n-2 < t < n-1, etc. A new term appears for each integer step to the left.
        //   Q = \frac{2}{n^n} \left( n - t \right)^n + \frac{2}{n^{n-1}} t \left( n - 1 - t \right)^{n-1} +
        //       \frac{n-1}{n^{n-1}} \left( t + 2 \right) \left( n - 2 - t \right)^{n-2} + \cdots
        // See Brown & Harvey, Journal of Statistical Software 26 (2008)
        // Durbin, "Distribution Theory for Tests Based on the Sample Distribution Function", 1973, p. 12, equation 2.4.8

        // All terms in this series are positive, so it does not suffer from cancelation. A corresponding recurrsive series that holds for t < n / 2
        // suffers from significant cancelation even for small n and is thus not useful for practical purposes.

        private double DurbinSeriesQ (double t) {

            if (t >= n) return (0.0);

            // first term
            double s = MoreMath.Pow((n - t) / n, n);

            // higher terms
            int jmax = (int) Math.Floor(n - t);
            double B = 1.0; // binomial coefficient ( n j )
            for (int j = 1; j <= jmax; j++) {
                B = B * (n - (j - 1)) / j;
                double C = Math.Pow((t + j) / n, j - 1) * Math.Pow((n - j - t) / n, n - j) * t / n;
                s += B * C;
            }

            return (2.0 * s);

        }

        private double DurbinSeriesQPrime (double t) {

            if (t >= n) return (0.0);

            // first term
            double s = 2.0 * Math.Pow((n - t) / n, n - 1);

            // higher terms
            int jmax = (int) Math.Floor(n - t);
            if (jmax == (n - t)) jmax--;
            for (int j = 1; j <= jmax; j++) {
                double B = Math.Exp(AdvancedIntegerMath.LogFactorial(n - 1) - AdvancedIntegerMath.LogFactorial(j) - AdvancedIntegerMath.LogFactorial(n - j));
                double C = Math.Pow((t + j) / n, j - 1) * Math.Pow((n - j - t) / n, n - j);
                double D = 1.0 + (j - 1) * t / (t + j) - (n - j) * t / (n - j - t);
                s += -2.0 * B * C * D;
            }
            return (s);
        }

        // Durbin derived a matrix formulation, later programmed by Marsaglia, that gives P by taking matrix powers.
        // If t = n D = k - h, where k is an integer and h is a fraction, the matrix is (2k - 1) X (2k - 1) and takes the form:
        //        { (1-h)/1!    1           0           0           0        }
        //        { (1-h^2)/2!  1/1!        1           0           0        }
        //    H = { (1-h^3)/3!  1/2!        1/1!        1           0        }
        //        { (1-h^4)/4!  1/3!        1/2!        1/1!        1        }
        //        { *           (1-h^4)/4!  (1-h^3)/3!  (1-h^2)/2!  (1-h)/1! }
        // The lower-left element is special and depends on the size of h:
        //    0 < h < 1/2:  * = (1 - 2 h^m) / m!
        //    1/2 < h < 1:  * = (1 - 2 h^m + (2h-1)^m)/ m!
        // Note the factor 2 in the first case and the additional term in the second case. The left probability is then given by
        // the middle element of H to the nth power:
        //    P = \frac{n!}{n^n} \left( H^n \right)_{k,k}
        // Note that, for a given value of t, the matrix H does not depend on n; n dependence only enters as the power to which H is taken.

        // While astounding, the matrix method is clearly quite onerous. Because all the entries of H are positive, though, it is not subject to
        // the cancelation errors that makes the series solution practically useless in the region t < n/2.

        // See Durbin, "Distribution THeory for Tests Based on the Sample Distribution Function", 1973 and
        // Marsaglia, Tsand, & Wang, "Evaluating Kolmogorov's Distribution", Journal of Statistical Software 8  (2003) 

        private double DurbinMatrixP (double t) {

            int k; double h;
            DecomposeInteger(t, out k, out h);

            SquareMatrix H = GetDurbinMatrix(k, h);

            SquareMatrix Hn = H.Power(n);
            //SquareMatrix Hn = NewMatrixPower(H, n);

            double f = AdvancedIntegerMath.Factorial(n) / MoreMath.Pow(n, n);
            return (f * Hn[k - 1, k - 1]);

        }

        private static void DecomposeInteger (double t, out int k, out double h) {

            // decompose t into an integer minus a fraction
            k = (int) Math.Ceiling(t);
            h = k - t;

        }

        private static SquareMatrix GetDurbinMatrix (int k, double h) {

            // dimension of matrix
            int m = 2 * k - 1;
            SquareMatrix H = new SquareMatrix(m);

            // populate the matrix along diagonals, since they share factorial factors

            // superdiagonal is all 1s
            for (int j = 1; j < m; j++) {
                H[j - 1, j] = 1.0;
            }

            // diagonal and subdiagonals
            double Fi = 1.0; // variable for 1/i!
            double hi = h; // variable for h^i
            for (int i = 1; i < m; i++) {
                // elements and ends of diagonal
                H[i - 1, 0] = Fi * (1.0 - hi);
                H[m - 1, m - i] = H[i - 1, 0];
                // elements in between
                for (int j = i + 1; j < m; j++) {
                    H[j - 1, j - i] = Fi;
                }
                // prepare for the next recurrsion
                hi = hi * h;
                Fi = Fi / (i + 1);
            }

            // lower left element
            double e = 1.0 - 2.0 * hi;
            if (h > 0.5) e += MoreMath.Pow(2.0 * h - 1.0, m);
            H[m - 1, 0] = Fi * e;

            return (H);

        }

        private static SquareMatrix GetDurbinMatrixPrime (int k, double h) {

            // dimension of matrix
            int m = 2 * k - 1;
            SquareMatrix H = new SquareMatrix(m);

            // left column and bottom row
            double Fi = 1.0;
            double hi = 1.0;
            for (int i = 1; i < m; i++) {
                H[i - 1, 0] = Fi * hi;
                H[m - 1, m - i] = H[i - 1, 0];
                Fi /= i;
                hi *= h;
            }

            // lower left element
            double e = 2.0 * hi;
            if (h > 0.5) e -= 2.0 * MoreMath.Pow(2.0 * h - 1.0, m - 1);
            H[m - 1, 0] = Fi * e;

            return (H);

        }

        private double DurbinMatrixPPrime (double t) {

            int k; double h;
            DecomposeInteger(t, out k, out h);

            // compute derivative of H
            SquareMatrix DH = GetDurbinMatrixPrime(k, h);

            //PrintMatrix(DH);

            // compute powers of H
            SquareMatrix[] PowerH = new SquareMatrix[n];
            PowerH[1] = GetDurbinMatrix(k, h);
            for (int i = 2; i < n; i++) {
                PowerH[i] = PowerH[1] * PowerH[i - 1];
            }

            // use D(H^n) = (DH) H^(n-1) + H (DH) H^(n-2) + H^2 (DH) H^(n-3) + \cdots + H^(n-2) (DH) H + H^(n-1) (DH)
            SquareMatrix HnP = DH * PowerH[n - 1];
            for (int i = 1; i < (n - 1); i++) {
                HnP += PowerH[i] * DH * PowerH[n - 1 - i];
            }
            HnP += PowerH[n - 1] * DH;

            double f = AdvancedIntegerMath.Factorial(n) / MoreMath.Pow(n, n);
            return (f * HnP[k - 1, k - 1]);
        }

        // We have pre-computed the mean, variance, and third central moment of t = n d for n < 8
        // by symbolic integration of the piece-wise polynomial expressions

        ///<inheritdoc />
        public override double Mean {
            get {
                switch (n) {
                    case 1:
                        // 3/4 = 0.750000 * sqrt(1)
                        return (3.0 / 4.0);
                    case 2:
                        // 13/12 = 0.766932 * sqrt(2)
                        return (13.0 / 12.0);
                    case 3:
                        // 293/216 = 0.783165 * sqrt(3)
                        return (293.0 / 216.0);
                    case 4:
                        // 813/512 = 0.793945 * sqrt(4)
                        return (813.0 / 512.0);
                    case 5:
                        // 134377/75000 = 0.801270 * sqrt(5)
                        return (134377.0 / 75000.0);
                    case 6:
                        // 1290643/653184 = 0.806668 * sqrt(6)
                        return (1290643.0 / 653184.0);
                    case 7:
                        // 7067335/3294172 = 0.810887 * sqrt(7)
                        return (7067335.0 / 3294172.0);
                    default:
                        // value tends toward \sqrt{\pi/2} \log(2) = 0.868731
                        // this will perform numerical integration, which will be very slow
                        return (base.Moment(1));
                }
            }
        }

        ///<inheritdoc />
        public override double Variance {
            get {
                switch (n) {
                    case 1:
                        // 0.020833 * 1
                        return (1.0 / 48.0);
                    case 2:
                        // 0.047611 * 2
                        return (7.0 / 72.0);
                    case 3:
                        // 0.055480 * 3
                        return (38827.0 / 233280.0);
                    case 4:
                        // 0.058290 * 4
                        return (183365.0 / 786432.0);
                    case 5:
                        // 11809828097/39375000000 = 0.059986 * 5
                        return (11809828097.0 / 39375000000.0);
                    case 6:
                        // 156623190071/426649337856 = 0.061183 * 6
                        return (156623190071.0 / 426649337856.0);
                    case 7:
                        // 42440671868125/97664122490256 = 0.062080 * 7
                        return (42440671868125.0 / 97664122490256.0);
                    default:
                        // value tends toward \pi/2 (\pi/6 - (log(2))^2) = 0.0677732
                        // this will perform numerical integration twice
                        // once to find the mean and again to find the variance
                        // this will be very slow
                        return (base.MomentAboutMean(2));
                }
            }
        }

        ///<inheritdoc />
        public override double Moment (int m) {
            if (m < 0) {
                throw new ArgumentOutOfRangeException("m");
            } else if (m == 0) {
                return (1.0);
            } else if (m == 1) {
                return (Mean);
            } else {
                return (base.Moment(m));
            }
        }

        ///<inheritdoc />
        public override double MomentAboutMean (int m) {
            if (m < 0) {
                throw new ArgumentOutOfRangeException("m");
            } else if (m == 0) {
                return (1.0);
            } else if (m == 1) {
                return (0.0);
            } else if (m == 2) {
                return (Variance);
            } else {
                return (base.MomentAboutMean(m));
            }
        }

    }

}
