using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Analysis;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        // The Gammow factor is the coefficient of the leading power of rho in the expansion of the CWF near the origin
        // It sets the order of magnitude of the function near the origin. Basically F ~ C, G ~ 1/C

        private static double CoulombFactorZero (double eta) {
            // \sqrt{\frac{x}{e^x - 1}}
            return (Math.Sqrt(1.0 / MoreMath.ReducedExpMinusOne(2.0 * Math.PI * eta)));
        }

        private static double CoulombFactor (int L, double eta) {

            // From the definition (A&S 14.1.7)
            //   C_L = 2^L e^{-\pi \eta / 2} | \Gamma(L + 1 + i \eta) | / \Gamma(2L + 2)
            // It follows that
            //   C_L / C_{L - 1} = \frac{ 2 | L + i \eta | }{ (2L + 1) (2L) }
            //                   = \frac{\sqrt{ L^2 + \eta^2}}{L (2L + 1)}
            // We can use this to recurse upward from C_0 to any desired C_L

            double C = CoulombFactorZero(eta);

            for (int k = 1; k <= L; k++) {
                C *= MoreMath.Hypot(k, eta) / (k * (2 * k + 1));
            }

            return (C);

            // For small L, the few flops for each recursion steps probably
            // add up to less than the cost of a complex Gamma function evaluation.
            // For large L, it might be better from a flop standpoint to evaluate
            // directly, but then we would need to consider that the two \Gamma functions
            // grow rapidly, so we would need some way to compute their ratio directly.

        }

        // This is an implementation of NIST 33.6.1:
        //    F_{\ell}(\eta, \rho) = C_{\ell}(\eta) \sum_{k = \ell + 1}^{\infty} A_{k} \rho^k
        //    F'_{\ell}(\eta, rho) = C_{\ell}(\eta) \sum_{k = \ell + 1}^{\infty} k A_{k} \rho^{k-1}
        // with the coefficients defined by the recurrence
        //    A_{k + 1} = 1 \qquad A_{k + 2} = \frac{\eta \rho}{L + 1}
        //    (k + \ell) (k - \ell - 1) A_k = 2 \eta A_{k-1} - A_{k-2}
        // We have used a form in which we track A_{k+1} \rho^{k - (\ell + 1)} so that we are not tracking
        // two quantities, A_k decreasing and \rho^k increasing, which might overflow before their product does.

        // Each term adds factors of order \rho^2 / (L + 1) and 2 \rho \rho / (L + 1). So for numerical
        // convergence we need \rho < \sqrt{X(L + 1)} and 2 \eta \rho < X (L + 1); X ~ 16 gets convergence
        // within ~ 30 terms.

        // There are circumstances where one term can be zero or negligible but the next term is not.
        // For example L=1, \eta=-1 causes the third term to vanish, which caused a bug. So we demand
        // that the series not change for two terms before returning.

        // Demanding that the 2nd term not overwhelm the first requires \eta \rho < (L + 1).
        // Demanding that the 3rd term not overwhelm the first when \eta = 0 requires \rho^2 < 2 (2L + 3).

        private static void CoulombF_Series (int L, double eta, double rho, out double F, out double FP) {

            double termPrevious = 1.0;
            double termCurrent = eta * rho / (L + 1);
            double fPrevious = termPrevious;
            double fCurrent = fPrevious + termCurrent;
            double fPrimePrevious = (L + 1) * termPrevious;
            double fPrimeCurrent = fPrimePrevious + (L + 2) * termCurrent;

            for (int j = 2; j < Global.SeriesMax; j++) {

                double termNext = rho * (2.0 * eta * termCurrent - rho * termPrevious) / (j * (2 * L + 1 + j));
                double fNext = fCurrent + termNext;
                double fPrimeNext = fPrimeCurrent + (L + 1 + j) * termNext;

                if (fNext == fPrevious && fPrimeNext == fPrimePrevious) {
                    double C = 1.0; // CoulombFactor(L, eta);
                    double rhoToPowerL = MoreMath.Pow(rho, L);
                    F = C * rhoToPowerL * rho * fNext;
                    FP = C * rhoToPowerL * fPrimeNext;
                    return;
                }

                termPrevious = termCurrent;
                fPrevious = fCurrent;
                fPrimePrevious = fPrimeCurrent;
                termCurrent = termNext;
                fCurrent = fNext;
                fPrimeCurrent = fPrimeNext;

            }
            throw new NonconvergenceException();

        }

        // series for L=0 for both F and G
        // this has the same convergence properties as the L != 0 series for F above

        private static void Coulomb_Zero_Series (double eta, double rho, out double F, out double FP, out double G, out double GP) {

            if (rho == 0.0) {
                double C = CoulombFactorZero(eta);
                F = 0.0;
                FP = C;
                G = 1.0 / C;
                GP = Double.NegativeInfinity;
                return;
            }

            double eta_rho = eta * rho;
            double rho_2 = rho * rho;

            double u0 = 0.0;
            double u1 = rho;
            double u = u0 + u1;
            double up = u1;

            double v0 = 1.0;
            double v1 = 0.0;
            double v = v0 + v1;

            for (int n = 2; n <= Global.SeriesMax; n++) {

                double u2 = (2.0 * eta_rho * u1 - rho_2 * u0) / (n * (n - 1));
                double v2 = (2.0 * eta_rho * v1 - rho_2 * v0 - 2.0 * eta * (2 * n - 1) * u2) / (n * (n - 1));

                double u_old = u;
                u += u2;
                up += n * u2;

                double v_old = v;
                v += v2;

                if ((u == u_old) && (v == v_old)) {

                    double C =  CoulombFactorZero(eta);
                    F = C * u;

                    FP = C * up / rho;

                    double r = AdvancedComplexMath.Psi(new Complex(1.0, eta)).Re + 2.0 * AdvancedMath.EulerGamma - 1.0;
                    G = (v + 2.0 * eta * u * (Math.Log(2.0 * rho) + r)) / C;

                    GP = (FP * G - 1.0) / F;

                    return;
                }

                u0 = u1; u1 = u2; v0 = v1; v1 = v2;

            }

            throw new NonconvergenceException();

        }

        // use Steed's method to compute F and G for a given L
        // the method uses a real continued fraction (1 constraint), an imaginary continued fraction (2 constraints)
        // and the Wronskian (4 constraints) to compute the 4 quantities F, F', G, G'
        // it is reliable past the truning point, but becomes slow if used far past the turning point

        // Solution (a la Barnett)
        //   F = s \left[ (f - p)^2 / q + q \right]^{-1/2}
        //   F' = f F
        //   G = g F \qquad g = (f - p) / q
        //   G' = (p g - q) F

        private static SolutionPair Coulomb_Steed (double L, double eta, double rho) {

            // compute CF1 (F'/F)
            double f = Coulomb_CF1(L, eta, rho, out int sign);

            // compute CF2 ((G' + iF')/(G + i F))
            Complex z = Coulomb_CF2(L, eta, rho);
            double p = z.Re;
            double q = z.Im;

            // use CF1, CF2, and Wronskian (FG' - GF' = 1) to solve for F, F', G, G' 
            double g = (f - p) / q;

            double F = sign / Math.Sqrt(g * g * q + q);
            double FP = f * F;
            double G = g * F;
            double GP = (p * g - q) * F;
            SolutionPair result = new SolutionPair(F, FP, G, GP);
            return (result);

        }

        // gives F'/F and sgn(F)
        // converges rapidly for rho < turning point; slowly for rho > turning point, but still converges

        private static double Coulomb_CF1 (double L, double eta, double rho, out int sign) {

            // maximum iterations
            int nmax = 2 * Global.SeriesMax;
            double rho0 = CoulombTurningPoint(L, eta);
            if (rho > rho0) nmax += (int) Math.Floor(2.0 * (rho - rho0));

            // use Wallis method of continued fraction evaluation

            double f = (L + 1.0) / rho + eta / (L + 1.0);
            sign = 1;

            double A0 = 1.0;
            double A1 = f;
            double B0 = 0.0;
            double B1 = 1.0;

            for (int n = 1; n < nmax; n++) {

                double f_old = f;

                // compute next term   
                double k = L + n;
                double t = eta / k;
                double a = -(1.0 + t * t);
                double b = (2.0 * k + 1.0) * (1.0 / rho + t / (k + 1.0));

                // apply it
                double A2 = b * A1 + a * A0;
                double B2 = b * B1 + a * B0;

                if (B2 != 0.0) {

                    // note B1 = 1 always, A1 = f always
                    f = A2 / B2;
                    if (B2 < 0) sign = -sign;

                    // check for convergence
                    if (f == f_old) {
                        return (f);
                    }

                    // renormalize by dividing by B2 and prepare for the next cycle
                    A1 = A1 / B2;
                    A2 = f;
                    B1 = B1 / B2;
                    B2 = 1.0;

                }

                A0 = A1;
                B0 = B1;
                A1 = A2;
                B1 = B2;
            }

            throw new NonconvergenceException();

        }

        // Steed and Barnett derived the following continued fraction for
        //   p + i q = \frac{G' + iF'}{G + i F}
        //   p + i q = i ( 1 - \eta / \rho) + 
        //     \frac{i}{\rho} \frac{i\eta - \ell)(i\eta + \ell + 1)}{2(\rho - \eta + i) +}
        //       \frac{(i \eta - \ell + 1)(i \eta + \ell + 2)}{2(rho - \eta + 2 i) +} \cdots
        // They derived it from the asymptotic series, but noted that this ratio
        // has much better convergence properties than the asypmtotic series
        // for the individual components from which it is derived.

        // NR gives the impresses that this converges for \rho greater than the turning point,
        // but it's actually more subtle. I don't know what the mathematical basin of convergence
        // is, but (a) it clearly exhibits numerical convergence to the correct value even far
        // inside the turning point and (b) it's easy to find arguments for which the number of
        // iterations required for numerical convergence is quite high even inside the turning
        // point and other arguments for which the required number of iterations is quite low
        // even inside the turning point. Very approximately, it looks empirically like the number of
        // iterations scales like ~20 (\eta / \rho)^{1/2}, with only a very weak dependence on \ell.
        // I would love to understand this in terms of the continued fraction expression,
        // but I haven't been able to.

        // The turning point is still relevent, though, because inside the turning point,
        // where G >> F, it becomes dominated by G'/G and the extraction of F and F'
        // looses accuracy. 

        private static Complex Coulomb_CF2 (double L, double eta, double rho) {

            Complex a = new Complex(1.0 + L, eta);
            Complex c = new Complex(-L, eta);

            double eta_squared = eta * eta;
            double two_rho_minus_eta = 2.0 * (rho - eta);
            Complex p0 = new Complex(-(eta_squared + L * (L + 1)), eta);
            Complex q0 = new Complex(two_rho_minus_eta, 2.0);

            Complex D = 1.0 / q0;
            Complex Df = p0 * D;
            Complex f = Df;

            // Predicting the iteration limit has proved difficult.
            int nmax = 100000;
            //int nmax = Global.SeriesMax;
            //nmax += (int) Math.Min(Math.Round(40.0 * Math.Sqrt(Math.Abs(eta) / rho)), 10000.0);
            //if (eta < 0) nmax += (int) Math.Floor(-9.0 * eta);

            for (int n = 1; n < nmax; n++) {

                Complex f_old = f;

                Complex p = new Complex((n + L + 1) * (n - L) - eta_squared, eta * (2 * n + 1));
                Complex q = new Complex(two_rho_minus_eta, 2 * (n + 1));

                D = 1.0 / (q + p * D);
                Df = (q * D - 1.0) * Df;
                f += Df;

                if (f == f_old) {
                    return (ComplexMath.I * f / rho + new Complex(0.0, 1.0 - eta / rho));
                }

            }

            throw new NonconvergenceException();

            // don't use Wallis algorithm for this continued fraction! it appears that sometimes noise in what
            // should be the last few iterations prevents convergence; Steed's method appears to do better
            // an example is L=0, eta=0.1, rho=14.0

        }

        // asymptotic region
        private static SolutionPair Coulomb_Asymptotic (double L, double eta, double rho) {

            // Abromowitz & Stegun 14.5 describes this asymptotic expansion

            
            RangeReduction.ReduceByPiHalves(rho, out long t00, out double t01);

            RangeReduction.ReduceByPiHalves(AdvancedComplexMath.LogGamma(new Complex(L + 1.0, eta)).Im - eta * Math.Log(2.0 * rho), out long t10, out double t11);

            RangeReduction.ReduceByOnes(t01 + t11, out long t20, out double t1);

            long t0 = t00 + t10 + t20 - (long) L;

            double s = RangeReduction.Sin(t0, t1);
            double c = RangeReduction.Cos(t0, t1);
            
            // determine the weights of sin and cos
            double f0 = 1.0;
            double g0 = 0.0;
            double fp0 = 0.0;
            double gp0 = 1.0 - eta / rho;

            double f = f0;
            double g = g0;
            double fp = fp0;
            double gp = gp0;

            for (int k = 0; k < Global.SeriesMax; k++) {

                // Remember the old values for comparison
                double f_old = f;
                double g_old = g;
                double fp_old = fp;
                //double gp_old = gp;

                // Compute the next corrections
                double q = 2 * (k + 1) * rho;
                double a = (2 * k + 1) * eta / q;
                double b = ((L * (L + 1) - k * (k + 1)) + eta * eta) / q;

                double f1 = a * f0 - b * g0;
                double g1 = a * g0 + b * f0;
                double fp1 = a * fp0 - b * gp0 - f1 / rho;
                double gp1 = a * gp0 + b * fp0 - g1 / rho;

                // Apply them
                f += f1;
                g += g1;
                fp += fp1;
                gp += gp1;

                // Test for convergence
                if ((f == f_old) && (g == g_old) && (fp == fp_old)) {
                    return (new SolutionPair(
                        g * c + f * s,
                        gp * c + fp * s,
                        f * c - g * s,
                        fp * c - gp * s)
                    );
                }

                // Prepare for the next iteration
                f0 = f1;
                g0 = g1;
                fp0 = fp1;
                gp0 = gp1;

            }

            throw new NonconvergenceException();

        }

        private static double Coulomb_Asymptotic_Limit (double L, double eta) {
            return (32.0 + 0.5 * (L * L + eta * eta));
        }

        private static double Coulomb_Series_Limit (int L, double eta) {

            // First correction term is \rho \eta / (L + 1). It should not overwhelm 1.
            // If \eta < 0, it's negative, so avoid catastrophic cancelation.
            // If \eta > 0, it can be larger without cancelation.
            double a = (L + 1) / Math.Abs(eta);
            if (eta < 0.0) {
                // I would prefer to make this 0.5, but to keep series more useful accept 0.875 for now.
                a *= 0.875;
            } else {
                // I would prefer to make this 1.5, but to keep series more useful accept 3.0 for now. 
                a *= 3.0;
            }

            // Second correction term, which is first is \eta = 0, is \rho^2 / 2 (2 L + 3).
            double b = Math.Sqrt(2 * (2 * L + 3));

            return (Math.Min(a, b));

        }

        // for rho < turning point, CWF are exponential; for rho > turning point, CWF are oscillatory
        // we use this in several branching calculations

        private static double CoulombTurningPoint (double L, double eta) {

            double p = L * (L + 1);
            double r = Math.Sqrt(p + eta * eta);

            if (eta >= 0.0) {
                return (r + eta);
            } else {
                return (p / (r - eta));
            }

        }

        /// <summary>
        /// Computes the regular and irregular Coulomb wave functions and their derivatives.
        /// </summary>
        /// <param name="L">The angular momentum number, which must be non-negative.</param>
        /// <param name="eta">The charge parameter, which can be positive or negative.</param>
        /// <param name="rho">The radial distance parameter, which must be non-negative.</param>
        /// <returns>The values of F, F', G, and G' for the given parameters.</returns>
        /// <remarks>
        /// <para>The Coulomb wave functions are the radial wave functions of a non-relativistic particle in a Coulomb
        /// potential.</para>
        /// <para>They satisfy the differential equation:</para>
        /// <img src="../images/CoulombODE.png" />
        /// <para>A repulsive potential is represented by &#x3B7; &gt; 0, an attractive potential by &#x3B7; &lt; 0.</para>
        /// <para>F is oscillatory in the region beyond the classical turning point. In the quantum tunneling region inside
        /// the classical turning point, F is exponentially suppressed and vanishes at the origin, while G grows exponentially and
        /// diverges at the origin.</para>
        /// <para>Many numerical libraries compute Coulomb wave functions in the quantum tunneling region using a WKB approximation,
        /// which accurately determine only the first few decimal digits; our library computes Coulomb wave functions even in this
        /// computationally difficult region to nearly full precision -- all but the last 4-5 decimal digits can be trusted.</para>
        /// <para>The irregular Coulomb wave functions G<sub>L</sub>(&#x3B7;,&#x3C1;) are the complementary independent solutions
        /// of the same differential equation.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="L"/> or <paramref name="rho"/> is negative.</exception>
        /// <seealso cref="CoulombF"/>
        /// <seealso cref="CoulombG"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Coulomb_wave_function" />
        /// <seealso href="http://mathworld.wolfram.com/CoulombWaveFunction.html" />
        public static SolutionPair Coulomb (int L, double eta, double rho) {

            if (L < 0) throw new ArgumentOutOfRangeException(nameof(L));
            if (rho < 0) throw new ArgumentOutOfRangeException(nameof(rho));

            if (rho == 0.0) {
                if (L == 0) {
                    double C = CoulombFactor(L, eta);
                    return (new SolutionPair(0.0, C, 1.0 / C, Double.NegativeInfinity));
                } else {
                    return (new SolutionPair(0.0, 0.0, Double.PositiveInfinity, Double.NegativeInfinity));
                }
            } else if (rho <= Coulomb_Series_Limit(0, eta)) {
                // Below the safe series radius for L=0, compute using the series.
                Coulomb_Zero_Series(eta, rho, out double F, out double FP, out double G, out double GP);

                // For higher L, recurse G upward, but compute F via the direct series.
                // G is safe to compute via recursion and F is not because G is increasing
                // rapidly and F is decreasing rapidly with increasing L. Since the series
                // for F was good for L = 0, it's certainly good for L > 0.
                if (L > 0) {
                    CoulombF_Series(L, eta, rho, out F, out FP);
                    double C = CoulombFactor(L, eta);
                    F *= C;
                    FP *= C;
                    Coulomb_Recurse_Upward(0, L, eta, rho, ref G, ref GP);
                }
                return (new SolutionPair(F, FP, G, GP));
            } else if (rho >= Coulomb_Asymptotic_Limit(L, eta)) {
                return (Coulomb_Asymptotic(L, eta, rho));
            } else if (rho >= CoulombTurningPoint(L, eta)) {
                return (Coulomb_Steed(L, eta, rho));
            } else {

                // This code is copied from CoulombG; factor it into a seperate method.
                double G, GP;
                if (rho >= Coulomb_Asymptotic_Limit(0, eta)) {
                    SolutionPair s = Coulomb_Asymptotic(0, eta, rho);
                    G = s.SecondSolutionValue;
                    GP = s.SecondSolutionDerivative;
                } else {
                    double rho1 = CoulombTurningPoint(0, eta);
                    if (rho >= rho1) {
                        SolutionPair s = Coulomb_Steed(0, eta, rho);
                        G = s.SecondSolutionValue;
                        GP = s.SecondSolutionDerivative;
                    } else {
                        SolutionPair s = Coulomb_Steed(0, eta, rho1);
                        G = s.SecondSolutionValue;
                        GP = s.SecondSolutionDerivative;

                        OdeResult r = FunctionMath.IntegrateConservativeOde(
                            (double x, double y) => (2.0 * eta / x - 1.0) * y,
                            rho1, G, GP, rho,
                            new OdeSettings() {
                                RelativePrecision = 2.5E-13,
                                AbsolutePrecision = 0.0,
                                EvaluationBudget = 25000
                            }
                        );

                        G = r.Y;
                        GP = r.YPrime;
                    }
                }

                Coulomb_Recurse_Upward(0, L, eta, rho, ref G, ref GP);

                // We can determine F and FP via CF1 and Wronskian, but if we are within
                // series limit, it's likely a little faster and more accurate.
                double F, FP;
                if (rho < Coulomb_Series_Limit(L, eta)) {
                    CoulombF_Series(L, eta, rho, out F, out FP);
                    double C = CoulombFactor(L, eta);
                    F *= C;
                    FP *= C;
                } else {
                    double f = Coulomb_CF1(L, eta, rho, out int _);
                    F = 1.0 / (f * G - GP);
                    FP = f * F;
                }
             
                return (new SolutionPair(F, FP, G, GP));

            }

        }

        /// <summary>
        /// Computes the regular Coulomb wave function.
        /// </summary>
        /// <param name="L">The angular momentum number, which must be non-negative.</param>
        /// <param name="eta">The charge parameter, which can be positive or negative.</param>
        /// <param name="rho">The radial distance parameter, which must be non-negative.</param>
        /// <returns>The value of F<sub>L</sub>(&#x3B7;,&#x3C1;).</returns>
        /// <remarks>
        /// <para>For information on the Coulomb wave functions, see the remarks on <see cref="Coulomb" />.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="L"/> or <paramref name="rho"/> is negative.</exception>
        /// <seealso cref="Coulomb"/>
        /// <seealso cref="CoulombG"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Coulomb_wave_function" />
        /// <seealso href="http://mathworld.wolfram.com/CoulombWaveFunction.html" />
        public static double CoulombF (int L, double eta, double rho) {

            if (L < 0) throw new ArgumentOutOfRangeException(nameof(L));
            if (rho < 0) throw new ArgumentOutOfRangeException(nameof(rho));

            double rho0 = Coulomb_Series_Limit(L, eta);
            if (rho <= rho0) {
                //if ((rho < 4.0 + 2.0 * Math.Sqrt(L)) && (Math.Abs(rho * eta) < 8.0  + 4.0 * L)) {
                // If rho and rho * eta are small enough, use the series expansion at the origin.
                CoulombF_Series(L, eta, rho, out double F, out double FP);
                double C = CoulombFactor(L, eta);
                return (C * F);
            } else if (rho >= Coulomb_Asymptotic_Limit(L, eta)) {
                // If rho is large enough, use the asymptotic expansion.
                SolutionPair s = Coulomb_Asymptotic(L, eta, rho);
                return (s.FirstSolutionValue);
            } else if (rho >= CoulombTurningPoint(L, eta)) {
                // Beyond the turning point, use the Barnett/Steed method.
                SolutionPair result = Coulomb_Steed(L, eta, rho);
                return (result.FirstSolutionValue);
            } else {
                // We are below the transition point but above the series limit.
                // Choose the approach with the least cost.
                ISolutionStrategy<double> above = new CoulombFFromSeriesAbove(L, eta, rho);
                ISolutionStrategy<double> left = new CoulombFFromOutwardIntegration(L, eta, rho);

                if (left.Cost < above.Cost) {
                    return left.Evaluate();
                } else {
                    return above.Evaluate();
                }

            }

        }

        private struct CoulombFFromSeriesAbove : ISolutionStrategy<double> {

            public CoulombFFromSeriesAbove (int L, double eta, double rho) {

                double L_above_a = rho * Math.Abs(eta) / (eta < 0.0 ? 0.875 : 2.0);
                double L_above_b = 0.5 * (0.5 * rho * rho - 3.0);
                L0 = (int) Math.Ceiling(Math.Min(Math.Max(L_above_a, L_above_b), Int32.MaxValue));
                Debug.Assert(L < L0);

                this.L = L;
                this.eta = eta;
                this.rho = rho;

                this.C0 = CoulombFactor(L0, eta);
            }

            private int L0, L;

            private double eta, rho;

            private double C0;

            public double Cost {
                get {
                    if (C0 == 0.0) return (Double.MaxValue);
                    // Assume series costs 100, each recursion step costs 10.
                    return (L0 == Int32.MaxValue ? Double.MaxValue : 100.0 + 10.0 * (L0 - L));
                }
            }

            public double Evaluate () {
                CoulombF_Series(L0, eta, rho, out double F, out double FP);
                double C = CoulombFactor(L0, eta);
                F *= C;
                FP *= C;
                Coulomb_Recurse_Downward(L0, L, eta, rho, ref F, ref FP);
                return (F);
            }

        }

        private class CoulombFFromOutwardIntegration : ISolutionStrategy<double> {

            public CoulombFFromOutwardIntegration (int L, double eta, double rho) {

                rho0 = Coulomb_Series_Limit(L, eta);
                Debug.Assert(rho0 < rho);

                this.L = L;
                this.eta = eta;
                this.rho = rho;
            }

            private int L;

            private double eta;

            private double rho0, rho;

            public double Cost {
                get {
                    // Assume series costs 100, each integration step costs 1000.
                    return (100.0 + 1000.0 * (rho - rho0));
                }
            }

            public double Evaluate () {

                CoulombF_Series(L, eta, rho0, out double F, out double FP);

                if ((F == 0.0) && (FP == 0.0)) return (0.0);

                OdeResult r = FunctionMath.IntegrateConservativeOde(
                    (double x, double y) => ((L * (L + 1) / x + 2.0 * eta) / x - 1.0) * y,
                    rho0, F, FP, rho,
                    new OdeSettings() {
                        RelativePrecision = 2.5E-13,
                        AbsolutePrecision = 0.0,
                        EvaluationBudget = 25000
                    }
                );

                double C = CoulombFactor(L, eta);

                return (C * r.Y);

            }

        }

        /// <summary>
        /// Computes the irregular Coulomb wave function.
        /// </summary>
        /// <param name="L">The angular momentum number, which must be non-negative.</param>
        /// <param name="eta">The charge parameter, which can be positive or negative.</param>
        /// <param name="rho">The radial distance parameter, which must be non-negative.</param>
        /// <returns>The value of G<sub>L</sub>(&#x3B7;,&#x3C1;).</returns>
        /// <remarks>
        /// <para>For information on the Coulomb wave functions, see the remarks on <see cref="Coulomb" />.</para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="L"/> or <paramref name="rho"/> is negative.</exception>
        /// <seealso cref="Coulomb"/>
        /// <seealso cref="CoulombF"/>
        /// <seealso href="http://en.wikipedia.org/wiki/Coulomb_wave_function" />
        /// <seealso href="http://mathworld.wolfram.com/CoulombWaveFunction.html" />
        public static double CoulombG (int L, double eta, double rho) {

            if (L < 0) throw new ArgumentOutOfRangeException(nameof(L));
            if (rho < 0) throw new ArgumentOutOfRangeException(nameof(rho));

            if (rho <= Coulomb_Series_Limit(0, eta)) {
                // For small enough rho, use the power series for L=0, then recurse upward to desired L.
                Coulomb_Zero_Series(eta, rho, out double F, out double FP, out double G, out double GP);
                Coulomb_Recurse_Upward(0, L, eta, rho, ref G, ref GP);
                return (G);
            } else if (rho >= Coulomb_Asymptotic_Limit(L, eta)) {
                // For large enough rho, use the asymptotic series.
                SolutionPair s = Coulomb_Asymptotic(L, eta, rho);
                return (s.SecondSolutionValue);
            } else if (rho >= CoulombTurningPoint(L, eta)) {
                // Beyond the turning point, use Steed's method.
                SolutionPair result = Coulomb_Steed(L, eta, rho);
                return (result.SecondSolutionValue);
            } else {
                // We will start at L=0 (which has a smaller turning point radius) and recurse up to the desired L.
                // This is okay because G increases with increasing L. We already know that we are beyond the L=0
                // series limit; otherwise we would have taken the branch for it above. We might still be
                // in the asymptotic region, in the Steed region, below the L=0 turning point (if \eta > 0).
                double G, GP;
                if (rho >= Coulomb_Asymptotic_Limit(0, eta)) {
                    SolutionPair s = Coulomb_Asymptotic(0, eta, rho);
                    G = s.SecondSolutionValue;
                    GP = s.SecondSolutionDerivative;
                } else {
                    double rho1 = CoulombTurningPoint(0, eta);
                    if (rho >= rho1) {
                        SolutionPair s = Coulomb_Steed(0, eta, rho);
                        G = s.SecondSolutionValue;
                        GP = s.SecondSolutionDerivative;
                    } else {
                        SolutionPair s = Coulomb_Steed(0, eta, rho1);
                        G = s.SecondSolutionValue;
                        GP = s.SecondSolutionDerivative;

                        OdeResult r = FunctionMath.IntegrateConservativeOde(
                            (double x, double y) => (2.0 * eta / x - 1.0) * y,
                            rho1, G, GP, rho,
                            new OdeSettings() {
                                RelativePrecision = 2.5E-13,
                                AbsolutePrecision = 0.0,
                                EvaluationBudget = 25000
                            }
                        );

                        G = r.Y;
                        GP = r.YPrime;
                    }
                }

                Coulomb_Recurse_Upward(0, L, eta, rho, ref G, ref GP);

                return (G);
            }

        }

        private static void Coulomb_Recurse_Upward (int L1, int L2, double eta, double rho, ref double U, ref double UP) {

            Debug.Assert(L2 >= L1);

            for (int K = L1 + 1; K <= L2; K++) {
                
                // compute some factors
                double S = Math.Sqrt(K * K + eta * eta);
                double T = K * K / rho + eta;

                // compute next higher function and derivative
                double U2 = (T * U - K * UP) / S;
                double UP2 = (S * U - T * U2) / K;

                // prepare for next iteration
                U = U2;
                UP = UP2;

                // This is a hacky fix to deal with the fact that when U = Double.PositiveInfinity and
                // U2 = Double.NegativeInfinity, UP2 = Double.PositiveInfinity - Double.PositiveInfinity = Double.NaN
                if (Double.IsPositiveInfinity(U)) {
                    UP = Double.NegativeInfinity;
                    break;
                }

            }

        }

        // The following recursion formula follow from A&S 14.2:
        //   \sqrt{L^2 + \eta^2} u_{L-1} = L u_L' + (L^2 / \rho + \eta) u_L
        //   L u_{L-1}' = (L^2  / \rho + \eta) u_{L-1} - \sqrt{L^2 + \eta^2} u_L
        // We use them for downward recursion.

        private static void Coulomb_Recurse_Downward (int L1, int L2, double eta, double rho, ref double U, ref double UP) {
            Debug.Assert(L2 <= L1);
            for (int K = L1; K > L2; K--) {
                double S = MoreMath.Hypot(K, eta);
                double T = K * K / rho + eta;

                double V = (T * U + K * UP) / S;
                UP = (T * V - S * U) / K;
                U = V;
            }

        }

    }

    internal interface ISolutionStrategy<T> {

        double Cost { get; }

        T Evaluate ();

    }

}