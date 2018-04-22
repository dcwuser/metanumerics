using System;
using System.Collections.Generic;
using System.Diagnostics;

using Meta.Numerics.Analysis;

namespace Meta.Numerics.Functions {

    public static partial class AdvancedMath {

        // The Gammow factor is the coefficient of the leading power of rho in the expansion of the CWF near the origin
        // It sets the order of magnitude of the function near the origin. Basically F ~ C, G ~ 1/C

        private static double CoulombFactorZero (double eta) {

            double x = Global.TwoPI * eta;

            if (Math.Abs(x) < 0.25) {
                // For small arguments, use the Taylor expansion of x/(e^x - 1) = \sum_k B_k / k! x^k,
                // which is the generating function of the Bernoulli numbers. Note B_1 = -1/2 is the only odd Bernoulli number.
                double xx = 1.0;
                double s = 1.0 - x / 2.0;
                for (int k = 1; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                    double s_old = s;
                    xx *= (x * x) / (2 * k) / (2 * k - 1);
                    s += AdvancedIntegerMath.Bernoulli[k] * xx;
                    if (s == s_old) return (Math.Sqrt(s));
                }
                throw new NonconvergenceException();
            } else {
                return (Math.Sqrt(x / (Math.Exp(x) - 1.0)));
            }

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
                C *= MoreMath.Hypot(k, eta) / k / (2 * k + 1);
            }

            return (C);

            // For small L, the few flops for each recursion steps probably
            // add up to less than the cost of a Gamma function evaluation.
            // For large L, it might be better from a flop standpoint to evaluate
            // directly, but then we would need to consider that the two \Gamma functions
            // grow rapidly, so we would need some way to compute their ratio directly.

        }

        /*
        private static double CoulombPhaseShiftZero (double eta) {

            if (eta < 0) return (CoulombPhaseShiftZero(-eta));

            if (eta < 0.125) {
                double t = -AdvancedMath.EulerGamma * eta;
                double meta2 = - eta * eta;
                for (int k = 0; k < Zeta.Length; k++) {
                    double t_old = t;
                    eta *= meta2;
                    double dt = Zeta[k] * eta / (2 * k + 3);
                    t -= dt;
                    if (t == t_old) return (t);
                }
                throw new NonconvergenceException();
            } else if (eta < 8.0) {
                return (AdvancedComplexMath.LogGamma(new Complex(1.0, eta)).Im);
            } else {
                double t = (Math.Log(eta) - 1.0) * eta;
                t = AdvancedMath.Reduce(t, 0.0);

                double meta2 = -eta * eta;

                for (int k = 1; k < AdvancedIntegerMath.Bernoulli.Length; k++) {
                    double t_old = t;
                    double dt = AdvancedIntegerMath.Bernoulli[k] / (2 * k) / (2 * k - 1) / eta;
                    t -= dt;
                    if (t == t_old) {
                        return (t + Math.PI / 4.0);
                    }
                    eta *= meta2;
                }
                throw new NonconvergenceException();
            }

        }

        private static readonly double[] Zeta = new double [] {
            1.20205690315959428540, // Zeta(3)
            1.03692775514336992633, // Zeta(5)
            1.00834927738192282684, // Zeta(7)
            1.00200839282608221442, // Zeta(9)
            1.00049418860411946456, // Zeta(11)
            1.00012271334757848915, // Zeta(13)
            1.00003058823630702049, // Zeta(15)
            1.0000076371976378998,  // Zeta(17)
            1.0000019082127165539   // Zeta(19)
        };

        private static double CoulombPhaseShift (int L, double eta) {

            double s = CoulombPhaseShiftZero(eta);
            for (int k = 1; k <= L; k++) {
                s += Math.Atan(eta / k);
            }
            return (s);

        }
        */

        // each new term introduces factors of rho^2 / (L+1) and 2 eta rho / (L+1), so for this to converge we need
        // rho < sqrt(X) (1 + sqrt(L)) and 2 eta rho < X (1 + L); X ~ 16 gets convergence within 30 terms

        private static void CoulombF_Series (int L, double eta, double rho, out double F, out double FP) {

            double eta_rho = eta * rho;
            double rho_2 = rho * rho;

            double u0 = 1.0;
            double u1 = eta_rho / (L + 1);
            double u = u0 + u1;
            double v = (L + 1) * u0 + (L + 2) * u1;

            for (int k = 2; k < Global.SeriesMax; k++) {

                double u2 = (2.0 * eta_rho * u1 - rho_2 * u0) / k / (2 * L + k + 1);
                double v2 = (L + 1 + k) * u2;

                double u_old = u;
                u += u2;
                v += v2;

                if ((k % 2 == 0) && (u == u_old)) {
                    double C = CoulombFactor(L, eta);
                    F = C * Math.Pow(rho, L + 1) * u;
                    FP = C * Math.Pow(rho, L) * v;
                    return;
                }

                u0 = u1; u1 = u2;

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

                double u2 = (2.0 * eta_rho * u1 - rho_2 * u0) / n / (n - 1);
                double v2 = (2.0 * eta_rho * v1 - rho_2 * v0 - 2.0 * eta * (2 * n - 1) * u2) / n / (n - 1);

                double u_old = u; u += u2; up += n * u2;
                double v_old = v; v += v2;

                if ((u == u_old) && (v == v_old)) {

                    double C =  CoulombFactorZero(eta);
                    F = C * u;

                    FP = C * up / rho;

                    double r = AdvancedComplexMath.Psi(new Complex(1.0, eta)).Re + 2.0 * AdvancedMath.EulerGamma - 1;
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

        private static SolutionPair Coulomb_Steed (double L, double eta, double rho) {

            // compute CF1 (F'/F)
            int sign;
            double f = Coulomb_CF1(L, eta, rho, out sign);

            // compute CF2 ((G' + iF')/(G + i F))
            Complex z = Coulomb_CF2(L, eta, rho);
            double p = z.Re;
            double q = z.Im;

            // use CF1, CF2, and Wronskian (FG' - GF' = 1) to solve for F, F', G, G' 
            double g = (f - p) / q;

            SolutionPair result = new SolutionPair();
            result.FirstSolutionValue = sign / Math.Sqrt(g * g * q + q);
            result.FirstSolutionDerivative = f * result.FirstSolutionValue;
            result.SecondSolutionValue = g * result.FirstSolutionValue;
            result.SecondSolutionDerivative = (p * g - q) * result.FirstSolutionValue;
            return (result);

        }

        // gives F'/F and sgn(F)
        // converges rapidly for rho < turning point; slowly for rho > turning point, but still converges

        private static double Coulomb_CF1 (double L, double eta, double rho, out int sign) {

            // maximum iterations
            int nmax = Global.SeriesMax;
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

        // computes (G' + iF')/(G + i F)
        // converges quickly for rho > turning point; does not converge at all below it

        private static Complex Coulomb_CF2 (double L, double eta, double rho) {

            Complex a = new Complex(1.0 + L, eta);
            Complex c = new Complex(-L, eta);

            Complex D = 1.0 / (new Complex(2.0 * (rho - eta), 2.0));
            Complex Df = a * c * D;
            Complex f = Df;

            int nmax = Global.SeriesMax;
            if (eta < 0) nmax += (int) Math.Floor(-8.0 * eta);

            for (int n = 1; n < nmax; n++) {

                Complex f_old = f;

                Complex p = (a + n) * (c + n);
                Complex q = new Complex(2.0 * (rho - eta), 2.0 * (n + 1));

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

            
            long t00; double t01;
            RangeReduction.ReduceByPiHalves(rho, out t00, out t01);

            long t10; double t11;
            RangeReduction.ReduceByPiHalves(AdvancedComplexMath.LogGamma(new Complex(L + 1.0, eta)).Im - eta * Math.Log(2.0 * rho), out t10, out t11);

            long t20; double t1;
            RangeReduction.ReduceByOnes(t01 + t11, out t20, out t1);

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

        // for rho < turning point, CWF are exponential; for rho > turning point, CWF are oscillatory
        // we use this in several branching calculations

        private static double CoulombTurningPoint (double L, double eta) {

            double p = L * (L + 1);
            double q = Math.Sqrt(p + eta * eta);

            if (eta >= 0.0) {
                return (q + eta);
            } else {
                return (p / (q - eta));
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
            } else if ((rho < 4.0) && Math.Abs(rho * eta) < 8.0) {

                // Below the safe series radius for L=0, compute using the series
                double F, FP, G, GP;
                Coulomb_Zero_Series(eta, rho, out F, out FP, out G, out GP);

                // For higher L, recurse G upward, but compute F via the direct series.
                // G is safe to compute via recursion and F is not because G is increasing
                // rapidly and F is decreasing rapidly with increasing L.
                if (L > 0) {
                    CoulombF_Series(L, eta, rho, out F, out FP);
                    Coulomb_Recurse_Upward(0, L, eta, rho, ref G, ref GP);
                }
                return (new SolutionPair(F, FP, G, GP));
            } else if (rho > 32.0 + (L * L + eta * eta) / 2.0) {
                return (Coulomb_Asymptotic(L, eta, rho));
            } else {
                double rho0 = CoulombTurningPoint(L, eta);
                if (rho > rho0) {
                    return (Coulomb_Steed(L, eta, rho));
                } else {

                    // First F
                    double F, FP;

                    double rho1 = Math.Min(
                        4.0 + 2.0 * Math.Sqrt(L),
                        (8.0 + 4.0 * L) / Math.Abs(eta)
                    );

                    if (rho < rho1) {
                        CoulombF_Series(L, eta, rho, out F, out FP);
                    } else {
                        CoulombF_Series(L, eta, rho1, out F, out FP);

                        OdeResult r = FunctionMath.IntegrateConservativeOde(
                            (double x, double y) => ((L * (L + 1) / x + 2.0 * eta ) / x - 1.0) * y,
                            rho1, F, FP, rho,
                            new OdeSettings() {
                                RelativePrecision = 2.5E-13,
                                AbsolutePrecision = 0.0,
                                EvaluationBudget = 8192 * 2
                            }
                        );

                        F = r.Y;
                        FP = r.YPrime;
                    }
                     
                    // Then G
                    double G, GP;

                    // For L = 0 the transition region is smaller, so we will determine G 
                    // at L = 0, where non-integration methods apply over a larger region,
                    // and then recurse upward.

                    double rho2 = CoulombTurningPoint(0, eta);

                    if (rho > 32.0 + eta * eta / 2.0) {
                        SolutionPair s = Coulomb_Asymptotic(0, eta, rho);
                        G = s.SecondSolutionValue;
                        GP = s.SecondSolutionDerivative;
                    } else if (rho > rho2) {
                        SolutionPair s = Coulomb_Steed(0, eta, rho);
                        G = s.SecondSolutionValue;
                        GP = s.SecondSolutionDerivative;
                    } else {
                        SolutionPair s = Coulomb_Steed(0, eta, rho2);
                        G = s.SecondSolutionValue;
                        GP = s.SecondSolutionDerivative;
                        
                        // Integrate inward from turning point.
                        // G increases and F decreases in this direction, so this is stable.
                        OdeResult r = FunctionMath.IntegrateConservativeOde(
                            (double x, double y) => (2.0 * eta / x - 1.0) * y,
                            rho2, G, GP, rho,
                            new OdeSettings() {
                                RelativePrecision = 2.5E-13,
                                AbsolutePrecision = 0.0,
                                EvaluationBudget = 8192 * 2
                            }
                        );

                        G = r.Y;
                        GP = r.YPrime;
                    }

                    Coulomb_Recurse_Upward(0, L, eta, rho, ref G, ref GP);

                    return (new SolutionPair(F, FP, G, GP));
                }
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

            if ((rho < 4.0 + 2.0 * Math.Sqrt(L)) && (Math.Abs(rho * eta) < 8.0  + 4.0 * L)) {
                // if rho and rho * eta are small enough, use the series expansion at the origin
                double F, FP;
                CoulombF_Series(L, eta, rho, out F, out FP);
                return (F);
            } else if (rho > 32.0 + (L * L + eta * eta) / 2.0) {
                // if rho is large enough, use the asymptotic expansion
                SolutionPair s = Coulomb_Asymptotic(L, eta, rho);
                return (s.FirstSolutionValue);
                //double F, G;
                //Coulomb_Asymptotic(L, eta, rho, out F, out G);
                //return (F);
            } else {
                // transition region
                if (rho >= CoulombTurningPoint(L, eta)) {
                    // beyond the turning point, use Steed's method
                    SolutionPair result = Coulomb_Steed(L, eta, rho);
                    return (result.FirstSolutionValue);
                } else {
                    // inside the turning point, integrate out from the series limit
                    return (CoulombF_Integrate(L , eta, rho));
                }
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

            if ((rho < 4.0) && Math.Abs(rho * eta) < 8.0) {
                // For small enough rho, use the power series for L=0, then recurse upward to desired L.
                double F, FP, G, GP;
                Coulomb_Zero_Series(eta, rho, out F, out FP, out G, out GP);
                Coulomb_Recurse_Upward(0, L, eta, rho, ref G, ref GP);
                return (G);
            } else if (rho > 32.0 + (L * L + eta * eta) / 2.0) {
                // For large enough rho, use the asymptotic series.
                SolutionPair s = Coulomb_Asymptotic(L, eta, rho);
                return (s.SecondSolutionValue);
            } else {
                // Transition region
                if (rho >= CoulombTurningPoint(L, eta)) {
                    // Beyond the turning point, use Steed's method.
                    SolutionPair result = Coulomb_Steed(L, eta, rho);
                    return (result.SecondSolutionValue);
                } else {
                    
                    // we will start at L=0 (which has a smaller turning point radius) and recurse up to the desired L
                    // this is okay because G increases with increasing L

                    double G, GP;

                    double rho0 = 2.0 * eta;
                    if (rho < rho0) {

                        // if inside the turning point even for L=0, start at the turning point and integrate in
                        // this is okay because G increases with decreasing rho

                        // use Steed's method at the turning point
                        // for large enough eta, we could use the turning point expansion at L=0, but it contributes
                        // a lot of code for little overall performance increase so we have chosen not to
                        SolutionPair result = Coulomb_Steed(0, eta, 2.0 * eta);
                        G = result.SecondSolutionValue;
                        GP = result.SecondSolutionDerivative;

                        OdeResult r = FunctionMath.IntegrateConservativeOde(
                            (double x, double y) => (2.0 * eta / x - 1.0) * y,
                            rho0, G, GP, rho,
                            new OdeSettings() {
                                RelativePrecision = 2.5E-13,
                                AbsolutePrecision = 0.0,
                                EvaluationBudget = 8192 * 2
                            }
                        );

                        G = r.Y;
                        GP = r.YPrime;

                    } else {

                        // if beyond the turning point for L=0, just use Steeds method

                        SolutionPair result = Coulomb_Steed(0, eta, rho);
                        G = result.SecondSolutionValue;
                        GP = result.SecondSolutionDerivative;

                    }


                    // Recurse up to the desired L.
                    Coulomb_Recurse_Upward(0, L, eta, rho, ref G, ref GP);

                    return (G);
                    
                }
            }

        }

        // Abromowitz & Rabinowitz asymptotic expansion for L = 0 at rho = 2 eta; coefficients from Isacsson
        // this is accurate at double precision down to about eta ~ 15

        /*
        private static BesselResult Coulomb_Zero_Turning_Expansion (double eta) {

            double beta = Math.Pow(2.0 * eta / 3.0, 1.0 / 3.0);
            double beta_squared = beta * beta;
            double[] b = new double[12];
            b[0] = 1.0;
            for (int i = 1; i < b.Length; i++) {
                b[i] = b[i - 1] * beta_squared;
            }
            double sqrtbeta = Math.Sqrt(beta);

            double F = AbromowitzA[0] * sqrtbeta * ( 1.0 - AbromowitzA[1] / b[2] - AbromowitzA[2] / b[3]
                - AbromowitzA[3] / b[5] - AbromowitzA[4] / b[6] - AbromowitzA[5] / b[8] - AbromowitzA[6] / b[9]
                - AbromowitzA[7] / b[11] );

            double G = AbromowitzA[0] * sqrtbeta * Math.Sqrt(3.0) * ( 1.0 + AbromowitzA[1] / b[2] - AbromowitzA[2] / b[3]
                + AbromowitzA[3] / b[5] - AbromowitzA[4] / b[6] + AbromowitzA[5] / b[8] - AbromowitzA[6] /b[9]
                + AbromowitzA[7] / b[11] );

            double FP = AbromowitzB[0] / sqrtbeta * ( 1.0 + AbromowitzB[1] / b[1] + AbromowitzB[2] / b[3]
                + AbromowitzB[3] / b[4] + AbromowitzB[4] / b[6] + AbromowitzB[5] / b[7] + AbromowitzB[6] / b[9]
                + AbromowitzB[7] / b[10] );

            double GP = AbromowitzB[0] / sqrtbeta * Math.Sqrt(3.0) * ( -1.0 + AbromowitzB[1] / b[1] - AbromowitzB[2] / b[3]
                + AbromowitzB[3] / b[4] - AbromowitzB[4] / b[6] + AbromowitzB[5] / b[7] - AbromowitzB[6] / b[9]
                + AbromowitzB[7] / b[10] );

            BesselResult result = new BesselResult();
            result.Regular = F;
            result.RegularPrime = FP;
            result.Irregular = G;
            result.IrregularPrime = GP;
            return (result);

        }


        // calculate and store the coefficients used in the Abromowitz expansion

        private static double[] AbromowitzA, AbromowitzB;

        static AdvancedMath () {

            double g1 = AdvancedMath.Gamma(1.0 / 3.0);
            double g2 = AdvancedMath.Gamma(2.0 / 3.0);
            double r12 = g1 / g2;

            AbromowitzA = new double[] {
                g1 / 2.0 / Math.Sqrt(Math.PI),
                2.0 / 35.0 / r12,
                32.0 / 8100.0,
                92672.0 / 73710000.0 / r12,
                6363008.0 / 35363790000.0,
                1176384772096.0 / 6114399291000000.0 / r12,
                525441777664.0 / 14608781649000000.0,
                181821895706607616.0 / 2525858347112100000000.0 / r12
            };
            AbromowitzB = new double[] {
                g2 / 2.0 / Math.Sqrt(Math.PI),
                1.0 / 15.0 * r12,
                8.0 / 56700.0,
                11488.0 / 18711000.0 * r12,
                25739264 / 417935700000.0,
                1246983424.0 / 18035532900000.0 * r12,
                277651871485952.0 / 14857990277130000000.0,
                1351563588999258112.0 / 64209977981849700000000.0 * r12
            };

        }
        */

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

            }

        }

        private static double CoulombF_Integrate (int L, double eta, double rho) {

            // start at the series limit
            double rho1 = Math.Min(
                4.0 + 2.0 * Math.Sqrt(L),
                (8.0 + 4.0 * L) / Math.Abs(eta)
            );

            double F, FP;
            CoulombF_Series(L, eta, rho1, out F, out FP);

            // TODO: switch so we integrate w/o the C factor, then apply it afterward
            if ((F == 0.0) && (FP == 0.0)) return (0.0);

            OdeResult r = FunctionMath.IntegrateConservativeOde(
                (double x, double y) => ((L * (L + 1) / x + 2.0 * eta) / x - 1.0) * y,
                rho1, F, FP, rho,
                new OdeSettings() {
                    RelativePrecision = 2.5E-13,
                    AbsolutePrecision = 0.0,
                    EvaluationBudget = 8192 * 2
                }
            );

            // return the result
            return (r.Y);

        }

    }

}