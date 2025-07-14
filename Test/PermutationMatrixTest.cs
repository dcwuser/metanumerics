using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Emit;
using System.Text;
using System.Threading.Tasks;
using FluentAssertions;

using Meta.Numerics;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Test {

    [TestClass]
    public class PermutationMatrixTest {

        [TestMethod]
        public void PermutationMatrixMultiplication () {

            Random rng = new Random(2);

            foreach (int n in TestUtilities.GenerateIntegerValues(2, 16).Take(4)) {


                Permutation a = Permutation.GetRandomPermutation(n, rng);

                PermutationMatrix ma = a.ToMatrix();
                ma.Dimension.Should().Be(n);

                SquareMatrix mac = Copy(ma);
                mac.Dimension.Should().Be(ma.Dimension);
                mac.FrobeniusNorm().Should().BeNearly(ma.FrobeniusNorm());
                mac.InfinityNorm().Should().BeNearly(ma.InfinityNorm());
                mac.OneNorm().Should().BeNearly(ma.OneNorm());
                mac.MaxNorm().Should().BeNearly(ma.MaxNorm());

                TestUtilities.IsNearlyEqual(mac.Inverse(), ma.Inverse()).Should().BeTrue();
                TestUtilities.IsNearlyEqual(mac.Transpose, ma.Transpose).Should().BeTrue();

                // Test that eigenvalues agree. To signal that an eigenvalue has appeared we set it to zero.
                Complex[] maEigenvalues = ma.Eigenvalues();
                Complex[] macEigenvalues = mac.Eigenvalues();
                for (int i = 0; i < maEigenvalues.Length; i++) {
                    Assert.IsTrue(maEigenvalues[i] != 0.0);
                    for (int j = 0; j < macEigenvalues.Length; j++) {
                        if (TestUtilities.IsNearlyEqual(macEigenvalues[i], macEigenvalues[j])) {
                            macEigenvalues[j] = 0.0;
                            break;
                        }
                    }
                }
                for (int j = 0; j < maEigenvalues.Length; j++) {
                    Assert.IsTrue(macEigenvalues[j] == 0.0);
                }

                Permutation ai = a.Inverse();
                PermutationMatrix mai = ai.ToMatrix();

                (ma * mai).Equals(Permutation.Identity(n).ToMatrix()).Should().BeTrue();          
                (mai * ma).Equals(Permutation.Identity(n).ToMatrix()).Should().BeTrue();


                Permutation b = Permutation.GetRandomPermutation(n, rng);
                PermutationMatrix mb = b.ToMatrix();

                (ma * mb).Equals((a * b).ToMatrix()).Should().BeTrue();
                (mb * ma).Equals((b * a).ToMatrix()).Should().BeTrue();

                SquareMatrix mbc = Copy(mb);
                (mac * mbc).Equals(Copy((a * b).ToMatrix())).Should().BeTrue();

            }

        }

        private static SquareMatrix Copy (AnySquareMatrix source) {
            SquareMatrix m = new SquareMatrix(source.Dimension);
            for (int r = 0; r < source.Dimension; r++) {
                for (int c = 0; c < source.Dimension; c++) {
                    m[r, c] = source[r, c];
                }
            }
            return m;
        }

        // This computation of the complex roots of unity is careful to ensure
        //   +1 and -1 (when it appears) are always exact
        //   conjugate roots are always exactly conjugate
        // Is use of CosPi and SinPi is suprufluous since the arguments are never reducable.
        // I notice for n=8 not all components are exactly equal.

        private static Complex[] RootsOfUnity(int n) {
            Complex[] roots = new Complex[n];
            roots[0] = Complex.One;
            double delta = 2.0 / n;
            int half_n = (n + 1) / 2;
            for (int i = 1; i < half_n; i++) {
                double theta = delta * i;
                roots[i] = new Complex(MoreMath.CosPi(theta), MoreMath.SinPi(theta));
                roots[n - i] = roots[i].Conjugate;
            }
            if (n % 2 == 0) roots[half_n] = -Complex.One;
            return roots;
        } 

    }
}
