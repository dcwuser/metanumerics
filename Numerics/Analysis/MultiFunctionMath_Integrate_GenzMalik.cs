using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Threading;

using Meta.Numerics.Functions;

namespace Meta.Numerics.Analysis {

    public static partial class MultiFunctionMath {

        private static UncertainValue Integrate_Adaptive (MultiFunctor f, CoordinateTransform[] map, IntegrationRegion r, EvaluationSettings settings) {

            // Create an evaluation rule
            GenzMalik75Rule rule = new GenzMalik75Rule(r.Dimension);

            // Use it on the whole region and put the answer in a linked list
            rule.Evaluate(f, map, r);
            LinkedList<IntegrationRegion> regionList = new LinkedList<IntegrationRegion>();
            regionList.AddFirst(r);

            // Iterate until convergence
            while (f.EvaluationCount < settings.EvaluationBudget) {

                // Add up value and errors in all regions.
                // While we are at it, take note of the region with the largest error.
                double value = 0.0;
                double error = 0.0;
                LinkedListNode<IntegrationRegion> regionNode = regionList.First;
                double maxError = 0.0;
                LinkedListNode<IntegrationRegion> maxErrorNode = null;
                while (regionNode != null) {
                    IntegrationRegion region = regionNode.Value;
                    value += region.Value;
                    error += region.Error;
                    //error += MoreMath.Sqr(region.Error);
                    if (region.Error > maxError) {
                        maxError = region.Error;
                        maxErrorNode = regionNode;
                    }
                    regionNode = regionNode.Next;
                }

                // Check for convergence.
                if ((error <= settings.AbsolutePrecision) || (error <= settings.RelativePrecision * Math.Abs(value))) {
                    return (new UncertainValue(value, error));
                }

                // Split the region with the largest error, and evaluate each subregion.
                IntegrationRegion maxErrorRegion = maxErrorNode.Value;
                regionList.Remove(maxErrorNode);
                IList<IntegrationRegion> subRegions = maxErrorRegion.Split(maxErrorRegion.SplitIndex);
                /*
                Countdown cnt = new Countdown(2);
                ThreadPool.QueueUserWorkItem((object state) => { rule.Evaluate(f, subRegions[0]); cnt.Signal(); });
                ThreadPool.QueueUserWorkItem((object state) => { rule.Evaluate(f, subRegions[1]); cnt.Signal(); });
                cnt.Wait();
                foreach (IntegrationRegion subRegion in subRegions) {
                    regionList.AddLast(subRegion);
                }
                */
                
                foreach (IntegrationRegion subRegion in subRegions) {
                    rule.Evaluate(f, map, subRegion);
                    regionList.AddLast(subRegion);
                }

                if (f.EvaluationCount >= settings.EvaluationBudget) {
                    throw new NonconvergenceException();
                }
                
            }

            throw new NonconvergenceException();

        }

    }

#if FUTURE
    public class Countdown {
        object _locker = new object();
        int _value;

        public Countdown () { }
        public Countdown (int initialCount) { _value = initialCount; }

        public void Signal () { AddCount(-1); }

        public void AddCount (int amount) {
            lock (_locker) {
                _value += amount;
                if (_value <= 0) Monitor.PulseAll(_locker);
            }
        }

        public void Wait () {
            lock (_locker)
                while (_value > 0)
                    Monitor.Wait(_locker);
        }
    }
#endif

    // Genz and Malik developed a 7th order cubature rule with an embedded 5th order rule for arbitrary dimensions.
    // For integration over the symmetric unit cube
    //   I = \integrate_{-1}^{+1} \! dx^{d} \, f(x_1, x_2, \cdots, x_d)
    // their estimate is
    //   I = w_0 f(0, 0, 0, \cdots, 0)                             {1 term}
    //        + w_1 [ f(a_1, 0, 0 , \cdots, 0) + \cdots ]          {2d terms}
    //        + w_2 [ f(a_2, 0, 0, \cdots, 0) + \cdots ]           {2d terms}
    //        + w_3 [ f(a_3, a_3, 0, \cdots, 0) + \cdots ]         {4d(d-1)/2 terms}
    //        + w_4 [ f(a_4, a_4, a_4, \cdots, a_4) + \cdots ]     {2^d terms}
    // where [ ] indicates summation over all positive and negative permutations of arguments, the argument values are
    //   a_1 = \sqrt{9/70}    a_2 = a_3 = \sqrt{9/10}    a_5 = \sqrt{9/19}
    // and the weights are
    //   w_0 = 2^d (12824 - 9120 d + 400 d^2) / 19683
    //   w_1 = 2^d (980 / 6561)
    //   w_2 = 2^d (1820 - 400 d) / 19683
    //   w_3 = 2^d (200 / 19683)
    //   w_4 = 6859 / 19683
    // This estimate is 7th order, i.e. exact for all linear combinations of terms x_1^{p_1} + x_2^{p_2} + \cdots x_d^{p_d} such that p_1 + p_2 + \cdots p_d \le 7.

    // Note as a check, in order to integrate a constant correctly, each weight times the number of terms of that weight must add to 2^d for all d.

    // The total number of evaluations as a function of dimension is n = 2^d + 2 d^2 + 2 d + 1.
    //   d     2     3     4     5     6     7     8
    //   n     17    33    57    93    149   241   401
    
    // The Genz-Malik is nice not only because it is a rule for arbitrary dimension, but also because there is an embedded (i.e. using the same points)
    // 5th order rule, which allows us to obtain an error estimate without any additional evaluations. The weights for the 5th order estimate are
    //   w_0 = 2^d (729 - 950 d + 50 d^2) / 729
    //   w_1 = 2^d 245 / 468
    //   w_2 = 2^d (265 - 100 d) / 1458
    //   w_3 = 2^d 25 / 729
    // and w_4 = 0.

    // It's amazingly difficult to get access to the original papers on cubature formulas. These results are quoted in
    // Huang & Oosterlee, "Adaptive integration for multi-factor portfolio credit loss models", available at
    // http://ta.twi.tudelft.nl/mf/users/oosterlee/oosterlee/huang_fs.pdf and other URLs.

    internal class GenzMalik75Rule {

        public GenzMalik75Rule (int d) {
            this.d = d;

            // compute and store 2^d
            td = 1 << d;

            // compute degree-7 weights
            w70 = (12824 - 9120 * d + 400 * d * d) / 19683.0;
            w71 = 980.0 / 6561.0;
            w72 = (1820 - 400 * d) / 19683.0;
            w73 = 200.0 / 19683.0;
            w74 = 6859.0 / 19683.0 / td;

            // compute degree-5 weights
            w50 = (729 - 950 * d + 50 * d * d) / 729.0;
            w51 = 245.0 / 486.0;
            w52 = (265 - 100 * d) / 1458.0;
            w53 = 25.0 / 729.0;
        }

        private readonly int d, td;

        private static readonly double a1 = Math.Sqrt(9.0 / 70.0), a2 = Math.Sqrt(9.0 / 10.0), a3 = a2, a4 = Math.Sqrt(9.0 / 19.0);

        private readonly double w70, w71, w72, w73, w74;
        private readonly double w50, w51, w52, w53;


        public void Evaluate (MultiFunctor f, IntegrationRegion r) {

            // evaluation at origin
            // Keep a vector of origin coordinates, we will often re-set components to them.
            double[] x0 = new double[d];
            double[] x = new double[d];
            for (int i = 0; i < x.Length; i++) {
                x0[i] = r.MapCoordinateFromSymmetricUnitInterval(i, 0.0);
                x[i] = x0[i];
            }
            double f0 = f.Evaluate(x);
            double I7 = w70 * f0;
            double I5 = w50 * f0;

            // near off-one-axis evaluations
            double f1 = 0.0;
            for (int i = 0; i < d; i++) {
                x[i] = r.MapCoordinateFromSymmetricUnitInterval(i, a1);
                f1 += f.Evaluate(x);
                x[i] = r.MapCoordinateFromSymmetricUnitInterval(i, -a1);
                f1 += f.Evaluate(x);
                x[i] = x0[i];
            }
            I7 += w71 * f1;
            I5 += w51 * f1;

            // far off-one-axis evaluations
            // while doing this, determine along which direction we will split, if asked
            int sIndex = 0; double sMax = 0.0;
            double f2 = 0.0;
            for (int i = 0; i < d; i++) {
                x[i] = r.MapCoordinateFromSymmetricUnitInterval(i, a2);
                double fp2 = f.Evaluate(x);
                f2 += fp2;
                x[i] = r.MapCoordinateFromSymmetricUnitInterval(i, -a2);
                double fm2 = f.Evaluate(x);
                f2 += fm2;
                x[i] = x0[i];
                // One way to view fm2 + fp2 - 2 f0 is as numerical 2nd derivative
                // Another way is as difference of f0 from value predicted by linear interpolation from fm2 and fp2
                double s = Math.Abs(fm2 + fp2 - 2.0 * f0);
                if (s > sMax) {
                    // Use the direction of largest difference for future splits.
                    sMax = s;
                    sIndex = i;
                } else if ((s == sMax) && (r.CoordinateWidth(i) > r.CoordinateWidth(sIndex))) {
                    // If two directions both have the same difference, pick the one with the larger width to split.
                    // This may seem like an unimportant corner case, but it turns out to be critical to convergence.
                    // Without this clause we get into cycles in which we split along the same direction forever.
                    // For example, given the 3D watson integral of [1 - \cos(x) \cos(y) \cos(z)]^{-1} over [0,\pi]^3,
                    // all directions have zero deficit from starting point because one of the cosines is always zero,
                    // and without this branch we split along the x-axis forever.
                    sIndex = i;
                }
            }
            I7 += w72 * f2;
            I5 += w52 * f2;
            r.SplitIndex = sIndex;

            // far off-two-axis evaluations
            double f3 = 0.0;
            for (int i = 0; i < d; i++) {
                for (int j = 0; j < i; j++) {
                    // ++
                    x[i] = r.MapCoordinateFromSymmetricUnitInterval(i, a3);
                    x[j] = r.MapCoordinateFromSymmetricUnitInterval(j, a3);
                    f3 += f.Evaluate(x);
                    // +-
                    x[j] = r.MapCoordinateFromSymmetricUnitInterval(j, -a3);
                    f3 += f.Evaluate(x);
                    // --
                    x[i] = r.MapCoordinateFromSymmetricUnitInterval(i, -a3);
                    f3 += f.Evaluate(x);
                    // -+
                    x[j] = r.MapCoordinateFromSymmetricUnitInterval(j, a3);
                    f3 += f.Evaluate(x);
                    x[i] = x0[i];
                    x[j] = x0[j];
                }
            }
            I7 += w73 * f3;
            I5 += w53 * f3;

            // mid off-all-axis evaluations
            // We need all 2^d permutations of + and - in each position.
            // So that we only need to change one component each time, we proceed in Gray code order.
            // We use a bit-vector to keep track of which component is + (0) and which is - (1).
            int state = 0;
            for (int j = 0; j < d; j++) {
                x[j] = r.MapCoordinateFromSymmetricUnitInterval(j, a4);
            }
            double f4 = f.Evaluate(x);
            for (int i = 0; i < (td - 1); i++) {
                int j = GrayFlipIndex(i);
                int mask = 1 << j;
                state = state ^ mask;
                x[j] = r.MapCoordinateFromSymmetricUnitInterval(j, ((state & mask) > 0) ? -a4 : a4);
                f4 += f.Evaluate(x);
            }
            I7 += w74 * f4;

            double V = r.Volume();

            r.Value = V * I7;
            r.Error = V * Math.Abs(I7 - I5);

        }

        private static int GrayFlipIndex (int j) {
            int i = 0;
            while ((j & 1) != 0) {
                j = j >> 1;
                i++;
            }
            return (i);
        }

        public void Evaluate (MultiFunctor f, CoordinateTransform[] map, IntegrationRegion r) {

            if (map == null) {
                Evaluate(f, r);
                return;
            }

            // Evaluation at origin. Keep a vector of origin coordinates and jacobian values.
            double[] x0 = new double[d];
            double[] x = new double[d];
            double[] v0 = new double[d];
            double[] v = new double[d];
            for (int i = 0; i < x.Length; i++) {
                x0[i] = map[i].Transform(r.MapCoordinateFromSymmetricUnitInterval(i, 0.0), out v0[i]);
                x[i] = x0[i];
                v[i] = v0[i];
            }
            double f0 = Jacobian(v) * f.Evaluate(x);
            double I7 = w70 * f0;
            double I5 = w50 * f0;

            // near off-one-axis evaluations
            double f1 = 0.0;
            for (int i = 0; i < d; i++) {
                x[i] = map[i].Transform(r.MapCoordinateFromSymmetricUnitInterval(i, a1), out v[i]);
                f1 += Jacobian(v) * f.Evaluate(x);
                x[i] = map[i].Transform(r.MapCoordinateFromSymmetricUnitInterval(i, -a1), out v[i]);
                f1 += Jacobian(v) * f.Evaluate(x);
                x[i] = x0[i];
                v[i] = v0[i];
            }
            I7 += w71 * f1;
            I5 += w51 * f1;

            // far off-one-axis evaluations
            int sIndex = 0; double sMax = 0.0;
            double f2 = 0.0;
            for (int i = 0; i < d; i++) {
                x[i] = map[i].Transform(r.MapCoordinateFromSymmetricUnitInterval(i, a2), out v[i]);
                double fp2 = Jacobian(v) * f.Evaluate(x);
                f2 += fp2;
                x[i] = map[i].Transform(r.MapCoordinateFromSymmetricUnitInterval(i, -a2), out v[i]);
                double fm2 = Jacobian(v) * f.Evaluate(x);
                f2 += fm2;
                x[i] = x0[i];
                v[i] = v0[i];

                double s = Math.Abs(fm2 + fp2 - 2.0 * f0);
                if (s > sMax) {
                    sMax = s;
                    sIndex = i;
                } else if ((s == sMax) && (r.CoordinateWidth(i) > r.CoordinateWidth(sIndex))) {
                    sIndex = i;
                }
            }
            I7 += w72 * f2;
            I5 += w52 * f2;

            // far off-two-axis evaluations
            double f3 = 0.0;
            for (int i = 0; i < d; i++) {
                for (int j = 0; j < i; j++) {
                    // ++
                    x[i] = map[i].Transform(r.MapCoordinateFromSymmetricUnitInterval(i, a3), out v[i]);
                    x[j] = map[j].Transform(r.MapCoordinateFromSymmetricUnitInterval(j, a3), out v[j]);
                    f3 += Jacobian(v) * f.Evaluate(x);
                    // +-
                    x[j] = map[j].Transform(r.MapCoordinateFromSymmetricUnitInterval(j, -a3), out v[j]);
                    f3 += Jacobian(v) * f.Evaluate(x);
                    // --
                    x[i] = map[i].Transform(r.MapCoordinateFromSymmetricUnitInterval(i, -a3), out v[i]);
                    f3 += Jacobian(v) * f.Evaluate(x);
                    // -+
                    x[j] = map[j].Transform(r.MapCoordinateFromSymmetricUnitInterval(j, a3), out v[j]);
                    f3 += Jacobian(v) * f.Evaluate(x);
                    x[i] = x0[i];
                    x[j] = x0[j];
                    v[i] = v0[i];
                    v[j] = v0[j];
                }
            }
            I7 += w73 * f3;
            I5 += w53 * f3;

            // mid off-all-axis evaluations
            // We need all 2^d permutations of + and - in each position.
            // So that we only need to change one component each time, we proceed in Gray code order.
            // We use a bit-vector to keep track of which component is + (0) and which is - (1).
            int state = 0;
            for (int j = 0; j < d; j++) {
                x[j] = map[j].Transform(r.MapCoordinateFromSymmetricUnitInterval(j, a4), out v[j]);
            }
            double f4 = Jacobian(v) * f.Evaluate(x);
            for (int i = 0; i < (td - 1); i++) {
                int j = GrayFlipIndex(i);
                int mask = 1 << j;
                state = state ^ mask;
                x[j] = map[j].Transform(r.MapCoordinateFromSymmetricUnitInterval(j, ((state & mask) > 0) ? -a4 : a4), out v[j]);
                f4 += Jacobian(v) * f.Evaluate(x);
            }
            I7 += w74 * f4;

            double V = r.Volume();

            r.Value = V * I7;
            r.Error = V * Math.Abs(I7 - I5);
            r.SplitIndex = sIndex;
        }

        private static double Jacobian (double[] j) {
            double result = 1.0;
            for (int i = 0; i < j.Length; i++) {
                result *= j[i];
            }
            return (result);
        }
    }

    internal class IntegrationRegion {

        public IntegrationRegion (IList<Interval> box) {
            this.box = box;
        }

        private readonly IList<Interval> box;

        public int Dimension {
            get {
                return(box.Count);
            }
        }

        public double Volume () {
            double V = 1.0;
            for (int i = 0; i < box.Count; i++) {
                V *= box[i].Width;
            }
            return (V);
        }

        public double MapCoordinateFromUnitInterval (int i, double t) {
            return (box[i].LeftEndpoint + t * box[i].Width);
        }

        public double MapCoordinateFromSymmetricUnitInterval (int i, double t) {
            return (box[i].Midpoint + t * box[i].Width / 2.0);
        }

        public double CoordinateWidth (int i) {
            return (box[i].Width);
        }

        public IList<IntegrationRegion> Split (int i) {
            Interval[] box1 = new Interval[box.Count];
            Interval[] box2 = new Interval[box.Count];
            for (int j = 0; j < box.Count; j++) {
                if (j == i) {
                    box1[j] = Interval.FromEndpoints(box[j].LeftEndpoint, box[j].Midpoint);
                    box2[j] = Interval.FromEndpoints(box[j].Midpoint, box[j].RightEndpoint);
                } else {
                    box1[j] = box[j];
                    box2[j] = box[j];
                }
            }
            return (new IntegrationRegion[] { new IntegrationRegion(box1), new IntegrationRegion(box2) });
        }

        public int SplitIndex { get; set; }

        public double Value { get; set; }

        public double Error { get; set; }

    }

}
