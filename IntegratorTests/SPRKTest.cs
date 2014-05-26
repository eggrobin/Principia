﻿using Integrators;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace IntegratorTests {

  [TestClass]
  public class SPRKTest {
    const int repeat_ = 100;
    const double benchmark_tmax_ = 1000;

    [TestMethod]
    public void HarmonicOscillator() {
      SymplecticPartitionedRungeKutta.Solution solution =
          SymplecticPartitionedRungeKutta.IncrementSPRK(
            computeHarmonicOscillatorForce,
            computeHarmonicOscillatorVelocity,
            q0: new double[] { 1 },
            p0: new double[] { 0 },
            t0: 0, tmax: 1000, Δt: .0001,
            coefficients: SymplecticPartitionedRungeKutta.Order5Optimal,
            samplingPeriod: 1);
      double qError = 0;
      double pError = 0;
      for (int i = 0; i < solution.time.Length; ++i) {
        qError = Math.Max(qError,
          Math.Abs(solution.position[i][0] - Math.Cos(solution.time[i])));
        pError = Math.Max(pError,
          Math.Abs(solution.momentum[i][0] + Math.Sin(solution.time[i])));
      }
      Console.WriteLine("qError = " + qError + ";");
      Console.WriteLine("pError = " + pError + ";");
      Assert.AreEqual(0, qError, 1E-12);
      Assert.AreEqual(0, pError, 1E-12);
    }

    [TestMethod]
    public void CsHarmonicOscillatorBenchmark() {
      Console.WriteLine(DateTime.Now);
      for (int i = 1; i <= repeat_; ++i) {
        SymplecticPartitionedRungeKutta.Solution solution =
            SymplecticPartitionedRungeKutta.IncrementSPRK(
                computeHarmonicOscillatorForce,
                computeHarmonicOscillatorVelocity,
                q0: new double[] { 1 },
                p0: new double[] { 0 },
                t0: 0, tmax: benchmark_tmax_, Δt: .0001,
                coefficients: SymplecticPartitionedRungeKutta.Order5Optimal,
                samplingPeriod: 0);
      }
      Console.WriteLine(DateTime.Now);

    }

    [TestMethod]
    public void CppHarmonicOscillatorBenchmark() {
      Console.WriteLine(DateTime.Now);
      for (int i = 1; i <= repeat_; ++i) {
        principia.integrators.SPRKTests.Run(benchmark_tmax_);
      }
      Console.WriteLine(DateTime.Now);
    }

    private void computeHarmonicOscillatorForce(double[] q, double t, ref  double[] result) {
      result[0] = -q[0];
    }

    private void computeHarmonicOscillatorVelocity(double[] p, ref  double[] result) {
      result[0] = p[0];
    }
  }
}