import 'dart:math';

import 'package:pious_squid/src/operations/operations_base.dart';

double _rootMeanSquareError(final Vector x) {
  var total = 0.0;
  for (var i = 0; i < x.length; i++) {
    total += x[i] * x[i];
  }
  return sqrt(total / x.length);
}

/// Gauss-Newton differential corrector.
class GaussNewton {
  GaussNewton._(); // disable constructor

  /// Perform Gauss-Newton differential correction for an initial guess [x0],
  /// residual function [f], and Jacobian function [fPrime].
  static Vector solve(final Vector x0, final Vector Function(Vector) f,
      final Matrix Function(Vector) fPrime,
      {final double tolerance = 1e-6,
      final double alpha = 1.0,
      final int maxIter = 1000,
      final bool safe = true,
      final bool printIter = false}) {
    var xn = x0;
    var rmsLast = double.infinity;
    var h = alpha;
    for (var i = 0; i < maxIter; i++) {
      final fx = f(xn);
      final jx = fPrime(xn);
      Vector dx;
      if (safe) {
        dx = jx.pseudoinverse().multiplyVector(fx).scale(h);
      } else {
        dx = jx.solve(fx).scale(h);
      }
      final dxm = dx.magnitude();
      final x1 = xn.subtract(dx);
      final rms = _rootMeanSquareError(f(x1));
      if (printIter) {
        print('${i + 1}: rms=$rms dx=$dxm h=$h x=$xn');
      }
      if (dxm <= tolerance) {
        return x1;
      }
      if (rms > rmsLast) {
        h *= 0.5;
        continue;
      }
      rmsLast = rms;
      xn = x1;
    }
    return xn;
  }
}
