import 'package:pious_squid/src/operations/operations_base.dart';

/// Gauss-Newton differential corrector.
class GaussNewton {
  GaussNewton._(); // disable constructor

  /// Perform Gauss-Newton differential correction for an initial guess [x0],
  /// residual function [f], and Jacobian function [fPrime].
  static Vector solve(final Vector x0, final Vector Function(Vector) f,
      final Matrix Function(Vector) fPrime,
      {final double tolerance = 1e-10,
      final double alpha = 1.0,
      final int maxIter = 1000,
      final bool printIter = false}) {
    var xn = x0;
    var rmsLast = double.infinity;
    var alphaT = alpha;
    for (var i = 0; i < maxIter; i++) {
      final fx = f(xn);
      final j = fPrime(xn);
      final jt = j.transpose();
      final js = jt.multiply(j).inverse().multiply(jt);
      final x1 = xn.subtract(js.multiplyVector(fx).scale(alphaT));
      final rms = f(x1).magnitude();
      if (printIter) {
        print('${i + 1}: $rms $xn');
      }
      if ((rmsLast - rms).abs() <= tolerance) {
        return x1;
      }
      if (rms > rmsLast) {
        alphaT *= 0.5;
        continue;
      }
      rmsLast = rms;
      xn = x1;
    }
    return xn;
  }
}
