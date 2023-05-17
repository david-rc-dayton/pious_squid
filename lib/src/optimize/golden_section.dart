import 'dart:math';

import 'package:pious_squid/src/operations/functions.dart';

/// Golden Secton bounded single value optimizer.
class GoldenSection {
  static final double _invPhi = 0.5 * (sqrt(5.0) - 1.0);
  static final double _invPhi2 = 0.5 * (3.0 - sqrt(5.0));

  /// Search for an optimal input value for function [f] that minimizes the
  /// output value.
  ///
  /// Takes [lower] and [upper] input search bounds, and an optional
  /// search [tolerance].
  static double search(
      final DifferentiableFunction f, final double lower, final double upper,
      {final double tolerance = 1e-5}) {
    var a = min(lower, upper);
    var b = max(lower, upper);
    var h = b - a;
    if (h <= tolerance) {
      return 0.5 * (a + b);
    }
    final n = (log(tolerance / h) / log(_invPhi)).ceil();
    var c = a + _invPhi2 * h;
    var d = a + _invPhi * h;
    var yc = f(c);
    var yd = f(d);
    for (var k = 0; k < n - 1; k++) {
      if (yc < yd) {
        b = d;
        d = c;
        yd = yc;
        h = _invPhi * h;
        c = a + _invPhi2 * h;
        yc = f(c);
      } else {
        a = c;
        c = d;
        yc = yd;
        h = _invPhi * h;
        d = a + _invPhi2 * h;
        yd = f(d);
      }
    }
    if (yc < yd) {
      return 0.5 * (a + d);
    }
    return 0.5 * (c + b);
  }
}
