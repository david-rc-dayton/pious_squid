import 'dart:math';

import 'package:pious_squid/src/operations/functions.dart';

/// Golden Secton bounded single value optimizer.
class GoldenSection {
  GoldenSection._(); // disable constructor

  static final double _grInv = 1.0 / (0.5 * (sqrt(5) + 1));

  static bool _check(final double fc, final double fd, final bool solveMax) =>
      solveMax ? fc > fd : fc < fd;

  /// Search for an optimal input value for function [f] that minimizes the
  /// output value.
  ///
  /// Takes [lower] and [upper] input search bounds, and an optional
  /// search [tolerance].
  static double search(
      final DifferentiableFunction f, final double lower, final double upper,
      {final double tolerance = 1e-5, final bool solveMax = false}) {
    var a = lower;
    var b = upper;
    var c = b - (b - a) * _grInv;
    var d = a + (b - a) * _grInv;
    while ((b - a).abs() > tolerance) {
      if (_check(f(c), f(d), solveMax)) {
        b = d;
      } else {
        a = c;
      }
      c = b - (b - a) * _grInv;
      d = a + (b - a) * _grInv;
    }
    return 0.5 * (b + a);
  }
}
