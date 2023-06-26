import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Ephemeris compressor.
class ChebyshevCompressor {
  /// Create a new [ChebyshevCompressor] object from an [Interpolator].
  ChebyshevCompressor(this._interpolator);

  /// Ephemeris interpolator to be compressed.
  final StateInterpolator _interpolator;

  /// Return the cosine of π times [x].
  static double _cosPi(final double x) => cos(pi * x);

  Vector3D _fitCoefficient(
      final int j, final int n, final double a, final double b) {
    var sumX = 0.0;
    var sumY = 0.0;
    var sumZ = 0.0;
    final h = 0.5;
    for (var i = 0; i < n; i++) {
      final x = _cosPi((i + h) / n);
      final state = _interpolator
          .interpolate(EpochUTC(x * (h * (b - a)) + (h * (b + a))))!;
      final fx = state.position.x;
      final fy = state.position.y;
      final fz = state.position.z;
      final nFac = _cosPi(j * (i + h) / n);
      sumX += fx * nFac;
      sumY += fy * nFac;
      sumZ += fz * nFac;
    }
    return Vector3D(sumX * (2 / n), sumY * (2 / n), sumZ * (2 / n));
  }

  ChebyshevCoefficients _fitWindow(
      final int coeffs, final double a, final double b) {
    final cx = Float64List(coeffs);
    final cy = Float64List(coeffs);
    final cz = Float64List(coeffs);
    for (var j = 0; j < coeffs; j++) {
      final result = _fitCoefficient(j, coeffs, a, b);
      cx[j] = result.x;
      cy[j] = result.y;
      cz[j] = result.z;
    }
    return ChebyshevCoefficients(a, b, cx, cy, cz);
  }

  /// Compress this object's interpolater, using the provided coefficients
  /// per revolution [cpr].
  ChebyshevInterpolator compress({final int cpr = 21}) {
    final (start, stop) = _interpolator.window();
    final period = _interpolator.interpolate(start)!.period();
    final coefficients = <ChebyshevCoefficients>[];
    var current = start;
    while (current < stop) {
      final step = min(period, stop.posix - current.posix);
      final segment = current.roll(step);
      coefficients.add(_fitWindow(cpr, current.posix, segment.posix));
      current = segment;
    }
    return ChebyshevInterpolator(coefficients);
  }
}
