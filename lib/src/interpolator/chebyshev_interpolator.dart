import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Chebyshev compressed ephemeris coefficients.
class ChebyshevCoefficients {
  /// Create a new [ChebyshevCoefficients] object.
  ChebyshevCoefficients(this.a, this.b, this._cx, this._cy, this._cz) {
    _cxd = _derivative(a, b, _cx);
    _cyd = _derivative(a, b, _cy);
    _czd = _derivative(a, b, _cz);
  }

  /// Coefficient start epoch _(POSIX seconds)_.
  final double a;

  /// Coefficient end epoch _(POSIX seconds)_.
  final double b;

  /// Position x-component coefficients.
  final Float64List _cx;

  /// Position y-component coefficients.
  final Float64List _cy;

  /// Position z-component coefficients.
  final Float64List _cz;

  /// Velocity x-component coefficients.
  late final Float64List _cxd;

  /// Velocity y-component coefficients.
  late final Float64List _cyd;

  /// Velocity z-component coefficients.
  late final Float64List _czd;

  /// Create a new array containing the coefficients needed to compute the
  /// derivative of the given Chebyshev coefficient set, using start time [a],
  /// end time [b], and coefficients [c].
  static Float64List _derivative(
      final double a, final double b, final Float64List c) {
    final n = c.length;
    final d = Float64List(n);
    d[n - 1] = 0;
    d[n - 2] = 2 * (n - 1) * c[n - 1];
    for (var k = n - 3; k >= 0; k--) {
      d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1];
    }
    for (var k = 0; k < n; k++) {
      d[k] *= 2 / (b - a);
    }
    return d;
  }

  /// Return the size _(bytes)_ of this coefficient set's cached data.
  int get sizeBytes => ((64 * 2) + (64 * 3 * _cx.length)) ~/ 8;

  /// Return the interpolated coefficient set [c] result at the provided
  /// epoch [t] _(POSIX seconds)_.
  double evaluate(final Float64List c, final double t) {
    final n = c.length;
    final x = (t - (0.5 * (b + a))) / (0.5 * (b - a));
    final alpha = 2 * x;
    final beta = -1;
    var y1 = 0.0;
    var y2 = 0.0;
    for (var k = n - 1; k >= 1; k--) {
      final tmp = y1;
      y1 = alpha * y1 + beta * y2 + c[k];
      y2 = tmp;
    }
    return x * y1 - y2 + 0.5 * c[0];
  }

  /// Return the interpolated position _(km)_ and velocity _(km/s)_ vectors at
  /// he provided time [t] _(POSIX seconds)_.
  (Vector, Vector) interpolate(final double t) {
    final x = evaluate(_cx, t);
    final y = evaluate(_cy, t);
    final z = evaluate(_cz, t);
    final xd = evaluate(_cxd, t);
    final yd = evaluate(_cyd, t);
    final zd = evaluate(_czd, t);
    final pos = Float64List(3);
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    final vel = Float64List(3);
    vel[0] = xd;
    vel[1] = yd;
    vel[2] = zd;
    return (Vector(pos), Vector(vel));
  }
}

/// Compressed Chebyshev ephemeris interpolator.
///
/// The [ChebyshevInterpolator] sacrifices state accuracy in order to gain
/// speed and memory requirements, so you can store many ephemerides in RAM
/// and interpolate through them quickly.
///
/// Using more [_coefficients] per revolution during lossy compression results
/// in increased accuracy and decreased performance.
class ChebyshevInterpolator extends StateInterpolator {
  /// Create a ne [ChebyshevInterpolator] object given the
  /// ephemeris [_coefficients].
  ChebyshevInterpolator(this._coefficients);
  final List<ChebyshevCoefficients> _coefficients;

  int _calcSizeBytes() {
    var output = 0;
    for (final coeffs in _coefficients) {
      output += coeffs.sizeBytes;
    }
    return output;
  }

  @override
  int get sizeBytes => _calcSizeBytes();

  @override
  J2000? interpolate(final EpochUTC epoch) {
    if (!inWindow(epoch)) {
      return null;
    }
    final coeffs = _matchCoefficients(epoch.posix);
    final (pos, vel) = coeffs.interpolate(epoch.posix);
    return J2000(epoch, pos, vel);
  }

  @override
  (EpochUTC, EpochUTC) window() =>
      (EpochUTC(_coefficients.first.a), EpochUTC(_coefficients.last.b));

  ChebyshevCoefficients _matchCoefficients(final double posix) {
    var left = 0;
    var right = _coefficients.length;
    while (left < right) {
      final middle = (left + right) >> 1;
      if (_coefficients[middle].b < posix) {
        left = middle + 1;
      } else {
        right = middle;
      }
    }
    return _coefficients[left];
  }
}
