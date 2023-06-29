import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Lagrange ephemeris interpolator.
///
/// The [LagrangeInterpolator] is capable of being an exceptionally accurate
/// interpolator at the expense of speed. Since only the position data is
/// needed for interpolation, memory requirements are lighter than interpolation
/// methods that require both position and velocity data.
///
/// The accuracy and speed of the interpolator can be adjusted using the
/// interpolator order parameter. A higher interpolation order results in more
/// accurate interpolation, but takes more time to interpolate. Some reasonable
/// interpolation order values are 8 or 10. This interpolator is fairly
/// tolerant of sparse ephemerides.
class LagrangeInterpolator extends StateInterpolator {
  /// Create a new [LagrangeInterpolator] object from the ephemeris
  /// time/position component arrays and interpolation [order]
  /// _(positive integer)_.
  LagrangeInterpolator(this._t, this._x, this._y, this._z, {this.order = 10});

  /// Create a new [LagrangeInterpolator] object from an array of inertial
  /// state vectors, and an optional interpolation [order] _(positive integer)_.
  factory LagrangeInterpolator.fromEphemeris(final List<J2000> ephemeris,
      {final int order = 10}) {
    final k = ephemeris.length;
    final t = Float64List(k);
    final x = Float64List(k);
    final y = Float64List(k);
    final z = Float64List(k);

    for (var i = 0; i < k; i++) {
      final state = ephemeris[i];
      t[i] = state.epoch.posix;
      x[i] = state.position.x;
      y[i] = state.position.y;
      z[i] = state.position.z;
    }
    return LagrangeInterpolator(t, x, y, z, order: order);
  }

  /// Epochs _(POSIX seconds)_.
  final Float64List _t;

  /// Position x-axis components _(km)_.
  final Float64List _x;

  /// Position y-axis components _(km)_.
  final Float64List _y;

  /// Position z-axis components _(km)_.
  final Float64List _z;

  /// Interpolation order _(positive integer)_.
  final int order;

  @override
  int get sizeBytes => (64 * 4 * _t.length) ~/ 8;

  @override
  J2000? interpolate(final EpochUTC epoch) {
    if (!inWindow(epoch)) {
      return null;
    }
    final posix = epoch.posix;
    final subDex = _slice(posix);
    final start = subDex.left;
    final stop = subDex.right;
    final ts = _t.sublist(start, stop);
    final xs = _x.sublist(start, stop);
    final ys = _y.sublist(start, stop);
    final zs = _z.sublist(start, stop);
    final position = Vector3D(_position(ts, xs, posix),
        _position(ts, ys, posix), _position(ts, zs, posix));
    final velocity = Vector3D(_velocity(ts, xs, posix),
        _velocity(ts, ys, posix), _velocity(ts, zs, posix));
    return J2000(epoch, position, velocity);
  }

  static double _position(
      final Float64List xs, final Float64List ys, final double x) {
    final k = xs.length - 1;
    var result = 0.0;
    for (var j = 0; j < k; j++) {
      var product = ys[j];
      for (var m = 0; m < k; m++) {
        if (j == m) {
          continue;
        }
        product *= (x - xs[m]) / (xs[j] - xs[m]);
      }
      result += product;
    }
    return result;
  }

  static double _velocity(
      final Float64List xs, final Float64List ys, final double x) {
    final k = xs.length;
    var result = 0.0;
    for (var j = 0; j < k; j++) {
      var total = 0.0;
      for (var i = 0; i < k; i++) {
        if (i == j) {
          continue;
        }
        var product = 1 / (xs[j] - xs[i]);
        for (var m = 0; m < k; m++) {
          if (m == i || m == j) {
            continue;
          }
          product *= (x - xs[m]) / (xs[j] - xs[m]);
        }
        total += product;
      }
      result += ys[j] * total;
    }
    return result;
  }

  static int _getClosest(final double target, final double t1, final int d1,
          final double t2, final int d2) =>
      ((target - t1) >= (t2 - target)) ? d2 : d1;

  ({int left, int right}) _slice(final double posix) {
    final n = _t.length;
    if (posix <= _t.first) {
      return (left: 0, right: order);
    }
    if (posix >= _t.last) {
      return (left: n - order, right: n);
    }

    var i = 0;
    var j = _t.length;
    var mid = 0;
    while (i < j) {
      mid = (i + j) >> 1;
      if (_t[mid] == posix) {
        break;
      }
      if (posix < _t[mid]) {
        if ((mid > 0) && (posix > _t[mid - 1])) {
          mid = _getClosest(posix, _t[mid - 1], mid - 1, _t[mid], mid);
          break;
        }
        j = mid;
      } else {
        if ((mid < _t.length - 1) && (posix < _t[mid + 1])) {
          mid = _getClosest(posix, _t[mid], mid, _t[mid + 1], mid + 1);
          break;
        }
        i = mid + 1;
      }
    }
    final offset = order ~/ 2;
    final left = mid - offset;
    final right = mid + offset - (order.isOdd ? 1 : 0);
    if (left < 0) {
      return (left: 0, right: order);
    }
    if (right > n) {
      return (left: n - order, right: n);
    }
    return (left: left, right: right);
  }

  @override
  EpochWindow window() => (EpochUTC(_t.first), EpochUTC(_t.last));
}
