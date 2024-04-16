import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Container for cubic spline data.
class CubicSpline {
  /// Create a new [CubicSpline] object.
  CubicSpline(this.t0, this.t1, this.a, this.b, this.c, this.d) : dt = t1 - t0;

  /// Create a new spline from a pair of states.
  factory CubicSpline.fromStates(final J2000 first, final J2000 last) {
    final t0 = first.epoch.posix;
    final p0 = first.position;
    final v0 = first.velocity;
    final t1 = last.epoch.posix;
    final p1 = last.position;
    final v1 = last.velocity;
    final dt = t1 - t0;

    final d = p0;
    final c = v0.scale(dt);
    final scaledV1 = v1.scale(dt);
    final a = scaledV1.subtract(p1.scale(2)).add(p0.scale(2)).add(v0.scale(dt));
    final b = p1.subtract(p0).subtract(v0.scale(dt)).subtract(a);
    return CubicSpline(t0, t1, a, b, c, d);
  }

  /// Sample start time _(POSIX seconds)_.
  final double t0;

  /// Sample end time _(POSIX seconds)_.
  final double t1;

  /// Spline duration _(seconds)_.
  final double dt;

  /// Spline coefficient A.
  final Vector3D a;

  /// Spline coefficient B.
  final Vector3D b;

  /// Spline coefficient C.
  final Vector3D c;

  /// Spline coefficient D.
  final Vector3D d;

  /// Interpolate position at the provided normalized time [tn].
  Vector3D _position(final double tn) =>
      d.add(c.scale(tn)).add(b.scale(tn * tn)).add(a.scale(tn * tn * tn));

  /// Interpolate velocity at the provided normalized time [tn].
  Vector3D _velocity(final double tn) =>
      c.scale(1 / dt).add(b.scale(2 * tn / dt)).add(a.scale(3 * tn * tn / dt));

  /// Interpolate position _(km)_ and velocity _(km/s)_ vectors at the
  /// provided time [t] _(POSIX seconds)_.
  PositionVelocity interpolate(final double t) {
    final tn = (t - t0) / (t1 - t0);
    return (position: _position(tn), velocity: _velocity(tn));
  }
}

/// Cubic spline ephemeris interpolator.
///
/// The [CubicSplineInterpolator] is a very fast and accurate interpolator
/// at the expense of memory due to the cached spline pairs used in the
/// interpolation operation. Accuracy is significantly impacted when using
/// sparse ephemerides.
class CubicSplineInterpolator extends StateInterpolator {
  /// Create a new [CubicSplineInterpolator] object from a set of [_splines].
  CubicSplineInterpolator(this._splines);

  /// Create a new [CubicSplineInterpolator] from an [ephemeris] list.
  factory CubicSplineInterpolator.fromEphemeris(final List<J2000> ephemeris) {
    final splines = <CubicSpline>[];
    for (var i = 0; i < ephemeris.length - 1; i++) {
      final e0 = ephemeris[i];
      final e1 = ephemeris[i + 1];
      splines.add(CubicSpline.fromStates(e0, e1));
    }
    return CubicSplineInterpolator(splines);
  }

  /// Cached splines.
  final List<CubicSpline> _splines;

  @override
  int get sizeBytes => (64 * 14 * _splines.length) ~/ 8;

  CubicSpline _matchSpline(final double posix) {
    var low = 0;
    var high = _splines.length - 1;

    while (low <= high) {
      final mid = low + ((high - low) >> 1);
      final midSpline = _splines[mid];

      if (posix < midSpline.t0) {
        high = mid - 1;
      } else if (posix > midSpline.t1) {
        low = mid + 1;
      } else {
        return _splines[mid];
      }
    }

    throw 'Corresponding spline not found for timestamp: $posix';
  }

  @override
  J2000? interpolate(final EpochUTC epoch) {
    if (!inWindow(epoch)) {
      return null;
    }
    final posix = epoch.posix;
    final splineVecs = _matchSpline(posix).interpolate(posix);
    return J2000(epoch, splineVecs.position, splineVecs.velocity);
  }

  @override
  EpochWindow window() =>
      (EpochUTC(_splines.first.t0), EpochUTC(_splines.last.t1));
}
