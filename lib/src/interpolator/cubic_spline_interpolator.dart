import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Container for cubic spline data.
class CubicSpline {
  /// Create a new [CubicSpline] object.
  CubicSpline(this.t0, this.p0, this.m0, this.t1, this.p1, this.m1);

  /// Sample start time _(POSIX seconds)_.
  final double t0;

  /// Sample start position vector _(km)_.
  final Vector3D p0;

  /// Sample start velocity vector _(km)_.
  final Vector3D m0;

  /// Sample end time _(POSIX seconds)_.
  final double t1;

  /// Sample end position vector _(km)_.
  final Vector3D p1;

  /// Sample end velocity vector _(km)_.
  final Vector3D m1;

  /// Interpolate position at the provided time [t] _(POSIX seconds)_.
  Vector3D _position(final double t) {
    final t2 = t * t;
    final t3 = t2 * t;
    final r0 = p0.scale(2 * t3 - 3 * t2 + 1);
    final v0 = m0.scale((t3 - 2 * t2 + t) * (t1 - t0));
    final r1 = p1.scale(-2 * t3 + 3 * t2);
    final v1 = m1.scale((t3 - t2) * (t1 - t0));
    return r0.add(v0).add(r1).add(v1);
  }

  /// Interpolate velocity at the provided time [t] _(POSIX seconds)_.
  Vector3D _velocity(final double t) {
    final t2 = t * t;
    final r0 = p0.scale(6 * t2 - 6 * t);
    final v0 = m0.scale((3 * t2 - 4 * t + 1) * (t1 - t0));
    final r1 = p1.scale(-6 * t2 + 6 * t);
    final v1 = m1.scale((3 * t2 - 2 * t) * (t1 - t0));
    return r0.add(v0).add(r1).add(v1).scale(1 / (t1 - t0));
  }

  /// Interpolate position _(km)_ and velocity _(km/s)_ vectors at the
  /// provided time [t] _(POSIX seconds)_.
  PositionVelocity interpolate(final double t) {
    final n = (t - t0) / (t1 - t0);
    return (position: _position(n), velocity: _velocity(n));
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
      final t0 = e0.epoch.posix;
      final p0 = e0.position;
      final m0 = e0.velocity;
      final e1 = ephemeris[i + 1];
      final t1 = e1.epoch.posix;
      final p1 = e1.position;
      final m1 = e1.velocity;
      splines.add(CubicSpline(t0, p0, m0, t1, p1, m1));
    }
    return CubicSplineInterpolator(splines);
  }

  /// Cached splines.
  final List<CubicSpline> _splines;

  @override
  int get sizeBytes => (64 * 14 * _splines.length) ~/ 8;

  CubicSpline _matchSpline(final double posix) {
    var left = 0;
    var right = _splines.length - 1;

    while (left <= right) {
      final mid = (left + right) >> 1;
      if (_splines[mid].t0 <= posix && posix <= _splines[mid].t1) {
        return _splines[mid];
      } else if (posix < _splines[mid].t1) {
        right = mid - 1;
      } else {
        left = mid + 1;
      }
    }

    return _splines[left];
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
