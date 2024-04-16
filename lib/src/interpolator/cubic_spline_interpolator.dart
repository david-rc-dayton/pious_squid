import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Container for cubic spline data.
class CubicSpline {
  /// Create a new [CubicSpline] object.
  CubicSpline(this.t0, final Vector3D p0, final Vector3D m0, this.t1,
      final Vector3D p1, final Vector3D m1) {
    final dx = p1.x - p0.x;
    final dy = p1.y - p0.y;
    final dz = p1.z - p0.z;
    final dvx = m1.x - m0.x;
    final dvy = m1.y - m0.y;
    final dvz = m1.z - m0.z;
    final dx2 = dx * dx;
    final dy2 = dy * dy;
    final dz2 = dz * dz;
    final dx3 = dx2 * dx;
    final dy3 = dy2 * dy;
    final dz3 = dz2 * dz;

    aPos = p0.x;
    bPos = m0.x;
    cPos = (3 * dx - 2 * m0.x - m1.x) / dx;
    dPos = (-2 * dx + m0.x + m1.x) / dx2;
    aVel = m0.x;
    bVel = (2 * dvx - 2 * m0.x - m1.x) / dx;
    cVel = (3 * dx * m0.x - 2 * dx * m1.x - dvx * dx - dvx * dx) / dx2;
    dVel = (-2 * dx * m0.x + dx * m1.x + dvx * dx) / dx3;
    ePos = p0.y;
    fPos = m0.y;
    gPos = (3 * dy - 2 * m0.y - m1.y) / dy;
    hPos = (-2 * dy + m0.y + m1.y) / dy2;
    eVel = m0.y;
    fVel = (2 * dvy - 2 * m0.y - m1.y) / dy;
    gVel = (3 * dy * m0.y - 2 * dy * m1.y - dvy * dy - dvy * dy) / dy2;
    hVel = (-2 * dy * m0.y + dy * m1.y + dvy * dy) / dy3;
    iPos = p0.z;
    jPos = m0.z;
    kPos = (3 * dz - 2 * m0.z - m1.z) / dz;
    lPos = (-2 * dz + m0.z + m1.z) / dz2;
    iVel = m0.z;
    jVel = (2 * dvz - 2 * m0.z - m1.z) / dz;
    kVel = (3 * dz * m0.z - 2 * dz * m1.z - dvz * dz - dvz * dz) / dz2;
    lVel = (-2 * dz * m0.z + dz * m1.z + dvz * dz) / dz3;
  }

  /// Sample start time _(POSIX seconds)_.
  final double t0;

  /// Sample end time _(POSIX seconds)_.
  final double t1;

  late final double aPos;
  late final double bPos;
  late final double cPos;
  late final double dPos;
  late final double aVel;
  late final double bVel;
  late final double cVel;
  late final double dVel;
  late final double ePos;
  late final double fPos;
  late final double gPos;
  late final double hPos;
  late final double eVel;
  late final double fVel;
  late final double gVel;
  late final double hVel;
  late final double iPos;
  late final double jPos;
  late final double kPos;
  late final double lPos;
  late final double iVel;
  late final double jVel;
  late final double kVel;
  late final double lVel;

  /// Interpolate position _(km)_ and velocity _(km/s)_ vectors at the
  /// provided time [t] _(POSIX seconds)_.
  PositionVelocity interpolate(final double t) {
    final n = (t - t0) / (t1 - t0);
    final n2 = n * n;
    final n3 = n2 * n;

    final xPos = aPos + bPos * n + cPos * n2 + dPos * n3;
    final yPos = ePos + fPos * n + gPos * n2 + hPos * n3;
    final zPos = iPos + jPos * n + kPos * n2 + lPos * n3;
    final xVel = aVel + bVel * n + cVel * n2 + dVel * n3;
    final yVel = eVel + fVel * n + gVel * n2 + hVel * n3;
    final zVel = iVel + jVel * n + kVel * n2 + lVel * n3;

    return (
      position: Vector3D(xPos, yPos, zPos),
      velocity: Vector3D(xVel, yVel, zVel)
    );
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
  int get sizeBytes => (64 * 25 * _splines.length) ~/ 8;

  CubicSpline _matchSpline(final double posix) {
    var left = 0;
    var right = _splines.length - 1;

    while (left <= right) {
      final mid = (left + right) ~/ 2;
      if (_splines[mid].t0 <= posix && posix <= _splines[mid].t1) {
        return _splines[mid];
      } else if (posix < _splines[mid].t0) {
        right = mid - 1;
      } else {
        left = mid + 1;
      }
    }

    if (right < 0) {
      return _splines.first;
    } else if (left >= _splines.length) {
      return _splines[-1];
    } else {
      if ((posix - _splines[right].t1).abs() <
          (posix - _splines[left].t1).abs()) {
        return _splines[right];
      } else {
        return _splines[left];
      }
    }
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
