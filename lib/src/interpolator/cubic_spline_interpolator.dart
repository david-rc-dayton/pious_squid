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

    _aPos = p0.x;
    _bPos = m0.x;
    _cPos = (3 * dx - 2 * m0.x - m1.x) / dx;
    _dPos = (-2 * dx + m0.x + m1.x) / dx2;
    _aVel = m0.x;
    _bVel = (2 * dvx - 2 * m0.x - m1.x) / dx;
    _cVel = (3 * dx * m0.x - 2 * dx * m1.x - dvx * dx - dvx * dx) / dx2;
    _dVel = (-2 * dx * m0.x + dx * m1.x + dvx * dx) / dx3;
    _ePos = p0.y;
    _fPos = m0.y;
    _gPos = (3 * dy - 2 * m0.y - m1.y) / dy;
    _hPos = (-2 * dy + m0.y + m1.y) / dy2;
    _eVel = m0.y;
    _fVel = (2 * dvy - 2 * m0.y - m1.y) / dy;
    _gVel = (3 * dy * m0.y - 2 * dy * m1.y - dvy * dy - dvy * dy) / dy2;
    _hVel = (-2 * dy * m0.y + dy * m1.y + dvy * dy) / dy3;
    _iPos = p0.z;
    _jPos = m0.z;
    _kPos = (3 * dz - 2 * m0.z - m1.z) / dz;
    _lPos = (-2 * dz + m0.z + m1.z) / dz2;
    _iVel = m0.z;
    _jVel = (2 * dvz - 2 * m0.z - m1.z) / dz;
    _kVel = (3 * dz * m0.z - 2 * dz * m1.z - dvz * dz - dvz * dz) / dz2;
    _lVel = (-2 * dz * m0.z + dz * m1.z + dvz * dz) / dz3;
  }

  /// Sample start time _(POSIX seconds)_.
  final double t0;

  /// Sample end time _(POSIX seconds)_.
  final double t1;

  late final double _aPos;
  late final double _bPos;
  late final double _cPos;
  late final double _dPos;
  late final double _aVel;
  late final double _bVel;
  late final double _cVel;
  late final double _dVel;
  late final double _ePos;
  late final double _fPos;
  late final double _gPos;
  late final double _hPos;
  late final double _eVel;
  late final double _fVel;
  late final double _gVel;
  late final double _hVel;
  late final double _iPos;
  late final double _jPos;
  late final double _kPos;
  late final double _lPos;
  late final double _iVel;
  late final double _jVel;
  late final double _kVel;
  late final double _lVel;

  /// Interpolate position _(km)_ and velocity _(km/s)_ vectors at the
  /// provided time [t] _(POSIX seconds)_.
  PositionVelocity interpolate(final double t) {
    final n = (t - t0) / (t1 - t0);
    final n2 = n * n;
    final n3 = n2 * n;

    final xPos = _aPos + _bPos * n + _cPos * n2 + _dPos * n3;
    final yPos = _ePos + _fPos * n + _gPos * n2 + _hPos * n3;
    final zPos = _iPos + _jPos * n + _kPos * n2 + _lPos * n3;
    final xVel = _aVel + _bVel * n + _cVel * n2 + _dVel * n3;
    final yVel = _eVel + _fVel * n + _gVel * n2 + _hVel * n3;
    final zVel = _iVel + _jVel * n + _kVel * n2 + _lVel * n3;

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
      final mid = (left + right) >> 1;
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
