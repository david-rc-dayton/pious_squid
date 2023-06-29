import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Two-body Velocity Verlet Blend interpolator.
///
/// The [VerletBlendInterpolator] retains the original ephemerides, so the
/// original _"truth"_ states can be retrieved if needed without imparting any
/// additional error, so this can be used to build other interpolator types.
/// The implementation is simple and very tolerant when working with sparse
/// ephemerides.
class VerletBlendInterpolator extends StateInterpolator {
  /// Create a new [VerletBlendInterpolator] interpolator object from an
  /// ephemeris array.
  VerletBlendInterpolator(this.ephemeris);

  /// Inertal state array.
  final List<J2000> ephemeris;

  @override
  int get sizeBytes => (64 * 7 * ephemeris.length) ~/ 8;

  @override
  EpochWindow window() => (ephemeris.first.epoch, ephemeris.last.epoch);

  static J2000 _getClosest(
          final double target, final J2000 s1, final J2000 s2) =>
      ((target - s1.epoch.posix) >= (s2.epoch.posix - target)) ? s2 : s1;

  J2000 _matchState(final EpochUTC epoch) {
    final target = epoch.posix;
    if (target <= ephemeris.first.epoch.posix) {
      return ephemeris.first;
    }
    if (target >= ephemeris.last.epoch.posix) {
      return ephemeris.last;
    }

    var i = 0;
    var j = ephemeris.length;
    var mid = 0;
    while (i < j) {
      mid = (i + j) >> 1;
      if (ephemeris[mid].epoch.posix == target) {
        return ephemeris[mid];
      }
      if (target < ephemeris[mid].epoch.posix) {
        if ((mid > 0) && (target > ephemeris[mid - 1].epoch.posix)) {
          return _getClosest(target, ephemeris[mid - 1], ephemeris[mid]);
        }
        j = mid;
      } else {
        if ((mid < ephemeris.length - 1) &&
            (target < ephemeris[mid + 1].epoch.posix)) {
          return _getClosest(target, ephemeris[mid], ephemeris[mid + 1]);
        }
        i = mid + 1;
      }
    }
    return ephemeris[mid];
  }

  static Vector3D _gravity(final Vector3D position) {
    final r = position.magnitude();
    return position.scale(-Earth.mu / (r * r * r));
  }

  static J2000 _integrate(final J2000 state, final double step) {
    final x0 = state.position;
    final a0 = _gravity(x0);
    final v0 = state.velocity;
    final x1 = x0.add(v0.scale(step)).add(a0.scale(0.5 * step * step));
    final a1 = _gravity(x1);
    final v1 = v0.add(a0.add(a1).scale(0.5 * step));
    return J2000(state.epoch.roll(step), x1, v1);
  }

  @override
  J2000? interpolate(final EpochUTC epoch) {
    if (!inWindow(epoch)) {
      return null;
    }
    var state = _matchState(epoch);
    while (state.epoch.posix != epoch.posix) {
      final delta = epoch.posix - state.epoch.posix;
      final stepMag = min(5.0, delta.abs());
      final stepSize = copySign(stepMag, delta);
      state = _integrate(state, stepSize);
    }
    return state;
  }

  /// Get the closest cached _"truth"_ state to the provided UTC epoch.
  ///
  /// This method returns `null` if cached ephemeris does not cover the
  /// provided epoch.
  J2000? getCachedState(final EpochUTC epoch) {
    if (!inWindow(epoch)) {
      return null;
    }
    return _matchState(epoch);
  }

  /// Convert this into a [CubicSplineInterpolator] object.
  CubicSplineInterpolator toCubicSpline() =>
      CubicSplineInterpolator.fromEphemeris(ephemeris);

  /// Convert this into a [LagrangeInterpolator] object.
  LagrangeInterpolator toLagrange({final int order = 10}) =>
      LagrangeInterpolator.fromEphemeris(ephemeris, order: order);
}
