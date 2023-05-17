import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/observation/observation_utils.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Topocentric right-ascension and declination.
class RadecTopocentric {
  /// Create a new [RadecTopocentric] object.
  RadecTopocentric(this.epoch, this.rightAscension, this.declination,
      [this.range,
      this.rightAscensionRate,
      this.declinationRate,
      this.rangeRate]);

  /// Create a new [RadecTopocentric] object, using degrees for the
  /// angular values.
  factory RadecTopocentric.fromDegrees(final EpochUTC epoch,
      final double rightAscensionDegrees, final double declinationDegrees,
      [final double? range,
      final double? rightAscensionRateDegrees,
      final double? declinationRateDegrees,
      final double? rangeRate]) {
    final rightAscensionRate = (rightAscensionRateDegrees != null)
        ? rightAscensionRateDegrees * deg2rad
        : null;
    final declinationRate = (declinationRateDegrees != null)
        ? declinationRateDegrees * deg2rad
        : null;
    return RadecTopocentric(
      epoch,
      rightAscensionDegrees * deg2rad,
      declinationDegrees * deg2rad,
      range,
      rightAscensionRate,
      declinationRate,
      rangeRate,
    );
  }

  /// Create a [RadecTopocentric] object from an inertial [state] and
  /// [site] vector.
  factory RadecTopocentric.fromStateVectors(
      final J2000 state, final J2000 site) {
    final p = state.position.add(site.position.negate());
    final pI = p.x;
    final pJ = p.y;
    final pK = p.z;
    final pMag = p.magnitude();
    final declination = asin(pK / pMag);
    final pDot = state.velocity.add(site.velocity.negate());
    final pIDot = pDot.x;
    final pJDot = pDot.y;
    final pKDot = pDot.z;
    final pIJMag = sqrt(pI * pI + pJ * pJ);
    double rightAscension;
    if (pIJMag != 0) {
      rightAscension = atan2(pJ, pI);
    } else {
      rightAscension = atan2(pJDot, pIDot);
    }
    final rangeRate = p.dot(pDot) / pMag;
    final rightAscensionRate =
        (pIDot * pJ - pJDot * pI) / (-(pJ * pJ) - (pI * pI));
    final declinationRate = (pKDot - rangeRate * sin(declination)) / pIJMag;
    return RadecTopocentric(state.epoch, rightAscension % twoPi, declination,
        pMag, rightAscensionRate, declinationRate, rangeRate);
  }

  /// Observation epoch.
  final EpochUTC epoch;

  /// Slant range _(km)_.
  final double? range;

  /// Right-ascension _(rad)_.
  final double rightAscension;

  /// Declination _(rad)_.
  final double declination;

  /// Slant range rate _(km/s)_.
  final double? rangeRate;

  /// Right-ascension rate _(rad/s)_.
  final double? rightAscensionRate;

  /// Declination rate _(rad/s)_.
  final double? declinationRate;

  /// Right-ascension _(°)_.
  double get rightAscensionDegrees => rightAscension * rad2deg;

  /// Declination _(°)_.
  double get declinationDegrees => declination * rad2deg;

  /// Right-ascension rate _(°/s)_.
  double? get rightAscensionRateDegrees =>
      (rightAscensionRate != null) ? rightAscensionRate! * rad2deg : null;

  /// Declination rate _(°/s)_.
  double? get declinationRateDegrees =>
      (declinationRate != null) ? declinationRate! * rad2deg : null;

  /// Return the position relative to the observer [site].
  ///
  /// An optional [range] _(km)_ value can be passed to override the value
  /// contained in this observation.
  Vector position(final Geodetic site, [final double? range]) {
    final r = range ?? this.range ?? 1.0;
    final s0 = site.toITRF(epoch).toJ2000();
    return radecToPosition(rightAscension, declination, r).add(s0.position);
  }

  /// Return the velocity relative to the observer [site].
  ///
  /// An optional [range] _(km)_ and [rangeRate] _(km/s)_ value can be passed
  /// to override the values contained in this observation.
  Vector velocity(final Geodetic site,
      [final double? range, final double? rangeRate]) {
    if (rightAscensionRate == null || declinationRate == null) {
      throw 'Velocity unsolvable, missing ra/dec rates.';
    }
    final r = range ?? this.range ?? 1.0;
    final rd = rangeRate ?? this.rangeRate ?? 0.0;
    final s0 = site.toITRF(epoch).toJ2000();
    return radecToVelocity(rightAscension, declination, r, rightAscensionRate!,
            declinationRate!, rd)
        .add(s0.velocity);
  }

  /// Convert this observation into a line-of-sight vector.
  Vector lineOfSight() {
    final ca = cos(rightAscension);
    final cd = cos(declination);
    final sa = sin(rightAscension);
    final sd = sin(declination);
    final result = Float64List(3);
    result[0] = cd * ca;
    result[1] = cd * sa;
    result[2] = sd;
    return Vector(result);
  }

  /// Calculate the angular distance _(rad)_ between this and another
  /// [RadecTopocentric] object.
  double angle(final RadecTopocentric radec,
          [final AngularDistanceMethod method =
              AngularDistanceMethod.cosine]) =>
      angularDistance(
          rightAscension, declination, radec.rightAscension, radec.declination,
          method: method);

  /// Calculate the angular distance _(°)_ between this and another
  /// [RadecTopocentric] object.
  double angleDegrees(final RadecTopocentric radec,
          [final AngularDistanceMethod method =
              AngularDistanceMethod.cosine]) =>
      angle(radec, method) * rad2deg;
}
