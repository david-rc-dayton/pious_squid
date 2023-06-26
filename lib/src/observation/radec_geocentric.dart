import 'dart:math';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/observation/observation_utils.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Geocentric right-ascension and declination.
class RadecGeocentric {
  /// Create a new [RadecGeocentric] object.
  RadecGeocentric(this.epoch, this.rightAscension, this.declination,
      [this.range,
      this.rightAscensionRate,
      this.declinationRate,
      this.rangeRate]);

  /// Create a new [RadecGeocentric] object, using degrees for the
  /// angular values.
  factory RadecGeocentric.fromDegrees(final EpochUTC epoch,
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
    return RadecGeocentric(
      epoch,
      rightAscensionDegrees * deg2rad,
      declinationDegrees * deg2rad,
      range,
      rightAscensionRate,
      declinationRate,
      rangeRate,
    );
  }

  /// Create a [RadecGeocentric] object from an inertial [state] vector.
  factory RadecGeocentric.fromStateVector(final J2000 state) {
    final rI = state.position.x;
    final rJ = state.position.y;
    final rK = state.position.z;
    final vI = state.velocity.x;
    final vJ = state.velocity.y;
    final vK = state.velocity.z;
    final rMag = state.position.magnitude();
    final declination = asin(rK / rMag);
    final rIJMag = sqrt(rI * rI + rJ * rJ);
    double rightAscension;
    if (rIJMag != 0) {
      rightAscension = atan2(rJ, rI);
    } else {
      rightAscension = atan2(vJ, vI);
    }
    final rangeRate = state.position.dot(state.velocity) / rMag;
    final rightAscensionRate = (vI * rJ - vJ * rI) / (-(rJ * rJ) - (rI * rI));
    final declinationRate = (vK - rangeRate * (rK / rMag)) / rIJMag;
    return RadecGeocentric(state.epoch, rightAscension % twoPi, declination,
        rMag, rightAscensionRate, declinationRate, rangeRate);
  }

  /// Observation epoch.
  final EpochUTC epoch;

  /// Right-ascension _(rad)_.
  final double rightAscension;

  /// Declination _(rad)_.
  final double declination;

  /// Slant range _(km)_.
  final double? range;

  /// Right-ascension rate _(rad/s)_.
  final double? rightAscensionRate;

  /// Declination rate _(rad/s)_.
  final double? declinationRate;

  /// Slant range rate _(km/s)_.
  final double? rangeRate;

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

  /// Return the position relative to the center of the Earth.
  ///
  /// An optional [range] _(km)_ value can be passed to override the value
  /// contained in this observation.
  Vector3D position([final double? range]) {
    final r = range ?? this.range ?? 1.0;
    return radecToPosition(rightAscension, declination, r);
  }

  /// Return the velocity relative to the centar of the Earth.
  ///
  /// An optional [range] _(km)_ and [rangeRate] _(km/s)_ value can be passed
  /// to override the value contained in this observation.
  Vector3D velocity([final double? range, final double? rangeRate]) {
    if (rightAscensionRate == null || declinationRate == null) {
      throw 'Velocity unsolvable, missing ra/dec rates.';
    }
    final r = range ?? this.range ?? 1.0;
    final rd = rangeRate ?? this.rangeRate ?? 0.0;
    return radecToVelocity(rightAscension, declination, r, rightAscensionRate!,
        declinationRate!, rd);
  }

  /// Calculate the angular distance _(rad)_ between this and another
  /// [RadecGeocentric] object.
  double angle(final RadecGeocentric radec,
          [final AngularDistanceMethod method =
              AngularDistanceMethod.cosine]) =>
      angularDistance(
          rightAscension, declination, radec.rightAscension, radec.declination,
          method: method);

  /// Calculate the angular distance _(°)_ between this and another
  /// [RadecGeocentric] object.
  double angleDegrees(final RadecGeocentric radec,
          [final AngularDistanceMethod method =
              AngularDistanceMethod.cosine]) =>
      angle(radec, method) * rad2deg;
}
