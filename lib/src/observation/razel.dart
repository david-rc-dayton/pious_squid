import 'dart:math';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Range, azimuth, and elevation.
class Razel {
  /// Create a new [Razel] object.
  Razel(this.epoch, this.range, this.azimuth, this.elevation,
      [this.rangeRate, this.azimuthRate, this.elevationRate]);

  /// Create a new [Razel] object, using degrees for the angular values.
  factory Razel.fromDegrees(final EpochUTC epoch, final double range,
      final double azimuthDegrees, final double elevationDegrees,
      [final double? rangeRate,
      final double? azimuthRateDegrees,
      final double? elevationRateDegrees]) {
    final azimuthRate =
        (azimuthRateDegrees != null) ? azimuthRateDegrees * deg2rad : null;
    final elevationRate =
        (elevationRateDegrees != null) ? elevationRateDegrees * deg2rad : null;
    return Razel(epoch, range, azimuthDegrees * deg2rad,
        elevationDegrees * deg2rad, rangeRate, azimuthRate, elevationRate);
  }

  /// Create a [Razel] object from an inertial [state] and
  /// [site] vector.
  factory Razel.fromStateVectors(final J2000 state, final J2000 site) {
    final stateEcef = state.toITRF();
    final siteEcef = site.toITRF();
    final po2 = halfPi;
    final r = stateEcef.position.subtract(siteEcef.position);
    final rDot = stateEcef.velocity;
    final geo = siteEcef.toGeodetic();
    final p = r.rotZ(geo.longitude).rotY(po2 - geo.latitude);
    final pDot = rDot.rotZ(geo.longitude).rotY(po2 - geo.latitude);
    final pS = p.x;
    final pE = p.y;
    final pZ = p.z;
    final pSDot = pDot.x;
    final pEDot = pDot.y;
    final pZDot = pDot.z;
    final pMag = p.magnitude();
    final pSEMag = sqrt(pS * pS + pE * pE);
    final elevation = asin(pZ / pMag);
    double azimuth;
    if (elevation != po2) {
      azimuth = atan2(-pE, pS) + pi;
    } else {
      azimuth = atan2(-pEDot, pSDot) + pi;
    }
    final rangeRate = p.dot(pDot) / pMag;
    final azimuthRate = (pSDot * pE - pEDot * pS) / (pS * pS + pE * pE);
    final elevationRate = (pZDot - rangeRate * sin(elevation)) / pSEMag;
    return Razel(state.epoch, pMag, azimuth % twoPi, elevation, rangeRate,
        azimuthRate, elevationRate);
  }

  /// Observation epoch.
  final EpochUTC epoch;

  /// Slant range _(km)_.
  final double range;

  /// Azimuth _(rad)_.
  final double azimuth;

  /// Elevation _(rad)_.
  final double elevation;

  /// Slant range rate _(km/s)_.
  final double? rangeRate;

  /// Azimuth rate _(rad/s)_.
  final double? azimuthRate;

  /// Elevation rate _(rad/s)_.
  final double? elevationRate;

  /// Azimuth _(°)_.
  double get azimuthDegrees => azimuth * rad2deg;

  /// Elevation _(°)_.
  double get elevationDegrees => elevation * rad2deg;

  /// Azimuth rate _(°/s)_.
  double? get azimuthRateDegrees =>
      (azimuthRate != null) ? azimuthRate! * rad2deg : null;

  /// Elevation rate _(°/s)_.
  double? get elevationRateDegrees =>
      (elevationRate != null) ? elevationRate! * rad2deg : null;

  /// Return the position relative to the observer [site].
  ///
  /// An optional azimuth [az] _(rad)_ and elevation [el] _(rad)_ value can be
  /// passed to override the values contained in this observation.
  Vector3D position(final J2000 site, [final double? az, final double? el]) {
    final ecef = site.toITRF();
    final geo = ecef.toGeodetic();
    final po2 = halfPi;
    final newAz = az ?? azimuth;
    final newEl = el ?? elevation;
    final sAz = sin(newAz);
    final cAz = cos(newAz);
    final sEl = sin(newEl);
    final cEl = cos(newEl);
    final pSez = Vector3D(-range * cEl * cAz, range * cEl * sAz, range * sEl);
    final rEcef = pSez
        .rotY(-(po2 - geo.latitude))
        .rotZ(-geo.longitude)
        .add(ecef.position);
    return ITRF(epoch, rEcef, Vector3D.origin).toJ2000().position;
  }

  /// Convert this observation into a [J2000] state vector.
  ///
  /// This will throw an error if the [rangeRate], [elevationRate], or
  /// [azimuthRate] are not defined.
  J2000 toStateVector(final J2000 site) {
    if (rangeRate == null || elevationRate == null || azimuthRate == null) {
      throw 'Cannot create state, required values are undefined.';
    }
    final ecef = site.toITRF();
    final geo = ecef.toGeodetic();
    final po2 = halfPi;
    final sAz = sin(azimuth);
    final cAz = cos(azimuth);
    final sEl = sin(elevation);
    final cEl = cos(elevation);
    final pSez = Vector3D(-range * cEl * cAz, range * cEl * sAz, range * sEl);
    final pDotSez = Vector3D(
        -rangeRate! * cEl * cAz +
            range * sEl * cAz * elevationRate! +
            range * cEl * sAz * azimuthRate!,
        rangeRate! * cEl * sAz -
            range * sEl * sAz * elevationRate! +
            range * cEl * cAz * azimuthRate!,
        rangeRate! * sEl + range * cEl * elevationRate!);
    final pEcef = pSez.rotY(-(po2 - geo.latitude)).rotZ(-geo.longitude);
    final pDotEcef = pDotSez.rotY(-(po2 - geo.latitude)).rotZ(-geo.longitude);
    final rEcef = pEcef.add(ecef.position);
    return ITRF(epoch, rEcef, pDotEcef).toJ2000();
  }

  /// Calculate the angular distance _(rad)_ between this and another
  /// [Razel] object.
  double angle(final Razel razel,
          [final AngularDistanceMethod method =
              AngularDistanceMethod.cosine]) =>
      angularDistance(azimuth, elevation, razel.azimuth, razel.elevation,
          method: method);

  /// Calculate the angular distance _(°)_ between this and another
  /// [Razel] object.
  double angleDegrees(final Razel razel,
          [final AngularDistanceMethod method =
              AngularDistanceMethod.cosine]) =>
      angle(razel, method) * rad2deg;
}
