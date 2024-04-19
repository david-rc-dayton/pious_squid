import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Sun metrics and operations.
class Sun {
  Sun._(); // disable constructor

  /// Sun gravitational parameter _(km³/s²)_.
  static const double mu = 1.32712428e11;

  /// Mean solar flux _(W/m²)_.
  static const double solarFlux = 1367.0;

  /// Mean solar pressure _(N/m²)_.
  static const double solarPressure = solarFlux / (speedOfLight * 1000);

  /// Sun umbra angle _(rad)_.
  static const double umbraAngle = 0.26411888 * deg2rad;

  /// Sun penumbra angle _(rad)_.
  static const double penumbraAngle = 0.26900424 * deg2rad;

  /// Sun radius _(km)_.
  static const double radius = 695500.0;

  /// Compute the Sun's distance from Earth, in kilometers, at the
  /// provided [epoch].
  static double distance(final EpochUTC epoch) {
    final mSun = 357.5291092 + 35999.05034 * epoch.toJulianCenturies();
    return 1.000140612 -
        0.016708617 * cos(mSun * deg2rad) -
        0.000139589 * cos(2.0 * mSun * deg2rad);
  }

  /// Calculate the Sun's ECI position _(km)_ for a given UTC [epoch].
  static Vector3D position(final EpochUTC epoch) {
    final jc = epoch.toTT().toJulianCenturies();
    final dtr = deg2rad;
    final lamSun = 280.460 + 36000.771 * jc;
    final mSun = 357.5291092 + 35999.05034 * jc;
    final lamEc = lamSun +
        1.914666471 * sin(mSun * dtr) +
        0.019994643 * sin(2.0 * mSun * dtr);
    final rMag = 1.000140612 -
        0.016708617 * cos(mSun * dtr) -
        0.000139589 * cos(2.0 * mSun * dtr);
    final obliq = 23.439291 - 0.0130042 * jc;
    final r = Vector3D(
        rMag * cos(lamEc * dtr),
        rMag * cos(obliq * dtr) * sin(lamEc * dtr),
        rMag * sin(obliq * dtr) * sin(lamEc * dtr));
    final rMOD = r.scale(astronomicalUnit);
    final p = Earth.precession(epoch);
    return rMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);
  }

  /// Calculate the Sun's apparent ECI position _(km)_ from Earth for a given
  /// UTC [epoch].
  static Vector3D positionApparent(final EpochUTC epoch) {
    final dSec = Sun.distance(epoch) / speedOfLight;
    return position(epoch.roll(-dSec));
  }

  /// Return `true` if the ECI satellite position [posSat] is in eclipse at the
  /// given UTC [epoch].
  static bool shadow(final EpochUTC epoch, final Vector3D posSat) {
    final posSun = positionApparent(epoch);
    var shadow = false;
    if (posSun.dot(posSat) < 0) {
      final angle = posSun.angle(posSat);
      final r = posSat.magnitude();
      final satHoriz = r * cos(angle);
      final satVert = r * sin(angle);
      final penVert = Earth.radiusEquator + tan(penumbraAngle) * satHoriz;
      if (satVert <= penVert) {
        shadow = true;
      }
    }
    return shadow;
  }

  /// Calculate eclipse angles given a satellite ECI position [satPos] _(km)_
  /// and Sun ECI position [sunPos] _(km)_.
  ///
  /// Returns a tuple containing the following three values:
  ///   - central body angle _(rad)_
  ///   - central body apparant radius _(rad)_
  ///   - sun apparant radius _(rad)_
  static (double bodyAngle, double bodyRadius, double sunRadius) eclipseAngles(
      final Vector3D satPos, final Vector3D sunPos) {
    final satSun = sunPos.subtract(satPos);
    final r = satPos.magnitude();
    return (
      // central body angle
      satSun.angle(satPos.negate()),
      // central body apparent radius
      asin(Earth.radiusEquator / r),
      // sun apparent radius
      asin(radius / satSun.magnitude())
    );
  }

  /// Calculate the lighting ratio given a satellite ECI position [satPos]
  /// _(km)_ and Sun ECI position [sunPos] _(km)_.
  ///
  /// Returns `1.0` if the satellite is fully illuminated and `0.0` when
  /// fully eclipsed.
  static double lightingRatio(final Vector3D satPos, final Vector3D sunPos) {
    final (sunSatAngle, aCent, aSun) = eclipseAngles(satPos, sunPos);
    if (sunSatAngle - aCent + aSun <= 1e-10) {
      return 0.0;
    } else if (sunSatAngle - aCent - aSun < -1e-10) {
      final ssa2 = sunSatAngle * sunSatAngle;
      final ssaInv = 1.0 / (2.0 * sunSatAngle);
      final ac2 = aCent * aCent;
      final as2 = aSun * aSun;
      final acAsDiff = ac2 - as2;
      final a1 = (ssa2 - acAsDiff) * ssaInv;
      final a2 = (ssa2 + acAsDiff) * ssaInv;
      final asr1 = a1 / aSun;
      final asr2 = as2 - a1 * a1;
      final acr1 = a2 / aCent;
      final acr2 = ac2 - a2 * a2;
      final p1 = as2 * acos(asr1) - a1 * sqrt(asr2);
      final p2 = ac2 * acos(acr1) - a2 * sqrt(acr2);
      return 1.0 - (p1 + p2) / (pi * as2);
    }
    return 1.0;
  }

  /// Calculate the Sun's angular diameter _(rad)_ from an ECI observer
  /// position [obsPos] and Sun position [sunPos] _(km)_.
  static double diameter(final Vector3D obsPos, final Vector3D sunPos) =>
      angularDiameter(radius * 2, obsPos.subtract(sunPos).magnitude(),
          AngularDiameterMethod.sphere);
}
