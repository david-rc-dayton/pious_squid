import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Moon metrics and operations.
class Moon {
  Moon._(); // disable constructor

  /// Moon gravitational parameter _(km³/s²)_.
  static double mu = 4902.799;

  /// Moon equatorial radius _(km)_.
  static double radiusEquator = 1738.0;

  /// Calculate the Moon's ECI position _(km)_ for a given UTC [epoch].
  static Vector3D position(final EpochUTC epoch) {
    final jc = epoch.toTDB().toJulianCenturies();
    final dtr = deg2rad;
    final lamEcl = 218.32 +
        481267.8813 * jc +
        6.29 * sin((134.9 + 477198.85 * jc) * dtr) -
        1.27 * sin((259.2 - 413335.38 * jc) * dtr) +
        0.66 * sin((235.7 + 890534.23 * jc) * dtr) +
        0.21 * sin((269.9 + 954397.70 * jc) * dtr) -
        0.19 * sin((357.5 + 35999.05 * jc) * dtr) -
        0.11 * sin((186.6 + 966404.05 * jc) * dtr);
    final phiEcl = 5.13 * sin((93.3 + 483202.03 * jc) * dtr) +
        0.28 * sin((228.2 + 960400.87 * jc) * dtr) -
        0.28 * sin((318.3 + 6003.18 * jc) * dtr) -
        0.17 * sin((217.6 - 407332.20 * jc) * dtr);
    final pllx = 0.9508 +
        0.0518 * cos((134.9 + 477198.85 * jc) * dtr) +
        0.0095 * cos((259.2 - 413335.38 * jc) * dtr) +
        0.0078 * cos((235.7 + 890534.23 * jc) * dtr) +
        0.0028 * cos((269.9 + 954397.7 * jc) * dtr);
    final obq = 23.439291 - 0.0130042 * jc;
    final rMag = 1.0 / sin(pllx * dtr);
    final r = Vector3D(
        rMag * (cos(phiEcl * dtr) * cos(lamEcl * dtr)),
        rMag *
            (cos(obq * dtr) * cos(phiEcl * dtr) * sin(lamEcl * dtr) -
                sin(obq * dtr) * sin(phiEcl * dtr)),
        rMag *
            (sin(obq * dtr) * cos(phiEcl * dtr) * sin(lamEcl * dtr) +
                cos(obq * dtr) * sin(phiEcl * dtr)));
    return r.scale(Earth.radiusEquator);
  }

  /// Calculate lunar illumination at a given UTC [epoch] and optional ECI
  /// observer position vector _(km)_.
  ///
  /// Result is `1.0` when the Moon fully illuminated and `0.0` is fully
  /// in shadow.
  static double illumination(final EpochUTC epoch, [final Vector3D? origin]) {
    final orig = origin ?? Vector3D.origin;
    final sunPos = Sun.position(epoch).subtract(orig);
    final moonPos = position(epoch).subtract(orig);
    final phaseAngle = sunPos.angle(moonPos);
    return 0.5 * (1 - cos(phaseAngle));
  }

  /// Calculate the Moon's angular diameter _(rad)_ from an ECI observer
  /// position [obsPos] and Moon position [moonPos] _(km)_.
  static double diameter(final Vector3D obsPos, final Vector3D moonPos) =>
      angularDiameter(radiusEquator * 2, obsPos.subtract(moonPos).magnitude(),
          AngularDiameterMethod.sphere);
}
