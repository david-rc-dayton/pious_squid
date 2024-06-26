import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/data/data_handler.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Earth precession `zeta` polynomial coefficients.
final Float64List _zetaPoly = Float64List.fromList([
  0.017998 * asec2rad,
  0.30188 * asec2rad,
  2306.2181 * asec2rad,
  0 * asec2rad,
]);

/// Earth precession `theta` polynomial coefficients.
final Float64List _thetaPoly = Float64List.fromList([
  -0.041833 * asec2rad,
  -0.42665 * asec2rad,
  2004.3109 * asec2rad,
  0 * asec2rad,
]);

/// Earth precession `zed` polynomial coefficients.
final Float64List _zedPoly = Float64List.fromList([
  0.018203 * asec2rad,
  1.09468 * asec2rad,
  2306.2181 * asec2rad,
  0 * asec2rad,
]);

/// Earth nutation Moon anomaly polynomial coefficients.
final Float64List _moonAnomPoly = Float64List.fromList([
  1.78e-5 * deg2rad,
  0.0086972 * deg2rad,
  (1325 * 360 + 198.8673981) * deg2rad,
  134.96298139 * deg2rad,
]);

/// Earth nutation Sun anomaly polynomial coefficients.
final Float64List _sunAnomPoly = Float64List.fromList([
  -3.3e-6 * deg2rad,
  -0.0001603 * deg2rad,
  (99 * 360 + 359.0503400) * deg2rad,
  357.52772333 * deg2rad,
]);

/// Earth nutation Moon latitude polynomial coefficients.
final Float64List _moonLatPoly = Float64List.fromList([
  3.1e-6 * deg2rad,
  -0.0036825 * deg2rad,
  (1342 * 360 + 82.0175381) * deg2rad,
  93.27191028 * deg2rad,
]);

/// Earth nutation Sun elongation polynomial coefficients.
final Float64List _sunElongPoly = Float64List.fromList([
  5.3e-6 * deg2rad,
  -0.0019142 * deg2rad,
  (1236 * 360 + 307.1114800) * deg2rad,
  297.85036306 * deg2rad,
]);

/// Earth nutation Moon right ascension polynomial coefficients.
final Float64List _moonRaanPoly = Float64List.fromList([
  2.2e-6 * deg2rad,
  0.0020708 * deg2rad,
  -(5 * 360 + 134.1362608) * deg2rad,
  125.04452222 * deg2rad,
]);

/// Earth nutation mean epsilon polynomial coefficients.
final Float64List _meanEpsilonPoly = Float64List.fromList([
  5.04e-7 * deg2rad,
  -1.64e-7 * deg2rad,
  -0.0130042 * deg2rad,
  23.439291 * deg2rad,
]);

/// Earth precession angles _(rad)_.
typedef PrecessionAngles = ({double zeta, double theta, double zed});

/// Earth nutation angles _(rad)_.
typedef NutationAngles = ({
  double dPsi,
  double dEps,
  double mEps,
  double eps,
  double eqEq,
  double gast
});

/// Earth metrics and operations.
class Earth {
  Earth._(); // disable constructor

  /// Earth gravitational parameter _(km²/s³)_.
  static const double mu = 398600.4415;

  /// Earth equatorial radius _(km)_.
  static const double radiusEquator = 6378.137;

  /// Earth coefficient of flattening _(unitless)_.
  static const double flattening = 1.0 / 298.25642;

  /// Earth polar radius _(km)_.
  static const double radiusPolar = radiusEquator * (1.0 - flattening);

  /// Earth mean radius _(km)_.
  static const double radiusMean = (2.0 * radiusEquator + radiusPolar) / 3.0;

  /// Earth eccentricity squared _(unitless)_.
  static const double eccentricitySquared = flattening * (2.0 - flattening);

  /// Earth J2 effect coefficient _(unitless)_.
  static const double j2 = 1.08262668355315e-3;

  /// Earth J3 effect coefficient _(unitless)_.
  static const double j3 = -2.53265648533224e-6;

  /// Earth J4 effect coefficient _(unitless)_.
  static const double j4 = -1.619621591367e-6;

  /// Earth J5 effect coefficient _(unitless)_.
  static const double j5 = -2.27296082868698e-7;

  /// Earth J6 effect coefficient _(unitless)_.
  static const double j6 = 5.40681239107085e-7;

  /// Earth rotation vector _(rad/s)_.
  static final Vector3D rotation = Vector3D(0, 0, 7.292115146706979e-5);

  /// Earth rotation vector _(rad/s)_, factoring length of day changes if
  /// Earth Orientation Parameters are available.
  static Vector3D rotationLod(final EpochUTC epoch) {
    final eop = DataHandler().getEop(epoch);
    return rotation.scale(1.0 - (eop.lod / secondsPerDay));
  }

  /// Calculate mean motion _(rad/s)_ from a given [semimajorAxis] _(km)_.
  static double smaToMeanMotion(final double semimajorAxis) =>
      sqrt(mu / (semimajorAxis * semimajorAxis * semimajorAxis));

  /// Calculate semimajor-axis _(km)_ from mean motion _(rad/s)_.
  static double meanMotionToSma(final double n) =>
      pow(Earth.mu / (n * n), 1.0 / 3.0).toDouble();

  /// Calculate semimajor-axis _(km)_ from a given number of revolutions per
  /// day [rpd].
  static double revsPerDayToSma(final double rpd) =>
      pow(mu, 1.0 / 3.0) / pow((twoPi * rpd) / secondsPerDay, 2.0 / 3.0);

  /// Calculate Earth [PrecessionAngles] at a given UTC [epoch].
  static PrecessionAngles precession(final EpochUTC epoch) {
    final t = epoch.toTT().toJulianCenturies();
    final zeta = evalPoly(t, _zetaPoly);
    final theta = evalPoly(t, _thetaPoly);
    final zed = evalPoly(t, _zedPoly);
    return (zeta: zeta, theta: theta, zed: zed);
  }

  /// Calculate Earth [NutationAngles] for a given UTC [epoch].
  static NutationAngles nutation(final EpochUTC epoch,
      {final int coeffs = 4, final bool useEop = false}) {
    final t = epoch.toTT().toJulianCenturies();
    final moonAnom = evalPoly(t, _moonAnomPoly);
    final sunAnom = evalPoly(t, _sunAnomPoly);
    final moonLat = evalPoly(t, _moonLatPoly);
    final sunElong = evalPoly(t, _sunElongPoly);
    final moonRaan = evalPoly(t, _moonRaanPoly);
    var deltaPsi = 0.0;
    var deltaEpsilon = 0.0;
    final dh = DataHandler();
    for (var i = 0; i < coeffs; i++) {
      final (a1, a2, a3, a4, a5, ai, bi, ci, di) = dh.getIau1980Coeffs(i);
      final arg = a1 * moonAnom +
          a2 * sunAnom +
          a3 * moonLat +
          a4 * sunElong +
          a5 * moonRaan;
      final sinC = ai + bi * t;
      final cosC = ci + di * t;
      deltaPsi += sinC * sin(arg);
      deltaEpsilon += cosC * cos(arg);
    }
    deltaPsi *= ttasec2rad;
    deltaEpsilon *= ttasec2rad;
    if (useEop) {
      final eop = DataHandler().getEop(epoch);
      deltaPsi += eop.dpsi;
      deltaEpsilon += eop.deps;
    }
    final meanEpsilon = evalPoly(t, _meanEpsilonPoly);
    final epsilon = meanEpsilon + deltaEpsilon;
    final eqEq = deltaPsi * cos(meanEpsilon) +
        (0.00264 * asec2rad) * sin(moonRaan) +
        (0.000063 * asec2rad) * sin(2.0 * moonRaan);
    final gast = epoch.gmstAngle() + eqEq;
    return (
      dPsi: deltaPsi,
      dEps: deltaEpsilon,
      mEps: meanEpsilon,
      eps: epsilon,
      eqEq: eqEq,
      gast: gast
    );
  }

  /// Convert a [semimajorAxis] _(km)_ to an eastward drift rate _(rad/day)_.
  static double smaToDrift(final double semimajorAxis) {
    final t = twoPi * sqrt(pow(semimajorAxis, 3) / mu) / secondsPerSiderealDay;
    return (1.0 - t) * twoPi;
  }

  /// Convert a [semimajorAxis] _(km)_ to an eastward drift rate _(°/day)_.
  static double smaToDriftDegrees(final double semimajorAxis) =>
      smaToDrift(semimajorAxis) * rad2deg;

  /// Convert an eastward [driftRate] _(rad/day)_ to a semimajor-axis _(km)_.
  static double driftToSemimajorAxis(final double driftRate) {
    final t = (-driftRate / twoPi + 1) * secondsPerSiderealDay;
    return pow((mu * t * t) / (4 * pi * pi), 1 / 3) as double;
  }

  /// Convert an eastward [driftRate] _(°/day)_ to a semimajor-axis _(km)_.
  static double driftDegreesToSma(final double driftRate) =>
      driftToSemimajorAxis(deg2rad * driftRate);

  /// Calculate the Earth's angular diameter _(rad)_ from an ECI satellite
  /// position [satPos] _(km)_.
  static double diameter(final Vector3D satPos) => angularDiameter(
      radiusEquator * 2, satPos.magnitude(), AngularDiameterMethod.sphere);
}
