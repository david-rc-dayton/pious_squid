import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/data/data_handler.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Earth precession `zeta` polynomial coefficients.
final Float64List _zetaPoly = Float64List.fromList(
    [0.017998 * asec2rad, 0.30188 * asec2rad, 2306.2181 * asec2rad, 0.0]);

/// Earth precession `theta` polynomial coefficients.
final Float64List _thetaPoly = Float64List.fromList(
    [-0.041833 * asec2rad, -0.42665 * asec2rad, 2004.3109 * asec2rad, 0.0]);

/// Earth precession `zed` polynomial coefficients.
final Float64List _zedPoly = Float64List.fromList(
    [0.018203 * asec2rad, 1.09468 * asec2rad, 2306.2181 * asec2rad, 0]);

/// Earth nutation Moon anomaly polynomial coefficients.
final Float64List _moonAnomPoly = Float64List.fromList([
  1.4343e-5 * deg2rad,
  0.0088553 * deg2rad,
  (1325.0 * 360.0 + 198.8675605) * deg2rad,
  134.96340251 * deg2rad,
]);

/// Earth nutation Sun anomaly polynomial coefficients.
final Float64List _sunAnomPoly = Float64List.fromList([
  3.8e-8 * deg2rad,
  -0.0001537 * deg2rad,
  (99.0 * 360.0 + 359.0502911) * deg2rad,
  357.52910918 * deg2rad
]);

/// Earth nutation Moon latitude polynomial coefficients.
final Float64List _moonLatPoly = Float64List.fromList([
  -2.88e-7 * deg2rad,
  -0.0035420 * deg2rad,
  (1342.0 * 360.0 + 82.0174577) * deg2rad,
  93.27209062 * deg2rad
]);

/// Earth nutation Sun elongation polynomial coefficients.
final Float64List _sunElongPoly = Float64List.fromList([
  1.831e-6 * deg2rad,
  -0.0017696 * deg2rad,
  (1236.0 * 360.0 + 307.1114469) * deg2rad,
  297.85019547 * deg2rad
]);

/// Earth nutation Moon right-ascension polynomial coefficients.
final Float64List _moonRaanPoly = Float64List.fromList([
  2.139e-6 * deg2rad,
  0.0020756 * deg2rad,
  -(5.0 * 360.0 + 134.1361851) * deg2rad,
  125.04455501 * deg2rad
]);

/// Earth nutation mean epsilon polynomial coefficients.
final Float64List _meanEpsilonPoly = Float64List.fromList([
  0.001813 * asec2rad,
  -0.00059 * asec2rad,
  -46.8150 * asec2rad,
  84381.448 * asec2rad,
]);

/// Container for Earth precession values.
class PrecessionAngles {
  /// Create a new [PrecessionAngles] object, given the precession
  /// angles _(rad)_.
  PrecessionAngles(this.zeta, this.theta, this.zed);

  /// Zeta angle _(rad)_.
  final double zeta;

  /// Theta angle _(rad)_.
  final double theta;

  /// Zed angle _(rad)_.
  final double zed;
}

/// Container for Earth nutation values.
class NutationAngles {
  /// Create a new [NutationAngles] object, given the nutation angles _(rad)_.
  NutationAngles(
      this.dPsi, this.dEps, this.mEps, this.eps, this.eqEq, this.gast);

  /// Delta-psi angle _(rad)_.
  final double dPsi;

  /// Delta-epsilon angle _(rad)_.
  final double dEps;

  /// Mean-epsilon angle _(rad)_.
  final double mEps;

  /// Epsilon angle _(rad)_.
  final double eps;

  /// Equation of equinoxes angle _(rad)_.
  final double eqEq;

  /// Greenwich apparent sidereal time angle _(rad)_.
  final double gast;
}

/// Earth metrics and operations.
class Earth {
  Earth._(); // disable constructor

  /// Earth gravitational parameter _(km²/s³)_.
  static const double mu = 398600.4415;

  /// Earth equatorial radius _(km)_.
  static const double radiusEquator = 6378.1363;

  /// Earth coefficient of flattening _(unitless)_.
  static const double flattening = 1.0 / 298.257223563;

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
  static final Vector rotation =
      Vector(Float64List.fromList([0, 0, 7.292115146706979e-5]));

  /// Calculate mean motion _(rad/s)_ from a given [semimajorAxis] _(km)_.
  static double smaToMeanMotion(final double semimajorAxis) =>
      sqrt(mu / (semimajorAxis * semimajorAxis * semimajorAxis));

  /// Calculate semimajor-axis _(km)_ from a given number of revolutions per
  /// day [rpd].
  static double revsPerDayToSma(final double rpd) =>
      pow(mu, 1 / 3) / pow((twoPi * rpd) / secondsPerDay, 2 / 3);

  /// Calculate Earth [PrecessionAngles] at a given UTC [epoch].
  static PrecessionAngles precession(final EpochUTC epoch) {
    final t = epoch.toTT().toJulianCenturies();
    final zeta = evalPoly(t, _zetaPoly);
    final theta = evalPoly(t, _thetaPoly);
    final zed = evalPoly(t, _zedPoly);
    return PrecessionAngles(zeta, theta, zed);
  }

  /// Calculate Earth [NutationAngles] for a given UTC [epoch].
  static NutationAngles nutation(final EpochUTC epoch) {
    final t = epoch.toTT().toJulianCenturies();
    final moonAnom = evalPoly(t, _moonAnomPoly);
    final sunAnom = evalPoly(t, _sunAnomPoly);
    final moonLat = evalPoly(t, _moonLatPoly);
    final sunElong = evalPoly(t, _sunElongPoly);
    final moonRaan = evalPoly(t, _moonRaanPoly);
    var deltaPsi = 0.0;
    var deltaEpsilon = 0.0;
    final dh = DataHandler();
    for (var i = 0; i < 4; i++) {
      final coeffs = dh.getIau1980Coeffs(i);
      final arg = coeffs.a1 * moonAnom +
          coeffs.a2 * sunAnom +
          coeffs.a3 * moonLat +
          coeffs.a4 * sunElong +
          coeffs.a5 * moonRaan;
      final sinC = coeffs.ai + coeffs.bi * t;
      final cosC = coeffs.ci + coeffs.di * t;
      deltaPsi += sinC * sin(arg);
      deltaEpsilon += cosC * cos(arg);
    }
    deltaPsi *= ttasec2rad;
    deltaEpsilon *= ttasec2rad;
    final meanEpsilon = evalPoly(t, _meanEpsilonPoly);
    final epsilon = meanEpsilon + deltaEpsilon;
    final eqEq = deltaPsi * cos(meanEpsilon) +
        (0.00264 * asec2rad) * sin(moonRaan) +
        (0.000063 * asec2rad) * sin(2.0 * moonRaan);
    final gast = epoch.gmstAngle() + eqEq;
    return NutationAngles(
        deltaPsi, deltaEpsilon, meanEpsilon, epsilon, eqEq, gast);
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
  static double diameter(final Vector satPos) => angularDiameter(
      radiusEquator * 2, satPos.magnitude(), AngularDiameterMethod.sphere);
}
