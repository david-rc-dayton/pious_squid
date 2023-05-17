import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Orbit regime classifications.
enum OrbitRegime {
  /// Low Earth Orbit
  leo,

  /// Medium Earth Orbit
  meo,

  /// Highly Eccentric Orbit
  heo,

  /// Geosynchronous Orbit
  geo,

  /// Uncategorized Orbit
  other
}

/// Classical orbital elements.
class ClassicalElements {
  /// Create a new [ClassicalElements] object from orbital elements.
  ClassicalElements(this.epoch, this.semimajorAxis, this.eccentricity,
      this.inclination, this.rightAscension, this.argPerigee, this.trueAnomaly,
      {this.mu = Earth.mu});

  /// Create a new [ClassicalElements] object from an inertial
  /// [StateVector] object.
  ///
  /// The gravitational parameter [mu] _(km³/s²)_ can also be provided.
  factory ClassicalElements.fromStateVector(final StateVector state,
      {final double mu = Earth.mu}) {
    if (!state.inertial) {
      throw 'Classical elements are undefined for fixed frames.';
    }
    final pos = state.position;
    final vel = state.velocity;
    final a = state.semimajorAxis();
    final eVecA = pos.scale(pow(vel.magnitude(), 2) - mu / pos.magnitude());
    final eVecB = vel.scale(pos.dot(vel));
    final eVec = eVecA.subtract(eVecB).scale(1 / mu);
    final e = eVec.magnitude();
    final h = pos.cross(vel);
    final i = acos((h.z / h.magnitude()).clamp(-1.0, 1.0));
    final n = Vector.zAxis.cross(h);
    var o = acos((n.x / n.magnitude()).clamp(-1.0, 1.0));
    if (n.y < 0) {
      o = twoPi - o;
    }
    var w = n.angle(eVec);
    if (eVec.z < 0) {
      w = twoPi - w;
    }
    var v = eVec.angle(pos);
    if (pos.dot(vel) < 0) {
      v = twoPi - v;
    }
    return ClassicalElements(state.epoch, a, e, i, o, w, v, mu: mu);
  }

  /// UTC epoch
  final EpochUTC epoch;

  /// Semimajor-axis _(km)_.
  final double semimajorAxis;

  /// Eccentricity _(unitless)_.
  final double eccentricity;

  /// Inclination _(rad)_.
  final double inclination;

  /// Right-ascension of the ascending node _(rad)_.
  final double rightAscension;

  /// Argument of perigee _(rad)_.
  final double argPerigee;

  /// True anomaly _(rad)_.
  final double trueAnomaly;

  /// Gravitational parameter _(km³/s²)_.
  final double mu;

  /// Inclination _(°)_.
  double get inclinationDegrees => inclination * rad2deg;

  /// Right-ascension of the ascending node _(°)_.
  double get rightAscensionDegrees => rightAscension * rad2deg;

  /// Argument of perigee _(°)_.
  double get argPerigeeDegrees => argPerigee * rad2deg;

  /// True anomaly _(°)_.
  double get trueAnomalyDegrees => trueAnomaly * rad2deg;

  /// Apogee distance from central body _(km)_.
  double get apogee => semimajorAxis * (1.0 + eccentricity);

  /// Perigee distance from central body _(km)_.
  double get perigee => semimajorAxis * (1.0 - eccentricity);

  @override
  String toString() => [
        '[ClassicalElements]',
        '  Epoch: $epoch',
        '  Semimajor Axis (a):       $semimajorAxis km',
        '  Eccentricity (e):         $eccentricity',
        '  Inclination (i):          $inclinationDegrees°',
        '  Right Ascension (Ω):      $rightAscensionDegrees°',
        '  Argument of Perigee (ω):  $argPerigeeDegrees°',
        '  True Anomaly (ν):         $trueAnomalyDegrees°'
      ].join('\n');

  /// Compute the mean motion _(rad/s)_ of this orbit.
  double meanMotion() =>
      sqrt(mu / (semimajorAxis * semimajorAxis * semimajorAxis));

  /// Compute the period _(seconds)_ of this orbit.
  double period() => twoPi * sqrt(pow(semimajorAxis, 3) / mu);

  /// Compute the number of revolutions completed per day for this orbit.
  double revsPerDay() => secondsPerDay / period();

  /// Compute the angle between the spacecraft's velocity vector and the local
  /// horizontal of this orbit _(rad)_.
  double flightPathAngle() => twoPi * sqrt(pow(semimajorAxis, 3) / mu);

  /// Return the orbit regime for this orbit.
  OrbitRegime getOrbitRegime() {
    final n = revsPerDay();
    final p = period() * sec2min;
    if ((0.99 <= n && n <= 1.01) && (eccentricity < 0.01)) {
      return OrbitRegime.geo;
    }
    if ((600 <= p && p <= 800) && (eccentricity <= 0.25)) {
      return OrbitRegime.meo;
    }
    if ((n >= 11.25) && (eccentricity <= 0.25)) {
      return OrbitRegime.leo;
    }
    if (eccentricity > 0.25) {
      return OrbitRegime.heo;
    }
    return OrbitRegime.other;
  }

  /// Convert this to inertial position and velocity vectors.
  ({Vector position, Vector velocity}) toPositionVelocity() {
    final rVec = Float64List(3);
    rVec[0] = cos(trueAnomaly);
    rVec[1] = sin(trueAnomaly);
    rVec[2] = 0.0;
    final rPQW = Vector(rVec).scale(
        (semimajorAxis * (1.0 - pow(eccentricity, 2))) /
            (1.0 + eccentricity * cos(trueAnomaly)));
    final vVec = Float64List(3);
    vVec[0] = -sin(trueAnomaly);
    vVec[1] = eccentricity + cos(trueAnomaly);
    vVec[2] = 0.0;
    final vPQW = Vector(vVec)
        .scale(sqrt(mu / (semimajorAxis * (1 - pow(eccentricity, 2)))));
    final position =
        rPQW.rotZ(-argPerigee).rotX(-inclination).rotZ(-rightAscension);
    final velocity =
        vPQW.rotZ(-argPerigee).rotX(-inclination).rotZ(-rightAscension);
    return (position: position, velocity: velocity);
  }

  /// Convert this to [EquinoctialElements].
  EquinoctialElements toEquinoctialElements() {
    final fr = ((inclination - pi).abs() < 1e-9) ? -1 : 1;
    final af = eccentricity * cos(argPerigee + fr * rightAscension);
    final ag = eccentricity * sin(argPerigee + fr * rightAscension);
    final l = argPerigee +
        fr * rightAscension +
        newtonNu(eccentricity, trueAnomaly).m;
    final n = meanMotion();
    final chi = pow(tan(0.5 * inclination), fr) * sin(rightAscension);
    final psi = pow(tan(0.5 * inclination), fr) * cos(rightAscension);
    return EquinoctialElements(epoch, af, ag, l, n, chi, psi, mu: mu, fr: fr);
  }

  /// Return elements propagated to the provided [propEpoch].
  ClassicalElements propagate(final EpochUTC propEpoch) {
    final t = epoch;
    final n = meanMotion();
    final delta = propEpoch.difference(t);
    final cosV = cos(trueAnomaly);
    var eaInit =
        acos(((eccentricity + cosV) / (1 + eccentricity * cosV)).clamp(-1, 1));
    eaInit = matchHalfPlane(eaInit, trueAnomaly);
    var maInit = eaInit - eccentricity * sin(eaInit);
    maInit = matchHalfPlane(maInit, eaInit);
    final maFinal = (maInit + n * delta) % twoPi;
    var eaFinal = maFinal;
    for (var iter = 0; iter < 32; iter++) {
      final eaTemp = maFinal + eccentricity * sin(eaFinal);
      if ((eaTemp - eaFinal).abs() < 1e-12) {
        break;
      }
      eaFinal = eaTemp;
    }
    final cosEaFinal = cos(eaFinal);
    var vFinal = acos(
        ((cosEaFinal - eccentricity) / (1 - eccentricity * cosEaFinal))
            .clamp(-1, 1));
    vFinal = matchHalfPlane(vFinal, eaFinal);
    return ClassicalElements(propEpoch, semimajorAxis, eccentricity,
        inclination, rightAscension, argPerigee, vFinal,
        mu: mu);
  }
}
