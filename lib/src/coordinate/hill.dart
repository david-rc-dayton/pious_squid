import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/maneuver/maneuver_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Hill Modified Equidistant Cyllindrical _(EQCM)_ coordinates.
class Hill extends RelativeState {
  /// Create a new [Hill] object, from relative state.
  Hill(this.epoch, final Vector3D position, final Vector3D velocity,
      final double semimajorAxis)
      : _semimajorAxis = semimajorAxis,
        _meanMotion = Earth.smaToMeanMotion(semimajorAxis),
        super(position, velocity);

  /// Create a new [Hill] object from ECI [state] and its relative
  /// motion [origin].
  factory Hill.fromJ2000(final J2000 state, final J2000 origin) {
    final magrtgt = origin.position.magnitude();
    final magrint = state.position.magnitude();
    final rotEciRsw =
        RelativeState.createMatrix(origin.position, origin.velocity);
    final vtgtrsw = rotEciRsw.multiplyVector3D(origin.velocity);
    final rintrsw = rotEciRsw.multiplyVector3D(state.position);
    final vintrsw = rotEciRsw.multiplyVector3D(state.velocity);

    final sinphiint = rintrsw.z / magrint;
    final phiint = asin(sinphiint);
    final cosphiint = cos(phiint);
    final lambdaint = atan2(rintrsw.y, rintrsw.x);
    final sinlambdaint = sin(lambdaint);
    final coslambdaint = cos(lambdaint);
    final lambdadottgt = vtgtrsw.y / magrtgt;

    final rhill =
        Vector3D(magrint - magrtgt, lambdaint * magrtgt, phiint * magrtgt);

    final rotRswSez = Matrix([
      [sinphiint * coslambdaint, sinphiint * sinlambdaint, -cosphiint],
      [-sinlambdaint, coslambdaint, 0.0],
      [cosphiint * coslambdaint, cosphiint * sinlambdaint, sinphiint]
    ]);

    final vintsez = rotRswSez.multiplyVector3D(vintrsw);
    final phidotint = -vintsez.x / magrint;
    final lambdadotint = vintsez.y / (magrint * cosphiint);

    final vhill = Vector3D(vintsez.z - vtgtrsw.x,
        magrtgt * (lambdadotint - lambdadottgt), magrtgt * phidotint);

    return Hill(origin.epoch, rhill, vhill, origin.semimajorAxis());
  }

  /// Create a new [Hill] object in a linear drift relative to an origin state,
  /// given the [radialPosition] _(km)_, [intrackPosition] _(km)_,
  /// [nodeVelocity] _(km/s)_, and [nodeOffsetTime] _(seconds)_.
  factory Hill.fromLinearDrift(
      final J2000 origin,
      final double radialPosition,
      final double intrackPosition,
      final double nodeVelocity,
      final double nodeOffsetTime) {
    final a = origin.semimajorAxis();
    final n = Earth.smaToMeanMotion(a);
    final yDot = (-3.0 * radialPosition * n) * 0.5;
    final z = nodeVelocity / n * sin(n * -nodeOffsetTime);
    final zDot = nodeVelocity * cos(n * -nodeOffsetTime);
    final r = Vector3D(radialPosition, intrackPosition, z);
    final v = Vector3D(0.0, yDot, zDot);
    return Hill(origin.epoch, r, v, a);
  }

  /// Create a new [Hill] object in a Natural Motion Circumnavigation _(NMC)_
  /// relative to an origin state, given the [majorAxisRange] _(km)_,
  /// [nodeVelocity] _(km/s)_, [nodeOffsetTime] (seconds), and an optional
  /// intrack [translation] _(km)_.
  factory Hill.fromNmc(final J2000 origin, final double majorAxisRange,
      final double nodeVelocity, final double nodeOffsetTime,
      [final double translation = 0.0]) {
    final a = origin.semimajorAxis();
    final n = Earth.smaToMeanMotion(a);
    final xDot = (majorAxisRange * n) * 0.5;
    final z = nodeVelocity / n * sin(n * -nodeOffsetTime);
    final zDot = nodeVelocity * cos(n * -nodeOffsetTime);
    final r = Vector3D(0.0, majorAxisRange + translation, z);
    final v = Vector3D(xDot, 0.0, zDot);
    return Hill(origin.epoch, r, v, a);
  }

  /// Create a new [Hill] object in a V-Bar perch relative to an origin state,
  /// given the [perchRange] _(km)_, [nodeVelocity] _(km/s)_, and
  /// [nodeOffsetTime] _(seconds)_.
  factory Hill.fromPerch(final J2000 origin, final double perchRange,
      final double nodeVelocity, final double nodeOffsetTime) {
    final a = origin.semimajorAxis();
    final n = Earth.smaToMeanMotion(a);
    final z = nodeVelocity / n * sin(n * -nodeOffsetTime);
    final zDot = nodeVelocity * cos(n * -nodeOffsetTime);
    final r = Vector3D(0.0, perchRange, z);
    final v = Vector3D(0.0, 0.0, zDot);
    return Hill(origin.epoch, r, v, a);
  }

  /// UTC epoch.
  final EpochUTC epoch;

  /// Origin semimajor-axis _(km)_.
  double _semimajorAxis;

  /// Origin mean motion _(rad/s)_.
  double _meanMotion;

  @override
  String get name => 'Hill';

  /// Origin semimajor-axis _(km)_.
  double get semimajorAxis => _semimajorAxis;
  set semimajorAxis(final double sma) {
    _semimajorAxis = sma;
    _meanMotion = Earth.smaToMeanMotion(_semimajorAxis);
  }

  /// Origin mean motion _(rad/s)_.
  double get meanMotion => _meanMotion;

  @override
  J2000 toJ2000(final J2000 origin) {
    final magrtgt = origin.position.magnitude();
    final magrint = magrtgt + position.x;
    final rotEciRsw =
        RelativeState.createMatrix(origin.position, origin.velocity);
    final vtgtrsw = rotEciRsw.multiplyVector3D(origin.velocity);

    final lambdadottgt = vtgtrsw.y / magrtgt;
    final lambdaint = position.y / magrtgt;
    final phiint = position.z / magrtgt;
    final sinphiint = sin(phiint);
    final cosphiint = cos(phiint);
    final sinlambdaint = sin(lambdaint);
    final coslambdaint = cos(lambdaint);

    final rotRswSez = Matrix([
      [sinphiint * coslambdaint, sinphiint * sinlambdaint, -cosphiint],
      [-sinlambdaint, coslambdaint, 0],
      [cosphiint * coslambdaint, cosphiint * sinlambdaint, sinphiint]
    ]);

    final rdotint = velocity.x + vtgtrsw.x;
    final lambdadotint = velocity.y / magrtgt + lambdadottgt;
    final phidotint = velocity.z / magrtgt;
    final vintsez = Vector3D(
        -magrint * phidotint, magrint * lambdadotint * cosphiint, rdotint);
    final vintrsw = rotRswSez.transpose().multiplyVector3D(vintsez);
    final vinteci = rotEciRsw.transpose().multiplyVector3D(vintrsw);

    final rintrsw = Vector3D(cosphiint * magrint * coslambdaint,
        cosphiint * magrint * sinlambdaint, sinphiint * magrint);

    final rinteci = rotEciRsw.transpose().multiplyVector3D(rintrsw);

    return J2000(origin.epoch, rinteci, vinteci);
  }

  /// Return the Clohessy-Wiltshire relative motion state transition matrix for
  /// elapsed time [t] _(seconds)_.
  Matrix transitionMatrix(final double t) {
    final n = _meanMotion;
    final sn = sin(n * t);
    final cs = cos(n * t);
    return Matrix([
      [4.0 - 3.0 * cs, 0.0, 0.0, sn / n, 2.0 * (1.0 - cs) / n, 0.0],
      [
        6.0 * (sn - n * t),
        1.0,
        0.0,
        -2.0 * (1.0 - cs) / n,
        (4.0 * sn - 3.0 * n * t) / n,
        0.0
      ],
      [0.0, 0.0, cs, 0.0, 0.0, sn / n],
      [3.0 * n * sn, 0.0, 0.0, cs, 2.0 * sn, 0.0],
      [-6.0 * n * (1.0 - cs), 0.0, 0.0, -2.0 * sn, 4.0 * cs - 3.0, 0.0],
      [0.0, 0.0, -n * sn, 0.0, 0.0, cs]
    ]);
  }

  /// Return this state transitioned by elapsed time [t] _(seconds)_.
  Hill transition(final double t) {
    final sysMat = transitionMatrix(t);
    final output = sysMat.multiplyVector(position.join(velocity));
    return Hill(epoch.roll(t), output.toVector3D(0), output.toVector3D(3),
        _semimajorAxis);
  }

  /// Return this state propagated to the [newEpoch].
  Hill propagate(final EpochUTC newEpoch) =>
      transition(newEpoch.difference(epoch));

  /// Return a copy of this state with the [maneuver] applied.
  Hill maneuver(final Thrust maneuver) {
    final state = propagate(maneuver.center);
    return Hill(state.epoch, state.position,
        state.velocity.add(maneuver.deltaV), state._semimajorAxis);
  }

  /// Return ephemeris over the provided [start] and [stop] epochs.
  ///
  /// Takes an optional ephemeris [step] size _(seconds)_.
  List<Hill> ephemeris(final EpochUTC start, final EpochUTC stop,
      [final double step = 60.0]) {
    final output = <Hill>[];
    var current = start;
    while (stop >= current) {
      output.add(propagate(current));
      current = current.roll(step);
    }
    return output;
  }

  /// Return the relative motion origin's orbital period _(seconds)_.
  double get period => twoPi / _meanMotion;

  /// Return a new [Hill] state at the next radial tangent relative to origin.
  Hill nextRadialTangent() {
    // save key components from current initial state
    final x = position.x;
    final xDot = velocity.x;
    final yDot = velocity.y;
    // solve for the time of nearest xDot == 0 occurence
    var t = atan(-xDot / (3.0 * _meanMotion * x + 2.0 * yDot)) / _meanMotion;
    // correct t by half period for any negative or zero solutions
    if (t <= 0) {
      t += 0.5 * period;
    } else if (t.isNaN) {
      t = 0.5 * period;
    }
    return propagate(epoch.roll(t));
  }

  /// Solve for a maneuver at the current epoch to reach the provided
  /// [waypoint].
  ///
  /// The crosstrack component of the maneuver will be set to zero if
  /// [ignoreCrosstrack] is set to `true`.
  Thrust solveManeuver(final Waypoint waypoint,
      {final bool ignoreCrosstrack = false}) {
    final t = waypoint.epoch.difference(epoch);
    final w = waypoint.relativePosition;
    // get matrix of CW equations for the given time
    final sysMat = transitionMatrix(t);
    // solve constant side of linear system using the first three equations in
    // the CW matrix; this is possible because starting and stopping positions
    // are known over the given interval; we are left with three unknown
    // velocities and three equations
    final posEquationMat = Matrix([
      [sysMat[0][0], sysMat[0][1], sysMat[0][2]],
      [sysMat[1][0], sysMat[1][1], sysMat[1][2]],
      [sysMat[2][0], sysMat[2][1], sysMat[2][2]],
    ]);
    final solnVector = w.subtract(posEquationMat.multiplyVector3D(position));
    // isolate velocity portions of the computed CW matrix
    final velEquationMat = Matrix([
      [sysMat[0][3], sysMat[0][4], sysMat[0][5]],
      [sysMat[1][3], sysMat[1][4], sysMat[1][5]],
      [sysMat[2][3], sysMat[2][4], sysMat[2][5]],
    ]);
    // save difference of desired velocity and current velocity to represent
    // the required burn magnitude
    var result = velEquationMat
        .inverse()
        .multiplyVector3D(solnVector)
        .subtract(velocity);
    if (ignoreCrosstrack) {
      result = Vector3D(result.x * 1000, result.y * 1000, 0.0);
    }
    return Thrust(epoch, result.x * 1000, result.y * 1000, result.z * 1000);
  }

  /// Calculate the maneuver sequence required to arrive at the provided
  /// [waypoints] given a UTC [pivot] epoch for the first maneuver, and an
  /// optional list of pre/post waypoint maneuvers.
  List<Thrust> maneuverSequence(
      final EpochUTC pivot, final List<Waypoint> waypoints,
      {List<Thrust>? preManeuvers, List<Thrust>? postManeuvers}) {
    var state = Hill(epoch, position, velocity, _semimajorAxis);
    preManeuvers ??= [];
    postManeuvers ??= [];
    final output = preManeuvers.sublist(0);
    output.sort((final a, final b) => a.center.compareTo(b.center));
    output
        .retainWhere((final mvr) => mvr.center >= epoch && mvr.center >= pivot);
    for (final mvr in output) {
      state = state.maneuver(mvr);
    }
    state = state.propagate(pivot);
    for (final wpt in waypoints) {
      final mvr = state.solveManeuver(wpt);
      state = state.maneuver(mvr);
      output.add(mvr);
    }
    output.addAll(postManeuvers);
    return output;
  }

  /// Return an update state reflecting the relative motion after the provided
  /// origin spacecraft [maneuver].
  Hill maneuverOrigin(final Thrust maneuver) {
    final state = propagate(maneuver.center);
    final vInit = sqrt(Earth.mu / _semimajorAxis);
    final vFinal = vInit - (maneuver.intrack * 1e-3);
    final aFinal = Earth.mu / (vFinal * vFinal);
    return Hill(state.epoch, state.position,
        state.velocity.subtract(maneuver.deltaV), aFinal);
  }
}
