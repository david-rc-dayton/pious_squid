import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/pious_squid_base.dart';

/// Relative state type.
enum RelativeStateType {
  /// Radial-Intrack-Crosstrack
  ric('RIC'),

  /// Modified Equidistance Cylindrical
  eqcm('EQCM');

  const RelativeStateType(this.name);

  /// Relative state type name.
  final String name;
}

Matrix _createMatrix3(final Vector3D position, final Vector3D velocity) {
  final ru = position.normalize();
  final cu = position.cross(velocity).normalize();
  final iu = cu.cross(ru).normalize();
  final output = Float64List(9);
  output[0] = ru.x;
  output[1] = ru.y;
  output[2] = ru.z;
  output[3] = iu.x;
  output[4] = iu.y;
  output[5] = iu.z;
  output[6] = cu.x;
  output[7] = cu.y;
  output[8] = cu.z;
  return Matrix(3, 3, output);
}

Matrix _createMatrix6(final Vector3D position, final Vector3D velocity) {
  final m = _createMatrix3(position, velocity);
  final output = Matrix(6, 6);
  output.setBlock(m, 0, 0);
  output.setBlock(m, 3, 3);
  return output;
}

RelativeState _ricFromJ2kMatrix(
    final J2000 state, final J2000 origin, final Matrix transform) {
  final w = origin.angularRateVector();
  final rRic =
      transform.multiplyVector3D(state.position.subtract(origin.position));
  final vRic = transform
      .multiplyVector3D(state.velocity.subtract(origin.velocity))
      .subtract(w.cross(rRic));
  return RelativeState(origin.epoch, rRic, vRic, origin.semimajorAxis(),
      type: RelativeStateType.ric);
}

J2000 _ricToJ2kMatrix(
    final RelativeState relative, final J2000 origin, final Matrix transform) {
  final w = origin.angularRateVector();
  final tt = transform.transpose();
  final tr = tt.multiplyVector3D(relative.position);
  final tv = tt.multiplyVector3D(relative.velocity.add(w.cross(tr)));
  return J2000(origin.epoch, origin.position.add(tr), origin.velocity.add(tv));
}

RelativeState _eqcmFromJ2kMatrix(
    final J2000 state, final J2000 origin, final Matrix transform) {
  final magrtgt = origin.position.magnitude();
  final magrint = state.position.magnitude();
  final vtgtrsw = transform.multiplyVector3D(origin.velocity);
  final rintrsw = transform.multiplyVector3D(state.position);
  final vintrsw = transform.multiplyVector3D(state.velocity);

  final sinphiint = rintrsw.z / magrint;
  final phiint = asin(sinphiint);
  final cosphiint = cos(phiint);
  final lambdaint = atan2(rintrsw.y, rintrsw.x);
  final sinlambdaint = sin(lambdaint);
  final coslambdaint = cos(lambdaint);
  final lambdadottgt = vtgtrsw.y / magrtgt;

  final rhill =
      Vector3D(magrint - magrtgt, lambdaint * magrtgt, phiint * magrtgt);

  final rotRswSezArr = Float64List(9);
  rotRswSezArr[0] = sinphiint * coslambdaint;
  rotRswSezArr[1] = sinphiint * sinlambdaint;
  rotRswSezArr[2] = -cosphiint;
  rotRswSezArr[3] = -sinlambdaint;
  rotRswSezArr[4] = coslambdaint;
  rotRswSezArr[6] = cosphiint * coslambdaint;
  rotRswSezArr[7] = cosphiint * sinlambdaint;
  rotRswSezArr[8] = sinphiint;
  final rotRswSez = Matrix(3, 3, rotRswSezArr);

  final vintsez = rotRswSez.multiplyVector3D(vintrsw);
  final phidotint = -vintsez.x / magrint;
  final lambdadotint = vintsez.y / (magrint * cosphiint);

  final vhill = Vector3D(vintsez.z - vtgtrsw.x,
      magrtgt * (lambdadotint - lambdadottgt), magrtgt * phidotint);

  return RelativeState(origin.epoch, rhill, vhill, origin.semimajorAxis(),
      type: RelativeStateType.eqcm);
}

J2000 _eqcmToJ2kMatrix(
    final RelativeState relative, final J2000 origin, final Matrix transform) {
  final tt = transform.transpose();
  final magrtgt = origin.position.magnitude();
  final magrint = magrtgt + relative.position.x;
  final vtgtrsw = transform.multiplyVector3D(origin.velocity);

  final lambdadottgt = vtgtrsw.y / magrtgt;
  final lambdaint = relative.position.y / magrtgt;
  final phiint = relative.position.z / magrtgt;
  final sinphiint = sin(phiint);
  final cosphiint = cos(phiint);
  final sinlambdaint = sin(lambdaint);
  final coslambdaint = cos(lambdaint);

  final rotRswSezArr = Float64List(9);
  rotRswSezArr[0] = sinphiint * coslambdaint;
  rotRswSezArr[1] = sinphiint * sinlambdaint;
  rotRswSezArr[2] = -cosphiint;
  rotRswSezArr[3] = -sinlambdaint;
  rotRswSezArr[4] = coslambdaint;
  rotRswSezArr[6] = cosphiint * coslambdaint;
  rotRswSezArr[7] = cosphiint * sinlambdaint;
  rotRswSezArr[8] = sinphiint;
  final rotRswSez = Matrix(3, 3, rotRswSezArr);

  final rdotint = relative.velocity.x + vtgtrsw.x;
  final lambdadotint = relative.velocity.y / magrtgt + lambdadottgt;
  final phidotint = relative.velocity.z / magrtgt;
  final vintsez = Vector3D(
      -magrint * phidotint, magrint * lambdadotint * cosphiint, rdotint);
  final vintrsw = rotRswSez.transpose().multiplyVector3D(vintsez);
  final vinteci = tt.multiplyVector3D(vintrsw);

  final rintrsw = Vector3D(cosphiint * magrint * coslambdaint,
      cosphiint * magrint * sinlambdaint, sinphiint * magrint);

  final rinteci = tt.multiplyVector3D(rintrsw);

  return J2000(origin.epoch, rinteci, vinteci);
}

/// Base class for relative states.
class RelativeState {
  /// Create a new [RelativeState] object.
  RelativeState(this.epoch, this.position, this.velocity, this._semimajorAxis,
      {this.type = RelativeStateType.ric})
      : _meanMotion = Earth.smaToMeanMotion(_semimajorAxis);

  /// Create a [RelativeState] from two J2000 states.
  factory RelativeState.fromJ2000(final J2000 state, final J2000 origin,
      {final RelativeStateType type = RelativeStateType.ric}) {
    final transform = _createMatrix3(origin.position, origin.velocity);
    switch (type) {
      case RelativeStateType.ric:
        return _ricFromJ2kMatrix(state, origin, transform);
      case RelativeStateType.eqcm:
        return _eqcmFromJ2kMatrix(state, origin, transform);
    }
  }

  /// Create a [RelativeState] for a linear drift relative to an origin state,
  /// given the [radialPosition] _(km)_, [intrackPosition] _(km)_,
  /// [nodeVelocity] _(km/s)_, and [nodeOffsetTime] _(seconds)_.
  factory RelativeState.fromLinearDrift(
      final J2000 origin,
      final double radialPosition,
      final double intrackPosition,
      final double nodeVelocity,
      final double nodeOffsetTime,
      {final RelativeStateType type = RelativeStateType.ric}) {
    final a = origin.semimajorAxis();
    final n = Earth.smaToMeanMotion(a);
    final yDot = (-3.0 * radialPosition * n) * 0.5;
    final z = nodeVelocity / n * sin(n * -nodeOffsetTime);
    final zDot = nodeVelocity * cos(n * -nodeOffsetTime);
    final r = Vector3D(radialPosition, intrackPosition, z);
    final v = Vector3D(0.0, yDot, zDot);
    return RelativeState(origin.epoch, r, v, a, type: type);
  }

  /// Create a [RelativeState] for a Natural Motion Circumnavigation _(NMC)_
  /// relative to an origin state, given the [majorAxisRange] _(km)_,
  /// [nodeVelocity] _(km/s)_, [nodeOffsetTime] (seconds), and an optional
  /// intrack [translation] _(km)_.
  factory RelativeState.fromNmc(final J2000 origin, final double majorAxisRange,
      final double nodeVelocity, final double nodeOffsetTime,
      {final double translation = 0.0,
      final RelativeStateType type = RelativeStateType.ric}) {
    final a = origin.semimajorAxis();
    final n = Earth.smaToMeanMotion(a);
    final xDot = (majorAxisRange * n) * 0.5;
    final z = nodeVelocity / n * sin(n * -nodeOffsetTime);
    final zDot = nodeVelocity * cos(n * -nodeOffsetTime);
    final r = Vector3D(0.0, majorAxisRange + translation, z);
    final v = Vector3D(xDot, 0.0, zDot);
    return RelativeState(origin.epoch, r, v, a, type: type);
  }

  /// Create [RelativeState] for a V-Bar perch relative to an origin state,
  /// given the [perchRange] _(km)_, [nodeVelocity] _(km/s)_, and
  /// [nodeOffsetTime] _(seconds)_.
  factory RelativeState.fromPerch(final J2000 origin, final double perchRange,
      final double nodeVelocity, final double nodeOffsetTime,
      {final RelativeStateType type = RelativeStateType.ric}) {
    final a = origin.semimajorAxis();
    final n = Earth.smaToMeanMotion(a);
    final z = nodeVelocity / n * sin(n * -nodeOffsetTime);
    final zDot = nodeVelocity * cos(n * -nodeOffsetTime);
    final r = Vector3D(0.0, perchRange, z);
    final v = Vector3D(0.0, 0.0, zDot);
    return RelativeState(origin.epoch, r, v, a, type: type);
  }

  /// Perform a pseudo-Lambert solve between relative positions [x0] and [xf]
  /// and time of flight [t0] and [t1].
  factory RelativeState.solveTransfer(final Vector3D x0, final Vector3D xf,
      final EpochUTC t0, final EpochUTC t1, final double sma,
      {final RelativeStateType type = RelativeStateType.ric}) {
    final dt = t1.difference(t0);
    final n = Earth.smaToMeanMotion(sma);
    final stm = RelativeState.transitionMatrix(dt, n);
    final rr = stm.getBlock(0, 0, 3, 3);
    final rv = stm.getBlock(0, 3, 3, 3);
    final result =
        rv.solve(xf.toVector().subtract(rr.multiplyVector(x0.toVector())));
    return RelativeState(t0, x0, result.toVector3D(0), sma, type: type);
  }

  /// Relative state type.
  final RelativeStateType type;

  /// Relative state epoch _(utc)_.
  final EpochUTC epoch;

  /// Relative position vector _(km)_.
  final Vector3D position;

  /// Relative velocity vector _(km)_.
  final Vector3D velocity;

  /// Origin orbit semimajor-axis _(km)_.
  double _semimajorAxis;

  /// Origin orbit mean-motion _(rad/s)_.
  double _meanMotion;

  /// Return the name of this coordinate frame.
  String get name => type.name;

  @override
  String toString() => [
        '[$name]',
        '  Epoch:          $epoch',
        '  Position:       ${position.toString(6)} km',
        '  Velocity:       ${velocity.toString(6)} km/s',
        '  Semimajor-Axis: ${_semimajorAxis.toStringAsFixed(3)} km',
        '  Mean Motion:    ${_meanMotion.toStringAsFixed(7)} rad/s'
      ].join('\n');

  /// Origin semimajor-axis _(km)_.
  double get semimajorAxis => _semimajorAxis;
  set semimajorAxis(final double sma) {
    _semimajorAxis = sma;
    _meanMotion = Earth.smaToMeanMotion(_semimajorAxis);
  }

  /// Origin mean motion _(rad/s)_.
  double get meanMotion => _meanMotion;

  /// Convert this to a [J2000] state vector object.
  J2000 toJ2000(final J2000 origin) {
    final transform = _createMatrix3(origin.position, origin.velocity);
    switch (type) {
      case RelativeStateType.ric:
        return _ricToJ2kMatrix(this, origin, transform);
      case RelativeStateType.eqcm:
        return _eqcmToJ2kMatrix(this, origin, transform);
    }
  }

  /// Create a 3x3 relative state transform matrix.
  static Matrix createMatrix3(
          final Vector3D position, final Vector3D velocity) =>
      _createMatrix3(position, velocity);

  /// Create a 6x6 relative state transform matrix.
  static Matrix createMatrix6(
          final Vector3D position, final Vector3D velocity) =>
      _createMatrix6(position, velocity);

  /// Return the relative range _(km)_.
  double range() => position.magnitude();

  /// Return the relative range rate _(km/s)_.
  double rangeRate() {
    final r = range();
    if (r == 0) {
      return 0.0;
    }
    return position.dot(velocity) / r;
  }

  /// Return the relative motion origin's orbital period _(seconds)_.
  double get period => twoPi / _meanMotion;

  /// Return the Clohessy-Wiltshire relative motion state transition matrix for
  /// elapsed time [t] _(seconds)_ and [meanMotion] _(rad/s)_.
  static Matrix transitionMatrix(final double t, final double meanMotion) {
    final n = meanMotion;
    final s = sin(n * t);
    final c = cos(n * t);
    final output = Float64List(36);
    output[0] = 4 - 3 * c;
    output[3] = s / n;
    output[4] = 2 * (1 - c) / n;
    output[6] = 6 * (s - n * t);
    output[7] = 1;
    output[9] = 2 * (c - 1) / n;
    output[10] = (4 * s - 3 * n * t) / n;
    output[14] = c;
    output[17] = s / n;
    output[18] = 3 * n * s;
    output[21] = c;
    output[22] = 2 * s;
    output[24] = 6 * n * (c - 1);
    output[27] = -2 * s;
    output[28] = 4 * c - 3;
    output[32] = -n * s;
    output[35] = c;
    return Matrix(6, 6, output);
  }

  /// Return the Clohessy-Wiltshire relative motion state transition matrix for
  /// a [newEpoch] _(UTC)_.
  Matrix transitionMatrixEpoch(final EpochUTC newEpoch) {
    final t = newEpoch.difference(epoch);
    return transitionMatrix(t, meanMotion);
  }

  /// Return this state transitioned by elapsed time [t] _(seconds)_.
  RelativeState transition(final double t) {
    final sysMat = RelativeState.transitionMatrix(t, meanMotion);
    final output = sysMat.multiplyVector(position.join(velocity));
    return RelativeState(epoch.roll(t), output.toVector3D(0),
        output.toVector3D(3), semimajorAxis,
        type: type);
  }

  /// Return this state propagated to the [newEpoch].
  RelativeState propagate(final EpochUTC newEpoch) =>
      transition(newEpoch.difference(epoch));

  /// Return a copy of this state with the [maneuver] applied.
  RelativeState maneuver(final Thrust maneuver) {
    final state = propagate(maneuver.center);
    return RelativeState(state.epoch, state.position,
        state.velocity.add(maneuver.deltaV), state.semimajorAxis,
        type: type);
  }

  /// Return ephemeris over the provided [start] and [stop] epochs.
  ///
  /// Takes an optional ephemeris [step] size _(seconds)_.
  List<RelativeState> ephemeris(final EpochUTC start, final EpochUTC stop,
      [final double step = 60.0]) {
    final output = <RelativeState>[];
    var current = start;
    while (stop >= current) {
      output.add(propagate(current));
      current = current.roll(step);
    }
    return output;
  }

  /// Return a new [Hill] state at the next radial tangent relative to origin.
  RelativeState nextRadialTangent() {
    // save key components from current initial state
    final x = position.x;
    final xDot = velocity.x;
    final yDot = velocity.y;
    // solve for the time of nearest xDot == 0 occurence
    var t = atan(-xDot / (3.0 * meanMotion * x + 2.0 * yDot)) / meanMotion;
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
    final dt = waypoint.epoch.difference(epoch);
    final n = meanMotion;
    final rfRel = waypoint.relativePosition.toVector();

    final r0 = position.toVector();
    final v0 = velocity.toVector();

    final stm = RelativeState.transitionMatrix(dt, n);
    final rr = stm.getBlock(0, 0, 3, 3);
    final rv = stm.getBlock(0, 3, 3, 3);
    var result = rv
        .solve(rfRel.subtract(rr.multiplyVector(r0)))
        .subtract(v0)
        .toVector3D(0);
    if (ignoreCrosstrack) {
      result = Vector3D(result.x, result.y, 0.0);
    }

    return Thrust(epoch, result.x * 1000, result.y * 1000, result.z * 1000);
  }

  /// Calculate the maneuver sequence required to arrive at the provided
  /// [waypoints] given a UTC [pivot] epoch for the first maneuver, and an
  /// optional list of pre/post waypoint maneuvers.
  List<Thrust> maneuverSequence(
      final EpochUTC pivot, final List<Waypoint> waypoints,
      {List<Thrust>? preManeuvers, List<Thrust>? postManeuvers}) {
    var state =
        RelativeState(epoch, position, velocity, semimajorAxis, type: type);
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
  RelativeState maneuverOrigin(final Thrust maneuver) {
    final state = propagate(maneuver.center);
    final vInit = sqrt(Earth.mu / _semimajorAxis);
    final vFinal = vInit - (maneuver.intrack * 1e-3);
    final aFinal = Earth.mu / (vFinal * vFinal);
    return RelativeState(state.epoch, state.position,
        state.velocity.subtract(maneuver.deltaV), aFinal,
        type: type);
  }
}
