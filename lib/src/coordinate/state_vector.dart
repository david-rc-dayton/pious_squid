import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Position and velocity [Vector3D] container.
typedef PositionVelocity = ({Vector3D position, Vector3D velocity});

/// Base class for state vectors.
abstract class StateVector {
  /// Create a new [StateVector] object, from [position] and [velocity] _(km)_.
  StateVector(this.epoch, this.position, this.velocity);

  /// UTC epoch.
  final EpochUTC epoch;

  /// Position vector _(km)_.
  final Vector3D position;

  /// Velocity vector _(km)_.
  final Vector3D velocity;

  /// Return the name of this coordinate frame.
  String get name;

  /// Return `true` if this coordinate frame is inertial.
  bool get inertial;

  @override
  String toString() => [
        '[$name]',
        '  Epoch: $epoch',
        '  Position: ${position.toString(6)} km',
        '  Velocity: ${velocity.toString(9)} km/s'
      ].join('\n');

  /// Return the state position and velocity as a single vector _(km,km/s)_.
  Vector posvel() => position.join(velocity);

  /// Return the mechanical energy _(km²/s²)_ of this orbit.
  double mechanicalEnergy() {
    final r = position.magnitude();
    final v = velocity.magnitude();
    return ((v * v) * 0.5) - Earth.mu / r;
  }

  /// Return the semimajor-axis _(km)_ of this orbit.
  double semimajorAxis() {
    final energy = mechanicalEnergy();
    return -Earth.mu / (2.0 * energy);
  }

  /// Return the period _(seconds)_ of this orbit.
  double period() {
    final a = semimajorAxis();
    return twoPi * sqrt((a * a * a) / Earth.mu);
  }

  /// Return the angular rate _(rad/s)_ of this orbit.
  double angularRate() {
    final a = semimajorAxis();
    return sqrt(Earth.mu / (a * a * a));
  }

  /// Return the angular rate _(rad/s)_ of this orbit as a vector.
  Vector3D angularRateVector() => Vector3D(0.0, 0.0, angularRate());

  /// Convert this to a [ClassicalElements] object.
  ///
  /// An optional gravitational parameter [mu] _(km²/s³)_ can be provided.
  ///
  /// An error will be thrown if this coordinate frame is not inertial.
  ClassicalElements toClassicalElements({final double mu = Earth.mu}) {
    if (!inertial) {
      throw 'Classical elements are undefined for fixed frames.';
    }
    return ClassicalElements.fromStateVector(this, mu: mu);
  }
}
