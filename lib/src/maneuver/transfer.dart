import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Container for a two-burn orbit transfer.
class TwoBurnOrbitTransfer {
  /// Create a new [TwoBurnOrbitTransfer] object.
  TwoBurnOrbitTransfer(
      this.vInit, this.vFinal, this.vTransA, this.vTransB, this.tTrans);

  /// Solve a two-burn Hohmann transfer to change semimajor-axis, given an
  /// initial semimajor-axis [rInit] _(km)_ and final semimajor-axis
  /// [rFinal] _(km)_.
  factory TwoBurnOrbitTransfer.hohmannTransfer(
      final double rInit, final double rFinal) {
    final vInit = sqrt(Earth.mu / rInit);
    final vFinal = sqrt(Earth.mu / rFinal);
    final vTransA = vInit * (sqrt((2.0 * rFinal) / (rInit + rFinal)) - 1.0);
    final vTransB = vFinal * (1.0 - sqrt((2 * rInit) / (rInit + rFinal)));
    final tTrans = pi * sqrt(pow(rInit + rFinal, 3) / (8.0 * Earth.mu));
    return TwoBurnOrbitTransfer(vInit, vFinal, vTransA, vTransB, tTrans);
  }

  /// Initial velocity _(km/s)_.
  final double vInit;

  /// Final velocity _(km/s)_.
  final double vFinal;

  /// First burn delta-velocity _(km/s)_.
  final double vTransA;

  /// Second burn delta-velocity _(km/s)_.
  final double vTransB;

  /// Transfer time-of-flight _(seconds)_.
  final double tTrans;

  /// Return the total delta-velocity magnitude for both maneuvers _(km/s)_.
  double get deltaV => vTransA.abs() + vTransB.abs();

  /// Convert this burn sequence into [Thrust] objects at the provided UTC
  /// [epoch], and optional thruster [durationRate] _(s/m/s)_.
  (Thrust, Thrust) toManeuvers(final EpochUTC epoch,
      {final double durationRate = 0.0}) {
    final mA =
        Thrust(epoch, 0.0, vTransA * 1000.0, 0.0, durationRate: durationRate);
    final mB = Thrust(epoch.roll(tTrans), 0.0, vTransB * 1000.0, 0.0,
        durationRate: durationRate);
    return (mA, mB);
  }
}
