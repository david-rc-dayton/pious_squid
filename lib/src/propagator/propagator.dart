import 'dart:math';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/optimize/optimize_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Propagator base class.
abstract class Propagator {
  /// Propagate state to the provided [epoch].
  J2000 propagate(final EpochUTC epoch);

  /// Generate a [VerletBlendInterpolator] containing ephemeris over the
  /// [start] and [stop] propagation period, with an optional
  /// ephemeris [interval].
  VerletBlendInterpolator ephemeris(final EpochUTC start, final EpochUTC stop,
      [final double interval = 60.0]) {
    final output = <J2000>[propagate(start)];
    var tempEpoch = start;
    while (tempEpoch <= stop) {
      tempEpoch = tempEpoch.roll(interval);
      output.add(propagate(tempEpoch));
    }
    return VerletBlendInterpolator(output);
  }

  /// Reset cached propagator state to its initial value.
  ///
  /// This will clear any post-maneuver states in the propagator.
  void reset();

  /// Store a checkpoint of this propagator's current state, and return the
  /// checkpoint index that can be used to [restore] the propagator state
  /// later.
  int checkpoint();

  /// Restore a state [checkpoint] at the provided index.
  void restore(final int index);

  /// Remove any stored [checkpoint] values from this propagator.
  void clearCheckpoints();

  /// Return the last propagated state.
  J2000 get state;

  /// Generate a list of [J2000] states integrating over a maneuver.
  ///
  /// If the maneuver is impulsive, the list will only contain a single state.
  List<J2000> maneuver(final Thrust maneuver, [final double interval = 60.0]);

  /// Generate a [VerletBlendInterpolator] containing maneuver ephemeris over
  /// the [start] and [finish] interval, with an optional ephemeris [interval].
  VerletBlendInterpolator ephemerisManeuver(
      final EpochUTC start, final EpochUTC finish, final List<Thrust> maneuvers,
      [final double interval = 60.0]);

  /// Propagate to the ascending node after the [start] epoch, and return the
  /// propagated state.
  J2000 propagateAscendingNode(final EpochUTC start) {
    final period = state.period();
    final step = period / 8;
    var current = start;
    final stop = current.roll(period);
    propagate(current);
    var previous = state.position.z;
    while (current <= stop) {
      current = current.roll(step);
      propagate(current);
      if ((state.position.z.sign == -previous.sign) && (state.velocity.z > 0)) {
        break;
      }
      previous = state.position.z;
    }
    final t = GoldenSection.search((final x) {
      propagate(EpochUTC(x));
      return state.position.z.abs();
    }, current.posix - step, current.posix, tolerance: 1e-3);
    propagate(EpochUTC(t));
    return state;
  }

  /// Propagate to the descending node after the [start] epoch, and return the
  /// propagated state.
  J2000 propagateDescendingNode(final EpochUTC start) {
    final period = state.period();
    final step = period / 8;
    var current = start;
    final stop = current.roll(period);
    propagate(current);
    var previous = state.position.z;
    while (current <= stop) {
      current = current.roll(step);
      propagate(current);
      if ((state.position.z.sign == -previous.sign) && (state.velocity.z < 0)) {
        break;
      }
      previous = state.position.z;
    }
    final t = GoldenSection.search((final x) {
      propagate(EpochUTC(x));
      return state.position.z.abs();
    }, current.posix - step, current.posix, tolerance: 1e-3);
    propagate(EpochUTC(t));
    return state;
  }

  /// Propagate to apogee after the [start] epoch, and return the propagated
  /// state.
  J2000 propagateApogee(final EpochUTC start) {
    propagate(start);
    var step = state.period() / 8;
    var current = start;
    var ta = state.toClassicalElements().trueAnomaly;
    while (step > 1e-3) {
      current = current.roll(step);
      propagate(current);
      final taNew = state.toClassicalElements().trueAnomaly;
      if (ta < pi && taNew >= pi) {
        current = current.roll(-step);
        step *= 0.5;
        continue;
      }
      ta = taNew;
    }
    propagate(current);
    return state;
  }

  /// Propagate to perigee after the [start] epoch, and return the propagated
  /// state.
  J2000 propagatePerigee(final EpochUTC start) {
    propagate(start);
    var step = state.period() / 8;
    var current = start;
    var ta = wrapAngle(state.toClassicalElements().trueAnomaly);
    while (step > 1e-3) {
      current = current.roll(step);
      propagate(current);
      final taNew = wrapAngle(state.toClassicalElements().trueAnomaly);
      if (ta < 0 && taNew >= 0) {
        current = current.roll(-step);
        step *= 0.5;
        continue;
      }
      ta = taNew;
    }
    propagate(current);
    return state;
  }
}
