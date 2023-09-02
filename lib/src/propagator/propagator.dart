import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
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

  /// Return the epoch of the ascending node after the [start] epoch.
  EpochUTC ascendingNodeEpoch(final EpochUTC start) {
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
    return EpochUTC(t);
  }

  /// Return the epoch of the descending node after the [start] epoch.
  EpochUTC descendingNodeEpoch(final EpochUTC start) {
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
    return EpochUTC(t);
  }

  /// Return the epoch of apogee after the [start] epoch.
  EpochUTC apogeeEpoch(final EpochUTC start) {
    final slice = 8;
    final period = state.period();
    final step = period / slice;
    var current = start;
    propagate(current);
    var tCache = current;
    var rCache = state.position.magnitude();
    for (var i = 0; i < slice; i++) {
      current = current.roll(step);
      final t = EpochUTC(GoldenSection.search((final x) {
        propagate(EpochUTC(x));
        return state.position.magnitude();
      }, current.posix - step, current.posix, tolerance: 1e-3, solveMax: true));
      propagate(t);
      final r = state.position.magnitude();
      if (r > rCache) {
        tCache = t;
        rCache = r;
      }
    }
    return tCache;
  }

  /// Return the epoch of perigee after the [start] epoch.
  EpochUTC perigeeEpoch(final EpochUTC start) {
    final slice = 8;
    final period = state.period();
    final step = period / slice;
    var current = start;
    propagate(current);
    var tCache = current;
    var rCache = state.position.magnitude();
    for (var i = 0; i < slice; i++) {
      current = current.roll(step);
      final t = EpochUTC(GoldenSection.search((final x) {
        propagate(EpochUTC(x));
        return state.position.magnitude();
      }, current.posix - step, current.posix,
          tolerance: 1e-3, solveMax: false));
      propagate(t);
      final r = state.position.magnitude();
      if (r < rCache) {
        tCache = t;
        rCache = r;
      }
    }
    return tCache;
  }
}
