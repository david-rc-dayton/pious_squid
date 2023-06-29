import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
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

  /// Reset cached propagator state to it's initial value.
  ///
  /// This will clear any post-maneuver states in the propagator.
  void reset();

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

  /// Return the epoch of apogee after the [start] epoch.
  EpochUTC apogeeEpoch(final EpochUTC start) {
    var current = start;
    var hub = start;
    var step = 300.0;
    while (true) {
      final r0 = propagate(current).position.magnitude();
      current = current.roll(step);
      final r1 = propagate(current).position.magnitude();
      current = current.roll(step);
      final r2 = propagate(current).position.magnitude();
      if (r1 > r0 && r1 > r2) {
        if (step < 1e-1) {
          break;
        }
        current = hub;
        step /= 2;
        continue;
      }
      hub = hub.roll(step);
      current = hub;
    }
    return hub.roll(step);
  }

  /// Return the epoch of perigee after the [start] epoch.
  EpochUTC perigeeEpoch(final EpochUTC start) {
    var current = start;
    var hub = start;
    var step = 300.0;
    while (true) {
      final r0 = propagate(current).position.magnitude();
      current = current.roll(step);
      final r1 = propagate(current).position.magnitude();
      current = current.roll(step);
      final r2 = propagate(current).position.magnitude();
      if (r1 < r0 && r1 < r2) {
        if (step < 1e-1) {
          break;
        }
        current = hub;
        step /= 2;
        continue;
      }
      hub = hub.roll(step);
      current = hub;
    }
    return hub.roll(step);
  }
}
