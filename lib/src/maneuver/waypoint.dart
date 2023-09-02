import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/optimize/optimize_base.dart';
import 'package:pious_squid/src/orbit_determination/orbit_determination_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Relative waypoint targeting.
class Waypoint {
  /// Create a new [Waypoint] object.
  Waypoint(this.epoch, this.relativePosition);

  /// Waypoint arrival time.
  final EpochUTC epoch;

  /// Position relative to target _(km)_.
  final Vector3D relativePosition;

  /// Return the perturbed error in a [maneuver] when compared against the
  /// target [waypoint] given an initial [state], [forceModel],
  /// [target] interpolator, and speculative relative maneuver
  /// [components] _(m/s)_.
  static double _error(
      final Waypoint waypoint,
      final Thrust maneuver,
      final J2000 state,
      final ForceModel forceModel,
      final StateInterpolator target,
      final Float64List components) {
    final testManeuver = Thrust(maneuver.center, components[0], components[1],
        components[2], maneuver.durationRate);
    final propagator = RungeKutta89Propagator(state, forceModel);
    final maneuverSteps = propagator.maneuver(testManeuver);
    final postManeuver = maneuverSteps.last;
    final interceptor = RungeKutta89Propagator(postManeuver, forceModel)
        .propagate(waypoint.epoch);
    final targetState = target.interpolate(waypoint.epoch);
    if (targetState == null) {
      throw 'Error calculation failed; epoch is outside the target interpolator ephemeris window.';
    }
    final expected = waypoint.relativePosition;
    final actual = RIC.fromJ2000(interceptor, targetState);
    return actual.position.distance(expected);
  }

  /// Generate a score function for refining perturbed [waypoint] maneuvers.
  ///
  /// The score function takes an array of speculative radial, intrack, and
  /// crosstrack components _(m/s)_ and returns the propagated error from the
  /// desired waypoint target.
  static CostFunction _refineManeuverScore(
          final Waypoint waypoint,
          final Thrust maneuver,
          final J2000 state,
          final ForceModel forceModel,
          final StateInterpolator target) =>
      (final Float64List components) =>
          _error(waypoint, maneuver, state, forceModel, target, components);

  /// Convert an array of [waypoints] into a maneuver sequence given an
  /// [interceptor] state, [pivot] epoch for the first burn to arrive at the
  /// first waypoint, and [target] ephemeris interpolator.
  ///
  /// Optional arguments are as follows:
  ///   - `preManeuvers`: maneuvers to execute before the pivot burn
  ///   - `postManeuvers`: maneuvers to execute after the last pivot burn
  ///   - `durationRate`: thruster duration rate _(s/m/s)_
  ///   - `forceModel`: interceptor force model, defaults to two-body
  ///   - `refine`: refine maneuvers to account for perturbations if `true`
  ///   - `maxIter`: maximum refinement iterations per maneuver
  ///   - `printIter`: print debug information on each refinement iteration
  static List<Thrust> toManeuvers(
      final J2000 interceptor,
      final EpochUTC pivot,
      final List<Waypoint> waypoints,
      final StateInterpolator target,
      final List<Thrust>? preManeuvers,
      final List<Thrust>? postManeuvers,
      {final double durationRate = 0.0,
      final ForceModel? forceModel,
      final bool refine = false,
      final int maxIter = 500,
      final bool printIter = false}) {
    final preMnv = preManeuvers ?? <Thrust>[];
    final postMnv = preManeuvers ?? <Thrust>[];
    var state = interceptor;
    if (preMnv.isNotEmpty) {
      for (final maneuver in preMnv) {
        final mvrStep =
            RungeKutta89Propagator(state, forceModel).maneuver(maneuver);
        state = mvrStep.last;
      }
    }
    state = RungeKutta89Propagator(state, forceModel).propagate(pivot);
    final pivotState = state;
    var waypointManeuvers = <Thrust>[];
    for (final wp in waypoints) {
      final targetWp = RIC(wp.relativePosition, Vector3D.origin);
      final targetState = target.interpolate(wp.epoch);
      if (targetState == null) {
        throw 'Waypoint outside target interpolator window.';
      }
      final wpState = targetWp.toJ2000(targetState);
      final tof = wp.epoch.difference(state.epoch);
      final revs = (tof / state.period()).floor();
      final shortPath = LambertIOD.useShortPath(state, targetState);
      final lambert = LambertIOD();
      final components = lambert.estimate(
          state.position, wpState.position, state.epoch, wp.epoch,
          posigrade: shortPath, nRev: revs);
      if (components == null) {
        throw 'Lambert solve result is null.';
      }
      final componentsRel =
          RIC.fromJ2000(components, state).velocity.scale(1e3);
      final maneuver = Thrust(state.epoch, componentsRel.x, componentsRel.y,
          componentsRel.z, durationRate);
      final tempProp = RungeKutta89Propagator(state, forceModel);
      tempProp.maneuver(maneuver);
      state = tempProp.propagate(wp.epoch);
      waypointManeuvers.add(maneuver);
    }
    if (refine) {
      waypointManeuvers = _refineManeuvers(waypoints, waypointManeuvers,
          pivotState, forceModel ?? (ForceModel()..setGravity()), target,
          maxIter: maxIter, printIter: printIter);
    }
    final output = <Thrust>[];
    output.addAll(preMnv);
    output.addAll(waypointManeuvers);
    output.addAll(postMnv);
    return output;
  }

  static List<Thrust> _refineManeuvers(
      final List<Waypoint> waypoints,
      final List<Thrust> maneuvers,
      final J2000 interceptor,
      final ForceModel forceModel,
      final StateInterpolator target,
      {final int maxIter = 500,
      final bool printIter = false}) {
    var state = interceptor;
    final output = <Thrust>[];
    for (var i = 0; i < maneuvers.length; i++) {
      final maneuver = maneuvers[i];
      final waypoint = waypoints[i];
      final guess = maneuver.deltaV.scale(1e3).toArray();
      final simplex = DownhillSimplex.generateSimplex(guess, 1e-1);
      final scoreFn =
          _refineManeuverScore(waypoint, maneuver, state, forceModel, target);
      final results = DownhillSimplex.solveSimplex(scoreFn, simplex,
          maxIter: maxIter,
          xTolerance: 1e-6,
          fTolerance: 1e-6,
          printIter: printIter);
      final tR = results[0];
      final tI = results[1];
      final tC = results[2];
      final newManeuver =
          Thrust(maneuver.center, tR, tI, tC, maneuver.durationRate);
      output.add(newManeuver);
      final mvrStep =
          RungeKutta89Propagator(state, forceModel).maneuver(newManeuver);
      state = mvrStep.last;
    }
    return output;
  }
}
