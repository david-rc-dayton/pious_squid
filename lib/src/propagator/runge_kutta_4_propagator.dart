import 'dart:math';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/interpolator/verlet_blend_interpolator.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/epoch_utc.dart';

/// Runge-Kutta 4 fixed numerical propagator.
class RungeKutta4Propagator extends Propagator {
  /// Create a new [RungeKutta4Propagator] object from an initial state vector and
  /// along with an optional [ForceModel] and [stepSize] in seconds.
  RungeKutta4Propagator(this._initState,
      [final ForceModel? forceModel, final double stepSize = 15.0])
      : _cacheState = _initState,
        _forceModel = forceModel ?? (ForceModel()..setGravity()),
        _stepSize = stepSize.abs();

  /// Initial state vector.
  final J2000 _initState;

  /// Propagator perturbation model.
  ForceModel _forceModel;

  /// Cache of last propagated state.
  J2000 _cacheState;

  final List<J2000> _checkpoints = [];

  /// Integration step size _(seconds)_.
  double _stepSize;

  /// Set the integrator step size to the provided number of [seconds].
  void setStepSize(final double seconds) => _stepSize = seconds.abs();

  /// Set numerical integration force model.
  void setForceModel(final ForceModel forceModel) {
    _forceModel = forceModel;
  }

  @override
  VerletBlendInterpolator ephemerisManeuver(
      final EpochUTC start, final EpochUTC finish, final List<Thrust> maneuvers,
      [final double interval = 60.0]) {
    final tMvr = maneuvers.sublist(0);
    tMvr.retainWhere(
        (final mvr) => (mvr.start >= start) || (mvr.stop <= finish));
    final ephemeris = <J2000>[];
    if (tMvr[0].start > start) {
      ephemeris.add(propagate(start));
    }
    for (final mvr in tMvr) {
      while (_cacheState.epoch < mvr.start) {
        final step = min(mvr.start.difference(_cacheState.epoch), interval);
        propagate(_cacheState.epoch.roll(step));
        if (_cacheState.epoch.posix != mvr.start.posix) {
          ephemeris.add(_cacheState);
        }
      }
      ephemeris.addAll(maneuver(mvr, interval));
    }
    while (_cacheState.epoch.posix < finish.posix) {
      final step = min(finish.difference(_cacheState.epoch), interval);
      propagate(_cacheState.epoch.roll(step));
      ephemeris.add(_cacheState);
    }
    return VerletBlendInterpolator(ephemeris);
  }

  @override
  List<J2000> maneuver(final Thrust maneuver, [final double interval = 60.0]) {
    if (maneuver.isImpulsive) {
      final output = [_cacheState];
      _cacheState = maneuver.apply(propagate(maneuver.center));
      propagate(_cacheState.epoch.roll(1e-3));
      output.add(_cacheState);
      return output;
    }
    var tState = propagate(maneuver.start);
    _forceModel.loadManeuver(maneuver);
    final ephemeris = [tState];
    while (tState.epoch < maneuver.stop) {
      final step = min(maneuver.stop.difference(tState.epoch), interval);
      tState = propagate(tState.epoch.roll(step));
      ephemeris.add(tState);
    }
    _forceModel.clearManeuver();
    return ephemeris;
  }

  Vector _kFn(final J2000 state, final double hArg, final Vector kArg) {
    final epoch = state.epoch.roll(hArg);
    final posvel = state.position.join(state.velocity);
    final result = posvel.add(kArg);
    final sample = J2000(epoch, result.toVector3D(0), result.toVector3D(3));
    return _forceModel.derivative(sample);
  }

  J2000 _integrate(final J2000 state, final double step) {
    final k1 = _kFn(state, 0, Vector.zero(6)).scale(step);
    final k2 = _kFn(state, 0.5 * step, k1.scale(0.5)).scale(step);
    final k3 = _kFn(state, 0.5 * step, k2.scale(0.5)).scale(step);
    final k4 = _kFn(state, step, k3).scale(step);
    final v1 = k1;
    final v2 = v1.add(k2.scale(2));
    final v3 = v2.add(k3.scale(2));
    final v4 = v3.add(k4);
    final tNext = state.epoch.roll(step);
    final posvel = state.position.join(state.velocity);
    final result = posvel.add(v4.scale(1 / 6));
    return J2000(tNext, result.toVector3D(0), result.toVector3D(3));
  }

  @override
  J2000 propagate(final EpochUTC epoch) {
    var delta = epoch.difference(_cacheState.epoch);
    while (delta != 0) {
      final direction = delta >= 0 ? 1 : -1;
      final dt = min(delta.abs(), _stepSize) * direction;
      _cacheState = _integrate(_cacheState, dt);
      delta = epoch.difference(_cacheState.epoch);
    }
    return _cacheState;
  }

  @override
  void reset() {
    _cacheState = _initState;
  }

  @override
  J2000 get state => _cacheState;

  @override
  int checkpoint() {
    _checkpoints.add(_cacheState);
    return _checkpoints.length - 1;
  }

  @override
  void clearCheckpoints() {
    _checkpoints.clear();
  }

  @override
  void restore(final int index) {
    _cacheState = _checkpoints[index];
  }
}
