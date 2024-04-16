import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Result of adaptive numerical integration.
class RkResult {
  /// Create a new [RkResult] object.
  RkResult(this.state, this.error, this.newStep);

  /// State vector result.
  final J2000 state;

  /// Local integration error.
  final double error;

  /// Proposed step size.
  final double newStep;
}

/// Runge-Kutta adaptive state checkpoint.
class _RkCheckpoint {
  _RkCheckpoint(this.cacheState, this.stepSize);

  /// Last cached state.
  J2000 cacheState;

  /// Last used step size (seconds).
  double stepSize;
}

/// Adaptive Runge-Kutta propagator base class.
abstract class RungeKuttaAdaptive extends Propagator {
  /// Create a new [RungeKuttaAdaptive] object from an initial state vector
  /// along with an optional [ForceModel] and [tolerance].
  RungeKuttaAdaptive(this._initState,
      [final ForceModel? forceModel, final double tolerance = 1e-9])
      : _cacheState = _initState,
        _forceModel = forceModel ?? (ForceModel()..setGravity()),
        _tolerance = max(_minTolerance, tolerance.abs());

  /// Initial state vector.
  final J2000 _initState;

  /// Propagator perturbation model.
  ForceModel _forceModel;

  /// Cache of last propagated state.
  J2000 _cacheState;

  final List<_RkCheckpoint> _checkpoints = [];

  /// Integrator local error tolerance.
  final double _tolerance;

  /// Integration step size _(seconds)_.
  double _stepSize = 60.0;

  /// Minimum allowable local error tolerance.
  static final double _minTolerance = 1e-15;

  /// Butcher tableau `A` values.
  Float64List get a;

  /// Butcher tableau `B` values.
  List<Float64List> get b;

  /// Butcher tableau `CH` values.
  Float64List get ch;

  /// Butcher tableau `C` values
  Float64List get c;

  /// Integrator order.
  int get order;

  @override
  J2000 get state => _cacheState;

  @override
  void reset() {
    _cacheState = _initState;
    _stepSize = 60.0;
  }

  /// Set numerical integration force model.
  void setForceModel(final ForceModel forceModel) {
    _forceModel = forceModel;
  }

  Vector _kfn(final EpochUTC epoch, final Vector rv, final double hArg,
      final Vector kArg, final double step) {
    final t = epoch.roll(hArg * step);
    final rvNew = rv.add(kArg);
    final sample = J2000(t, rvNew.toVector3D(0), rvNew.toVector3D(3));
    return _forceModel.derivative(sample).scale(step);
  }

  RkResult _integrate(final J2000 state, final double step) {
    final k = List.filled(a.length, Vector.origin3, growable: false);
    final y = state.position.join(state.velocity);
    for (var i = 0; i < a.length; i++) {
      var kArg = Vector.origin6;
      if (i != 0) {
        for (var j = 0; j < i; j++) {
          kArg = kArg.add(k[j].scale(b[i][j]));
        }
      }
      k[i] = _kfn(state.epoch, y, a[i], kArg, step);
    }
    var y1 = y;
    var y2 = y;
    for (var i = 0; i < k.length; i++) {
      y1 = y1.add(k[i].scale(ch[i]));
      y2 = y2.add(k[i].scale(c[i]));
    }
    final teVal = y1.distance(y2);
    var hNew = (0.9 * step * pow(_tolerance / teVal, 1.0 / order)).abs();
    final hOld = step.abs();
    hNew = hNew.clamp(0.2 * hOld, 5.0 * hOld);
    hNew = hNew.clamp(1e-5, 1000.0);
    return RkResult(
        J2000(state.epoch.roll(step), y1.toVector3D(0), y1.toVector3D(3)),
        teVal,
        hNew);
  }

  @override
  J2000 propagate(final EpochUTC epoch) {
    var delta = epoch.difference(_cacheState.epoch);
    while (delta != 0) {
      final direction = delta >= 0 ? 1 : -1;
      final dt = min(delta.abs(), _stepSize) * direction;
      final result = _integrate(_cacheState, dt);
      _stepSize = result.newStep;
      if (result.error > _tolerance) {
        continue;
      }
      _cacheState = result.state;
      delta = epoch.difference(_cacheState.epoch);
    }
    return _cacheState;
  }

  @override
  List<J2000> maneuver(final Thrust maneuver, [final double interval = 60.0]) {
    if (maneuver.isImpulsive) {
      final output = [_cacheState];
      _cacheState = maneuver.apply(propagate(maneuver.center));
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
  int checkpoint() {
    _checkpoints.add(_RkCheckpoint(_cacheState, _stepSize));
    return _checkpoints.length - 1;
  }

  @override
  void clearCheckpoints() {
    _checkpoints.clear();
  }

  @override
  void restore(final int index) {
    final checkpoint = _checkpoints[index];
    _cacheState = checkpoint.cacheState;
    _stepSize = checkpoint.stepSize;
  }
}
