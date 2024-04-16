import 'dart:math';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Kepler analytical two-body propagator.
class KeplerPropagator extends Propagator {
  /// Create a new [KeplerPropagator] object from orbital elements.
  KeplerPropagator(this._initElements)
      : _elements = _initElements,
        _cacheState = J2000.fromClassicalElements(_initElements);

  final ClassicalElements _initElements;
  ClassicalElements _elements;
  J2000 _cacheState;
  final List<J2000> _checkpoints = [];

  @override
  J2000 get state => _cacheState;

  @override
  J2000 propagate(final EpochUTC epoch) {
    _cacheState = J2000.fromClassicalElements(_elements.propagate(epoch));
    return _cacheState;
  }

  @override
  void reset() {
    _elements = _initElements;
    _cacheState = J2000.fromClassicalElements(_elements);
  }

  @override
  List<J2000> maneuver(final Thrust maneuver, [final double interval = 60]) {
    final output = [_cacheState];
    _cacheState = maneuver.apply(propagate(maneuver.center));
    _elements = _cacheState.toClassicalElements();
    output.add(_cacheState);
    return output;
  }

  @override
  VerletBlendInterpolator ephemerisManeuver(
      final EpochUTC start, final EpochUTC finish, final List<Thrust> maneuvers,
      [final double interval = 60.0]) {
    final tMvr = maneuvers.sublist(0);
    tMvr.retainWhere(
        (final mvr) => (mvr.center >= start) || (mvr.center <= finish));
    final ephemeris = <J2000>[];
    if (tMvr[0].start > start) {
      ephemeris.add(propagate(start));
    }
    for (final mvr in tMvr) {
      while (_cacheState.epoch < mvr.center) {
        final step = min(mvr.center.difference(_cacheState.epoch), interval);
        propagate(_cacheState.epoch.roll(step));
        if (_cacheState.epoch.posix != mvr.center.posix) {
          ephemeris.add(_cacheState);
        }
      }
      ephemeris.addAll(maneuver(mvr, interval));
    }
    while (_cacheState.epoch < finish) {
      final step = min(finish.difference(_cacheState.epoch), interval);
      propagate(_cacheState.epoch.roll(step));
      ephemeris.add(_cacheState);
    }
    return VerletBlendInterpolator(ephemeris);
  }

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
