import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/interpolator/interpolator_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// SGP4 propagator.
class Sgp4Propagator extends Propagator {
  /// Create a new [Sgp4Propagator] object from a [TLE].
  Sgp4Propagator(this.tle) : _cacheState = tle.state.toJ2000();

  /// [TLE] object used for propagation.
  final TLE tle;

  J2000 _cacheState;
  final List<J2000> _checkpoints = [];

  @override
  J2000 get state => _cacheState;

  @override
  VerletBlendInterpolator ephemerisManeuver(
      final EpochUTC start, final EpochUTC finish, final List<Thrust> maneuvers,
      [final double interval = 60.0]) {
    throw 'Maneuvers cannot be modelled with SGP4.';
  }

  @override
  List<J2000> maneuver(final Thrust maneuver, [final double interval = 60.0]) {
    throw 'Maneuvers cannot be modelled with SGP4.';
  }

  @override
  J2000 propagate(final EpochUTC epoch) {
    _cacheState = tle.propagate(epoch).toJ2000();
    return _cacheState;
  }

  @override
  void reset() {
    _cacheState = tle.state.toJ2000();
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
