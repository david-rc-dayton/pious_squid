import 'dart:typed_data';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/optimize/optimize_base.dart';
import 'package:pious_squid/src/orbit_determination/orbit_determination_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Modified Gooding angles-only initial orbit determination.
class ModifiedGoodingIOD {
  /// Create a new [ModifiedGoodingIOD] object given a list of
  /// [ObservationOptical] objects and an optional gravitational
  /// parameter [_mu].
  ModifiedGoodingIOD(this._observations, [this._mu = Earth.mu]);
  final List<ObservationOptical> _observations;

  /// Gravitational parameter _(km²/s³)_.
  final double _mu;

  J2000 _createInitial(
      final double r0, final double rN, final int nRev, final bool direction) {
    final iod = GoodingIOD(_observations.first,
        _observations[_observations.length ~/ 2], _observations.last, _mu);
    return iod.solve(r0, rN, nRev: nRev, direction: direction);
  }

  CostFunction _createErrorFunction(final EpochUTC aprioriEpoch) {
    final forceModel = ForceModel()..setGravity(_mu);
    double scoreFn(final Float64List x) {
      final position = Vector3D(x[0], x[1], x[2]);
      final velocity = Vector3D(x[3], x[4], x[5]);
      final state = J2000(aprioriEpoch, position, velocity);
      final propagator = RungeKutta89Propagator(state, forceModel);
      var total = 0.0;
      for (var i = 0; i < _observations.length; i++) {
        final oC = _observations[i];
        final sC = propagator.propagate(oC.epoch);
        final pC = oC.site;
        final expected = oC.observation.lineOfSight();
        final actual = RadecTopocentric.fromStateVectors(sC, pC).lineOfSight();
        final error = expected.angle(actual);
        total += error;
      }
      return total;
    }

    return scoreFn;
  }

  /// Attempt to solve a state estimate given a range estimate for the first
  /// and last observation in the observation list [r0] and [rN] _(km)_.
  J2000 solve(final double r0, final double rN,
      {final int nRev = 0,
      final bool direction = true,
      final double posSearch = 10.0,
      final double velSearch = 0.1,
      final double tolerance = 1e-6,
      final bool printIter = false}) {
    if (_observations.length < 3) {
      throw 'At least 3 observations required for Gooding IOD.';
    }
    final init = _createInitial(r0, rN, nRev, direction);
    final guess = init.position.join(init.velocity).toArray();
    final solveFn = _createErrorFunction(init.epoch);
    final simplex = <Float64List>[
      guess,
      Float64List.fromList([
        guess[0] + posSearch,
        guess[1],
        guess[2],
        guess[3],
        guess[4],
        guess[5]
      ]),
      Float64List.fromList([
        guess[0],
        guess[1] + posSearch,
        guess[2],
        guess[3],
        guess[4],
        guess[5]
      ]),
      Float64List.fromList([
        guess[0],
        guess[1],
        guess[2] + posSearch,
        guess[3],
        guess[4],
        guess[5]
      ]),
      Float64List.fromList([
        guess[0],
        guess[1],
        guess[2],
        guess[3] + velSearch,
        guess[4],
        guess[5]
      ]),
      Float64List.fromList([
        guess[0],
        guess[1],
        guess[2],
        guess[3],
        guess[4] + velSearch,
        guess[5]
      ]),
      Float64List.fromList([
        guess[0],
        guess[1],
        guess[2],
        guess[3],
        guess[4],
        guess[5] + velSearch
      ]),
    ];
    final result = DownhillSimplex.solveSimplex(solveFn, simplex,
        adaptive: true,
        xTolerance: tolerance,
        fTolerance: tolerance,
        printIter: printIter);
    return J2000(init.epoch, Vector3D(result[0], result[1], result[2]),
        Vector3D(result[3], result[4], result[5]));
  }
}
