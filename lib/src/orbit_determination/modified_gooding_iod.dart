import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/orbit_determination/orbit_determination_base.dart';

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

  /// Attempt to solve a state estimate given a range estimate for the first
  /// and last observation in the observation list [r0] and [rN] _(km)_.
  J2000 solve(final double r0, final double rN,
      {final int nRev = 0,
      final bool direction = true,
      final bool printIter = false}) {
    if (_observations.length < 3) {
      throw 'At least 3 observations required for Gooding IOD.';
    }
    final init = _createInitial(r0, rN, nRev, direction);
    final fm = ForceModel()..setGravity();
    final result = GaussNewtonOD.solveOptical(_observations, init, fm,
        printIter: printIter);
    return result.state;
  }
}
