import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Gibbs 3-position inital orbit determination.
class GibbsIOD {
  /// Create a new [GibbsIOD] object, with an optional gravitational
  /// parameter [mu].
  GibbsIOD([this.mu = Earth.mu]);

  /// Abort solve if position plane exceeds this value.
  static final double _coplanarThreshold = 5.0 * deg2rad;

  /// Gravitational parameter _(km²/s³)_.
  final double mu;

  /// Attempt to create a state estimate from three inertial position vectors.
  ///
  /// Throws an error if the positions are not coplanar.
  J2000 solve(final Vector3D r1, final Vector3D r2, final Vector3D r3,
      final EpochUTC t2, final EpochUTC t3) {
    final num = r1.normalize().dot(r2.normalize().cross(r3.normalize()));
    final alpha = halfPi - acos(num);
    if (alpha.abs() > _coplanarThreshold) {
      throw 'Orbits are not coplanar.';
    }

    final r1m = r1.magnitude();
    final r2m = r2.magnitude();
    final r3m = r3.magnitude();

    final d = r1.cross(r2).add(r2.cross(r3).add(r3.cross(r1)));
    final n = r2
        .cross(r3)
        .scale(r1m)
        .add(r3.cross(r1).scale(r2m))
        .add(r1.cross(r2).scale(r3m));
    final b = d.cross(r2);
    final s =
        r1.scale(r2m - r3m).add(r2.scale(r3m - r1m).add(r3.scale(r1m - r2m)));

    final nm = n.magnitude();
    final dm = d.magnitude();

    final vm = sqrt(mu / (nm * dm));
    final vlEci = b.scale(vm / r2m).add(s.scale(vm));

    final pv = J2000(t2, r2, vlEci);

    final forceModel = ForceModel()..setGravity(mu);
    final orbit = RungeKutta89Propagator(pv, forceModel);

    final pv2 = J2000(t2, r2, vlEci.negate());
    final orbit2 = RungeKutta89Propagator(pv2, forceModel);

    final estP3 = orbit.propagate(t3).position;
    final dist = estP3.subtract(r3).magnitude();
    final estP3_2 = orbit2.propagate(t3).position;
    final dist2 = estP3_2.subtract(r3).magnitude();

    if (dist <= dist2) {
      orbit.reset();
      return orbit.propagate(t2);
    }
    orbit2.reset();
    return orbit2.propagate(t2);
  }
}
