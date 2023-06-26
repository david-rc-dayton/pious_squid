import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Herrik-Gibbs 3-position initial orbit determination.
///
/// Possibly better than regular Gibbs IOD for closely spaced position
/// vectors _(less than 5°)_.
class HerrickGibbsIOD {
  /// Create a new [HerrickGibbsIOD] object with optional
  /// gravitational parameter [mu].
  HerrickGibbsIOD([this.mu = Earth.mu]);

  /// Gravitation parameter _(km²/s³)_.
  final double mu;

  /// Attempt to create a state estimate from three inertial position vectors.
  J2000 solve(final Vector3D r1, final EpochUTC t1, final Vector3D r2,
      final EpochUTC t2, final Vector3D r3, final EpochUTC t3) {
    final dt31 = t3.difference(t1);
    final dt32 = t3.difference(t2);
    final dt21 = t2.difference(t1);
    final r1m = r1.magnitude();
    final r2m = r2.magnitude();
    final r3m = r3.magnitude();
    final vA = r1.scale(
        -dt32 * ((1.0 / (dt21 * dt31)) + (mu / (12.0 * r1m * r1m * r1m))));
    final vB = r2.scale((dt32 - dt21) *
        ((1.0 / (dt21 * dt32)) + (mu / (12.0 * r2m * r2m * r2m))));
    final vC = r3.scale(
        dt21 * ((1.0 / (dt32 * dt31)) + (mu / (12.0 * r3m * r3m * r3m))));
    final v2 = vA.add(vB).add(vC);
    return J2000(t2, r2, v2);
  }
}
