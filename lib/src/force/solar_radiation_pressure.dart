import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Solar radiation pressure model.
class SolarRadiationPressure extends Force {
  /// Create a new [SolarRadiationPressure] object.
  SolarRadiationPressure(this.mass, this.area, this.reflectCoeff);

  /// Spacecraft mass _(kg)_.
  final double mass;

  /// Spacecraft cross-sectional area _(m²)_.
  final double area;

  /// Reflectivity coefficient (unitless).
  final double reflectCoeff;

  /// Solar pressure _(N/m²)_;
  static final double _kRef = 4.56e-6 * pow(astronomicalUnit, 2);

  @override
  Vector3D acceleration(final J2000 state) {
    final rSun = Sun.positionApparent(state.epoch);
    final r = state.position.subtract(rSun);
    final rMag = r.magnitude();
    final r2 = rMag * rMag;
    final ratio = Sun.lightingRatio(state.position, rSun);
    final p = ratio * _kRef / r2;
    final flux = r.scale(p / rMag);
    return flux.scale((area * reflectCoeff / mass) * 1e-3);
  }
}
