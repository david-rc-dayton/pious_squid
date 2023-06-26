import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Simple central-body gravity model.
class Gravity implements Force {
  /// Create a new [Gravity] object with optional gravitational
  /// parameter [mu] _(km²/s³)_.
  Gravity([this.mu = Earth.mu]);

  /// Gravitational parameter _(km²/s³)_.
  final double mu;

  /// Calculate the inertial acceleration vector _(km/s²)_ due to
  /// a spherical central body for the given [state] vector.
  Vector3D _spherical(final J2000 state) {
    final rMag = state.position.magnitude();
    return state.position.scale(-mu / (rMag * rMag * rMag));
  }

  @override
  Vector3D acceleration(final J2000 state) => _spherical(state);
}
