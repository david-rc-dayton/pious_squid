import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Radial-Intrack-Crosstrack relative state.
class RIC extends RelativeState {
  /// Create a new [RIC] object, from relative [position] and [velocity] _(km)_.
  RIC(super.position, super.velocity);

  /// Create a new [RIC] object given an inertial [state] vector, its relative
  /// motion origin, and its tranformation matrix.
  factory RIC.fromJ2000Matrix(
      final J2000 state, final J2000 origin, final Matrix transform) {
    final dr = state.position.subtract(origin.position);
    final dv = state.velocity.subtract(origin.velocity);
    return RIC(transform.multiplyVector3D(dr), transform.multiplyVector3D(dv));
  }

  /// Create a new [RIC] object given an inertial [state] vector and its
  /// relative motion origin.
  factory RIC.fromJ2000(final J2000 state, final J2000 origin) =>
      RIC.fromJ2000Matrix(state, origin,
          RelativeState.createMatrix(origin.position, origin.velocity));

  @override
  String get name => 'RIC';

  /// Convert this state to J2000 given an [origin] state and [tranform]
  /// matrix.
  J2000 toJ2000Matrix(final J2000 origin, final Matrix transform) {
    final tt = transform.transpose();
    final tr = tt.multiplyVector3D(position);
    final tv = tt.multiplyVector3D(velocity);
    return J2000(
        origin.epoch, origin.position.add(tr), origin.velocity.add(tv));
  }

  @override
  J2000 toJ2000(final J2000 origin) => toJ2000Matrix(
      origin, RelativeState.createMatrix(origin.position, origin.velocity));
}
