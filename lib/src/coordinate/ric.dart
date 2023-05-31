import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Radial-Intrack-Crosstrack relative state.
class RIC extends RelativeState {
  /// Create a new [RIC] object, from relative [position] and [velocity] _(km)_.
  RIC(final Vector position, final Vector velocity) : super(position, velocity);

  /// Create a new [RIC] object given an inertial [state] vector and its
  /// relative motion origin.
  factory RIC.fromJ2000(final J2000 state, final J2000 origin) {
    final q = RelativeState.createMatrix(origin.position, origin.velocity);
    final dr = state.position.subtract(origin.position);
    final dv = state.velocity.subtract(origin.velocity);
    return RIC(q.multiplyVector(dr), q.multiplyVector(dv));
  }

  @override
  String get name => 'RIC';

  @override
  J2000 toJ2000(final J2000 origin) {
    final qT = RelativeState.createMatrix(origin.position, origin.velocity)
        .transpose();
    final tr = qT.multiplyVector(position);
    final tv = qT.multiplyVector(velocity);
    return J2000(
        origin.epoch, origin.position.add(tr), origin.velocity.add(tv));
  }
}
