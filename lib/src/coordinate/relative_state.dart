import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Base class for relative states.
abstract class RelativeState {
  /// Create a new [RelativeState] object.
  RelativeState(this.position, this.velocity);

  /// Relative position vector _(km)_.
  final Vector position;

  /// Relative velocity vector _(km)_.
  final Vector velocity;

  /// Return the name of this coordinate frame.
  String get name;

  /// Convert this to a [J2000] state vector object.
  J2000 toJ2000(final J2000 origin);

  /// Create a relative frame transform matrix from an inertial [position] and
  /// [velocity] vector _(km)_.
  static Matrix createMatrix(final Vector position, final Vector velocity) {
    final ru = position.normalize();
    final cu = position.cross(velocity).normalize();
    final iu = cu.cross(ru).normalize();
    return Matrix([
      [ru.x, ru.y, ru.z],
      [iu.x, iu.y, iu.z],
      [cu.x, cu.y, cu.z]
    ]);
  }

  /// Return the relative range _(km)_.
  double range() => position.magnitude();

  /// Return the relative range rate _(km/s)_.
  double rangeRate() {
    final direction = (position.angle(velocity) > halfPi) ? -1.0 : 1.0;
    return velocity.magnitude() * direction;
  }
}
