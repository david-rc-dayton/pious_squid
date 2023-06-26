import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Third-body gravity model.
class ThirdBodyGravity implements Force {
  /// Create a new [ThirdBodyGravity] object with the selected bodies enabled.
  ThirdBodyGravity({this.moon = false, this.sun = false});

  /// Moon gravity enabled if `true`.
  final bool moon;

  /// Sun gravity enabled if `true`.
  final bool sun;

  static Vector3D _moonGravity(final J2000 state) {
    final rMoon = Moon.position(state.epoch);
    final aNum = rMoon.subtract(state.position);
    final aDen = pow(aNum.magnitude(), 3);
    final bNum = rMoon;
    final bDen = pow(rMoon.magnitude(), 3);
    final gravity = aNum.scale(1 / aDen).add(bNum.scale(-1 / bDen));
    return gravity.scale(Moon.mu);
  }

  static Vector3D _sunGravity(final J2000 state) {
    final rSun = Sun.positionApparent(state.epoch);
    final aNum = rSun.subtract(state.position);
    final aDen = pow(aNum.magnitude(), 3);
    final bNum = rSun;
    final bDen = pow(rSun.magnitude(), 3);
    final gravity = aNum.scale(1 / aDen).add(bNum.scale(-1 / bDen));
    return gravity.scale(Sun.mu);
  }

  @override
  Vector3D acceleration(final J2000 state) {
    var accVec = Vector3D.origin;
    if (moon) {
      accVec = accVec.add(_moonGravity(state));
    }
    if (sun) {
      accVec = accVec.add(_sunGravity(state));
    }
    return accVec;
  }
}
