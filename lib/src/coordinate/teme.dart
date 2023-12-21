import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';

/// True Equator Mean Equinox _(TEME)_ state vector
class TEME extends StateVector {
  /// Create a new [TEME] object given inertial [position] and [velocity].
  TEME(super.epoch, super.position, super.velocity);

  /// Create a new [TEME] object from a [ClassicalElements] object.
  factory TEME.fromClassicalElements(final ClassicalElements elements) {
    final rv = elements.toPositionVelocity();
    return TEME(elements.epoch, rv.position, rv.velocity);
  }

  @override
  String get name => 'TEME';

  @override
  bool get inertial => true;

  /// Convert this to a [J2000] state vector object.
  J2000 toJ2000() {
    final p = Earth.precession(epoch);
    final n = Earth.nutation(epoch);
    final eps = n.mEps + n.dEps;
    final dPsiCosEps = n.dPsi * cos(eps);
    final rMOD =
        position.rotZ(-dPsiCosEps).rotX(eps).rotZ(n.dPsi).rotX(-n.mEps);
    final vMOD =
        velocity.rotZ(-dPsiCosEps).rotX(eps).rotZ(n.dPsi).rotX(-n.mEps);
    final rJ2K = rMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);
    final vJ2K = vMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);
    return J2000(epoch, rJ2K, vJ2K);
  }
}
