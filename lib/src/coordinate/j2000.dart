import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';

/// J2000 state vector.
class J2000 extends StateVector {
  /// Create a new [J2000] state vector object.
  J2000(super.epoch, super.position, super.velocity);

  /// Create a new [J2000] object from a [ClassicalElements] object.
  factory J2000.fromClassicalElements(final ClassicalElements elements) {
    final rv = elements.toPositionVelocity();
    return J2000(elements.epoch, rv.position, rv.velocity);
  }

  @override
  String get name => 'J2000';

  @override
  bool get inertial => true;

  /// Convert this to an [ITRF] state vector object.
  ITRF toITRF() {
    final p = Earth.precession(epoch);
    final n = Earth.nutation(epoch);
    final ast = epoch.gmstAngle() + n.eqEq;
    final rMOD = position.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    final vMOD = velocity.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    final rTOD = rMOD.rotX(n.mEps).rotZ(-n.dPsi).rotX(-n.eps);
    final vTOD = vMOD.rotX(n.mEps).rotZ(-n.dPsi).rotX(-n.eps);
    final rPEF = rTOD.rotZ(ast);
    final vPEF = vTOD.rotZ(ast).add(Earth.rotation.negate().cross(rPEF));
    return ITRF(epoch, rPEF, vPEF);
  }

  /// Convert this to an [TEME] state vector object.
  TEME toTEME() {
    final p = Earth.precession(epoch);
    final n = Earth.nutation(epoch);
    final eps = n.mEps + n.dEps;
    final dPsiCosEps = n.dPsi * cos(eps);
    final rMOD = position.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    final vMOD = velocity.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    final rTEME = rMOD.rotX(n.mEps).rotZ(-n.dPsi).rotX(-eps).rotZ(dPsiCosEps);
    final vTEME = vMOD.rotX(n.mEps).rotZ(-n.dPsi).rotX(-eps).rotZ(dPsiCosEps);
    return TEME(epoch, rTEME, vTEME);
  }
}
