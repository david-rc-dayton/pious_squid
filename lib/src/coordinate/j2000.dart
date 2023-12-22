import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/data/data_base.dart';

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
    final w = DataHandler().getEop(epoch);
    final ast = epoch.gmstAngle() + n.eqEq;
    final rMOD = position.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    final vMOD = velocity.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    final rTOD = rMOD.rotX(n.mEps).rotZ(-n.dPsi).rotX(-n.eps);
    final vTOD = vMOD.rotX(n.mEps).rotZ(-n.dPsi).rotX(-n.eps);
    final rPEF = rTOD.rotZ(ast);
    final vPEF =
        vTOD.rotZ(ast).add(Earth.rotationLod(epoch).negate().cross(rPEF));
    final rITRF = rPEF.rotX(-w.y).rotY(-w.x);
    final vITRF = vPEF.rotX(-w.y).rotY(-w.x);
    return ITRF(epoch, rITRF, vITRF);
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

  /// Convert this to a [GCRF] state vector object.
  GCRF toGCRF() {
    final p = Earth.precession(epoch);

    final n0 = Earth.nutation(epoch, coeffs: 106);
    final n1 = Earth.nutation(epoch, coeffs: 106, useEop: true);

    var rMOD = position.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    var vMOD = velocity.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    final rTOD = rMOD.rotX(n0.mEps).rotZ(-n0.dPsi).rotX(-n0.eps);
    final vTOD = vMOD.rotX(n0.mEps).rotZ(-n0.dPsi).rotX(-n0.eps);

    rMOD = rTOD.rotX(n1.eps).rotZ(n1.dPsi).rotX(-n1.mEps);
    vMOD = vTOD.rotX(n1.eps).rotZ(n1.dPsi).rotX(-n1.mEps);
    final rGCRF = rMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);
    final vGCRF = vMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);

    return GCRF(epoch, rGCRF, vGCRF);
  }
}
