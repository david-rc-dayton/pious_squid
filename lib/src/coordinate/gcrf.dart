import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';

/// GCRF state vector.
class GCRF extends StateVector {
  /// Create a new [GCRF] state vector object.
  GCRF(super.epoch, super.position, super.velocity);

  /// Create a new [GCRF] object from a [ClassicalElements] object.
  factory GCRF.fromClassicalElements(final ClassicalElements elements) {
    final rv = elements.toPositionVelocity();
    return GCRF(elements.epoch, rv.position, rv.velocity);
  }

  @override
  String get name => 'GCRF';

  @override
  bool get inertial => true;

  /// Convert this to a [J2000] state vector object.
  J2000 toJ2000() {
    final p = Earth.precession(epoch);

    final n0 = Earth.nutation(epoch, coeffs: 106, useEop: true);
    final n1 = Earth.nutation(epoch, coeffs: 106);

    var rMOD = position.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    var vMOD = velocity.rotZ(-p.zeta).rotY(p.theta).rotZ(-p.zed);
    final rTOD = rMOD.rotX(n0.mEps).rotZ(-n0.dPsi).rotX(-n0.eps);
    final vTOD = vMOD.rotX(n0.mEps).rotZ(-n0.dPsi).rotX(-n0.eps);

    rMOD = rTOD.rotX(n1.eps).rotZ(n1.dPsi).rotX(-n1.mEps);
    vMOD = vTOD.rotX(n1.eps).rotZ(n1.dPsi).rotX(-n1.mEps);
    final rJ2000 = rMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);
    final vJ2000 = vMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);

    return J2000(epoch, rJ2000, vJ2000);
  }
}
