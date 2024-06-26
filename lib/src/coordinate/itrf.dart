import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/data/data_base.dart';

/// International Terrestrial Reference Frame _(ITRF)_
class ITRF extends StateVector {
  /// Create a new [ITRF] state vector object.
  ITRF(super.epoch, super.position, super.velocity);

  @override
  String get name => 'ITRF';

  @override
  bool get inertial => false;

  /// Return height above the Earth's surface _(km)_.
  double getHeight() {
    final a = Earth.radiusEquator;
    final e2 = Earth.eccentricitySquared;
    final r = position.magnitude();
    final sl = position.z / r;
    final cl2 = 1 - sl * sl;
    final coeff = sqrt((1 - e2) / (1 - e2 * cl2));
    return r - a * coeff;
  }

  /// Convert this to a [J2000] state vector object.
  J2000 toJ2000() {
    final p = Earth.precession(epoch);
    final n = Earth.nutation(epoch);
    final w = DataHandler().getEop(epoch);
    final ast = epoch.gmstAngle() + n.eqEq;
    final rPEF = position.rotY(w.x).rotX(w.y);
    final vPEF = velocity.rotY(w.x).rotX(w.y);
    final rTOD = rPEF.rotZ(-ast);
    final vTOD = vPEF.add(Earth.rotationLod(epoch).cross(rPEF)).rotZ(-ast);
    final rMOD = rTOD.rotX(n.eps).rotZ(n.dPsi).rotX(-n.mEps);
    final vMOD = vTOD.rotX(n.eps).rotZ(n.dPsi).rotX(-n.mEps);
    final rJ2000 = rMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);
    final vJ2000 = vMOD.rotZ(p.zed).rotY(-p.theta).rotZ(p.zeta);
    return J2000(epoch, rJ2000, vJ2000);
  }

  /// Convert this to a [Geodetic] coordinate object.
  Geodetic toGeodetic() {
    final sma = Earth.radiusEquator;
    final esq = Earth.eccentricitySquared;
    final x = position.x;
    final y = position.y;
    final z = position.z;
    final lon = atan2(y, x);
    final r = sqrt(x * x + y * y);
    final phi = atan(z / r);
    var lat = phi;
    var c = 0.0;
    for (var i = 0; i < 12; i++) {
      final slat = sin(lat);
      c = sma / sqrt(1 - esq * slat * slat);
      lat = atan((z + c * esq * slat) / r);
    }
    final alt = (pi * 0.5 - lat.abs() > (pi / 180))
        ? r / cos(lat) - c
        : z / sin(lat) - (c * (1 - esq));
    return Geodetic(lat, lon, alt);
  }
}
