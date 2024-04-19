import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/data/data_handler.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Harris-Priester atmospheric drag force model.
///
/// Atmospheric density model assumes mean solar flux.
class AtmosphericDrag implements Force {
  /// Create a new [AtmosphericDrag] object.
  ///
  /// The cosine exponent should be a number between `2` for low inclination
  /// orbits and `6` for polar orbits.
  AtmosphericDrag(this.dragCoeff, this.cosine);

  /// Drag coefficient (unitless).
  final double dragCoeff;

  /// Cosine exponent.
  final int cosine;

  /// Return Harris-Priester atmospheric density _(kg/m³)_ for a given [state]
  /// and cosine exponent [n].
  static double _getHPDensity(final ITRF state, final int n) {
    final hpa = DataHandler().getHpAtmosphere(state.epoch, state.getHeight());
    if (hpa == null) {
      return 0.0;
    }
    final sunPos = Sun.positionApparent(state.epoch);
    final sunVec = J2000(state.epoch, sunPos, Vector3D.origin)
        .toITRF()
        .position
        .normalize();
    final bulVec = sunVec.rotZ(-30.0 * deg2rad);
    final cosPsi = bulVec.normalize().dot(state.position.normalize());
    final c2Psi2 = 0.5 * (1.0 + cosPsi);
    final cPsi2 = sqrt(c2Psi2);
    final cosPow = (cPsi2 > 1e-12) ? c2Psi2 * pow(cPsi2, n - 2) : 0.0;
    final altitude = hpa.height;
    final (h: h0, minD: min0, maxD: max0) = hpa.hp0;
    final (h: h1, minD: min1, maxD: max1) = hpa.hp1;
    final dH = (h0 - altitude) / (h0 - h1);
    final rhoMin = min0 * pow(min1 / min0, dH);
    if (cosPow == 0) {
      return rhoMin;
    }
    final rhoMax = max0 * pow(max1 / max0, dH);
    return rhoMin + (rhoMax - rhoMin) * cosPow;
  }

  @override
  Vector3D acceleration(final J2000 state, [final double bcoeffInv = 0.0]) {
    if (bcoeffInv == 0.0) {
      return Vector3D.origin;
    }
    final itrfState = state.toITRF();
    final density = _getHPDensity(itrfState, cosine);
    if (density == 0) {
      return Vector3D.origin;
    }
    final rotation =
        ITRF(state.epoch, Earth.rotation, Vector3D.origin).toJ2000().position;
    final vRel =
        state.velocity.subtract(rotation.cross(state.position)).scale(1000.0);
    final vm = vRel.magnitude();
    final fScale = -0.5 * density * (dragCoeff * bcoeffInv) * vm;
    return vRel.scale(fScale / 1000.0);
  }
}
