import 'dart:typed_data';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/data/data_handler.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Complex Earth gravity model, accounting for EGM-96 zonal, sectoral, and
/// tesseral geopotential perturbations.
class EarthGravity implements Force {
  /// Create a new [EarthGravity] object with the given [degree] and [order].
  ///
  /// Degree and order should be a number between `0` and `36`. A degree of `0`
  /// will model a spherical Earth, and an order of `0` will only model zonal
  /// perturbations. Zonal, sectoral, and tesseral perturbations will be modeled
  /// if the degree and order are both non-zero.
  EarthGravity(final int degree, final int order)
      : degree = degree.clamp(0, 36),
        order = order.clamp(0, 36) {
    _asphericalFlag = degree >= 2;
  }

  /// Geopotential degree _(zonal)_.
  final int degree;

  /// Geopotential order _(sectoral)_.
  final int order;

  /// Use aspherical perturbations if `true`.
  late final bool _asphericalFlag;

  /// Calculate the inertial acceleration vector _(km/s²)_ due to
  /// a spherical Earth for the given [state] vector.
  Vector _spherical(final J2000 state) {
    final rMag = state.position.magnitude();
    return state.position.scale(-Earth.mu / (rMag * rMag * rMag));
  }

  /// Calculate the inertial acceleration vector delta _(km/s²)_ between a
  /// spherical Earth and aspherical Earth for the given [state] vector.
  Vector _aspherical(final J2000 state) {
    final posEcef = state.toITRF().position;
    final ri = 1.0 / posEcef.magnitude();
    final xor = posEcef.x * ri;
    final yor = posEcef.y * ri;
    final zor = posEcef.z * ri;

    final ep = zor;
    final reor = Earth.radiusEquator * ri;
    var reorn = reor;
    final muor2 = Earth.mu * ri * ri;

    var sumH = 0.0;
    var sumGm = 0.0;
    var sumJ = 0.0;
    var sumK = 0.0;

    final cTil = Float64List(order + 4);
    final sTil = Float64List(order + 4);

    final pN = Float64List(order + 4);
    final pNm1 = Float64List(order + 4);
    final pNm2 = Float64List(order + 4);

    pNm2[0] = 1.0;
    pNm1[0] = ep;
    pNm1[1] = 1.0;
    cTil[0] = 1.0;
    cTil[1] = xor;
    sTil[1] = yor;

    final dh = DataHandler();
    for (var nm2 = 0, nm1 = 1, n = 2, np1 = 3;
        n <= degree;
        nm2++, nm1++, n++, np1++) {
      final twonm1 = 2.0 * n - 1.0;

      reorn *= reor;
      final cN0 = dh.getEgm96Coeffs(n, 0).clm;

      pN[0] = (twonm1 * ep * pNm1[0] - nm1 * pNm2[0]) / n;
      pN[1] = pNm2[1] + twonm1 * pNm1[0];
      pN[2] = pNm2[2] + twonm1 * pNm1[1];

      var sumHn = pN[1] * cN0;
      var sumGmn = pN[0] * cN0 * np1;

      if (order > 0) {
        var sumJn = 0.0;
        var sumKn = 0.0;

        cTil[n] = cTil[1] * cTil[nm1] - sTil[1] * sTil[nm1];
        sTil[n] = sTil[1] * cTil[nm1] + cTil[1] * sTil[nm1];

        final lim = (n < order) ? n : order;

        for (var mm2 = -1, mm1 = 0, m = 1, mp1 = 2, mp2 = 3;
            m <= lim;
            mm2++, mm1++, m++, mp1++, mp2++) {
          pN[mp1] = pNm2[mp1] + twonm1 * pNm1[m];

          final dm = m;
          final npmp1 = n + mp1;

          final pNm = pN[m];
          final pNmp1 = pN[mp1];

          final coeffs = dh.getEgm96Coeffs(n, m);
          final cNm = coeffs.clm;
          final sNm = coeffs.slm;

          final mxPnm = dm * pNm;
          final bNmtil = cNm * cTil[m] + sNm * sTil[m];
          final pNmBnm = pNm * bNmtil;
          final bNmtm1 = cNm * cTil[mm1] + sNm * sTil[mm1];
          final aNmtm1 = cNm * sTil[mm1] - sNm * cTil[mm1];

          sumHn += pNmp1 * bNmtil;
          sumGmn += npmp1 * pNmBnm;
          sumJn += mxPnm * bNmtm1;
          sumKn -= mxPnm * aNmtm1;
        }

        sumJ += reorn * sumJn;
        sumK += reorn * sumKn;
      }

      sumH += reorn * sumHn;
      sumGm += reorn * sumGmn;

      if (n < degree) {
        for (var i = 0; i <= n; i++) {
          pNm2[i] = pNm1[i];
          pNm1[i] = pN[i];
        }
      }
    }

    final lambda = sumGm + ep * sumH;
    final g = Float64List(3);
    g[0] = -muor2 * (lambda * xor - sumJ);
    g[1] = -muor2 * (lambda * yor - sumK);
    g[2] = -muor2 * (lambda * zor - sumH);
    return ITRF(state.epoch, Vector(g), Vector.origin3).toJ2000().position;
  }

  @override
  Vector acceleration(final J2000 state) {
    var accVec = _spherical(state);
    if (_asphericalFlag) {
      accVec = accVec.add(_aspherical(state));
    }
    return accVec;
  }
}
