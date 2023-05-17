import 'dart:math';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/sgp4/pointer.dart';

/// ----------------------------------------------------------------------------
///
///                           procedure dscom
///
///  this procedure provides deep space common items used by both the secular
///    and periodics subroutines.  input is provided as shown. this routine
///    used to be called dpper, but the functions inside weren't well organized.
///
///  author        : david vallado                  719-573-2600   28 jun 2005
///
///  inputs        :
///    epoch       -
///    ep          - eccentricity
///    argpp       - argument of perigee
///    tc          -
///    inclp       - inclination
///    nodep       - right ascension of ascending node
///    np          - mean motion
///
///  outputs       :
///    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
///    day         -
///    e3          -
///    ee2         -
///    em          - eccentricity
///    emsq        - eccentricity squared
///    gam         -
///    peo         -
///    pgho        -
///    pho         -
///    pinco       -
///    plo         -
///    rtemsq      -
///    se2, se3         -
///    sgh2, sgh3, sgh4        -
///    sh2, sh3, si2, si3, sl2, sl3, sl4         -
///    s1, s2, s3, s4, s5, s6, s7          -
///    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
///    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
///    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
///    nm          - mean motion
///    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
///    zmol        -
///    zmos        -
///
///  locals        :
///    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
///    betasq      -
///    cc          -
///    ctem, stem        -
///    x1, x2, x3, x4, x5, x6, x7, x8          -
///    xnodce      -
///    xnoi        -
///    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
///    zcosi  , zsini  , zcosil , zsinil ,
///    zx          -
///    zy          -
///
///  coupling      :
///    none.
///
///  references    :
///    hoots, roehrich, norad spacetrack report #3 1980
///    hoots, norad spacetrack report #6 1986
///    hoots, schumacher and glover 2004
///    vallado, crawford, hujsak, kelso  2006
/// ----------------------------------------------------------------------------
void dscom(
    final double epoch,
    final double ep,
    final double argpp,
    final double tc,
    final double inclp,
    final double nodep,
    final double np,
    final Pointer<double> snodm,
    final Pointer<double> cnodm,
    final Pointer<double> sinim,
    final Pointer<double> cosim,
    final Pointer<double> sinomm,
    final Pointer<double> cosomm,
    final Pointer<double> day,
    final Pointer<double> e3,
    final Pointer<double> ee2,
    final Pointer<double> em,
    final Pointer<double> emsq,
    final Pointer<double> gam,
    final Pointer<double> peo,
    final Pointer<double> pgho,
    final Pointer<double> pho,
    final Pointer<double> pinco,
    final Pointer<double> plo,
    final Pointer<double> rtemsq,
    final Pointer<double> se2,
    final Pointer<double> se3,
    final Pointer<double> sgh2,
    final Pointer<double> sgh3,
    final Pointer<double> sgh4,
    final Pointer<double> sh2,
    final Pointer<double> sh3,
    final Pointer<double> si2,
    final Pointer<double> si3,
    final Pointer<double> sl2,
    final Pointer<double> sl3,
    final Pointer<double> sl4,
    final Pointer<double> s1,
    final Pointer<double> s2,
    final Pointer<double> s3,
    final Pointer<double> s4,
    final Pointer<double> s5,
    final Pointer<double> s6,
    final Pointer<double> s7,
    final Pointer<double> ss1,
    final Pointer<double> ss2,
    final Pointer<double> ss3,
    final Pointer<double> ss4,
    final Pointer<double> ss5,
    final Pointer<double> ss6,
    final Pointer<double> ss7,
    final Pointer<double> sz1,
    final Pointer<double> sz2,
    final Pointer<double> sz3,
    final Pointer<double> sz11,
    final Pointer<double> sz12,
    final Pointer<double> sz13,
    final Pointer<double> sz21,
    final Pointer<double> sz22,
    final Pointer<double> sz23,
    final Pointer<double> sz31,
    final Pointer<double> sz32,
    final Pointer<double> sz33,
    final Pointer<double> xgh2,
    final Pointer<double> xgh3,
    final Pointer<double> xgh4,
    final Pointer<double> xh2,
    final Pointer<double> xh3,
    final Pointer<double> xi2,
    final Pointer<double> xi3,
    final Pointer<double> xl2,
    final Pointer<double> xl3,
    final Pointer<double> xl4,
    final Pointer<double> nm,
    final Pointer<double> z1,
    final Pointer<double> z2,
    final Pointer<double> z3,
    final Pointer<double> z11,
    final Pointer<double> z12,
    final Pointer<double> z13,
    final Pointer<double> z21,
    final Pointer<double> z22,
    final Pointer<double> z23,
    final Pointer<double> z31,
    final Pointer<double> z32,
    final Pointer<double> z33,
    final Pointer<double> zmol,
    final Pointer<double> zmos) {
  /* -------------------------- constants ------------------------- */
  final zes = 0.01675;
  final zel = 0.05490;
  final c1ss = 2.9864797e-6;
  final c1l = 4.7968065e-7;
  final zsinis = 0.39785416;
  final zcosis = 0.91744867;
  final zcosgs = 0.1945905;
  final zsings = -0.98088458;

  /* --------------------- local variables ------------------------ */
  int lsflg;
  double a1,
      a2,
      a3,
      a4,
      a5,
      a6,
      a7,
      a8,
      a9,
      a10,
      betasq,
      cc,
      ctem,
      stem,
      x1,
      x2,
      x3,
      x4,
      x5,
      x6,
      x7,
      x8,
      xnodce,
      xnoi,
      zcosg,
      zcosgl,
      zcosh,
      zcoshl,
      zcosi,
      zcosil,
      zsing,
      zsingl,
      zsinh,
      zsinhl,
      zsini,
      zsinil,
      zx,
      zy;

  nm.value = np;
  em.value = ep;
  snodm.value = sin(nodep);
  cnodm.value = cos(nodep);
  sinomm.value = sin(argpp);
  cosomm.value = cos(argpp);
  sinim.value = sin(inclp);
  cosim.value = cos(inclp);
  emsq.value = em.value * em.value;
  betasq = 1.0 - emsq.value;
  rtemsq.value = sqrt(betasq);

  /* ----------------- initialize lunar solar terms --------------- */
  peo.value = 0.0;
  pinco.value = 0.0;
  plo.value = 0.0;
  pgho.value = 0.0;
  pho.value = 0.0;
  day.value = epoch + 18261.5 + tc / 1440.0;
  xnodce = (4.5236020 - 9.2422029e-4 * day.value) % twoPi;
  stem = sin(xnodce);
  ctem = cos(xnodce);
  zcosil = 0.91375164 - 0.03568096 * ctem;
  zsinil = sqrt(1.0 - zcosil * zcosil);
  zsinhl = 0.089683511 * stem / zsinil;
  zcoshl = sqrt(1.0 - zsinhl * zsinhl);
  gam.value = 5.8351514 + 0.0019443680 * day.value;
  zx = 0.39785416 * stem / zsinil;
  zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
  zx = atan2(zx, zy);
  zx = gam.value + zx - xnodce;
  zcosgl = cos(zx);
  zsingl = sin(zx);

  /* ------------------------- do solar terms --------------------- */
  zcosg = zcosgs;
  zsing = zsings;
  zcosi = zcosis;
  zsini = zsinis;
  zcosh = cnodm.value;
  zsinh = snodm.value;
  cc = c1ss;
  xnoi = 1.0 / nm.value;

  for (lsflg = 1; lsflg <= 2; lsflg++) {
    a1 = zcosg * zcosh + zsing * zcosi * zsinh;
    a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
    a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
    a8 = zsing * zsini;
    a9 = zsing * zsinh + zcosg * zcosi * zcosh;
    a10 = zcosg * zsini;
    a2 = cosim.value * a7 + sinim.value * a8;
    a4 = cosim.value * a9 + sinim.value * a10;
    a5 = -sinim.value * a7 + cosim.value * a8;
    a6 = -sinim.value * a9 + cosim.value * a10;

    x1 = a1 * cosomm.value + a2 * sinomm.value;
    x2 = a3 * cosomm.value + a4 * sinomm.value;
    x3 = -a1 * sinomm.value + a2 * cosomm.value;
    x4 = -a3 * sinomm.value + a4 * cosomm.value;
    x5 = a5 * sinomm.value;
    x6 = a6 * sinomm.value;
    x7 = a5 * cosomm.value;
    x8 = a6 * cosomm.value;

    z31.value = 12.0 * x1 * x1 - 3.0 * x3 * x3;
    z32.value = 24.0 * x1 * x2 - 6.0 * x3 * x4;
    z33.value = 12.0 * x2 * x2 - 3.0 * x4 * x4;
    z1.value = 3.0 * (a1 * a1 + a2 * a2) + z31.value * emsq.value;
    z2.value = 6.0 * (a1 * a3 + a2 * a4) + z32.value * emsq.value;
    z3.value = 3.0 * (a3 * a3 + a4 * a4) + z33.value * emsq.value;
    z11.value = -6.0 * a1 * a5 + emsq.value * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
    z12.value = -6.0 * (a1 * a6 + a3 * a5) +
        emsq.value * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
    z13.value = -6.0 * a3 * a6 + emsq.value * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
    z21.value = 6.0 * a2 * a5 + emsq.value * (24.0 * x1 * x5 - 6.0 * x3 * x7);
    z22.value = 6.0 * (a4 * a5 + a2 * a6) +
        emsq.value * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
    z23.value = 6.0 * a4 * a6 + emsq.value * (24.0 * x2 * x6 - 6.0 * x4 * x8);
    z1.value = z1.value + z1.value + betasq * z31.value;
    z2.value = z2.value + z2.value + betasq * z32.value;
    z3.value = z3.value + z3.value + betasq * z33.value;
    s3.value = cc * xnoi;
    s2.value = -0.5 * s3.value / rtemsq.value;
    s4.value = s3.value * rtemsq.value;
    s1.value = -15.0 * em.value * s4.value;
    s5.value = x1 * x3 + x2 * x4;
    s6.value = x2 * x3 + x1 * x4;
    s7.value = x2 * x4 - x1 * x3;

    /* ----------------------- do lunar terms ------------------- */
    if (lsflg == 1) {
      ss1.value = s1.value;
      ss2.value = s2.value;
      ss3.value = s3.value;
      ss4.value = s4.value;
      ss5.value = s5.value;
      ss6.value = s6.value;
      ss7.value = s7.value;
      sz1.value = z1.value;
      sz2.value = z2.value;
      sz3.value = z3.value;
      sz11.value = z11.value;
      sz12.value = z12.value;
      sz13.value = z13.value;
      sz21.value = z21.value;
      sz22.value = z22.value;
      sz23.value = z23.value;
      sz31.value = z31.value;
      sz32.value = z32.value;
      sz33.value = z33.value;
      zcosg = zcosgl;
      zsing = zsingl;
      zcosi = zcosil;
      zsini = zsinil;
      zcosh = zcoshl * cnodm.value + zsinhl * snodm.value;
      zsinh = snodm.value * zcoshl - cnodm.value * zsinhl;
      cc = c1l;
    }
  }

  zmol.value = (4.7199672 + 0.22997150 * day.value - gam.value) % twoPi;
  zmos.value = (6.2565837 + 0.017201977 * day.value) % twoPi;

  /* ------------------------ do solar terms ---------------------- */
  se2.value = 2.0 * ss1.value * ss6.value;
  se3.value = 2.0 * ss1.value * ss7.value;
  si2.value = 2.0 * ss2.value * sz12.value;
  si3.value = 2.0 * ss2.value * (sz13.value - sz11.value);
  sl2.value = -2.0 * ss3.value * sz2.value;
  sl3.value = -2.0 * ss3.value * (sz3.value - sz1.value);
  sl4.value = -2.0 * ss3.value * (-21.0 - 9.0 * emsq.value) * zes;
  sgh2.value = 2.0 * ss4.value * sz32.value;
  sgh3.value = 2.0 * ss4.value * (sz33.value - sz31.value);
  sgh4.value = -18.0 * ss4.value * zes;
  sh2.value = -2.0 * ss2.value * sz22.value;
  sh3.value = -2.0 * ss2.value * (sz23.value - sz21.value);

  /* ------------------------ do lunar terms ---------------------- */
  ee2.value = 2.0 * s1.value * s6.value;
  e3.value = 2.0 * s1.value * s7.value;
  xi2.value = 2.0 * s2.value * z12.value;
  xi3.value = 2.0 * s2.value * (z13.value - z11.value);
  xl2.value = -2.0 * s3.value * z2.value;
  xl3.value = -2.0 * s3.value * (z3.value - z1.value);
  xl4.value = -2.0 * s3.value * (-21.0 - 9.0 * emsq.value) * zel;
  xgh2.value = 2.0 * s4.value * z32.value;
  xgh3.value = 2.0 * s4.value * (z33.value - z31.value);
  xgh4.value = -18.0 * s4.value * zel;
  xh2.value = -2.0 * s2.value * z22.value;
  xh3.value = -2.0 * s2.value * (z23.value - z21.value);

  //#include "debug2.cpp"
} // dscom
