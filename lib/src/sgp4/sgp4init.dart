import 'dart:math';

import 'package:pious_squid/src/sgp4/dpper.dart';
import 'package:pious_squid/src/sgp4/dscom.dart';
import 'package:pious_squid/src/sgp4/dsinit.dart';
import 'package:pious_squid/src/sgp4/elsetrec.dart';
import 'package:pious_squid/src/sgp4/gravconst.dart';
import 'package:pious_squid/src/sgp4/initl.dart';
import 'package:pious_squid/src/sgp4/pointer.dart';
import 'package:pious_squid/src/sgp4/sgp4.dart';

/// ----------------------------------------------------------------------------
///
///                             procedure sgp4init
///
///  this procedure initializes variables for sgp4.
///
///  author        : david vallado                  719-573-2600   28 jun 2005
///
///  inputs        :
///    opsmode     - mode of operation afspc or improved 'a', 'i'
///    whichconst  - which set of constants to use  72, 84
///    satn        - satellite number
///    bstar       - sgp4 type drag coefficient              kg/m2er
///    ecco        - eccentricity
///    epoch       - epoch time in days from jan 0, 1950. 0 hr
///    argpo       - argument of perigee (output if ds)
///    inclo       - inclination
///    mo          - mean anomaly (output if ds)
///    no          - mean motion
///    nodeo       - right ascension of ascending node
///
///  outputs       :
///    satrec      - common values for subsequent calls
///    return code - non-zero on error.
///                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
///                   2 - mean motion less than 0.0
///                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
///                   4 - semi-latus rectum < 0.0
///                   5 - epoch elements are sub-orbital
///                   6 - satellite has decayed
///
///  locals        :
///    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
///    cc1sq  , cc2    , cc3
///    coef   , coef1
///    cosio4      -
///    day         -
///    dndt        -
///    em          - eccentricity
///    emsq        - eccentricity squared
///    eeta        -
///    etasq       -
///    gam         -
///    argpm       - argument of perigee
///    nodem       -
///    inclm       - inclination
///    mm          - mean anomaly
///    nm          - mean motion
///    perige      - perigee
///    pinvsq      -
///    psisq       -
///    qzms24      -
///    rtemsq      -
///    s1, s2, s3, s4, s5, s6, s7          -
///    sfour       -
///    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
///    sz1, sz2, sz3
///    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
///    tc          -
///    temp        -
///    temp1, temp2, temp3       -
///    tsi         -
///    xpidot      -
///    xhdot1      -
///    z1, z2, z3          -
///    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
///
///  coupling      :
///    getgravconst-
///    initl       -
///    dscom       -
///    dpper       -
///    dsinit      -
///    sgp4        -
///
///  references    :
///    hoots, roehrich, norad spacetrack report #3 1980
///    hoots, norad spacetrack report #6 1986
///    hoots, schumacher and glover 2004
///    vallado, crawford, hujsak, kelso  2006
/// ----------------------------------------------------------------------------
bool sgp4init(
    final Sgp4GravConst whichconst,
    final String opsmode,
    final String satn,
    final double epoch,
    final double xbstar,
    final double xndot,
    final double xnddot,
    final double xecco,
    final double xargpo,
    final double xinclo,
    final double xmo,
    final double xnoKozai,
    final double xnodeo,
    final ElsetRec satrec) {
  /* --------------------- local variables ------------------------ */
  double cc1sq,
      cc2,
      cc3,
      coef,
      coef1,
      cosio4,
      eeta,
      etasq,
      perige,
      pinvsq,
      psisq,
      qzms24,
      sfour,
      tc,
      temp,
      temp1,
      temp2,
      temp3,
      tsi,
      xpidot,
      xhdot1,
      qzms2t,
      ss,
      x2o3,
      delmotemp,
      qzms2ttemp,
      qzms24temp;
  final r = <double>[0, 0, 0], v = <double>[0, 0, 0];
  final ainv = Pointer(0.0),
      ao = Pointer(0.0),
      con42 = Pointer(0.0),
      cosio = Pointer(0.0),
      cosio2 = Pointer(0.0),
      eccsq = Pointer(0.0),
      omeosq = Pointer(0.0),
      posq = Pointer(0.0),
      rp = Pointer(0.0),
      rteosq = Pointer(0.0),
      sinio = Pointer(0.0),
      snodm = Pointer(0.0),
      cnodm = Pointer(0.0),
      sinim = Pointer(0.0),
      cosim = Pointer(0.0),
      sinomm = Pointer(0.0),
      cosomm = Pointer(0.0),
      emsq = Pointer(0.0),
      day = Pointer(0.0),
      em = Pointer(0.0),
      gam = Pointer(0.0),
      rtemsq = Pointer(0.0),
      nm = Pointer(0.0),
      s1 = Pointer(0.0),
      s2 = Pointer(0.0),
      s3 = Pointer(0.0),
      s4 = Pointer(0.0),
      s5 = Pointer(0.0),
      s6 = Pointer(0.0),
      s7 = Pointer(0.0),
      ss1 = Pointer(0.0),
      ss2 = Pointer(0.0),
      ss3 = Pointer(0.0),
      ss4 = Pointer(0.0),
      ss5 = Pointer(0.0),
      ss6 = Pointer(0.0),
      ss7 = Pointer(0.0),
      sz1 = Pointer(0.0),
      sz2 = Pointer(0.0),
      sz3 = Pointer(0.0),
      sz11 = Pointer(0.0),
      sz12 = Pointer(0.0),
      sz13 = Pointer(0.0),
      sz21 = Pointer(0.0),
      sz22 = Pointer(0.0),
      sz23 = Pointer(0.0),
      sz31 = Pointer(0.0),
      sz32 = Pointer(0.0),
      sz33 = Pointer(0.0),
      z1 = Pointer(0.0),
      z2 = Pointer(0.0),
      z3 = Pointer(0.0),
      z11 = Pointer(0.0),
      z12 = Pointer(0.0),
      z13 = Pointer(0.0),
      z21 = Pointer(0.0),
      z22 = Pointer(0.0),
      z23 = Pointer(0.0),
      z31 = Pointer(0.0),
      z32 = Pointer(0.0),
      z33 = Pointer(0.0),
      argpm = Pointer(0.0),
      inclm = Pointer(0.0),
      mm = Pointer(0.0),
      nodem = Pointer(0.0),
      dndt = Pointer(0.0);

  /* ------------------------ initialization --------------------- */
  // sgp4fix divisor for divide by zero check on inclination
  // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
  // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
  const temp4 = 1.5e-12;

  /* ----------- set all near earth variables to zero ------------ */
  satrec.isimp.value = 0;
  satrec.method.value = 'n';
  satrec.aycof.value = 0.0;
  satrec.con41.value = 0.0;
  satrec.cc1.value = 0.0;
  satrec.cc4.value = 0.0;
  satrec.cc5.value = 0.0;
  satrec.d2.value = 0.0;
  satrec.d3.value = 0.0;
  satrec.d4.value = 0.0;
  satrec.delmo.value = 0.0;
  satrec.eta.value = 0.0;
  satrec.argpdot.value = 0.0;
  satrec.omgcof.value = 0.0;
  satrec.sinmao.value = 0.0;
  satrec.t.value = 0.0;
  satrec.t2cof.value = 0.0;
  satrec.t3cof.value = 0.0;
  satrec.t4cof.value = 0.0;
  satrec.t5cof.value = 0.0;
  satrec.x1mth2.value = 0.0;
  satrec.x7thm1.value = 0.0;
  satrec.mdot.value = 0.0;
  satrec.nodedot.value = 0.0;
  satrec.xlcof.value = 0.0;
  satrec.xmcof.value = 0.0;
  satrec.nodecf.value = 0.0;

  /* ----------- set all deep space variables to zero ------------ */
  satrec.irez.value = 0;
  satrec.d2201.value = 0.0;
  satrec.d2211.value = 0.0;
  satrec.d3210.value = 0.0;
  satrec.d3222.value = 0.0;
  satrec.d4410.value = 0.0;
  satrec.d4422.value = 0.0;
  satrec.d5220.value = 0.0;
  satrec.d5232.value = 0.0;
  satrec.d5421.value = 0.0;
  satrec.d5433.value = 0.0;
  satrec.dedt.value = 0.0;
  satrec.del1.value = 0.0;
  satrec.del2.value = 0.0;
  satrec.del3.value = 0.0;
  satrec.didt.value = 0.0;
  satrec.dmdt.value = 0.0;
  satrec.dnodt.value = 0.0;
  satrec.domdt.value = 0.0;
  satrec.e3.value = 0.0;
  satrec.ee2.value = 0.0;
  satrec.peo.value = 0.0;
  satrec.pgho.value = 0.0;
  satrec.pho.value = 0.0;
  satrec.pinco.value = 0.0;
  satrec.plo.value = 0.0;
  satrec.se2.value = 0.0;
  satrec.se3.value = 0.0;
  satrec.sgh2.value = 0.0;
  satrec.sgh3.value = 0.0;
  satrec.sgh4.value = 0.0;
  satrec.sh2.value = 0.0;
  satrec.sh3.value = 0.0;
  satrec.si2.value = 0.0;
  satrec.si3.value = 0.0;
  satrec.sl2.value = 0.0;
  satrec.sl3.value = 0.0;
  satrec.sl4.value = 0.0;
  satrec.gsto.value = 0.0;
  satrec.xfact.value = 0.0;
  satrec.xgh2.value = 0.0;
  satrec.xgh3.value = 0.0;
  satrec.xgh4.value = 0.0;
  satrec.xh2.value = 0.0;
  satrec.xh3.value = 0.0;
  satrec.xi2.value = 0.0;
  satrec.xi3.value = 0.0;
  satrec.xl2.value = 0.0;
  satrec.xl3.value = 0.0;
  satrec.xl4.value = 0.0;
  satrec.xlamo.value = 0.0;
  satrec.zmol.value = 0.0;
  satrec.zmos.value = 0.0;
  satrec.atime.value = 0.0;
  satrec.xli.value = 0.0;
  satrec.xni.value = 0.0;

  /* ------------------------ earth constants ----------------------- */
  // sgp4fix identify constants and allow alternate values
  // this is now the only call for the constants
  getgravconst(whichconst, satrec.tumin, satrec.mus, satrec.radiusearthkm,
      satrec.xke, satrec.j2, satrec.j3, satrec.j4, satrec.j3oj2);

  //-------------------------------------------------------------------------

  satrec.error.value = 0;
  satrec.operationmode.value = opsmode;
  // new alpha5 or 9-digit number
  satrec.satnum.value = satn;

  // sgp4fix - note the following variables are also passed directly via satrec.
  // it is possible to streamline the sgp4init call by deleting the "x"
  // variables, but the user would need to set the satrec.* values first. we
  // include the additional assignments in case twoline2rv is not used.
  satrec.bstar.value = xbstar;
  // sgp4fix allow additional parameters in the struct
  satrec.ndot.value = xndot;
  satrec.nddot.value = xnddot;
  satrec.ecco.value = xecco;
  satrec.argpo.value = xargpo;
  satrec.inclo.value = xinclo;
  satrec.mo.value = xmo;
  // sgp4fix rename variables to clarify which mean motion is intended
  satrec.no_kozai.value = xnoKozai;
  satrec.nodeo.value = xnodeo;

  // single averaged mean elements
  satrec.am.value = satrec.em.value = satrec.im.value =
      satrec.Om.value = satrec.mm.value = satrec.nm.value = 0.0;

  /* ------------------------ earth constants ----------------------- */
  // sgp4fix identify constants and allow alternate values no longer needed
  // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
  ss = 78.0 / satrec.radiusearthkm.value + 1.0;
  // sgp4fix use multiply for speed instead of pow
  qzms2ttemp = (120.0 - 78.0) / satrec.radiusearthkm.value;
  qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
  x2o3 = 2.0 / 3.0;

  satrec.init.value = 'y';
  satrec.t.value = 0.0;

  // sgp4fix remove satn as it is not needed in initl
  initl(
      satrec.xke.value,
      satrec.j2.value,
      satrec.ecco.value,
      epoch,
      satrec.inclo.value,
      satrec.no_kozai.value,
      satrec.operationmode.value,
      satrec.method,
      ainv,
      ao,
      satrec.con41,
      con42,
      cosio,
      cosio2,
      eccsq,
      omeosq,
      posq,
      rp,
      rteosq,
      sinio,
      satrec.gsto,
      satrec.no_unkozai);
  satrec.a.value =
      pow(satrec.no_unkozai.value * satrec.tumin.value, -2.0 / 3.0) as double;
  satrec.alta.value = satrec.a.value * (1.0 + satrec.ecco.value) - 1.0;
  satrec.altp.value = satrec.a.value * (1.0 - satrec.ecco.value) - 1.0;
  satrec.error.value = 0;

  // sgp4fix remove this check as it is unnecessary
  // the mrt check in sgp4 handles decaying satellite cases even if the starting
  // condition is below the surface of te earth
  //     if (rp < 1.0)
  //       {
  //         printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
  //         satrec.error = 5;
  //       }

  if ((omeosq.value >= 0.0) || (satrec.no_unkozai.value >= 0.0)) {
    satrec.isimp.value = 0;
    if (rp.value < (220.0 / satrec.radiusearthkm.value + 1.0)) {
      satrec.isimp.value = 1;
    }
    sfour = ss;
    qzms24 = qzms2t;
    perige = (rp.value - 1.0) * satrec.radiusearthkm.value;

    /* - for perigees below 156 km, s and qoms2t are altered - */
    if (perige < 156.0) {
      sfour = perige - 78.0;
      if (perige < 98.0) {
        sfour = 20.0;
      }
      // sgp4fix use multiply for speed instead of pow
      qzms24temp = (120.0 - sfour) / satrec.radiusearthkm.value;
      qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
      sfour = sfour / satrec.radiusearthkm.value + 1.0;
    }
    pinvsq = 1.0 / posq.value;

    tsi = 1.0 / (ao.value - sfour);
    satrec.eta.value = ao.value * satrec.ecco.value * tsi;
    etasq = satrec.eta.value * satrec.eta.value;
    eeta = satrec.ecco.value * satrec.eta.value;
    psisq = (1.0 - etasq).abs();
    coef = qzms24 * pow(tsi, 4.0);
    coef1 = coef / pow(psisq, 3.5);
    cc2 = coef1 *
        satrec.no_unkozai.value *
        (ao.value * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
            0.375 *
                satrec.j2.value *
                tsi /
                psisq *
                satrec.con41.value *
                (8.0 + 3.0 * etasq * (8.0 + etasq)));
    satrec.cc1.value = satrec.bstar.value * cc2;
    cc3 = 0.0;
    if (satrec.ecco.value > 1.0e-4) {
      cc3 = -2.0 *
          coef *
          tsi *
          satrec.j3oj2.value *
          satrec.no_unkozai.value *
          sinio.value /
          satrec.ecco.value;
    }
    satrec.x1mth2.value = 1.0 - cosio2.value;
    satrec.cc4.value = 2.0 *
        satrec.no_unkozai.value *
        coef1 *
        ao.value *
        omeosq.value *
        (satrec.eta.value * (2.0 + 0.5 * etasq) +
            satrec.ecco.value * (0.5 + 2.0 * etasq) -
            satrec.j2.value *
                tsi /
                (ao.value * psisq) *
                (-3.0 *
                        satrec.con41.value *
                        (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                    0.75 *
                        satrec.x1mth2.value *
                        (2.0 * etasq - eeta * (1.0 + etasq)) *
                        cos(2.0 * satrec.argpo.value)));
    satrec.cc5.value = 2.0 *
        coef1 *
        ao.value *
        omeosq.value *
        (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
    cosio4 = cosio2.value * cosio2.value;
    temp1 = 1.5 * satrec.j2.value * pinvsq * satrec.no_unkozai.value;
    temp2 = 0.5 * temp1 * satrec.j2.value * pinvsq;
    temp3 =
        -0.46875 * satrec.j4.value * pinvsq * pinvsq * satrec.no_unkozai.value;
    satrec.mdot.value = satrec.no_unkozai.value +
        0.5 * temp1 * rteosq.value * satrec.con41.value +
        0.0625 *
            temp2 *
            rteosq.value *
            (13.0 - 78.0 * cosio2.value + 137.0 * cosio4);
    satrec.argpdot.value = -0.5 * temp1 * con42.value +
        0.0625 * temp2 * (7.0 - 114.0 * cosio2.value + 395.0 * cosio4) +
        temp3 * (3.0 - 36.0 * cosio2.value + 49.0 * cosio4);
    xhdot1 = -temp1 * cosio.value;
    satrec.nodedot.value = xhdot1 +
        (0.5 * temp2 * (4.0 - 19.0 * cosio2.value) +
                2.0 * temp3 * (3.0 - 7.0 * cosio2.value)) *
            cosio.value;
    xpidot = satrec.argpdot.value + satrec.nodedot.value;
    satrec.omgcof.value = satrec.bstar.value * cc3 * cos(satrec.argpo.value);
    satrec.xmcof.value = 0.0;
    if (satrec.ecco.value > 1.0e-4) {
      satrec.xmcof.value = -x2o3 * coef * satrec.bstar.value / eeta;
    }
    satrec.nodecf.value = 3.5 * omeosq.value * xhdot1 * satrec.cc1.value;
    satrec.t2cof.value = 1.5 * satrec.cc1.value;
    // sgp4fix for divide by zero with xinco = 180 deg
    if ((cosio.value + 1.0).abs() > 1.5e-12) {
      satrec.xlcof.value = -0.25 *
          satrec.j3oj2.value *
          sinio.value *
          (3.0 + 5.0 * cosio.value) /
          (1.0 + cosio.value);
    } else {
      satrec.xlcof.value = -0.25 *
          satrec.j3oj2.value *
          sinio.value *
          (3.0 + 5.0 * cosio.value) /
          temp4;
    }
    satrec.aycof.value = -0.5 * satrec.j3oj2.value * sinio.value;
    // sgp4fix use multiply for speed instead of pow
    delmotemp = 1.0 + satrec.eta.value * cos(satrec.mo.value);
    satrec.delmo.value = delmotemp * delmotemp * delmotemp;
    satrec.sinmao.value = sin(satrec.mo.value);
    satrec.x7thm1.value = 7.0 * cosio2.value - 1.0;

    /* --------------- deep space initialization ------------- */
    if ((2 * pi / satrec.no_unkozai.value) >= 225.0) {
      satrec.method.value = 'd';
      satrec.isimp.value = 1;
      tc = 0.0;
      inclm.value = satrec.inclo.value;

      dscom(
          epoch,
          satrec.ecco.value,
          satrec.argpo.value,
          tc,
          satrec.inclo.value,
          satrec.nodeo.value,
          satrec.no_unkozai.value,
          snodm,
          cnodm,
          sinim,
          cosim,
          sinomm,
          cosomm,
          day,
          satrec.e3,
          satrec.ee2,
          em,
          emsq,
          gam,
          satrec.peo,
          satrec.pgho,
          satrec.pho,
          satrec.pinco,
          satrec.plo,
          rtemsq,
          satrec.se2,
          satrec.se3,
          satrec.sgh2,
          satrec.sgh3,
          satrec.sgh4,
          satrec.sh2,
          satrec.sh3,
          satrec.si2,
          satrec.si3,
          satrec.sl2,
          satrec.sl3,
          satrec.sl4,
          s1,
          s2,
          s3,
          s4,
          s5,
          s6,
          s7,
          ss1,
          ss2,
          ss3,
          ss4,
          ss5,
          ss6,
          ss7,
          sz1,
          sz2,
          sz3,
          sz11,
          sz12,
          sz13,
          sz21,
          sz22,
          sz23,
          sz31,
          sz32,
          sz33,
          satrec.xgh2,
          satrec.xgh3,
          satrec.xgh4,
          satrec.xh2,
          satrec.xh3,
          satrec.xi2,
          satrec.xi3,
          satrec.xl2,
          satrec.xl3,
          satrec.xl4,
          nm,
          z1,
          z2,
          z3,
          z11,
          z12,
          z13,
          z21,
          z22,
          z23,
          z31,
          z32,
          z33,
          satrec.zmol,
          satrec.zmos);
      dpper(
          satrec.e3.value,
          satrec.ee2.value,
          satrec.peo.value,
          satrec.pgho.value,
          satrec.pho.value,
          satrec.pinco.value,
          satrec.plo.value,
          satrec.se2.value,
          satrec.se3.value,
          satrec.sgh2.value,
          satrec.sgh3.value,
          satrec.sgh4.value,
          satrec.sh2.value,
          satrec.sh3.value,
          satrec.si2.value,
          satrec.si3.value,
          satrec.sl2.value,
          satrec.sl3.value,
          satrec.sl4.value,
          satrec.t.value,
          satrec.xgh2.value,
          satrec.xgh3.value,
          satrec.xgh4.value,
          satrec.xh2.value,
          satrec.xh3.value,
          satrec.xi2.value,
          satrec.xi3.value,
          satrec.xl2.value,
          satrec.xl3.value,
          satrec.xl4.value,
          satrec.zmol.value,
          satrec.zmos.value,
          inclm.value,
          satrec.init.value,
          satrec.ecco,
          satrec.inclo,
          satrec.nodeo,
          satrec.argpo,
          satrec.mo,
          satrec.operationmode.value);

      argpm.value = 0.0;
      nodem.value = 0.0;
      mm.value = 0.0;

      dsinit(
          satrec.xke.value,
          cosim.value,
          emsq.value,
          satrec.argpo.value,
          s1.value,
          s2.value,
          s3.value,
          s4.value,
          s5.value,
          sinim.value,
          ss1.value,
          ss2.value,
          ss3.value,
          ss4.value,
          ss5.value,
          sz1.value,
          sz3.value,
          sz11.value,
          sz13.value,
          sz21.value,
          sz23.value,
          sz31.value,
          sz33.value,
          satrec.t.value,
          tc,
          satrec.gsto.value,
          satrec.mo.value,
          satrec.mdot.value,
          satrec.no_unkozai.value,
          satrec.nodeo.value,
          satrec.nodedot.value,
          xpidot,
          z1.value,
          z3.value,
          z11.value,
          z13.value,
          z21.value,
          z23.value,
          z31.value,
          z33.value,
          satrec.ecco.value,
          eccsq.value,
          em,
          argpm,
          inclm,
          mm,
          nm,
          nodem,
          satrec.irez,
          satrec.atime,
          satrec.d2201,
          satrec.d2211,
          satrec.d3210,
          satrec.d3222,
          satrec.d4410,
          satrec.d4422,
          satrec.d5220,
          satrec.d5232,
          satrec.d5421,
          satrec.d5433,
          satrec.dedt,
          satrec.didt,
          satrec.dmdt,
          dndt,
          satrec.dnodt,
          satrec.domdt,
          satrec.del1,
          satrec.del2,
          satrec.del3,
          satrec.xfact,
          satrec.xlamo,
          satrec.xli,
          satrec.xni);
    }

    /* ----------- set variables if not deep space ----------- */
    if (satrec.isimp.value != 1) {
      cc1sq = satrec.cc1.value * satrec.cc1.value;
      satrec.d2.value = 4.0 * ao.value * tsi * cc1sq;
      temp = satrec.d2.value * tsi * satrec.cc1.value / 3.0;
      satrec.d3.value = (17.0 * ao.value + sfour) * temp;
      satrec.d4.value = 0.5 *
          temp *
          ao.value *
          tsi *
          (221.0 * ao.value + 31.0 * sfour) *
          satrec.cc1.value;
      satrec.t3cof.value = satrec.d2.value + 2.0 * cc1sq;
      satrec.t4cof.value = 0.25 *
          (3.0 * satrec.d3.value +
              satrec.cc1.value * (12.0 * satrec.d2.value + 10.0 * cc1sq));
      satrec.t5cof.value = 0.2 *
          (3.0 * satrec.d4.value +
              12.0 * satrec.cc1.value * satrec.d3.value +
              6.0 * satrec.d2.value * satrec.d2.value +
              15.0 * cc1sq * (2.0 * satrec.d2.value + cc1sq));
    }
  } // if omeosq = 0 ...

  /* finally propogate to zero epoch to initialize all others. */
  // sgp4fix take out check to let satellites process until they are actually below earth surface
  //       if(satrec.error == 0)
  sgp4(satrec, 0.0, r, v);

  satrec.init.value = 'n';

  //#include "debug6.cpp"
  //sgp4fix return boolean. satrec.error contains any error codes
  return true;
} // sgp4init
