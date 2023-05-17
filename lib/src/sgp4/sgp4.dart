import 'dart:math';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/sgp4/dpper.dart';
import 'package:pious_squid/src/sgp4/dspace.dart';
import 'package:pious_squid/src/sgp4/elsetrec.dart';
import 'package:pious_squid/src/sgp4/pointer.dart';

///-----------------------------------------------------------------------------
///
///                             procedure sgp4
///
///  this procedure is the sgp4 prediction model from space command. this is an
///  updated and combined version of sgp4 and sdp4, which were originally
///  published separately in spacetrack report #3. this version follows the
///  methodology from the aiaa paper (2006) describing the history and
///  development of the code.
///
///  author        : david vallado                  719-573-2600   28 jun 2005
///
///  inputs        :
///    satrec	 - initialised structure from sgp4init() call.
///    tsince	 - time since epoch (minutes)
///
///  outputs       :
///    r           - position vector                     km
///    v           - velocity                            km/sec
///  return code - non-zero on error.
///                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
///                   2 - mean motion less than 0.0
///                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
///                   4 - semi-latus rectum < 0.0
///                   5 - epoch elements are sub-orbital
///                   6 - satellite has decayed
///
///  locals        :
///    am          -
///    axnl, aynl        -
///    betal       -
///    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
///    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
///    cosisq  , cossu   , sinsu   , cosu    , sinu
///    delm        -
///    delomg      -
///    dndt        -
///    eccm        -
///    emsq        -
///    ecose       -
///    el2         -
///    eo1         -
///    eccp        -
///    esine       -
///    argpm       -
///    argpp       -
///    omgadf      -c
///    pl          -
///    r           -
///    rtemsq      -
///    rdotl       -
///    rl          -
///    rvdot       -
///    rvdotl      -
///    su          -
///    t2  , t3   , t4    , tc
///    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
///    u   , ux   , uy    , uz     , vx     , vy     , vz
///    inclm       - inclination
///    mm          - mean anomaly
///    nm          - mean motion
///    nodem       - right asc of ascending node
///    xinc        -
///    xincp       -
///    xl          -
///    xlm         -
///    mp          -
///    xmdf        -
///    xmx         -
///    xmy         -
///    nodedf      -
///    xnode       -
///    nodep       -
///    np          -
///
///  coupling      :
///    getgravconst- no longer used. Variables are conatined within satrec
///    dpper
///    dpspace
///
///  references    :
///    hoots, roehrich, norad spacetrack report #3 1980
///    hoots, norad spacetrack report #6 1986
///    hoots, schumacher and glover 2004
///    vallado, crawford, hujsak, kelso  2006
/// ----------------------------------------------------------------------------
bool sgp4(final ElsetRec satrec, final double tsince, final List<double> r,
    final List<double> v) {
  double am,
      axnl,
      aynl,
      betal,
      cosim,
      cnod,
      cos2u,
      coseo1 = 0.0,
      cosi,
      cosip,
      cosisq,
      cossu,
      cosu,
      delm,
      delomg,
      emsq,
      ecose,
      el2,
      eo1,
      esine,
      argpdf,
      pl,
      mrt = 0.0,
      mvt,
      rdotl,
      rl,
      rvdot,
      rvdotl,
      sinim,
      sin2u,
      sineo1 = 0.0,
      sini,
      sinip,
      sinsu,
      sinu,
      snod,
      su,
      t2,
      t3,
      t4,
      tem5,
      temp,
      temp1,
      temp2,
      tempa,
      tempe,
      templ,
      u,
      ux,
      uy,
      uz,
      vx,
      vy,
      vz,
      xinc,
      xl,
      xlm,
      xmdf,
      xmx,
      xmy,
      nodedf,
      xnode,
      tc,
      x2o3,
      vkmpersec,
      delmtemp;
  int ktr;
  final em = Pointer(0.0),
      argpm = Pointer(0.0),
      inclm = Pointer(0.0),
      mm = Pointer(0.0),
      nodem = Pointer(0.0),
      dndt = Pointer(0.0),
      nm = Pointer(0.0),
      ep = Pointer(0.0),
      xincp = Pointer(0.0),
      nodep = Pointer(0.0),
      argpp = Pointer(0.0),
      mp = Pointer(0.0);

  /* ------------------ set mathematical constants --------------- */
  // sgp4fix divisor for divide by zero check on inclination
  // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
  // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
  const temp4 = 1.5e-12;
  x2o3 = 2.0 / 3.0;
  // sgp4fix identify constants and allow alternate values
  // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
  vkmpersec = satrec.radiusearthkm.value * satrec.xke.value / 60.0;

  /* --------------------- clear sgp4 error flag ----------------- */
  satrec.t.value = tsince;
  satrec.error.value = 0;

  /* ------- update for secular gravity and atmospheric drag ----- */
  xmdf = satrec.mo.value + satrec.mdot.value * satrec.t.value;
  argpdf = satrec.argpo.value + satrec.argpdot.value * satrec.t.value;
  nodedf = satrec.nodeo.value + satrec.nodedot.value * satrec.t.value;
  argpm.value = argpdf;
  mm.value = xmdf;
  t2 = satrec.t.value * satrec.t.value;
  nodem.value = nodedf + satrec.nodecf.value * t2;
  tempa = 1.0 - satrec.cc1.value * satrec.t.value;
  tempe = satrec.bstar.value * satrec.cc4.value * satrec.t.value;
  templ = satrec.t2cof.value * t2;

  if (satrec.isimp.value != 1) {
    delomg = satrec.omgcof.value * satrec.t.value;
    // sgp4fix use mutliply for speed instead of pow
    delmtemp = 1.0 + satrec.eta.value * cos(xmdf);
    delm = satrec.xmcof.value *
        (delmtemp * delmtemp * delmtemp - satrec.delmo.value);
    temp = delomg + delm;
    mm.value = xmdf + temp;
    argpm.value = argpdf - temp;
    t3 = t2 * satrec.t.value;
    t4 = t3 * satrec.t.value;
    tempa = tempa -
        satrec.d2.value * t2 -
        satrec.d3.value * t3 -
        satrec.d4.value * t4;
    tempe = tempe +
        satrec.bstar.value *
            satrec.cc5.value *
            (sin(mm.value) - satrec.sinmao.value);
    templ = templ +
        satrec.t3cof.value * t3 +
        t4 * (satrec.t4cof.value + satrec.t.value * satrec.t5cof.value);
  }

  nm.value = satrec.no_unkozai.value;
  em.value = satrec.ecco.value;
  inclm.value = satrec.inclo.value;
  if (satrec.method.value == 'd') {
    tc = satrec.t.value;
    dspace(
        satrec.irez.value,
        satrec.d2201.value,
        satrec.d2211.value,
        satrec.d3210.value,
        satrec.d3222.value,
        satrec.d4410.value,
        satrec.d4422.value,
        satrec.d5220.value,
        satrec.d5232.value,
        satrec.d5421.value,
        satrec.d5433.value,
        satrec.dedt.value,
        satrec.del1.value,
        satrec.del2.value,
        satrec.del3.value,
        satrec.didt.value,
        satrec.dmdt.value,
        satrec.dnodt.value,
        satrec.domdt.value,
        satrec.argpo.value,
        satrec.argpdot.value,
        satrec.t.value,
        tc,
        satrec.gsto.value,
        satrec.xfact.value,
        satrec.xlamo.value,
        satrec.no_unkozai.value,
        satrec.atime,
        em,
        argpm,
        inclm,
        satrec.xli,
        mm,
        satrec.xni,
        nodem,
        dndt,
        nm);
  } // if method = d

  if (nm.value <= 0.0) {
    //         printf("# error nm %f\n", nm);
    satrec.error.value = 2;
    // sgp4fix add return
    return false;
  }
  am = pow(satrec.xke.value / nm.value, x2o3) * tempa * tempa;
  nm.value = satrec.xke.value / pow(am, 1.5);
  em.value = em.value - tempe;

  // fix tolerance for error recognition
  // sgp4fix am is fixed from the previous nm check
  if ((em.value >= 1.0) || (em.value < -0.001) /* || (am < 0.95)*/) {
    //         printf("# error em %f\n", em);
    satrec.error.value = 1;
    // sgp4fix to return if there is an error in eccentricity
    return false;
  }
  // sgp4fix fix tolerance to avoid a divide by zero
  if (em.value < 1.0e-6) {
    em.value = 1.0e-6;
  }
  mm.value = mm.value + satrec.no_unkozai.value * templ;
  xlm = mm.value + argpm.value + nodem.value;
  emsq = em.value * em.value;
  temp = 1.0 - emsq;

  nodem.value = nodem.value % twoPi;
  argpm.value = argpm.value % twoPi;
  xlm = xlm % twoPi;
  mm.value = (xlm - argpm.value - nodem.value) % twoPi;

  // sgp4fix recover singly averaged mean elements
  satrec.am.value = am;
  satrec.em.value = em.value;
  satrec.im.value = inclm.value;
  satrec.Om.value = nodem.value;
  satrec.om.value = argpm.value;
  satrec.mm.value = mm.value;
  satrec.nm.value = nm.value;

  /* ----------------- compute extra mean quantities ------------- */
  sinim = sin(inclm.value);
  cosim = cos(inclm.value);

  /* -------------------- add lunar-solar periodics -------------- */
  ep.value = em.value;
  xincp.value = inclm.value;
  argpp.value = argpm.value;
  nodep.value = nodem.value;
  mp.value = mm.value;
  sinip = sinim;
  cosip = cosim;
  if (satrec.method.value == 'd') {
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
        satrec.inclo.value,
        'n',
        ep,
        xincp,
        nodep,
        argpp,
        mp,
        satrec.operationmode.value);
    if (xincp.value < 0.0) {
      xincp.value = -xincp.value;
      nodep.value = nodep.value + pi;
      argpp.value = argpp.value - pi;
    }
    if ((ep.value < 0.0) || (ep.value > 1.0)) {
      //            printf("# error ep %f\n", ep);
      satrec.error.value = 3;
      // sgp4fix add return
      return false;
    }
  } // if method = d

  /* -------------------- long period periodics ------------------ */
  if (satrec.method.value == 'd') {
    sinip = sin(xincp.value);
    cosip = cos(xincp.value);
    satrec.aycof.value = -0.5 * satrec.j3oj2.value * sinip;
    // sgp4fix for divide by zero for xincp = 180 deg
    if ((cosip + 1.0).abs() > 1.5e-12) {
      satrec.xlcof.value = -0.25 *
          satrec.j3oj2.value *
          sinip *
          (3.0 + 5.0 * cosip) /
          (1.0 + cosip);
    } else {
      satrec.xlcof.value =
          -0.25 * satrec.j3oj2.value * sinip * (3.0 + 5.0 * cosip) / temp4;
    }
  }
  axnl = ep.value * cos(argpp.value);
  temp = 1.0 / (am * (1.0 - ep.value * ep.value));
  aynl = ep.value * sin(argpp.value) + temp * satrec.aycof.value;
  xl = mp.value + argpp.value + nodep.value + temp * satrec.xlcof.value * axnl;

  /* --------------------- solve kepler's equation --------------- */
  u = (xl - nodep.value) % twoPi;
  eo1 = u;
  tem5 = 9999.9;
  ktr = 1;
  //   sgp4fix for kepler iteration
  //   the following iteration needs better limits on corrections
  while ((tem5.abs() >= 1.0e-12) && (ktr <= 10)) {
    sineo1 = sin(eo1);
    coseo1 = cos(eo1);
    tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
    tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
    if (tem5.abs() >= 0.95) {
      tem5 = tem5 > 0.0 ? 0.95 : -0.95;
    }
    eo1 = eo1 + tem5;
    ktr = ktr + 1;
  }

  /* ------------- short period preliminary quantities ----------- */
  ecose = axnl * coseo1 + aynl * sineo1;
  esine = axnl * sineo1 - aynl * coseo1;
  el2 = axnl * axnl + aynl * aynl;
  pl = am * (1.0 - el2);
  if (pl < 0.0) {
    //         printf("# error pl %f\n", pl);
    satrec.error.value = 4;
    // sgp4fix add return
    return false;
  } else {
    rl = am * (1.0 - ecose);
    rdotl = sqrt(am) * esine / rl;
    rvdotl = sqrt(pl) / rl;
    betal = sqrt(1.0 - el2);
    temp = esine / (1.0 + betal);
    sinu = am / rl * (sineo1 - aynl - axnl * temp);
    cosu = am / rl * (coseo1 - axnl + aynl * temp);
    su = atan2(sinu, cosu);
    sin2u = (cosu + cosu) * sinu;
    cos2u = 1.0 - 2.0 * sinu * sinu;
    temp = 1.0 / pl;
    temp1 = 0.5 * satrec.j2.value * temp;
    temp2 = temp1 * temp;

    /* -------------- update for short period periodics ------------ */
    if (satrec.method.value == 'd') {
      cosisq = cosip * cosip;
      satrec.con41.value = 3.0 * cosisq - 1.0;
      satrec.x1mth2.value = 1.0 - cosisq;
      satrec.x7thm1.value = 7.0 * cosisq - 1.0;
    }
    mrt = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41.value) +
        0.5 * temp1 * satrec.x1mth2.value * cos2u;
    su = su - 0.25 * temp2 * satrec.x7thm1.value * sin2u;
    xnode = nodep.value + 1.5 * temp2 * cosip * sin2u;
    xinc = xincp.value + 1.5 * temp2 * cosip * sinip * cos2u;
    mvt = rdotl -
        nm.value * temp1 * satrec.x1mth2.value * sin2u / satrec.xke.value;
    rvdot = rvdotl +
        nm.value *
            temp1 *
            (satrec.x1mth2.value * cos2u + 1.5 * satrec.con41.value) /
            satrec.xke.value;

    /* --------------------- orientation vectors ------------------- */
    sinsu = sin(su);
    cossu = cos(su);
    snod = sin(xnode);
    cnod = cos(xnode);
    sini = sin(xinc);
    cosi = cos(xinc);
    xmx = -snod * cosi;
    xmy = cnod * cosi;
    ux = xmx * sinsu + cnod * cossu;
    uy = xmy * sinsu + snod * cossu;
    uz = sini * sinsu;
    vx = xmx * cossu - cnod * sinsu;
    vy = xmy * cossu - snod * sinsu;
    vz = sini * cossu;

    /* --------- position and velocity (in km and km/sec) ---------- */
    r[0] = (mrt * ux) * satrec.radiusearthkm.value;
    r[1] = (mrt * uy) * satrec.radiusearthkm.value;
    r[2] = (mrt * uz) * satrec.radiusearthkm.value;
    v[0] = (mvt * ux + rvdot * vx) * vkmpersec;
    v[1] = (mvt * uy + rvdot * vy) * vkmpersec;
    v[2] = (mvt * uz + rvdot * vz) * vkmpersec;
  } // if pl > 0

  // sgp4fix for decaying satellites
  if (mrt < 1.0) {
    //         printf("# decay condition %11.6f \n",mrt);
    satrec.error.value = 6;
    return false;
  }

  //#include "debug7.cpp"
  return true;
} // sgp4
