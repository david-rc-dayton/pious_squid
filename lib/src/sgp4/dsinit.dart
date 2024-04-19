// ignore_for_file: parameter_assignments

import 'dart:math';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/pointer.dart';

/// ----------------------------------------------------------------------------
///
///                           procedure dsinit
///
///  this procedure provides deep space contributions to mean motion dot due
///    to geopotential resonance with half day and one day orbits.
///
///  author        : david vallado                  719-573-2600   28 jun 2005
///
///  inputs        :
///    xke         - reciprocal of tumin
///    cosim, sinim-
///    emsq        - eccentricity squared
///    argpo       - argument of perigee
///    s1, s2, s3, s4, s5      -
///    ss1, ss2, ss3, ss4, ss5 -
///    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
///    t           - time
///    tc          -
///    gsto        - greenwich sidereal time                   rad
///    mo          - mean anomaly
///    mdot        - mean anomaly dot (rate)
///    no          - mean motion
///    nodeo       - right ascension of ascending node
///    nodedot     - right ascension of ascending node dot (rate)
///    xpidot      -
///    z1, z3, z11, z13, z21, z23, z31, z33 -
///    eccm        - eccentricity
///    argpm       - argument of perigee
///    inclm       - inclination
///    mm          - mean anomaly
///    xn          - mean motion
///    nodem       - right ascension of ascending node
///
///  outputs       :
///    em          - eccentricity
///    argpm       - argument of perigee
///    inclm       - inclination
///    mm          - mean anomaly
///    nm          - mean motion
///    nodem       - right ascension of ascending node
///    irez        - flag for resonance           0-none, 1-one day, 2-half day
///    atime       -
///    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
///    dedt        -
///    didt        -
///    dmdt        -
///    dndt        -
///    dnodt       -
///    domdt       -
///    del1, del2, del3        -
///    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
///    theta       -
///    xfact       -
///    xlamo       -
///    xli         -
///    xni
///
///  locals        :
///    ainv2       -
///    aonv        -
///    cosisq      -
///    eoc         -
///    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
///    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
///    sini2       -
///    temp        -
///    temp1       -
///    theta       -
///    xno2        -
///
///  coupling      :
///    getgravconst- no longer used
///
///  references    :
///    hoots, roehrich, norad spacetrack report #3 1980
///    hoots, norad spacetrack report #6 1986
///    hoots, schumacher and glover 2004
///    vallado, crawford, hujsak, kelso  2006
///	----------------------------------------------------------------------------
void dsinit(
    // sgp4fix just send in xke as a constant and eliminate getgravconst call
    // gravconsttype whichconst,
    final double xke,
    final double cosim,
    double emsq,
    final double argpo,
    final double s1,
    final double s2,
    final double s3,
    final double s4,
    final double s5,
    final double sinim,
    final double ss1,
    final double ss2,
    final double ss3,
    final double ss4,
    final double ss5,
    final double sz1,
    final double sz3,
    final double sz11,
    final double sz13,
    final double sz21,
    final double sz23,
    final double sz31,
    final double sz33,
    final double t,
    final double tc,
    final double gsto,
    final double mo,
    final double mdot,
    final double no,
    final double nodeo,
    final double nodedot,
    final double xpidot,
    final double z1,
    final double z3,
    final double z11,
    final double z13,
    final double z21,
    final double z23,
    final double z31,
    final double z33,
    final double ecco,
    final double eccsq,
    final Pointer<double> em,
    final Pointer<double> argpm,
    final Pointer<double> inclm,
    final Pointer<double> mm,
    final Pointer<double> nm,
    final Pointer<double> nodem,
    final Pointer<int> irez,
    final Pointer<double> atime,
    final Pointer<double> d2201,
    final Pointer<double> d2211,
    final Pointer<double> d3210,
    final Pointer<double> d3222,
    final Pointer<double> d4410,
    final Pointer<double> d4422,
    final Pointer<double> d5220,
    final Pointer<double> d5232,
    final Pointer<double> d5421,
    final Pointer<double> d5433,
    final Pointer<double> dedt,
    final Pointer<double> didt,
    final Pointer<double> dmdt,
    final Pointer<double> dndt,
    final Pointer<double> dnodt,
    final Pointer<double> domdt,
    final Pointer<double> del1,
    final Pointer<double> del2,
    final Pointer<double> del3,
    final Pointer<double> xfact,
    final Pointer<double> xlamo,
    final Pointer<double> xli,
    final Pointer<double> xni) {
  /* --------------------- local variables ------------------------ */
  double ainv2,
      aonv,
      cosisq,
      eoc,
      f220,
      f221,
      f311,
      f321,
      f322,
      f330,
      f441,
      f442,
      f522,
      f523,
      f542,
      f543,
      g200,
      g201,
      g211,
      g300,
      g310,
      g322,
      g410,
      g422,
      g520,
      g521,
      g532,
      g533,
      ses,
      sgs,
      sghl,
      sghs,
      shs,
      shll,
      sis,
      sini2,
      sls,
      temp,
      temp1,
      theta,
      xno2,
      emo,
      emsqo;

  final q22 = 1.7891679e-6;
  final q31 = 2.1460748e-6;
  final q33 = 2.2123015e-7;
  final root22 = 1.7891679e-6;
  final root44 = 7.3636953e-9;
  final root54 = 2.1765803e-9;
  final rptim =
      4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
  final root32 = 3.7393792e-7;
  final root52 = 1.1428639e-7;
  final x2o3 = 2.0 / 3.0;
  final znl = 1.5835218e-4;
  final zns = 1.19459e-5;

  // sgp4fix identify constants and allow alternate values
  // just xke is used here so pass it in rather than have multiple calls
  // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

  /* -------------------- deep space initialization ------------ */
  irez.value = 0;
  if ((nm.value < 0.0052359877) && (nm.value > 0.0034906585)) {
    irez.value = 1;
  }
  if ((nm.value >= 8.26e-3) && (nm.value <= 9.24e-3) && (em.value >= 0.5)) {
    irez.value = 2;
  }

  /* ------------------------ do solar terms ------------------- */
  ses = ss1 * zns * ss5;
  sis = ss2 * zns * (sz11 + sz13);
  sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
  sghs = ss4 * zns * (sz31 + sz33 - 6.0);
  shs = -zns * ss2 * (sz21 + sz23);
  // sgp4fix for 180 deg incl
  if ((inclm.value < 5.2359877e-2) || (inclm.value > pi - 5.2359877e-2)) {
    shs = 0.0;
  }
  if (sinim != 0.0) {
    shs = shs / sinim;
  }
  sgs = sghs - cosim * shs;

  /* ------------------------- do lunar terms ------------------ */
  dedt.value = ses + s1 * znl * s5;
  didt.value = sis + s2 * znl * (z11 + z13);
  dmdt.value = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
  sghl = s4 * znl * (z31 + z33 - 6.0);
  shll = -znl * s2 * (z21 + z23);
  // sgp4fix for 180 deg incl
  if ((inclm.value < 5.2359877e-2) || (inclm.value > pi - 5.2359877e-2)) {
    shll = 0.0;
  }
  domdt.value = sgs + sghl;
  dnodt.value = shs;
  if (sinim != 0.0) {
    domdt.value = domdt.value - cosim / sinim * shll;
    dnodt.value = dnodt.value + shll / sinim;
  }

  /* ----------- calculate deep space resonance effects -------- */
  dndt.value = 0.0;
  theta = (gsto + tc * rptim) % twoPi;
  em.value = em.value + dedt.value * t;
  inclm.value = inclm.value + didt.value * t;
  argpm.value = argpm.value + domdt.value * t;
  nodem.value = nodem.value + dnodt.value * t;
  mm.value = mm.value + dmdt.value * t;
  //   sgp4fix for negative inclinations
  //   the following if statement should be commented out
  //if (inclm < 0.0)
  //  {
  //    inclm  = -inclm;
  //    argpm  = argpm - pi;
  //    nodem = nodem + pi;
  //  }

  /* -------------- initialize the resonance terms ------------- */
  if (irez.value != 0) {
    aonv = pow(nm.value / xke, x2o3) as double;

    /* ---------- geopotential resonance for 12 hour orbits ------ */
    if (irez.value == 2) {
      cosisq = cosim * cosim;
      emo = em.value;
      em.value = ecco;
      emsqo = emsq;
      emsq = eccsq;
      eoc = em.value * emsq;
      g201 = -0.306 - (em.value - 0.64) * 0.440;

      if (em.value <= 0.65) {
        g211 = 3.616 - 13.2470 * em.value + 16.2900 * emsq;
        g310 = -19.302 + 117.3900 * em.value - 228.4190 * emsq + 156.5910 * eoc;
        g322 =
            -18.9068 + 109.7927 * em.value - 214.6334 * emsq + 146.5816 * eoc;
        g410 = -41.122 + 242.6940 * em.value - 471.0940 * emsq + 313.9530 * eoc;
        g422 =
            -146.407 + 841.8800 * em.value - 1629.014 * emsq + 1083.4350 * eoc;
        g520 =
            -532.114 + 3017.977 * em.value - 5740.032 * emsq + 3708.2760 * eoc;
      } else {
        g211 = -72.099 + 331.819 * em.value - 508.738 * emsq + 266.724 * eoc;
        g310 =
            -346.844 + 1582.851 * em.value - 2415.925 * emsq + 1246.113 * eoc;
        g322 =
            -342.585 + 1554.908 * em.value - 2366.899 * emsq + 1215.972 * eoc;
        g410 =
            -1052.797 + 4758.686 * em.value - 7193.992 * emsq + 3651.957 * eoc;
        g422 = -3581.690 +
            16178.110 * em.value -
            24462.770 * emsq +
            12422.520 * eoc;
        if (em.value > 0.715) {
          g520 =
              -5149.66 + 29936.92 * em.value - 54087.36 * emsq + 31324.56 * eoc;
        } else {
          g520 = 1464.74 - 4664.75 * em.value + 3763.64 * emsq;
        }
      }
      if (em.value < 0.7) {
        g533 = -919.22770 +
            4988.6100 * em.value -
            9064.7700 * emsq +
            5542.21 * eoc;
        g521 = -822.71072 +
            4568.6173 * em.value -
            8491.4146 * emsq +
            5337.524 * eoc;
        g532 =
            -853.66600 + 4690.2500 * em.value - 8624.7700 * emsq + 5341.4 * eoc;
      } else {
        g533 = -37995.780 +
            161616.52 * em.value -
            229838.20 * emsq +
            109377.94 * eoc;
        g521 = -51752.104 +
            218913.95 * em.value -
            309468.16 * emsq +
            146349.42 * eoc;
        g532 = -40023.880 +
            170470.89 * em.value -
            242699.48 * emsq +
            115605.82 * eoc;
      }

      sini2 = sinim * sinim;
      f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
      f221 = 1.5 * sini2;
      f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
      f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
      f441 = 35.0 * sini2 * f220;
      f442 = 39.3750 * sini2 * sini2;
      f522 = 9.84375 *
          sinim *
          (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) +
              0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
      f523 = sinim *
          (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq) +
              6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
      f542 = 29.53125 *
          sinim *
          (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
      f543 = 29.53125 *
          sinim *
          (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));
      xno2 = nm.value * nm.value;
      ainv2 = aonv * aonv;
      temp1 = 3.0 * xno2 * ainv2;
      temp = temp1 * root22;
      d2201.value = temp * f220 * g201;
      d2211.value = temp * f221 * g211;
      temp1 = temp1 * aonv;
      temp = temp1 * root32;
      d3210.value = temp * f321 * g310;
      d3222.value = temp * f322 * g322;
      temp1 = temp1 * aonv;
      temp = 2.0 * temp1 * root44;
      d4410.value = temp * f441 * g410;
      d4422.value = temp * f442 * g422;
      temp1 = temp1 * aonv;
      temp = temp1 * root52;
      d5220.value = temp * f522 * g520;
      d5232.value = temp * f523 * g532;
      temp = 2.0 * temp1 * root54;
      d5421.value = temp * f542 * g521;
      d5433.value = temp * f543 * g533;
      xlamo.value = (mo + nodeo + nodeo - theta - theta) % twoPi;
      xfact.value =
          mdot + dmdt.value + 2.0 * (nodedot + dnodt.value - rptim) - no;
      em.value = emo;
      emsq = emsqo;
    }

    /* ---------------- synchronous resonance terms -------------- */
    if (irez.value == 1) {
      g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
      g310 = 1.0 + 2.0 * emsq;
      g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
      f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
      f311 =
          0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
      f330 = 1.0 + cosim;
      f330 = 1.875 * f330 * f330 * f330;
      del1.value = 3.0 * nm.value * nm.value * aonv * aonv;
      del2.value = 2.0 * del1.value * f220 * g200 * q22;
      del3.value = 3.0 * del1.value * f330 * g300 * q33 * aonv;
      del1.value = del1.value * f311 * g310 * q31 * aonv;
      xlamo.value = (mo + nodeo + argpo - theta) % twoPi;
      xfact.value =
          mdot + xpidot - rptim + dmdt.value + domdt.value + dnodt.value - no;
    }

    /* ------------ for sgp4, initialize the integrator ---------- */
    xli.value = xlamo.value;
    xni.value = no;
    atime.value = 0.0;
    nm.value = no + dndt.value;
  }

  //#include "debug3.cpp"
} // dsinit
