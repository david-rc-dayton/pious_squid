import 'dart:math';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/pointer.dart';

/// ----------------------------------------------------------------------------
///
///                           procedure dspace
///
///  this procedure provides deep space contributions to mean elements for
///    perturbing third body.  these effects have been averaged over one
///    revolution of the sun and moon.  for earth resonance effects, the
///    effects have been averaged over no revolutions of the satellite.
///    (mean motion)
///
///  author        : david vallado                  719-573-2600   28 jun 2005
///
///  inputs        :
///    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
///    dedt        -
///    del1, del2, del3  -
///    didt        -
///    dmdt        -
///    dnodt       -
///    domdt       -
///    irez        - flag for resonance           0-none, 1-one day, 2-half day
///    argpo       - argument of perigee
///    argpdot     - argument of perigee dot (rate)
///    t           - time
///    tc          -
///    gsto        - gst
///    xfact       -
///    xlamo       -
///    no          - mean motion
///    atime       -
///    em          - eccentricity
///    ft          -
///    argpm       - argument of perigee
///    inclm       - inclination
///    xli         -
///    mm          - mean anomaly
///    xni         - mean motion
///    nodem       - right ascension of ascending node
///
///  outputs       :
///    atime       -
///    em          - eccentricity
///    argpm       - argument of perigee
///    inclm       - inclination
///    xli         -
///    mm          - mean anomaly
///    xni         -
///    nodem       - right ascension of ascending node
///    dndt        -
///    nm          - mean motion
///
///  locals        :
///    delt        -
///    ft          -
///    theta       -
///    x2li        -
///    x2omi       -
///    xl          -
///    xldot       -
///    xnddt       -
///    xndt        -
///    xomi        -
///
///  coupling      :
///    none        -
///
///  references    :
///    hoots, roehrich, norad spacetrack report #3 1980
///    hoots, norad spacetrack report #6 1986
///    hoots, schumacher and glover 2004
///    vallado, crawford, hujsak, kelso  2006
/// ----------------------------------------------------------------------------
void dspace(
    final int irez,
    final double d2201,
    final double d2211,
    final double d3210,
    final double d3222,
    final double d4410,
    final double d4422,
    final double d5220,
    final double d5232,
    final double d5421,
    final double d5433,
    final double dedt,
    final double del1,
    final double del2,
    final double del3,
    final double didt,
    final double dmdt,
    final double dnodt,
    final double domdt,
    final double argpo,
    final double argpdot,
    final double t,
    final double tc,
    final double gsto,
    final double xfact,
    final double xlamo,
    final double no,
    final Pointer<double> atime,
    final Pointer<double> em,
    final Pointer<double> argpm,
    final Pointer<double> inclm,
    final Pointer<double> xli,
    final Pointer<double> mm,
    final Pointer<double> xni,
    final Pointer<double> nodem,
    final Pointer<double> dndt,
    final Pointer<double> nm) {
  int iretn;
  double delt,
      ft,
      theta,
      x2li,
      x2omi,
      xl,
      xldot = 0.0,
      xnddt = 0.0,
      xndt = 0.0,
      xomi;

  final fasx2 = 0.13130908;
  final fasx4 = 2.8843198;
  final fasx6 = 0.37448087;
  final g22 = 5.7686396;
  final g32 = 0.95240898;
  final g44 = 1.8014998;
  final g52 = 1.0508330;
  final g54 = 4.4108898;
  final rptim =
      4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
  final stepp = 720.0;
  final stepn = -720.0;
  final step2 = 259200.0;

  /* ----------- calculate deep space resonance effects ----------- */
  dndt.value = 0.0;
  theta = (gsto + tc * rptim) % twoPi;
  em.value = em.value + dedt * t;

  inclm.value = inclm.value + didt * t;
  argpm.value = argpm.value + domdt * t;
  nodem.value = nodem.value + dnodt * t;
  mm.value = mm.value + dmdt * t;

  //   sgp4fix for negative inclinations
  //   the following if statement should be commented out
  //  if (inclm < 0.0)
  // {
  //    inclm = -inclm;
  //    argpm = argpm - pi;
  //    nodem = nodem + pi;
  //  }

  /* - update resonances : numerical (euler-maclaurin) integration - */
  /* ------------------------- epoch restart ----------------------  */
  //   sgp4fix for propagator problems
  //   the following integration works for negative time steps and periods
  //   the specific changes are unknown because the original code was so convoluted

  // sgp4fix take out atime = 0.0 and fix for faster operation
  ft = 0.0;
  if (irez != 0) {
    // sgp4fix streamline check
    if ((atime.value == 0.0) ||
        (t * atime.value <= 0.0) ||
        (t.abs() < atime.value.abs())) {
      atime.value = 0.0;
      xni.value = no;
      xli.value = xlamo;
    }
    // sgp4fix move check outside loop
    if (t > 0.0) {
      delt = stepp;
    } else {
      delt = stepn;
    }

    iretn = 381; // added for do loop
    while (iretn == 381) {
      /* ------------------- dot terms calculated ------------- */
      /* ----------- near - synchronous resonance terms ------- */
      if (irez != 2) {
        xndt = del1 * sin(xli.value - fasx2) +
            del2 * sin(2.0 * (xli.value - fasx4)) +
            del3 * sin(3.0 * (xli.value - fasx6));
        xldot = xni.value + xfact;
        xnddt = del1 * cos(xli.value - fasx2) +
            2.0 * del2 * cos(2.0 * (xli.value - fasx4)) +
            3.0 * del3 * cos(3.0 * (xli.value - fasx6));
        xnddt = xnddt * xldot;
      } else {
        /* --------- near - half-day resonance terms -------- */
        xomi = argpo + argpdot * atime.value;
        x2omi = xomi + xomi;
        x2li = xli.value + xli.value;
        xndt = d2201 * sin(x2omi + xli.value - g22) +
            d2211 * sin(xli.value - g22) +
            d3210 * sin(xomi + xli.value - g32) +
            d3222 * sin(-xomi + xli.value - g32) +
            d4410 * sin(x2omi + x2li - g44) +
            d4422 * sin(x2li - g44) +
            d5220 * sin(xomi + xli.value - g52) +
            d5232 * sin(-xomi + xli.value - g52) +
            d5421 * sin(xomi + x2li - g54) +
            d5433 * sin(-xomi + x2li - g54);
        xldot = xni.value + xfact;
        xnddt = d2201 * cos(x2omi + xli.value - g22) +
            d2211 * cos(xli.value - g22) +
            d3210 * cos(xomi + xli.value - g32) +
            d3222 * cos(-xomi + xli.value - g32) +
            d5220 * cos(xomi + xli.value - g52) +
            d5232 * cos(-xomi + xli.value - g52) +
            2.0 *
                (d4410 * cos(x2omi + x2li - g44) +
                    d4422 * cos(x2li - g44) +
                    d5421 * cos(xomi + x2li - g54) +
                    d5433 * cos(-xomi + x2li - g54));
        xnddt = xnddt * xldot;
      }

      /* ----------------------- integrator ------------------- */
      // sgp4fix move end checks to end of routine
      if ((t - atime.value).abs() >= stepp) {
        iretn = 381;
      } else // exit here
      {
        ft = t - atime.value;
        iretn = 0;
      }

      if (iretn == 381) {
        xli.value = xli.value + xldot * delt + xndt * step2;
        xni.value = xni.value + xndt * delt + xnddt * step2;
        atime.value = atime.value + delt;
      }
    } // while iretn = 381

    nm.value = xni.value + xndt * ft + xnddt * ft * ft * 0.5;
    xl = xli.value + xldot * ft + xndt * ft * ft * 0.5;
    if (irez != 1) {
      mm.value = xl - 2.0 * nodem.value + 2.0 * theta;
      dndt.value = nm.value - no;
    } else {
      mm.value = xl - nodem.value - argpm.value + theta;
      dndt.value = nm.value - no;
    }
    nm.value = no + dndt.value;
  }

  //#include "debug4.cpp"
} // dsspace
