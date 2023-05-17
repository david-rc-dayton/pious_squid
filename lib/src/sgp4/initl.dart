import 'dart:math';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/sgp4/pointer.dart';
import 'package:pious_squid/src/sgp4/utilities.dart';

/// ----------------------------------------------------------------------------
///
///                           procedure initl
///
///  this procedure initializes the spg4 propagator. all the initialization is
///    consolidated here instead of having multiple loops inside other routines.
///
///  author        : david vallado                  719-573-2600   28 jun 2005
///
///  inputs        :
///    satn        - satellite number - not needed, placed in satrec
///    xke         - reciprocal of tumin
///    j2          - j2 zonal harmonic
///    ecco        - eccentricity                           0.0 - 1.0
///    epoch       - epoch time in days from jan 0, 1950. 0 hr
///    inclo       - inclination of satellite
///    no          - mean motion of satellite
///
///  outputs       :
///    ainv        - 1.0 / a
///    ao          - semi major axis
///    con41       -
///    con42       - 1.0 - 5.0 cos(i)
///    cosio       - cosine of inclination
///    cosio2      - cosio squared
///    eccsq       - eccentricity squared
///    method      - flag for deep space                    'd', 'n'
///    omeosq      - 1.0 - ecco * ecco
///    posq        - semi-parameter squared
///    rp          - radius of perigee
///    rteosq      - square root of (1.0 - ecco*ecco)
///    sinio       - sine of inclination
///    gsto        - gst at time of observation               rad
///    no          - mean motion of satellite
///
///  locals        :
///    ak          -
///    d1          -
///    del         -
///    adel        -
///    po          -
///
///  coupling      :
///    getgravconst- no longer used
///    gstime      - find greenwich sidereal time from the julian date
///
///  references    :
///    hoots, roehrich, norad spacetrack report #3 1980
///    hoots, norad spacetrack report #6 1986
///    hoots, schumacher and glover 2004
///    vallado, crawford, hujsak, kelso  2006
///	----------------------------------------------------------------------------
void initl(
    // sgp4fix satn not needed. include in satrec in case needed later
    // int satn,
    // sgp4fix just pass in xke and j2
    // gravconsttype whichconst,
    final double xke,
    final double j2,
    final double ecco,
    final double epoch,
    final double inclo,
    final double noKozai,
    final String opsmode,
    final Pointer<String> method,
    final Pointer<double> ainv,
    final Pointer<double> ao,
    final Pointer<double> con41,
    final Pointer<double> con42,
    final Pointer<double> cosio,
    final Pointer<double> cosio2,
    final Pointer<double> eccsq,
    final Pointer<double> omeosq,
    final Pointer<double> posq,
    final Pointer<double> rp,
    final Pointer<double> rteosq,
    final Pointer<double> sinio,
    final Pointer<double> gsto,
    final Pointer<double> noUnkozai) {
  /* --------------------- local variables ------------------------ */
  double ak, d1, del, adel, po, x2o3;

  // sgp4fix use old way of finding gst
  double ds70;
  double ts70, tfrac, c1, thgr70, fk5r, c1p2p;

  /* ----------------------- earth constants ---------------------- */
  // sgp4fix identify constants and allow alternate values
  // only xke and j2 are used here so pass them in directly
  // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
  x2o3 = 2.0 / 3.0;

  /* ------------- calculate auxillary epoch quantities ---------- */
  eccsq.value = ecco * ecco;
  omeosq.value = 1.0 - eccsq.value;
  rteosq.value = sqrt(omeosq.value);
  cosio.value = cos(inclo);
  cosio2.value = cosio.value * cosio.value;

  /* ------------------ un-kozai the mean motion ----------------- */
  ak = pow(xke / noKozai, x2o3) as double;
  d1 = 0.75 * j2 * (3.0 * cosio2.value - 1.0) / (rteosq.value * omeosq.value);
  del = d1 / (ak * ak);
  adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
  del = d1 / (adel * adel);
  noUnkozai.value = noKozai / (1.0 + del);

  ao.value = pow(xke / (noUnkozai.value), x2o3) as double;
  sinio.value = sin(inclo);
  po = ao.value * omeosq.value;
  con42.value = 1.0 - 5.0 * cosio2.value;
  con41.value = -con42.value - cosio2.value - cosio2.value;
  ainv.value = 1.0 / ao.value;
  posq.value = po * po;
  rp.value = ao.value * (1.0 - ecco);
  method.value = 'n';

  // sgp4fix modern approach to finding sidereal time
  //   if (opsmode == 'a')
  //      {
  // sgp4fix use old way of finding gst
  // count integer number of days from 0 jan 1970
  ts70 = epoch - 7305.0;
  ds70 = (ts70 + 1.0e-8).floorToDouble();
  tfrac = ts70 - ds70;
  // find greenwich location at epoch
  c1 = 1.72027916940703639e-2;
  thgr70 = 1.7321343856509374;
  fk5r = 5.07551419432269442e-15;
  c1p2p = c1 + twoPi;
  var gsto1 = (thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r) % twoPi;
  if (gsto1 < 0.0) {
    gsto1 = gsto1 + twoPi;
  }
  //    }
  //    else
  gsto.value = gstime_SGP4(epoch + 2433281.5);

  //#include "debug5.cpp"
} // initl
