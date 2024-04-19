import 'dart:math';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/pointer.dart';

/// ----------------------------------------------------------------------------
///
///                           procedure dpper
///
///  this procedure provides deep space long period periodic contributions
///    to the mean elements.  by design, these periodics are zero at epoch.
///    this used to be dscom which included initialization, but it's really a
///    recurring function.
///
///  author        : david vallado                  719-573-2600   28 jun 2005
///
///  inputs        :
///    e3          -
///    ee2         -
///    peo         -
///    pgho        -
///    pho         -
///    pinco       -
///    plo         -
///    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
///    t           -
///    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
///    zmol        -
///    zmos        -
///    ep          - eccentricity                           0.0 - 1.0
///    inclo       - inclination - needed for lyddane modification
///    nodep       - right ascension of ascending node
///    argpp       - argument of perigee
///    mp          - mean anomaly
///
///  outputs       :
///    ep          - eccentricity                           0.0 - 1.0
///    inclp       - inclination
///    nodep        - right ascension of ascending node
///    argpp       - argument of perigee
///    mp          - mean anomaly
///
///  locals        :
///    alfdp       -
///    betdp       -
///    cosip  , sinip  , cosop  , sinop  ,
///    dalf        -
///    dbet        -
///    dls         -
///    f2, f3      -
///    pe          -
///    pgh         -
///    ph          -
///    pinc        -
///    pl          -
///    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
///    sll   , sls
///    xls         -
///    xnoh        -
///    zf          -
///    zm          -
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
void dpper(
    final double e3,
    final double ee2,
    final double peo,
    final double pgho,
    final double pho,
    final double pinco,
    final double plo,
    final double se2,
    final double se3,
    final double sgh2,
    final double sgh3,
    final double sgh4,
    final double sh2,
    final double sh3,
    final double si2,
    final double si3,
    final double sl2,
    final double sl3,
    final double sl4,
    final double t,
    final double xgh2,
    final double xgh3,
    final double xgh4,
    final double xh2,
    final double xh3,
    final double xi2,
    final double xi3,
    final double xl2,
    final double xl3,
    final double xl4,
    final double zmol,
    final double zmos,
    final double inclo,
    final String init,
    final Pointer<double> ep,
    final Pointer<double> inclp,
    final Pointer<double> nodep,
    final Pointer<double> argpp,
    final Pointer<double> mp,
    final String opsmode) {
  /* --------------------- local variables ------------------------ */
  double alfdp,
      betdp,
      cosip,
      cosop,
      dalf,
      dbet,
      dls,
      f2,
      f3,
      pe,
      pgh,
      ph,
      pinc,
      pl,
      sel,
      ses,
      sghl,
      sghs,
      shll,
      shs,
      sil,
      sinip,
      sinop,
      sinzf,
      sis,
      sll,
      sls,
      xls,
      xnoh,
      zf,
      zm;

  /* ---------------------- constants ----------------------------- */
  final zns = 1.19459e-5;
  final zes = 0.01675;
  final znl = 1.5835218e-4;
  final zel = 0.05490;

  /* --------------- calculate time varying periodics ----------- */
  zm = zmos + zns * t;
  // be sure that the initial call has time set to zero
  if (init == 'y') {
    zm = zmos;
  }
  zf = zm + 2.0 * zes * sin(zm);
  sinzf = sin(zf);
  f2 = 0.5 * sinzf * sinzf - 0.25;
  f3 = -0.5 * sinzf * cos(zf);
  ses = se2 * f2 + se3 * f3;
  sis = si2 * f2 + si3 * f3;
  sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
  sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
  shs = sh2 * f2 + sh3 * f3;
  zm = zmol + znl * t;
  if (init == 'y') {
    zm = zmol;
  }
  zf = zm + 2.0 * zel * sin(zm);
  sinzf = sin(zf);
  f2 = 0.5 * sinzf * sinzf - 0.25;
  f3 = -0.5 * sinzf * cos(zf);
  sel = ee2 * f2 + e3 * f3;
  sil = xi2 * f2 + xi3 * f3;
  sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
  sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
  shll = xh2 * f2 + xh3 * f3;
  pe = ses + sel;
  pinc = sis + sil;
  pl = sls + sll;
  pgh = sghs + sghl;
  ph = shs + shll;

  if (init == 'n') {
    pe = pe - peo;
    pinc = pinc - pinco;
    pl = pl - plo;
    pgh = pgh - pgho;
    ph = ph - pho;
    inclp.value = inclp.value + pinc;
    ep.value = ep.value + pe;
    sinip = sin(inclp.value);
    cosip = cos(inclp.value);

    /* ----------------- apply periodics directly ------------ */
    //  sgp4fix for lyddane choice
    //  strn3 used original inclination - this is technically feasible
    //  gsfc used perturbed inclination - also technically feasible
    //  probably best to readjust the 0.2 limit value and limit discontinuity
    //  0.2 rad = 11.45916 deg
    //  use next line for original strn3 approach and original inclination
    //  if (inclo >= 0.2)
    //  use next line for gsfc version and perturbed inclination
    if (inclp.value >= 0.2) {
      ph = ph / sinip;
      pgh = pgh - cosip * ph;
      argpp.value = argpp.value + pgh;
      nodep.value = nodep.value + ph;
      mp.value = mp.value + pl;
    } else {
      /* ---- apply periodics with lyddane modification ---- */
      sinop = sin(nodep.value);
      cosop = cos(nodep.value);
      alfdp = sinip * sinop;
      betdp = sinip * cosop;
      dalf = ph * cosop + pinc * cosip * sinop;
      dbet = -ph * sinop + pinc * cosip * cosop;
      alfdp = alfdp + dalf;
      betdp = betdp + dbet;
      nodep.value = nodep.value % twoPi;
      //  sgp4fix for afspc written intrinsic functions
      // nodep used without a trigonometric function ahead
      if ((nodep.value < 0.0) && (opsmode == 'a')) {
        nodep.value = nodep.value + twoPi;
      }
      xls = mp.value + argpp.value + cosip * nodep.value;
      dls = pl + pgh - pinc * nodep.value * sinip;
      xls = xls + dls;
      xnoh = nodep.value;
      nodep.value = atan2(alfdp, betdp);
      //  sgp4fix for afspc written intrinsic functions
      // nodep used without a trigonometric function ahead
      if ((nodep.value < 0.0) && (opsmode == 'a')) {
        nodep.value = nodep.value + twoPi;
      }
      if ((xnoh - nodep.value).abs() > pi) {
        if (nodep.value < xnoh) {
          nodep.value = nodep.value + twoPi;
        } else {
          nodep.value = nodep.value - twoPi;
        }
      }
      mp.value = mp.value + pl;
      argpp.value = xls - mp.value - cosip * nodep.value;
    }
  } // if init == 'n'

  //#include "debug1.cpp"
} // dpper
