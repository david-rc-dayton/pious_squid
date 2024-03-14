// ignore_for_file: non_constant_identifier_names, public_member_api_docs

import 'package:pious_squid/src/operations/pointer.dart';

/// Element set record.
class ElsetRec {
  final satnum = Pointer('');
  final epochyr = Pointer(0), epochtynumrev = Pointer(0);
  final error = Pointer(0);
  final operationmode = Pointer('');
  final init = Pointer(''), method = Pointer('');

  /* Near Earth */
  final isimp = Pointer(0);
  final aycof = Pointer(0.0),
      con41 = Pointer(0.0),
      cc1 = Pointer(0.0),
      cc4 = Pointer(0.0),
      cc5 = Pointer(0.0),
      d2 = Pointer(0.0),
      d3 = Pointer(0.0),
      d4 = Pointer(0.0),
      delmo = Pointer(0.0),
      eta = Pointer(0.0),
      argpdot = Pointer(0.0),
      omgcof = Pointer(0.0),
      sinmao = Pointer(0.0),
      t = Pointer(0.0),
      t2cof = Pointer(0.0),
      t3cof = Pointer(0.0),
      t4cof = Pointer(0.0),
      t5cof = Pointer(0.0),
      x1mth2 = Pointer(0.0),
      x7thm1 = Pointer(0.0),
      mdot = Pointer(0.0),
      nodedot = Pointer(0.0),
      xlcof = Pointer(0.0),
      xmcof = Pointer(0.0),
      nodecf = Pointer(0.0);

  /* Deep Space */
  final irez = Pointer(0);
  final d2201 = Pointer(0.0),
      d2211 = Pointer(0.0),
      d3210 = Pointer(0.0),
      d3222 = Pointer(0.0),
      d4410 = Pointer(0.0),
      d4422 = Pointer(0.0),
      d5220 = Pointer(0.0),
      d5232 = Pointer(0.0),
      d5421 = Pointer(0.0),
      d5433 = Pointer(0.0),
      dedt = Pointer(0.0),
      del1 = Pointer(0.0),
      del2 = Pointer(0.0),
      del3 = Pointer(0.0),
      didt = Pointer(0.0),
      dmdt = Pointer(0.0),
      dnodt = Pointer(0.0),
      domdt = Pointer(0.0),
      e3 = Pointer(0.0),
      ee2 = Pointer(0.0),
      peo = Pointer(0.0),
      pgho = Pointer(0.0),
      pho = Pointer(0.0),
      pinco = Pointer(0.0),
      plo = Pointer(0.0),
      se2 = Pointer(0.0),
      se3 = Pointer(0.0),
      sgh2 = Pointer(0.0),
      sgh3 = Pointer(0.0),
      sgh4 = Pointer(0.0),
      sh2 = Pointer(0.0),
      sh3 = Pointer(0.0),
      si2 = Pointer(0.0),
      si3 = Pointer(0.0),
      sl2 = Pointer(0.0),
      sl3 = Pointer(0.0),
      sl4 = Pointer(0.0),
      gsto = Pointer(0.0),
      xfact = Pointer(0.0),
      xgh2 = Pointer(0.0),
      xgh3 = Pointer(0.0),
      xgh4 = Pointer(0.0),
      xh2 = Pointer(0.0),
      xh3 = Pointer(0.0),
      xi2 = Pointer(0.0),
      xi3 = Pointer(0.0),
      xl2 = Pointer(0.0),
      xl3 = Pointer(0.0),
      xl4 = Pointer(0.0),
      xlamo = Pointer(0.0),
      zmol = Pointer(0.0),
      zmos = Pointer(0.0),
      atime = Pointer(0.0),
      xli = Pointer(0.0),
      xni = Pointer(0.0);

  final a = Pointer(0.0),
      altp = Pointer(0.0),
      alta = Pointer(0.0),
      epochdays = Pointer(0.0),
      jdsatepoch = Pointer(0.0),
      jdsatepochF = Pointer(0.0),
      nddot = Pointer(0.0),
      ndot = Pointer(0.0),
      bstar = Pointer(0.0),
      rcse = Pointer(0.0),
      inclo = Pointer(0.0),
      nodeo = Pointer(0.0),
      ecco = Pointer(0.0),
      argpo = Pointer(0.0),
      mo = Pointer(0.0),
      no_kozai = Pointer(0.0);

  // sgp4fix add new variables from tle
  final classification = Pointer(''), intldesg = Pointer('');
  final ephtype = Pointer(0);
  final elnum = Pointer(0), revnum = Pointer(0);

  // sgp4fix add unkozai'd variable
  final no_unkozai = Pointer(0.0);

  // sgp4fix add singly averaged variables
  final am = Pointer(0.0),
      em = Pointer(0.0),
      im = Pointer(0.0),
      Om = Pointer(0.0),
      om = Pointer(0.0),
      mm = Pointer(0.0),
      nm = Pointer(0.0);

  // sgp4fix add constant parameters to eliminate mutliple calls during execution
  final tumin = Pointer(0.0),
      mus = Pointer(0.0),
      radiusearthkm = Pointer(0.0),
      xke = Pointer(0.0),
      j2 = Pointer(0.0),
      j3 = Pointer(0.0),
      j4 = Pointer(0.0),
      j3oj2 = Pointer(0.0);

  //       Additional elements to capture relevant TLE and object information:
  final dia_mm = Pointer(0); // RSO dia in mm
  final period_sec = Pointer(0.0); // Period in seconds
  final active = Pointer(0); // "Active S/C" flag (0=n, 1=y)
  final not_orbital = Pointer(0); // "Orbiting S/C" flag (0=n, 1=y)
  final rcs_m2 = Pointer(0.0); // "RCS (m^2)" storage
}
