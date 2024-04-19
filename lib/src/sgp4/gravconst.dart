import 'dart:math';

import 'package:pious_squid/src/operations/pointer.dart';

/// SGP4 Gravitational Constants
enum Sgp4GravConst {
  /// WGS-72 _(old)_
  wgs72old,

  /// WGS-72
  wgs72,

  /// WGS-84
  wgs84
}

/// ----------------------------------------------------------------------------
///
///                           function getgravconst
///
///  this function gets constants for the propagator. note that mu is identified
///  to facilitiate comparisons with newer models. the common useage is wgs72.
///
///  author        : david vallado                  719-573-2600   21 jul 2006
///
///  inputs        :
///    whichconst  - which set of constants to use  wgs72old, wgs72, wgs84
///
///  outputs       :
///    tumin       - minutes in one time unit
///    mu          - earth gravitational parameter
///    radiusearthkm - radius of the earth in km
///    xke         - reciprocal of tumin
///    j2, j3, j4  - un-normalized zonal harmonic values
///    j3oj2       - j3 divided by j2
///
///  locals        :
///
///  coupling      :
///    none
///
///  references    :
///    norad spacetrack report #3
///    vallado, crawford, hujsak, kelso  2006
///	----------------------------------------------------------------------------
void getgravconst(
    final Sgp4GravConst whichconst,
    final Pointer<double> tumin,
    final Pointer<double> mus,
    final Pointer<double> radiusearthkm,
    final Pointer<double> xke,
    final Pointer<double> j2,
    final Pointer<double> j3,
    final Pointer<double> j4,
    final Pointer<double> j3oj2) {
  switch (whichconst) {
    // -- wgs-72 low precision str#3 constants --
    case Sgp4GravConst.wgs72old:
      mus.value = 398600.79964; // in km3 / s2
      radiusearthkm.value = 6378.135; // km
      xke.value = 0.0743669161; // reciprocal of tumin
      tumin.value = 1.0 / xke.value;
      j2.value = 0.001082616;
      j3.value = -0.00000253881;
      j4.value = -0.00000165597;
      j3oj2.value = j3.value / j2.value;
      break;
    // ------------ wgs-72 constants ------------
    case Sgp4GravConst.wgs72:
      mus.value = 398600.8; // in km3 / s2
      radiusearthkm.value = 6378.135; // km
      xke.value = 60.0 /
          sqrt(radiusearthkm.value *
              radiusearthkm.value *
              radiusearthkm.value /
              mus.value);
      tumin.value = 1.0 / xke.value;
      j2.value = 0.001082616;
      j3.value = -0.00000253881;
      j4.value = -0.00000165597;
      j3oj2.value = j3.value / j2.value;
      break;
    case Sgp4GravConst.wgs84:
      // ------------ wgs-84 constants ------------
      mus.value = 398600.5; // in km3 / s2
      radiusearthkm.value = 6378.137; // km
      xke.value = 60.0 /
          sqrt(radiusearthkm.value *
              radiusearthkm.value *
              radiusearthkm.value /
              mus.value);
      tumin.value = 1.0 / xke.value;
      j2.value = 0.00108262998905;
      j3.value = -0.00000253215306;
      j4.value = -0.00000161098761;
      j3oj2.value = j3.value / j2.value;
      break;
  }
} // getgravconst
