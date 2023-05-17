// ignore_for_file: parameter_assignments

import 'dart:math';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/sgp4/elsetrec.dart';
import 'package:pious_squid/src/sgp4/gravconst.dart';
import 'package:pious_squid/src/sgp4/pointer.dart';
import 'package:pious_squid/src/sgp4/sgp4init.dart';
import 'package:pious_squid/src/sgp4/utilities.dart';

final _blankLine = ''.padLeft(130, ' ').split('');

/// ----------------------------------------------------------------------------
///
///                           function twoline2rv
///
///  this function converts the two line element set character string data to
///    variables and initializes the sgp4 variables. several intermediate varaibles
///    and quantities are determined. note that the result is a structure so multiple
///    satellites can be processed simaltaneously without having to reinitialize. the
///    verification mode is an important option that permits quick checks of any
///    changes to the underlying technical theory. this option works using a
///    modified tle file in which the start, stop, and delta time values are
///    included at the end of the second line of data. this only works with the
///    verification mode. the catalog mode simply propagates from -1440 to 1440 min
///    from epoch and is useful when performing entire catalog runs.
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs        :
///    longstr1    - first line of the tle
///    longstr2    - second line of the tle
///    typerun     - type of run                    verification 'v', catalog 'c',
///                                                 manual 'm'
///    typeinput   - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
///    opsmode     - mode of operation afspc or improved 'a', 'i'
///    whichconst  - which set of constants to use  72, 84
///
///  outputs       :
///    satrec      - structure containing all the sgp4 satellite information
///
///  coupling      :
///    getgravconst-
///    days2mdhms  - conversion of days to month, day, hour, minute, second
///    jday        - convert day month year hour minute second into julian date
///    sgp4init    - initialize the sgp4 variables
///
///  references    :
///    norad spacetrack report #3
///    vallado, crawford, hujsak, kelso  2006
/// ----------------------------------------------------------------------------
void twoline2rv(String line1, String line2, final String opsmode,
    final Sgp4GravConst whichconst, final ElsetRec satrec) {
  final longstr1 = [...line1.split(''), ..._blankLine].sublist(0, 130);
  final longstr2 = [...line2.split(''), ..._blankLine].sublist(0, 130);

  const xpdotp = 1440.0 / (2.0 * pi); // 229.1831180523293

  int j;
  final cardnumb = Pointer(0);
  // sgp4fix include in satrec
  // long revnum = 0, elnum = 0;
  // char classification, intldesg[11];
  var year = 0;
  final nexp = Pointer(0),
      ibexp = Pointer(0),
      mon = Pointer(0),
      day = Pointer(0),
      hr = Pointer(0),
      minute = Pointer(0);
  final sec = Pointer(0.0);

  // sgp4fix no longer needed
  // getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );

  satrec.error.value = 0;

  // set the implied decimal points since doing a formated read
  // fixes for bad input data values (missing, ...)
  for (j = 10; j <= 15; j++) {
    if (longstr1[j] == ' ') {
      longstr1[j] = '_';
    }
  }

  if (longstr1[44] != ' ') {
    longstr1[43] = longstr1[44];
  }
  longstr1[44] = '.';
  if (longstr1[7] == ' ') {
    longstr1[7] = 'U';
  }
  if (longstr1[9] == ' ') {
    longstr1[9] = '.';
  }
  for (j = 45; j <= 49; j++) {
    if (longstr1[j] == ' ') {
      longstr1[j] = '0';
    }
  }
  if (longstr1[51] == ' ') {
    longstr1[51] = '0';
  }
  if (longstr1[53] != ' ') {
    longstr1[52] = longstr1[53];
  }
  longstr1[53] = '.';
  longstr2[25] = '.';
  for (j = 26; j <= 32; j++) {
    if (longstr2[j] == ' ') {
      longstr2[j] = '0';
    }
  }
  if (longstr1[62] == ' ') {
    longstr1[62] = '0';
  }
  if (longstr1[68] == ' ') {
    longstr1[68] = '0';
  }

  // ---- parse line1 ----
  // 0         1         2         3         4         5         6
  // 012345678901234567890123456789012345678901234567890123456789012345678
  // 1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927 (pre)
  // 1 25544U 98067A_  08264.51782528 -.00002182 .00000-0-.11606-4 0  2927 (post)
  line1 = longstr1.join('');
  cardnumb.value = int.parse(line1.substring(0, 1));
  satrec.satnum.value = line1.substring(2, 7).trim();
  satrec.classification.value = line1.substring(7, 8);
  satrec.intldesg.value = line1.substring(9, 17).replaceAll('_', '').trim();
  satrec.epochyr.value = int.parse(line1.substring(18, 20));
  satrec.epochdays.value = double.parse(line1.substring(20, 32));
  satrec.ndot.value = double.parse(line1.substring(33, 43));
  satrec.nddot.value = double.parse(line1.substring(43, 50));
  nexp.value = int.parse(line1.substring(50, 52));
  satrec.bstar.value = double.parse(line1.substring(52, 59));
  ibexp.value = int.parse(line1.substring(59, 61));
  satrec.ephtype.value = int.parse(line1.substring(62, 63));
  satrec.elnum.value = int.parse(line1.substring(64, 68));

  // ---- parse line2 ----
  // 0         1         2         3         4         5         6
  // 012345678901234567890123456789012345678901234567890123456789012345678
  // 2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537 (pre)
  // 2 25544  51.6416 247.4627.0006703 130.5360 325.0288 15.72125391563537 (post)
  line2 = longstr2.join('');
  cardnumb.value = int.parse(line2.substring(0, 1));
  satrec.satnum.value = line2.substring(2, 7).trim();
  satrec.inclo.value = double.parse(line2.substring(8, 16));
  satrec.nodeo.value = double.parse(line2.substring(17, 25));
  satrec.ecco.value = double.parse(line2.substring(25, 33));
  satrec.argpo.value = double.parse(line2.substring(34, 42));
  satrec.mo.value = double.parse(line2.substring(43, 51));
  satrec.no_kozai.value = double.parse(line2.substring(51, 63));
  satrec.revnum.value = int.parse(line2.substring(63, 68));

  // ---- find no, ndot, nddot ----
  satrec.no_kozai.value = satrec.no_kozai.value / xpdotp; //* rad/min
  satrec.nddot.value = satrec.nddot.value * pow(10.0, nexp.value);
  satrec.bstar.value = satrec.bstar.value * pow(10.0, ibexp.value);

  // ---- convert to sgp4 units ----
  // satrec.a    = pow( satrec.no_kozai*tumin , (-2.0/3.0) );
  satrec.ndot.value = satrec.ndot.value / (xpdotp * 1440.0); //* ? * minperday
  satrec.nddot.value = satrec.nddot.value / (xpdotp * 1440.0 * 1440);

  // ---- find standard orbital elements ----
  satrec.inclo.value = satrec.inclo.value * deg2rad;
  satrec.nodeo.value = satrec.nodeo.value * deg2rad;
  satrec.argpo.value = satrec.argpo.value * deg2rad;
  satrec.mo.value = satrec.mo.value * deg2rad;

  // sgp4fix not needed here
  // satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
  // satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;

  // ----------------------------------------------------------------
  // find sgp4epoch time of element set
  // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
  // and minutes from the epoch (time)
  // ----------------------------------------------------------------

  // ---------------- temp fix for years from 1957-2056 -------------------
  // --------- correct fix will occur when year is 4-digit in tle ---------
  if (satrec.epochyr.value < 57) {
    year = satrec.epochyr.value + 2000;
  } else {
    year = satrec.epochyr.value + 1900;
  }

  days2mdhms_SGP4(year, satrec.epochdays.value, mon, day, hr, minute, sec);
  jday_SGP4(year, mon.value, day.value, hr.value, minute.value, sec.value,
      satrec.jdsatepoch, satrec.jdsatepochF);

  // ---------------- initialize the orbit at sgp4epoch -------------------
  sgp4init(
      whichconst,
      opsmode,
      satrec.satnum.value,
      (satrec.jdsatepoch.value + satrec.jdsatepochF.value) - 2433281.5,
      satrec.bstar.value,
      satrec.ndot.value,
      satrec.nddot.value,
      satrec.ecco.value,
      satrec.argpo.value,
      satrec.inclo.value,
      satrec.mo.value,
      satrec.no_kozai.value,
      satrec.nodeo.value,
      satrec);
} // twoline2rv
