// ignore_for_file: parameter_assignments, non_constant_identifier_names

import 'dart:math';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/pointer.dart';

/// ----------------------------------------------------------------------------
///
///                           function gstime_SGP4
///
///  this function finds the greenwich sidereal time.
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    jdut1       - julian date in ut1             days from 4713 bc
///
///  outputs       :
///    gstime      - greenwich sidereal time        0 to 2pi rad
///
///  locals        :
///    temp        - temporary variable for doubles   rad
///    tut1        - julian centuries from the
///                  jan 1, 2000 12 h epoch (ut1)
///
///  coupling      :
///    none
///
///  references    :
///    vallado       2013, 187, eq 3-45
/// --------------------------------------------------------------------------
double gstime_SGP4(final double jdut1) {
  double temp, tut1;

  tut1 = (jdut1 - 2451545.0) / 36525.0;
  temp = -6.2e-6 * tut1 * tut1 * tut1 +
      0.093104 * tut1 * tut1 +
      (876600.0 * 3600 + 8640184.812866) * tut1 +
      67310.54841; // sec
  temp = (temp * deg2rad / 240.0) % twoPi; //360/86400 = 1/240, to deg, to rad

  // ------------------------ check quadrants ---------------------
  if (temp < 0.0) {
    temp += twoPi;
  }

  return temp;
} // gstime

/// Calculate the sign of input variable [x].
double sgn_SGP4(final double x) {
  if (x < 0.0) {
    return -1.0;
  } else {
    return 1.0;
  }
} // sgn

/// ----------------------------------------------------------------------------
///
///                           function mag_SGP4
///
///  this procedure finds the magnitude of a vector.
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    vec         - vector
///
///  outputs       :
///    mag         - answer
///
///  locals        :
///    none.
///
///  coupling      :
///    none.
/// ----------------------------------------------------------------------------
double mag_SGP4(final List<double> x) =>
    sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
// mag

/// ----------------------------------------------------------------------------
///
///                           procedure cross_SGP4
///
///  this procedure crosses two vectors.
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    vec1        - vector number 1
///    vec2        - vector number 2
///
///  outputs       :
///    outvec      - vector result of a x b
///
///  locals        :
///    none.
///
///  coupling      :
///    mag           magnitude of a vector
/// ----------------------------------------------------------------------------
void cross_SGP4(final List<double> vec1, final List<double> vec2,
    final List<double> outvec) {
  outvec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  outvec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  outvec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
} // end cross

/// ----------------------------------------------------------------------------
///
///                           function dot_SGP4
///
///  this function finds the dot product of two vectors.
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    vec1        - vector number 1
///    vec2        - vector number 2
///
///  outputs       :
///    dot         - result
///
///  locals        :
///    none.
///
///  coupling      :
///    none.
/// ----------------------------------------------------------------------------
double dot_SGP4(final List<double> x, final List<double> y) =>
    x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; // dot

/// ----------------------------------------------------------------------------
///
///                           procedure angle_SGP4
///
///  this procedure calculates the angle between two vectors.  the output is
///    set to 999999.1 to indicate an undefined value.  be sure to check for
///    this at the output phase.
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    vec1        - vector number 1
///    vec2        - vector number 2
///
///  outputs       :
///    theta       - angle between the two vectors  -pi to pi
///
///  locals        :
///    temp        - temporary real variable
///
///  coupling      :
///    dot           dot product of two vectors
/// ----------------------------------------------------------------------------
double angle_SGP4(final List<double> vec1, final List<double> vec2) {
  double small, undefined, magv1, magv2, temp;
  small = 0.00000001;
  undefined = 999999.1;

  magv1 = mag_SGP4(vec1);
  magv2 = mag_SGP4(vec2);

  if (magv1 * magv2 > small * small) {
    temp = dot_SGP4(vec1, vec2) / (magv1 * magv2);
    if (temp.abs() > 1.0) {
      temp = sgn_SGP4(temp) * 1.0;
    }
    return acos(temp);
  } else {
    return undefined;
  }
} // angle

/// ----------------------------------------------------------------------------
///
///                           function asinh_SGP4
///
///  this function evaluates the inverse hyperbolic sine function.
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    xval        - angle value                                  any real
///
///  outputs       :
///    arcsinh     - result                                       any real
///
///  locals        :
///    none.
///
///  coupling      :
///    none.
/// ----------------------------------------------------------------------------
double asinh_SGP4(final double xval) => log(xval + sqrt(xval * xval + 1.0));
// asinh

/// ----------------------------------------------------------------------------
///
///                           function newtonnu_SGP4
///
///  this function solves keplers equation when the true anomaly is known.
///    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
///    the parabolic limit at 168\F8 is arbitrary. the hyperbolic anomaly is also
///    limited. the hyperbolic sine is used because it's not double valued.
///
///  author        : david vallado                  719-573-2600   27 may 2002
///
///  revisions
///    vallado     - fix small                                     24 sep 2002
///
///  inputs          description                    range / units
///    ecc         - eccentricity                   0.0  to
///    nu          - true anomaly                   -2pi to 2pi rad
///
///  outputs       :
///    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 \F8
///    m           - mean anomaly                   0.0  to 2pi rad       151.7425 \F8
///
///  locals        :
///    e1          - eccentric anomaly, next value  rad
///    sine        - sine of e
///    cose        - cosine of e
///    ktr         - index
///
///  coupling      :
///    asinh       - arc hyperbolic sine
///
///  references    :
///    vallado       2013, 77, alg 5
/// ----------------------------------------------------------------------------
void newtonnu_SGP4(final double ecc, final double nu, final Pointer<double> e0,
    final Pointer<double> m) {
  double small, sine, cose;

  // ---------------------  implementation   ---------------------
  e0.value = 999999.9;
  m.value = 999999.9;
  small = 0.00000001;

  // --------------------------- circular ------------------------
  if (ecc.abs() < small) {
    m.value = nu;
    e0.value = nu;
  } else
  // ---------------------- elliptical -----------------------
  if (ecc < 1.0 - small) {
    sine = (sqrt(1.0 - ecc * ecc) * sin(nu)) / (1.0 + ecc * cos(nu));
    cose = (ecc + cos(nu)) / (1.0 + ecc * cos(nu));
    e0.value = atan2(sine, cose);
    m.value = e0.value - ecc * sin(e0.value);
  } else
  // -------------------- hyperbolic  --------------------
  if (ecc > 1.0 + small) {
    if ((ecc > 1.0) && (nu.abs() + 0.00001 < pi - acos(1.0 / ecc))) {
      sine = (sqrt(ecc * ecc - 1.0) * sin(nu)) / (1.0 + ecc * cos(nu));
      e0.value = asinh_SGP4(sine);
      m.value = ecc * sinh(e0.value) - e0.value;
    }
  } else
  // ----------------- parabolic ---------------------
  if (nu.abs() < 168.0 * pi / 180.0) {
    e0.value = tan(nu * 0.5);
    m.value = e0.value + (e0.value * e0.value * e0.value) / 3.0;
  }

  if (ecc < 1.0) {
    m.value = m.value % 2.0 * pi;
    if (m.value < 0.0) {
      m.value = m.value + 2.0 * pi;
    }
    e0.value = e0.value % 2.0 * pi;
  }
} // newtonnu

/// ----------------------------------------------------------------------------
///
///                           function rv2coe_SGP4
///
///  this function finds the classical orbital elements given the geocentric
///    equatorial position and velocity vectors.
///
///  author        : david vallado                  719-573-2600   21 jun 2002
///
///  revisions
///    vallado     - fix special cases                              5 sep 2002
///    vallado     - delete extra check in inclination code        16 oct 2002
///    vallado     - add constant file use                         29 jun 2003
///    vallado     - add mu                                         2 apr 2007
///
///  inputs          description                    range / units
///    r           - ijk position vector            km
///    v           - ijk velocity vector            km / s
///    mu          - gravitational parameter        km3 / s2
///
///  outputs       :
///    p           - semilatus rectum               km
///    a           - semimajor axis                 km
///    ecc         - eccentricity
///    incl        - inclination                    0.0  to pi rad
///    omega       - right ascension of ascending node    0.0  to 2pi rad
///    argp        - argument of perigee            0.0  to 2pi rad
///    nu          - true anomaly                   0.0  to 2pi rad
///    m           - mean anomaly                   0.0  to 2pi rad
///    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
///    truelon     - true longitude            (ce) 0.0  to 2pi rad
///    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
///
///  locals        :
///    hbar        - angular momentum h vector      km2 / s
///    ebar        - eccentricity     e vector
///    nbar        - line of nodes    n vector
///    c1          - v**2 - u/r
///    rdotv       - r dot v
///    hk          - hk unit vector
///    sme         - specfic mechanical energy      km2 / s2
///    i           - index
///    e           - eccentric, parabolic,
///                  hyperbolic anomaly             rad
///    temp        - temporary variable
///    typeorbit   - type of orbit                  ee, ei, ce, ci
///
///  coupling      :
///    mag         - magnitude of a vector
///    cross       - cross product of two vectors
///    angle       - find the angle between two vectors
///    newtonnu    - find the mean anomaly
///
///  references    :
///    vallado       2013, 113, alg 9, ex 2-5
/// ----------------------------------------------------------------------------
void rv2coe_SGP4(
    final List<double> r,
    final List<double> v,
    final double mus,
    final Pointer<double> p,
    final Pointer<double> a,
    final Pointer<double> ecc,
    final Pointer<double> incl,
    final Pointer<double> omega,
    final Pointer<double> argp,
    final Pointer<double> nu,
    final Pointer<double> m,
    final Pointer<double> arglat,
    final Pointer<double> truelon,
    final Pointer<double> lonper) {
  double undefined,
      small,
      magr,
      magv,
      magn,
      sme,
      rdotv,
      infinite,
      temp,
      c1,
      hk,
      magh;
  final hbar = <double>[0, 0, 0],
      nbar = <double>[0, 0, 0],
      ebar = <double>[0, 0, 0];
  final e = Pointer<double>(0.0);

  int i;
  // switch this to an integer msvs seems to have probelms with this and strncpy_s
  //char typeorbit[2];
  int typeorbit;
  // here
  // typeorbit = 1 = 'ei'
  // typeorbit = 2 = 'ce'
  // typeorbit = 3 = 'ci'
  // typeorbit = 4 = 'ee'

  small = 0.00000001;
  undefined = 999999.1;
  infinite = 999999.9;

  // -------------------------  implementation   -----------------
  magr = mag_SGP4(r);
  magv = mag_SGP4(v);

  // ------------------  find h n and e vectors   ----------------
  cross_SGP4(r, v, hbar);
  magh = mag_SGP4(hbar);
  if (magh > small) {
    nbar[0] = -hbar[1];
    nbar[1] = hbar[0];
    nbar[2] = 0.0;
    magn = mag_SGP4(nbar);
    c1 = magv * magv - mus / magr;
    rdotv = dot_SGP4(r, v);
    for (i = 0; i <= 2; i++) {
      ebar[i] = (c1 * r[i] - rdotv * v[i]) / mus;
    }
    ecc.value = mag_SGP4(ebar);

    // ------------  find a e and semi-latus rectum   ----------
    sme = (magv * magv * 0.5) - (mus / magr);
    if (sme.abs() > small) {
      a.value = -mus / (2.0 * sme);
    } else {
      a.value = infinite;
    }
    p.value = magh * magh / mus;

    // -----------------  find inclination   -------------------
    hk = hbar[2] / magh;
    incl.value = acos(hk);

    // --------  determine type of orbit for later use  --------
    // ------ elliptical, parabolic, hyperbolic inclined -------
    //#ifdef _MSC_VER  // chk if compiling under MSVS
    //		   strcpy_s(typeorbit, 2 * sizeof(char), "ei");
    //#else
    //		   strcpy(typeorbit, "ei");
    //#endif
    typeorbit = 1;

    if (ecc.value < small) {
      // ----------------  circular equatorial ---------------
      if ((incl.value < small) | ((incl.value - pi).abs() < small)) {
        //#ifdef _MSC_VER
        //				   strcpy_s(typeorbit, sizeof(typeorbit), "ce");
        //#else
        //				   strcpy(typeorbit, "ce");
        //#endif
        typeorbit = 2;
      } else {
        // --------------  circular inclined ---------------
        //#ifdef _MSC_VER
        //				   strcpy_s(typeorbit, sizeof(typeorbit), "ci");
        //#else
        //				   strcpy(typeorbit, "ci");
        //#endif
        typeorbit = 3;
      }
    } else {
      // - elliptical, parabolic, hyperbolic equatorial --
      if ((incl.value < small) | ((incl.value - pi).abs() < small)) {
        //#ifdef _MSC_VER
        //				   strcpy_s(typeorbit, sizeof(typeorbit), "ee");
        //#else
        //				   strcpy(typeorbit, "ee");
        //#endif
        typeorbit = 4;
      }
    }

    // ----------  find right ascension of the ascending node ------------
    if (magn > small) {
      temp = nbar[0] / magn;
      if (temp.abs() > 1.0) {
        temp = sgn_SGP4(temp);
      }
      omega.value = acos(temp);
      if (nbar[1] < 0.0) {
        omega.value = twoPi - omega.value;
      }
    } else {
      omega.value = undefined;
    }

    // ---------------- find argument of perigee ---------------
    //if (strcmp(typeorbit, "ei") == 0)
    if (typeorbit == 1) {
      argp.value = angle_SGP4(nbar, ebar);
      if (ebar[2] < 0.0) {
        argp.value = twoPi - argp.value;
      }
    } else {
      argp.value = undefined;
    }

    // ------------  find true anomaly at epoch    -------------
    //if (typeorbit[0] == 'e')
    if ((typeorbit == 1) || (typeorbit == 4)) {
      nu.value = angle_SGP4(ebar, r);
      if (rdotv < 0.0) {
        nu.value = twoPi - nu.value;
      }
    } else {
      nu.value = undefined;
    }

    // ----  find argument of latitude - circular inclined -----
    //if (strcmp(typeorbit, "ci") == 0)
    if (typeorbit == 3) {
      arglat.value = angle_SGP4(nbar, r);
      if (r[2] < 0.0) {
        arglat.value = twoPi - arglat.value;
      }
      m.value = arglat.value;
    } else {
      arglat.value = undefined;
    }

    // -- find longitude of perigee - elliptical equatorial ----
    //if ((ecc>small) && (strcmp(typeorbit, "ee") == 0))
    if ((ecc.value > small) && (typeorbit == 4)) {
      temp = ebar[0] / ecc.value;
      if (temp.abs() > 1.0) {
        temp = sgn_SGP4(temp);
      }
      lonper.value = acos(temp);
      if (ebar[1] < 0.0) {
        lonper.value = twoPi - lonper.value;
      }
      if (incl.value > halfPi) {
        lonper.value = twoPi - lonper.value;
      }
    } else {
      lonper.value = undefined;
    }

    // -------- find true longitude - circular equatorial ------
    //if ((magr>small) && (strcmp(typeorbit, "ce") == 0))
    if ((magr > small) && (typeorbit == 2)) {
      temp = r[0] / magr;
      if (temp.abs() > 1.0) {
        temp = sgn_SGP4(temp);
      }
      truelon.value = acos(temp);
      if (r[1] < 0.0) {
        truelon.value = twoPi - truelon.value;
      }
      if (incl.value > halfPi) {
        truelon.value = twoPi - truelon.value;
      }
      m.value = truelon.value;
    } else {
      truelon.value = undefined;
    }

    // ------------ find mean anomaly for all orbits -----------
    //if (typeorbit[0] == 'e')
    if ((typeorbit == 1) || (typeorbit == 4)) {
      newtonnu_SGP4(ecc.value, nu.value, e, m);
    }
  } else {
    p.value = undefined;
    a.value = undefined;
    ecc.value = undefined;
    incl.value = undefined;
    omega.value = undefined;
    argp.value = undefined;
    nu.value = undefined;
    m.value = undefined;
    arglat.value = undefined;
    truelon.value = undefined;
    lonper.value = undefined;
  }
} // rv2coe

/// ----------------------------------------------------------------------------
///
///                           procedure jday_SGP4
///
///  this procedure finds the julian date given the year, month, day, and time.
///  the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
///
///  algorithm     : calculate the answer in one step for efficiency
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    year        - year                           1900 .. 2100
///    mon         - month                          1 .. 12
///    day         - day                            1 .. 28,29,30,31
///    hr          - universal time hour            0 .. 23
///    min         - universal time min             0 .. 59
///    sec         - universal time sec             0.0 .. 59.999
///
///  outputs       :
///    jd          - julian date                    days from 4713 bc
///    jdfrac      - julian date fraction into day  days from 4713 bc
///
///  locals        :
///    none.
///
///  coupling      :
///    none.
///
///  references    :
///    vallado       2013, 183, alg 14, ex 3-4
/// ----------------------------------------------------------------------------
void jday_SGP4(
    final int year,
    final int mon,
    final int day,
    final int hr,
    final int minute,
    final double sec,
    final Pointer<double> jd,
    final Pointer<double> jdFrac) {
  jd.value = 367.0 * year -
      ((7 * (year + ((mon + 9) / 12.0).floor())) * 0.25).floor() +
      (275 * mon / 9.0).floor() +
      day +
      1721013.5; // use - 678987.0 to go to mjd directly
  jdFrac.value = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;

  // check that the day and fractional day are correct
  if (jdFrac.value.abs() > 1.0) {
    final dtt = jdFrac.value.floorToDouble();
    jd.value = jd.value + dtt;
    jdFrac.value = jdFrac.value - dtt;
  }

  // - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
} // jday

/// ----------------------------------------------------------------------------
///
///                           procedure days2mdhms_SGP4
///
///  this procedure converts the day of the year, days, to the equivalent month
///    day, hour, minute and second.
///
///  algorithm     : set up array for the number of days per month
///                  find leap year - use 1900 because 2000 is a leap year
///                  loop through a temp value while the value is < the days
///                  perform int conversions to the correct day and month
///                  convert remainder into h m s using type conversions
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    year        - year                           1900 .. 2100
///    days        - julian day of the year         1.0  .. 366.0
///
///  outputs       :
///    mon         - month                          1 .. 12
///    day         - day                            1 .. 28,29,30,31
///    hr          - hour                           0 .. 23
///    min         - minute                         0 .. 59
///    sec         - second                         0.0 .. 59.999
///
///  locals        :
///    dayofyr     - day of year
///    temp        - temporary extended values
///    inttemp     - temporary int value
///    i           - index
///    lmonth[13]  - int array containing the number of days per month
///
///  coupling      :
///    none.
/// ----------------------------------------------------------------------------
void days2mdhms_SGP4(
    final int year,
    final double days,
    final Pointer<int> mon,
    final Pointer<int> day,
    final Pointer<int> hr,
    final Pointer<int> minute,
    final Pointer<double> sec) {
  int i, inttemp, dayofyr;
  double temp;
  final lmonth = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

  dayofyr = days.floor();
  /* ----------------- find month and day of month ---------------- */
  if ((year % 4) == 0) {
    lmonth[2] = 29;
  }

  i = 1;
  inttemp = 0;
  while ((dayofyr > inttemp + lmonth[i]) && (i < 12)) {
    inttemp = inttemp + lmonth[i];
    i++;
  }
  mon.value = i;
  day.value = dayofyr - inttemp;

  /* ----------------- find hours minutes and seconds ------------- */
  temp = (days - dayofyr) * 24.0;
  hr.value = temp.floor();
  temp = (temp - hr.value) * 60.0;
  minute.value = temp.floor();
  sec.value = (temp - minute.value) * 60.0;
} // days2mdhms

/// ----------------------------------------------------------------------------
///
///                           procedure invjday_SGP4
///
///  this procedure finds the year, month, day, hour, minute and second
///  given the julian date. tu can be ut1, tdt, tdb, etc.
///
///  algorithm     : set up starting values
///                  find leap year - use 1900 because 2000 is a leap year
///                  find the elapsed days through the year in a loop
///                  call routine to find each individual value
///
///  author        : david vallado                  719-573-2600    1 mar 2001
///
///  inputs          description                    range / units
///    jd          - julian date                    days from 4713 bc
///    jdfrac      - julian date fraction into day  days from 4713 bc
///
///  outputs       :
///    year        - year                           1900 .. 2100
///    mon         - month                          1 .. 12
///    day         - day                            1 .. 28,29,30,31
///    hr          - hour                           0 .. 23
///    min         - minute                         0 .. 59
///    sec         - second                         0.0 .. 59.999
///
///  locals        :
///    days        - day of year plus fractional
///                  portion of a day               days
///    tu          - julian centuries from 0 h
///                  jan 0, 1900
///    temp        - temporary double values
///    leapyrs     - number of leap years from 1900
///
///  coupling      :
///    days2mdhms  - finds month, day, hour, minute and second given days
///                  and year
///
///  references    :
///    vallado       2013, 203, alg 22, ex 3-13
/// ----------------------------------------------------------------------------
void invjday_SGP4(
    double jd,
    double jdfrac,
    final Pointer<int> year,
    final Pointer<int> mon,
    final Pointer<int> day,
    final Pointer<int> hr,
    final Pointer<int> minute,
    final Pointer<double> sec) {
  int leapyrs;
  double dt, days, tu, temp;

  // check jdfrac for multiple days
  if (jdfrac.abs() >= 1.0) {
    jd = jd + jdfrac.floor();
    jdfrac = jdfrac - jdfrac.floor();
  }

  // check for fraction of a day included in the jd
  dt = jd - jd.floor() - 0.5;
  if (dt.abs() > 0.00000001) {
    jd = jd - dt;
    jdfrac = jdfrac + dt;
  }

  /* --------------- find year and days of the year --------------- */
  temp = jd - 2415019.5;
  tu = temp / 365.25;
  year.value = 1900 + tu.floor();
  leapyrs = ((year.value - 1901) * 0.25).floor();

  days = (temp - ((year.value - 1900) * 365.0 + leapyrs)).floorToDouble();

  /* ------------ check for case of beginning of a year ----------- */
  if (days + jdfrac < 1.0) {
    year.value = year.value - 1;
    leapyrs = ((year.value - 1901) * 0.25).floor();
    days = (temp - ((year.value - 1900) * 365.0 + leapyrs)).floorToDouble();
  }

  /* ----------------- find remaining data  ------------------------- */
  days2mdhms_SGP4(year.value, days + jdfrac, mon, day, hr, minute, sec);
} // invjday
