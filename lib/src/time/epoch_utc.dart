import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/data/data_handler.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/time/epoch.dart';
import 'package:pious_squid/src/time/time_base.dart';
import 'package:pious_squid/src/time/time_scales.dart';

final Float64List _gmstPoly = Float64List.fromList(
    [-6.2e-6, 0.093104, 876600 * 3600 + 8640184.812866, 67310.54841]);

const List<List<int>> _dayOfYearLookup = [
  [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],
  [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
];

bool _isLeapYear(final int year) =>
    (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);

int _dayOfYear(final int year, final int month, final int day) {
  final dex = _isLeapYear(year) ? 1 : 0;
  return _dayOfYearLookup[dex][month - 1] + day - 1;
}

double _dateToPosix(final int year, final int month, final int day,
    final int hour, final int minute, final double second) {
  final days = _dayOfYear(year, month, day);
  final yearMod = year - 1900;
  return (minute * 60 +
          hour * 3600 +
          days * 86400 +
          (yearMod - 70) * 31536000 +
          ((yearMod - 69) ~/ 4) * 86400 -
          ((yearMod - 1) ~/ 100) * 86400 +
          ((yearMod + 299) ~/ 400) * 86400) +
      second;
}

/// Universal Coordinated Time _(UTC)_.
class EpochUTC extends Epoch {
  /// Create a new [EpochUTC] epoch given a [posix] timestamp _(s)_.
  EpochUTC(final double posix) : super(posix);

  /// Create a new [EpochUTC] object from an ISO formatted string of the form
  /// `yyyy-mm-ddThh:mm:ss.sssZ`, e.g.: `2022-11-26T12:14:11.362Z`.
  factory EpochUTC.fromDateTimeString(final String dateTimeString) {
    final dts = dateTimeString.trim().toUpperCase().endsWith('Z')
        ? dateTimeString
        : '${dateTimeString}Z';
    return EpochUTC(
        DateTime.parse(dts).toUtc().millisecondsSinceEpoch.toDouble() / 1000);
  }

  /// Create a new [EpochUTC] from the number of seconds from the J2000 epoch
  /// in Terrestrial Time: `2000-01-01T12:00:00.000 TT`.
  factory EpochUTC.fromJ2000TTSeconds(final double seconds) {
    final tInit = EpochUTC(seconds + 946728000);
    final ls = DataHandler().getLeapSeconds(tInit.toJulianDate());
    return tInit.roll(-32.184 - ls);
  }

  /// Create a new [EpochUTC] object from a definitive ephemeris epoch string of
  /// the form `ddd/yyyy hh:mm:ss.sss`, e.g: `156/2021 00:00:00.000`.
  factory EpochUTC.fromDefinitiveString(final String definitiveString) {
    final fields = definitiveString.trim().split(' ');
    final dateFields = fields[0].split('/');
    final day = int.parse(dateFields[0]);
    final year = int.parse(dateFields[1]);
    final timeField = fields[1];
    final dts = DateTime.parse('$year-01-01T${timeField}Z')
        .add(Duration(days: day - 1));
    return EpochUTC(dts.toUtc().millisecondsSinceEpoch.toDouble() / 1000);
  }

  /// Create a new [EpochUTC] object representing the current system time.
  EpochUTC.now()
      : super(DateTime.now().toUtc().millisecondsSinceEpoch.toDouble() / 1000);

  /// Create a new [EpochUTC] object for a given [year], [month], [day], [hour],
  /// [minute], and [second].
  EpochUTC.fromDate(final int year, final int month, final int day,
      final int hour, final int minute, final double second)
      : super(_dateToPosix(year, month, day, hour, minute, second));

  /// Create a new [EpochUTC] from from a [DateTime] object.
  EpochUTC.fromDateTime(final DateTime dt)
      : super(dt.toUtc().millisecondsSinceEpoch / 1000);

  /// Create a new [EpochUTC] object, offset by the provided number of
  /// [seconds].
  EpochUTC roll(final double seconds) => EpochUTC(posix + seconds);

  /// Convert to Modified Julian Date _(MJD)_.
  double toMjd() => toJulianDate() - 2400000.5;

  /// Convert to Modified Julian Date _(MJD)_ as used by the Goddard Space
  /// Flight Center _(GSFC)_.
  double toMjdGsfc() => toMjd() - 29999.5;

  /// Convert to International Atomic Time _(TAI)_.
  EpochTAI toTAI() {
    final ls = DataHandler().getLeapSeconds(toJulianDate());
    return EpochTAI(posix + ls);
  }

  /// Convert to Terrestrial Time _(TT)_.
  EpochTT toTT() => EpochTT(toTAI().posix + 32.184);

  /// Convert to Barycentric Dynamical Time _(TDB)_.
  EpochTDB toTDB() {
    final tt = toTT();
    final tTT = tt.toJulianCenturies();
    final mEarth = (357.5277233 + 35999.05034 * tTT) * deg2rad;
    final seconds = 0.001658 * sin(mEarth) + 0.00001385 * sin(2 * mEarth);
    return EpochTDB(tt.posix + seconds);
  }

  /// Convert to a Global Positioning System _(GPS)_ formatted epoch.
  EpochGPS toGPS() {
    final ls = DataHandler().getLeapSeconds(toJulianDate());
    final delta = roll(ls - EpochGPS.offset).difference(EpochGPS.reference);
    final week = delta / secondsPerWeek;
    final weekFloor = week.floor();
    final seconds = (week - weekFloor) * secondsPerWeek;
    return EpochGPS(weekFloor, seconds);
  }

  /// Calculate the Greenwich Mean Sidereal Time _(GMST)_ angle _(rad)_.
  double gmstAngle() {
    final t = toJulianCenturies();
    final seconds = evalPoly(t, _gmstPoly);
    var result = ((seconds / 240) * deg2rad) % twoPi;
    if (result < 0) {
      result += twoPi;
    }
    return result;
  }

  /// Calculate the Greenwich Mean Sidereal Time _(GMST)_ angle _(°)_.
  double gmstAngleDegrees() => gmstAngle() * rad2deg;
}
