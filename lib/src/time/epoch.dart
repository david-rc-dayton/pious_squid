import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Year, month, day, hour, minute, and second.
typedef GregorianFields = (int, int, int, int, int, double);

/// Time stamped value container.
class TimeStamped<T> {
  /// Create a new time stamped [value] container at the provided [epoch].
  TimeStamped(this.epoch, this.value);

  /// Timestamp epoch.
  final EpochUTC epoch;

  /// Timestamped value.
  final T value;
}

/// Base class for [Epoch] data.
class Epoch implements Comparable {
  /// Create a new [Epoch] object given the number of seconds elapsed since the
  /// [posix] epoch _(`1970-01-01T00:00:00.000`)_ in the [Epoch] time scale.
  Epoch(this.posix);

  /// Number of seconds since the POSIX epoch _(`1970-01-01T00:00:00.000`)_ in
  /// this object's time scale.
  final double posix;

  @override
  String toString() => toDateTime().toIso8601String();

  /// Convert this to an Excel spreadsheet string.
  String toExcelString() => toString().substring(0, 19);

  /// Return the difference _(s)_ between this and another [epoch]/
  double difference(final Epoch epoch) => posix - epoch.posix;

  /// Check if this has the same timestamp as the provided [epoch].
  bool equals(final Epoch epoch) => posix == epoch.posix;

  /// Convert to a [DateTime] object.
  DateTime toDateTime() =>
      DateTime.fromMillisecondsSinceEpoch((posix * 1000).toInt(), isUtc: true);

  /// Convert to Julian date.
  double toJulianDate() => (posix / secondsPerDay) + 2440587.5;

  /// Convert to Julian centuries.
  double toJulianCenturies() => (toJulianDate() - 2451545) / 36525;

  /// Convert to Gregorian date fields for year, month, day, hour, minute, and
  /// second.
  GregorianFields toGregorianFields() {
    final jdfull = toJulianDate();
    final jdfrac = jdfull - jdfull.floor();
    final jd = jdfull.floor();
    var temp = jd - 2415019.5;
    final tu = temp / 365.25;
    var year = 1900 + tu.floor();
    var leapyrs = ((year - 1901) * 0.25).floor();
    var days = (temp - ((year - 1900) * 365 + leapyrs)).floor();
    if (days + jdfrac < 1) {
      year = year - 1;
      leapyrs = ((year - 1901) * 0.25).floor();
      days = (temp - ((year - 1900) * 365 + leapyrs)).floor();
    }
    final lmonth = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    if (year % 4 == 0) {
      lmonth[2] = 29;
    }
    var i = 1;
    var inttemp = 0;
    while ((days > inttemp + lmonth[i]) && (i < 12)) {
      inttemp += lmonth[i];
      i++;
    }
    final mon = i;
    final day = days - inttemp;
    temp = ((days + jdfrac) - days) * 24.0;
    final hr = temp.floor();
    temp = (temp - hr) * 60.0;
    final minute = temp.floor();
    final sec = (temp - minute) * 60.0;
    return (year, mon, day, hr, minute, sec);
  }

  /// Return the year of this [Epoch] object.
  int year() {
    final jdfull = toJulianDate();
    final jdfrac = jdfull - jdfull.floor();
    final jd = jdfull.floor();
    final temp = jd - 2415019.5;
    final tu = temp / 365.25;
    var year = 1900 + tu.floor();
    final leapyrs = ((year - 1901) * 0.25).floor();
    final days = (temp - ((year - 1900) * 365 + leapyrs)).floor();
    if (days + jdfrac < 1) {
      year = year - 1;
    }
    return year;
  }

  @override
  int compareTo(covariant final Epoch other) => posix.compareTo(other.posix);

  /// Check if this is later than the [other] epoch.
  bool operator >(final Epoch other) => posix > other.posix;

  /// Check if this is later or the same as the [other] epoch.
  bool operator >=(final Epoch other) => posix >= other.posix;

  /// Check if this is earlier than the [other] epoch.
  bool operator <(final Epoch other) => posix < other.posix;

  /// Check if this is earlier or the same as the [other] epoch.
  bool operator <=(final Epoch other) => posix <= other.posix;
}
