import 'package:pious_squid/src/operations/constants.dart';

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
