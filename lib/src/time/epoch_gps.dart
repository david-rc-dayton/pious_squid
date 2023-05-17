import 'dart:math';

import 'package:pious_squid/src/data/data_handler.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Global Positioning System _(GPS)_ formatted epoch.
class EpochGPS {
  /// Create a new GPS epoch given the [week] since reference epoch, and number
  /// of [seconds] into the [week].
  EpochGPS(this.week, this.seconds);

  /// Number of weeks since the GPS reference epoch.
  final int week;

  /// Number of seconds into the week.
  final double seconds;

  /// GPS reference epoch.
  static final EpochUTC reference =
      EpochUTC.fromDateTimeString('1980-01-06T00:00:00.000Z');

  /// GPS leap second difference from TAI/UTC offsets.
  static final double offset = 19;

  /// Get GPS week accounting for 10-bit rollover.
  int get week10Bit => week % pow(2, 10).toInt();

  /// Get GPS week accounting for 13-bit rollover.
  int get week13Bit => week % pow(2, 13).toInt();

  @override
  String toString() => '$week:${seconds.toStringAsFixed(3)}';

  /// Convert this to a UTC epoch.
  EpochUTC toUTC() {
    final init = reference.roll(week * secondsPerWeek + seconds);
    final ls = DataHandler().getLeapSeconds(init.toJulianDate());
    return init.roll(-(ls - offset));
  }
}
