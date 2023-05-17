/// Leap second value tuples.
final List<(double, double)> leapSeconds = [
  (2441317.5, 10),
  (2441499.5, 11),
  (2441683.5, 12),
  (2442048.5, 13),
  (2442413.5, 14),
  (2442778.5, 15),
  (2443144.5, 16),
  (2443509.5, 17),
  (2443874.5, 18),
  (2444239.5, 19),
  (2444786.5, 20),
  (2445151.5, 21),
  (2445516.5, 22),
  (2446247.5, 23),
  (2447161.5, 24),
  (2447892.5, 25),
  (2448257.5, 26),
  (2448804.5, 27),
  (2449169.5, 28),
  (2449534.5, 29),
  (2450083.5, 30),
  (2450630.5, 31),
  (2451179.5, 32),
  (2453736.5, 33),
  (2454832.5, 34),
  (2456109.5, 35),
  (2457204.5, 36),
  (2457754.5, 37),
];

/// Leap second data.
class LeapSecond {
  /// Create a new [LeapSecond] object.
  LeapSecond(this.jd, this.offset);

  /// Julian date.
  final double jd;

  /// Offset seconds.
  final double offset;
}

/// Leap second data container.
class LeapSecondData {
  /// Create a new [LeapSecondData] container object given an array
  /// of [offsets].
  LeapSecondData(this._offsets) {
    final o0 = _offsets.first;
    final oN = _offsets.last;
    _jdFirst = o0.jd;
    _jdLast = oN.jd;
    _offsetFirst = o0.offset;
    _offsetLast = oN.offset;
  }

  /// Create a new [LeapSecondData] object from an array of leap second
  /// value tuples [vals].
  factory LeapSecondData.fromVals(final List<(double, double)> vals) {
    final output = <LeapSecond>[];
    for (final v in vals) {
      final (jd, offset) = v;
      output.add(LeapSecond(jd, offset));
    }
    return LeapSecondData(output);
  }

  /// Leap second offsets.
  final List<LeapSecond> _offsets;

  /// First Julian date.
  late final double _jdFirst;

  /// Last Julian date.
  late final double _jdLast;

  /// First offset seconds.
  late final double _offsetFirst;

  /// Last offset seconds.
  late final double _offsetLast;

  /// Return the number of leap seconds for a given Julian date [jd].
  double getLeapSeconds(final double jd) {
    if (jd >= _jdLast) {
      return _offsetLast;
    }
    if (jd <= _jdFirst) {
      return _offsetFirst;
    }
    for (var i = 0; i < _offsets.length - 2; i++) {
      if (jd >= _offsets[i].jd && jd < _offsets[i + 1].jd) {
        return _offsets[i].offset;
      }
    }
    return 0;
  }
}

/// Leap second data container.
final leapSecondData = LeapSecondData.fromVals(leapSeconds);
