import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/sgp4/elsetrec.dart';
import 'package:pious_squid/src/sgp4/gravconst.dart';
import 'package:pious_squid/src/sgp4/sgp4.dart';
import 'package:pious_squid/src/sgp4/twoline2rv.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// SGP4 Operation Mode.
enum Sgp4OpsMode {
  /// Air Force Space Command
  afspc('a'),

  /// Improved
  improved('i');

  const Sgp4OpsMode(this.value);

  /// Return the string value of this enumeration.
  final String value;
}

/// Two-line element set.
class TLE {
  /// Create a new [TLE] object from [line1] and [line2] of the element set.
  TLE(this.line1, this.line2,
      {final Sgp4OpsMode opsMode = Sgp4OpsMode.afspc,
      final Sgp4GravConst gravConst = Sgp4GravConst.wgs72})
      : epoch = _parseEpoch(line1.substring(18, 32)),
        satnum = _parseSatnum(line1.substring(2, 7)) {
    twoline2rv(line1, line2, opsMode.value, gravConst, _satrec);
  }

  /// Line 1 of the element set.
  final String line1;

  /// Line 2 of the element set.
  final String line2;

  /// Element set epoch.
  final EpochUTC epoch;

  /// NORAD satellite number.
  final int satnum;
  final ElsetRec _satrec = ElsetRec();
  static final Map<String, String> _alpha5 = {
    'A': '10',
    'B': '11',
    'C': '12',
    'D': '13',
    'E': '14',
    'F': '15',
    'G': '16',
    'H': '17',
    'J': '18',
    'K': '19',
    'L': '20',
    'M': '21',
    'N': '22',
    'P': '23',
    'Q': '24',
    'R': '25',
    'S': '26',
    'T': '27',
    'U': '28',
    'V': '29',
    'W': '30',
    'X': '31',
    'Y': '32',
    'Z': '33'
  };

  @override
  String toString() => '$line1\n$line2';

  static EpochUTC _parseEpoch(final String epochStr) {
    var year = int.parse(epochStr.substring(0, 2));
    if (year >= 57) {
      year += 1900;
    } else {
      year += 2000;
    }
    final days = double.parse(epochStr.substring(2, 14)) - 1;
    return EpochUTC.fromDateTimeString('$year-01-01T00:00:00.000Z')
        .roll(days * secondsPerDay);
  }

  static int _parseSatnum(final String satnumStr) {
    final values = satnumStr.toUpperCase().split('');
    for (var i = 0; i < values.length; i++) {
      final value = values[i];
      if (_alpha5.containsKey(value)) {
        values[i] = _alpha5[value]!;
      }
    }
    return int.parse(values.join(''));
  }

  /// Propagate this [TLE] to the provided [epoch].
  TEME propagate(final EpochUTC epoch) {
    final r = Float64List(3);
    final v = Float64List(3);
    sgp4(_satrec, epoch.difference(this.epoch) / 60.0, r, v);
    return TEME(epoch, Vector3D(r[0], r[1], r[2]), Vector3D(v[0], v[1], v[2]));
  }

  TEME _currentState() {
    final r = Float64List(3);
    final v = Float64List(3);
    sgp4(_satrec, 0.0, r, v);
    return TEME(epoch, Vector3D(r[0], r[1], r[2]), Vector3D(v[0], v[1], v[2]));
  }

  /// Return the state of this [TLE] at [epoch].
  TEME get state => _currentState();
}
