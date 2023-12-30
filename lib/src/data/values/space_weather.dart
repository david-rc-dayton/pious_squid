import 'dart:collection';

import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/time/time_base.dart';

final _swCsvMatch = RegExp('^[0-9]{4}-[0-9]{2}-[0-9]{2},');

/// Container for space weather data.
class SpaceWeather {
  /// Store space weather data on a given Modified Julian Date.
  SpaceWeather(this.mjd, this.f107Obs, this.f107Adj, this.f107ObsCtr81,
      this.f107ObsLst81, this.f107AdjCtr81, this.f107AdjLst81,
      {this.isEmpty = false});

  /// Create a new [SpaceWeather] object from a CelesTrack SW CSV [line].
  factory SpaceWeather.fromCsvLine(final String line) {
    final fields = line.split(',');
    final epochFields = fields[0].split('-');
    final mjd = EpochUTC.fromDate(int.parse(epochFields[0]),
            int.parse(epochFields[1]), int.parse(epochFields[2]))
        .toMjd();
    final f107Obs = double.parse(fields[24]);
    final f107Adj = double.parse(fields[25]);
    final f107ObsCtr81 = double.parse(fields[27]);
    final f107ObsLst81 = double.parse(fields[28]);
    final f107AdjCtr81 = double.parse(fields[29]);
    final f107AdjLst81 = double.parse(fields[30]);
    return SpaceWeather(mjd, f107Obs, f107Adj, f107ObsCtr81, f107ObsLst81,
        f107AdjCtr81, f107AdjLst81);
  }

  /// Create an empty [SpaceWeather] object with the provided Modified
  /// Julian Date.
  factory SpaceWeather.empty(final double mjd) =>
      SpaceWeather(mjd, 0, 0, 0, 0, 0, 0, isEmpty: true);

  /// Modified Julian Date
  final double mjd;

  /// Observed 10.7-cm Solar Radio Flux _(F10.7)_.
  final double f107Obs;

  /// 10.7-cm Solar Radio Flux _(F10.7)_ adjusted to 1 AU.
  final double f107Adj;

  /// Centered 81-day arithmetic average of F10.7 _(observed)_.
  final double f107ObsCtr81;

  /// Last 81-day arithmetic average of F10.7 _(observed)_.
  final double f107ObsLst81;

  /// Centered 81-day arithmetic average of F10.7 _(adjusted)_.
  final double f107AdjCtr81;

  /// Last 81-day arithmetic average of F10.7 _(adjusted)_.
  final double f107AdjLst81;

  /// `true` if this space weather entry is not based on actual data
  final bool isEmpty;

  /// Linearly interpolate between two [SpaceWeather] objects, for a given
  /// Modified Julian Date.
  SpaceWeather interpolate(final SpaceWeather sw, final double tMjd) =>
      SpaceWeather(
          tMjd,
          linearInterpolate(tMjd, mjd, f107Obs, sw.mjd, sw.f107Obs),
          linearInterpolate(tMjd, mjd, f107Adj, sw.mjd, sw.f107Adj),
          linearInterpolate(tMjd, mjd, f107ObsCtr81, sw.mjd, sw.f107ObsCtr81),
          linearInterpolate(tMjd, mjd, f107ObsLst81, sw.mjd, sw.f107ObsLst81),
          linearInterpolate(tMjd, mjd, f107AdjCtr81, sw.mjd, sw.f107AdjCtr81),
          linearInterpolate(tMjd, mjd, f107AdjLst81, sw.mjd, sw.f107AdjLst81));
}

/// Container for a set of [SpaceWeather] entries;
class SpaceWeatherData {
  /// Create a new [SpaceWeatherData] object.
  SpaceWeatherData();

  final List<SpaceWeather> _entries = [];

  /// Latest accessed [SpaceWeather] cache, indexed by Modified Julian Date.
  final HashMap<int, SpaceWeather> _cache = HashMap();

  /// Modified Julian Date of the earliest stored entry.
  double _mjdMin = 0;

  /// Modified Julian Date of the latest stored entry.
  double _mjdMax = 0;

  /// [SpaceWeather] cache limit.
  static const int _cacheLimit = 30;

  /// Clear any stored [SpaceWeather] entries.
  void clearEntries() {
    _entries.clear();
    _cache.clear();
    _mjdMin = 0;
    _mjdMax = 0;
  }

  /// Update [SpaceWeather] data from a CelesTrak space weather CSV string.
  void updateFromCsv(final String csv) {
    clearEntries();
    final lines = csv.split('\n');
    lines.removeWhere((final line) => !_swCsvMatch.hasMatch(line));
    for (final line in lines) {
      _entries.add(SpaceWeather.fromCsvLine(line));
    }
    _mjdMin = _entries.first.mjd;
    _mjdMax = _entries.last.mjd;
  }

  /// Lookup and cache [SpaceWeather] data for the provided Modified
  /// Julian Date.
  SpaceWeather _lookupSw(final int mjd) {
    for (var i = 0; i < _entries.length; i++) {
      final f = _entries[i];
      if (mjd == f.mjd) {
        _cache[mjd] = f;
        return f;
      }
    }
    final fEmpty = SpaceWeather.empty(mjd.toDouble());
    _cache[mjd] = fEmpty;
    return fEmpty;
  }

  /// Find and interpolate [SpaceWeather] data for the provided Modified
  /// Julian Date.
  SpaceWeather _matchSw(final double mjd) {
    if (_cache.length > _cacheLimit) {
      _cache.clear();
    }
    final m0 = mjd.toInt();
    final m1 = m0 + 1;
    final f0 = _cache[m0] ?? _lookupSw(m0);
    final f1 = _cache[m1] ?? _lookupSw(m1);
    return f0.interpolate(f1, mjd);
  }

  /// Get interpolated [SpaceWeather] data for the provided Modified Julian Date.
  ///
  /// An empty [SpaceWeather] object will be returned if there isn't data
  /// available for the provided epoch.
  SpaceWeather getSpaceWeather(final double mjd) {
    if (_entries.isEmpty || mjd < _mjdMin || mjd > _mjdMax - 1) {
      return SpaceWeather.empty(mjd);
    }
    return _matchSw(mjd);
  }
}

/// Global space weather data container.
final swData = SpaceWeatherData();
