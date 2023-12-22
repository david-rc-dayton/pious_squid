import 'dart:collection';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';

final _eopCsvMatch = RegExp('^[0-9]{4}-[0-9]{2}-[0-9]{2},');

/// Earth Orientation Parameter data.
class EOP {
  /// Create a new [EOP] object, given a Modified Julian Date, polar motion x/y
  /// axis _(rad)_, UT1-UTC _(seconds)_, Length of Day _(seconds)_, delta-psi
  /// _(rad)_, and delta-epsilon _(rad)_.
  EOP(this.mjd, this.x, this.y, this.dut1, this.lod, this.dpsi, this.deps,
      {this.isEmpty = false});

  /// Create a new [EOP] object from a CelesTrack EOP CSV [line].
  factory EOP.fromCsvLine(final String line) {
    final fields = line.split(',');
    final mjd = double.parse(fields[1]);
    final x = double.parse(fields[2]) * asec2rad;
    final y = double.parse(fields[3]) * asec2rad;
    final dut1 = double.parse(fields[4]);
    final lod = double.parse(fields[5]);
    final dpsi = double.parse(fields[6]) * asec2rad;
    final deps = double.parse(fields[7]) * asec2rad;
    return EOP(mjd, x, y, dut1, lod, dpsi, deps);
  }

  /// Create an empty [EOP] object with the provided Modified Julian Date.
  factory EOP.empty(final double mjd) =>
      EOP(mjd, 0, 0, 0, 0, 0, 0, isEmpty: true);

  /// Modified Julian Date
  final double mjd;

  /// Polar motion x-axis _(rad)_
  final double x;

  /// Polar motion y-axis _(rad)_
  final double y;

  /// UT1-UTC _(seconds)_
  final double dut1;

  /// Length of Day _(seconds)_
  final double lod;

  /// Delta-Psi _(rad)_
  final double dpsi;

  /// Delta-Epsilon _(rad)_
  final double deps;

  /// `true` if this EOP entry is not based on actual data
  final bool isEmpty;

  @override
  String toString() => [
        'MJD: ${mjd.toInt()}',
        'X: ${(x * rad2asec).toStringAsFixed(6)}"',
        'Y: ${(y * rad2asec).toStringAsFixed(6)}"',
        'DUT1: ${dut1.toStringAsFixed(7)}s',
        'LOD: ${lod.toStringAsFixed(7)}s',
        'DPSI: ${(dpsi * rad2asec).toStringAsFixed(6)}"',
        'DEPS: ${(deps * rad2asec).toStringAsFixed(6)}"'
      ].join('  ');

  /// Linearly interpolate between two [EOD] objects, for a given Modified
  /// Julian Date.
  EOP interpolate(final EOP eop, final double tMjd) => EOP(
      tMjd,
      linearInterpolate(tMjd, mjd, x, eop.mjd, eop.x),
      linearInterpolate(tMjd, mjd, y, eop.mjd, eop.y),
      linearInterpolate(tMjd, mjd, dut1, eop.mjd, eop.dut1),
      linearInterpolate(tMjd, mjd, lod, eop.mjd, eop.lod),
      linearInterpolate(tMjd, mjd, dpsi, eop.mjd, eop.dpsi),
      linearInterpolate(tMjd, mjd, deps, eop.mjd, eop.deps));
}

/// Earth Orientation Parameter data container.
class EOPData {
  /// Create a new [EOPData] object.
  EOPData();

  /// [EOP] data entries.
  final List<EOP> _entries = [];

  /// Latest accessed [EOP] cache, indexed by Modified Julian Date.
  final HashMap<int, EOP> _cache = HashMap();

  /// Modified Julian Date of the earliest stored entry.
  double _mjdMin = 0;

  /// Modified Julian Date of the latest stored entry.
  double _mjdMax = 0;

  /// [EOP] cache limit.
  static const int _cacheLimit = 30;

  /// Clear any stored [EOP] entries.
  void clearEntries() {
    _entries.clear();
    _cache.clear();
    _mjdMin = 0;
    _mjdMax = 0;
  }

  /// Update [EOP] data from a CelesTrak Earth Orientation Parameter CSV string.
  void updateFromCsv(final String csv) {
    final lines = csv.split('\n');
    lines.removeWhere((final line) => !_eopCsvMatch.hasMatch(line));
    clearEntries();
    for (final line in lines) {
      _entries.add(EOP.fromCsvLine(line));
    }
    _mjdMin = _entries.first.mjd;
    _mjdMax = _entries.last.mjd;
  }

  /// Lookup and cache [EOP] data for the provided Modified Julian Date.
  EOP _lookupEop(final int mjd) {
    for (var i = 0; i < _entries.length; i++) {
      final f = _entries[i];
      if (mjd == f.mjd) {
        _cache[mjd] = f;
        return f;
      }
    }
    final fEmpty = EOP.empty(mjd.toDouble());
    _cache[mjd] = fEmpty;
    return fEmpty;
  }

  /// Find and interpolate [EOP] data for the provided Modified Julian Date.
  EOP _matchEop(final double mjd) {
    if (_cache.length > _cacheLimit) {
      _cache.clear();
    }
    final m0 = mjd.toInt();
    final m1 = m0 + 1;
    final f0 = _cache[m0] ?? _lookupEop(m0);
    final f1 = _cache[m1] ?? _lookupEop(m1);
    return f0.interpolate(f1, mjd);
  }

  /// Get interpolated [EOP] data for the provided Modified Julian Date.
  ///
  /// An empty [EOP] object will be returned if there isn't data available
  /// for the provided epoch.
  EOP getEOP(final double mjd) {
    if (_entries.isEmpty || mjd < _mjdMin || mjd > _mjdMax - 1) {
      return EOP.empty(mjd);
    }
    return _matchEop(mjd);
  }
}

/// Global Earth Orientation Parameter data container.
final eopData = EOPData();
