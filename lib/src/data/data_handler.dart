import 'package:pious_squid/pious_squid.dart';

/// Astrodynamic data management singleton.
class DataHandler {
  /// Create a new [DataHandler] object.
  factory DataHandler() => _instance;

  /// Private constructor for singleton.
  DataHandler._internal();

  /// [DataHandler] singleton instance.
  static final DataHandler _instance = DataHandler._internal();

  /// Return de-normalized [Egm96] coefficients for a given [l] and [m] index.
  Egm96Entry getEgm96Coeffs(final int l, final int m) =>
      egm96Data.getCoeffs(l, m);

  /// Return [Iau1980] nutation coefficients for a given [row].
  Iau1980Entry getIau1980Coeffs(final int row) => iau1980Data.getCoeffs(row);

  /// Return the leap second offset for a given Julian date [jd].
  double getLeapSeconds(final double jd) => leapSecondData.getLeapSeconds(jd);

  /// Return atmospheric parameters for a given [height] above Earth's
  /// surface _(km)_
  HpAtmosphereResult? getHpAtmosphere(final double height) =>
      hpAtmosphereData.getAtmosphere(height);

  /// Update Earth Orientation parameters from Celestrak CSV data.
  void updateEopFromCsv(final String csv) {
    eopData.updateFromCsv(csv);
  }

  /// Clear any loaded Earth Orentation Parameter data.
  void clearEop() {
    eopData.clearEntries();
  }

  /// Get Earth Orientation Parameter data for the provided epoch.
  ///
  /// An empty [EOP] object will be returned if there isn't data available
  /// for the provided epoch.
  EOP getEop(final EpochUTC epoch) => eopData.getEOP(epoch.toMjd());
}
