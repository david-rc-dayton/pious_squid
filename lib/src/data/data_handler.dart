import 'package:pious_squid/src/data/data_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

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
  HpAtmosphereResult? getHpAtmosphere(
      final EpochUTC epoch, final double height) {
    final mjd = epoch.toMjd();
    final sw = swData.getSpaceWeather(mjd);
    if (sw.isEmpty) {
      return hpAtmosphereData.getAtmosphere(height);
    }
    return hpAtmosphereFluxData.getAtmosphere(height, sw.f107AdjCtr81);
  }

  /// Update Earth Orientation parameters from Celestrak CSV data.
  void updateEarthOrientationParametersFromCsv(final String csv) {
    eopData.updateFromCsv(csv);
  }

  /// Clear any loaded Earth orentation parameter data.
  void clearEarthOrientationParameters() {
    eopData.clearEntries();
  }

  /// Get Earth orientation parameter data for the provided epoch.
  ///
  /// An empty [EOP] object will be returned if there isn't data available
  /// for the provided epoch.
  EOP getEop(final EpochUTC epoch) => eopData.getEOP(epoch.toMjd());

  /// Update space weather parameters from Celestrak CSV data.
  void updateSpaceWeatherFromCsv(final String csv) {
    swData.updateFromCsv(csv);
  }

  /// Clear any loaded space weather data.
  void clearSpaceWeather() {
    swData.clearEntries();
  }

  /// Get space weather data for the provided epoch.
  ///
  /// An empty [SpaceWeather] object will be returned if there isn't data
  /// available for the provided epoch.
  SpaceWeather getSpaceWeather(final EpochUTC epoch) =>
      swData.getSpaceWeather(epoch.toMjd());
}
