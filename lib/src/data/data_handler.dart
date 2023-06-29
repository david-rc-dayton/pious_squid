import 'package:pious_squid/src/data/values/egm96.dart';
import 'package:pious_squid/src/data/values/hp_atmosphere.dart';
import 'package:pious_squid/src/data/values/iau1980.dart';
import 'package:pious_squid/src/data/values/leap_second.dart';

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
}
