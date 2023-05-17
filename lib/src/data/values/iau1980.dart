/// Array of the first 4 IAU-1980 coefficients.
final List<(int, int, int, int, int, double, double, double, double)> iau1980 =
    [
  (0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9),
  (0, 0, 2, -2, 2, -13187, -1.6, 5736, -3.1),
  (0, 0, 2, 0, 2, -2274, -0.2, 977, -0.5),
  (0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5),
];

/// IAU-1980 nutation data
class Iau1980 {
  /// Create a new [Iau1980] object from a tuple of coefficients.
  Iau1980(this.a1, this.a2, this.a3, this.a4, this.a5, this.ai, this.bi,
      this.ci, this.di);

  /// a1 coefficient
  final int a1;

  /// a2 coefficient
  final int a2;

  /// a3 coefficient
  final int a3;

  /// a4 coefficient
  final int a4;

  /// a5 coefficient
  final int a5;

  /// Ai coefficient
  final double ai;

  /// Bi coefficient
  final double bi;

  /// Ci coefficient
  final double ci;

  /// Di coefficient
  final double di;
}

/// Container for IAU-1980 data.
class Iau1980Data {
  /// Create a new [Iau1980Data] object.
  Iau1980Data(this._coeffs);

  /// Create a new [Iau1980Data] container object from an array of IAU-1980
  /// coefficient tuples [coeffs].
  factory Iau1980Data.fromCoeffs(
      final List<(int, int, int, int, int, double, double, double, double)>
          coeffs) {
    final output = <Iau1980>[];
    for (final c in coeffs) {
      final (a1, a2, a3, a4, a5, ai, bi, ci, di) = c;
      output.add(Iau1980(a1, a2, a3, a4, a5, ai, bi, ci, di));
    }
    return Iau1980Data(output);
  }

  /// IAU-1980 coefficients.
  final List<Iau1980> _coeffs;

  /// Get IAU-1980 coefficients for a given row number.
  Iau1980 getCoeffs(final int row) => _coeffs[row];
}

/// IAU-1980 data container.
final iau1980Data = Iau1980Data.fromCoeffs(iau1980);
