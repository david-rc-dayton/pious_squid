/// Harris-Priester atmosphere entry for height, min, and max density values.
typedef HpAtmosphereEntry = ({double h, double minD, double maxD});

/// Harris-Priester atmosphere data, assuming mean solar flux.
final List<HpAtmosphereEntry> hpAtmosphere = [
  (h: 100, minD: 4.974e-7, maxD: 4.974e-7),
  (h: 120, minD: 2.49e-8, maxD: 2.49e-8),
  (h: 130, minD: 8.377e-9, maxD: 8.71e-9),
  (h: 140, minD: 3.899e-9, maxD: 4.059e-9),
  (h: 150, minD: 2.122e-9, maxD: 2.215e-9),
  (h: 160, minD: 1.263e-9, maxD: 1.344e-9),
  (h: 170, minD: 8.008e-10, maxD: 8.758e-10),
  (h: 180, minD: 5.283e-10, maxD: 6.01e-10),
  (h: 190, minD: 3.617e-10, maxD: 4.297e-10),
  (h: 200, minD: 2.557e-10, maxD: 3.162e-10),
  (h: 210, minD: 1.839e-10, maxD: 2.396e-10),
  (h: 220, minD: 1.341e-10, maxD: 1.853e-10),
  (h: 230, minD: 9.949e-11, maxD: 1.455e-10),
  (h: 240, minD: 7.488e-11, maxD: 1.157e-10),
  (h: 250, minD: 5.709e-11, maxD: 9.308e-11),
  (h: 260, minD: 4.403e-11, maxD: 7.555e-11),
  (h: 270, minD: 3.43e-11, maxD: 6.182e-11),
  (h: 280, minD: 2.697e-11, maxD: 5.095e-11),
  (h: 290, minD: 2.139e-11, maxD: 4.226e-11),
  (h: 300, minD: 1.708e-11, maxD: 3.526e-11),
  (h: 320, minD: 1.099e-11, maxD: 2.511e-11),
  (h: 340, minD: 7.214e-12, maxD: 1.819e-11),
  (h: 360, minD: 4.824e-12, maxD: 1.337e-11),
  (h: 380, minD: 3.274e-12, maxD: 9.955e-12),
  (h: 400, minD: 2.249e-12, maxD: 7.492e-12),
  (h: 420, minD: 1.558e-12, maxD: 5.684e-12),
  (h: 440, minD: 1.091e-12, maxD: 4.355e-12),
  (h: 460, minD: 7.701e-13, maxD: 3.362e-12),
  (h: 480, minD: 5.474e-13, maxD: 2.612e-12),
  (h: 500, minD: 3.916e-13, maxD: 2.042e-12),
  (h: 520, minD: 2.819e-13, maxD: 1.605e-12),
  (h: 540, minD: 2.042e-13, maxD: 1.267e-12),
  (h: 560, minD: 1.488e-13, maxD: 1.005e-12),
  (h: 580, minD: 1.092e-13, maxD: 7.997e-13),
  (h: 600, minD: 8.07e-14, maxD: 6.39e-13),
  (h: 620, minD: 6.012e-14, maxD: 5.123e-13),
  (h: 640, minD: 4.519e-14, maxD: 4.121e-13),
  (h: 660, minD: 3.43e-14, maxD: 3.325e-13),
  (h: 680, minD: 2.632e-14, maxD: 2.691e-13),
  (h: 700, minD: 2.043e-14, maxD: 2.185e-13),
  (h: 720, minD: 1.607e-14, maxD: 1.779e-13),
  (h: 740, minD: 1.281e-14, maxD: 1.452e-13),
  (h: 760, minD: 1.036e-14, maxD: 1.19e-13),
  (h: 780, minD: 8.496e-15, maxD: 9.776e-14),
  (h: 800, minD: 7.069e-15, maxD: 8.059e-14),
  (h: 840, minD: 4.68e-15, maxD: 5.741e-14),
  (h: 880, minD: 3.2e-15, maxD: 4.21e-14),
  (h: 920, minD: 2.21e-15, maxD: 3.13e-14),
  (h: 960, minD: 1.56e-15, maxD: 2.36e-14),
  (h: 1000, minD: 1.15e-15, maxD: 1.81e-14),
];

/// Harris-Priester atmospheric density bracket.
class HpAtmosphereResult {
  /// Create a new [HpAtmosphereResult] object.
  HpAtmosphereResult(this.height, this.hp0, this.hp1);

  /// Height above Earth's surface _(km)_.
  final double height;

  /// Lower bound for atmospheric parameters.
  final HpAtmosphereEntry hp0;

  /// Upper bound for atmospheric parameters.
  final HpAtmosphereEntry hp1;
}

/// Container for Harris-Priester atmosphere data.
class HpAtmosphereData {
  /// Create a new [HpAtmosphereData] object given an array of
  /// [HpAtmosphere] objects.
  HpAtmosphereData(this._table)
      : _hMin = _table.first.h,
        _hMax = _table.last.h;

  /// Array of Harris-Priester atmospheric parameters
  final List<HpAtmosphereEntry> _table;

  /// Minimum atmosphere height _(km)_
  final double _hMin;

  /// Maximum atmosphere height _(km)_
  final double _hMax;

  /// Return atmospheric parameters for a given [height] above Earth's
  /// surface _(km)_.
  HpAtmosphereResult? getAtmosphere(final double height) {
    if (height < _hMin || height > _hMax) {
      return null;
    }
    var index = 0;
    while (index < _table.length - 2 && height > _table[index + 1].h) {
      index++;
    }
    return HpAtmosphereResult(height, _table[index], _table[index + 1]);
  }
}

/// Harris-Priester atmosphere data container.
final hpAtmosphereData = HpAtmosphereData(hpAtmosphere);
