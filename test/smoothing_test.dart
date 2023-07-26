import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

final data = <(EpochUTC, double)>[
  (EpochUTC(17.0), 94),
  (EpochUTC(13.1), 73),
  (EpochUTC(12.2), 59),
  (EpochUTC(15.3), 80),
  (EpochUTC(16.4), 93),
  (EpochUTC(14.5), 85),
  (EpochUTC(16.6), 66),
  (EpochUTC(16.7), 79),
  (EpochUTC(18.8), 77),
  (EpochUTC(19.9), 81),
]..sort((final a, final b) => a.$1.compareTo(b.$1));

final epochs = data.map((final e) => e.$1).toList();
final values = data.map((final e) => e.$2).toList();

void main() {
  group('Smoothing', () {
    test('Exponential Simple', () {
      final smoothed = ExponentialSmoothing.smooth(values, 0.2);
      final expected = [
        59.000,
        61.800,
        66.440,
        69.152,
        73.922,
        72.337,
        73.670,
        77.736,
        77.589,
        78.271
      ];
      for (var i = 0; i < smoothed.length; i++) {
        expect(smoothed[i], closeTo(expected[i], 1e-2));
      }
    });

    test('Exponential Double', () {
      final smoothed = ExponentialSmoothing.smoothDouble(values, 0.2, 0.9);
      final expected = [
        59.000,
        73.000,
        86.600,
        96.192,
        103.551,
        102.139,
        97.104,
        92.818,
        86.202,
        80.052
      ];
      for (var i = 0; i < smoothed.length; i++) {
        expect(smoothed[i], closeTo(expected[i], 1e-2));
      }
    });

    test('Exponential Time', () {
      final smoothed = ExponentialSmoothing.smoothTime(epochs, values, 2.0);
      final expected = [
        59.000,
        64.073,
        74.608,
        76.386,
        83.414,
        81.757,
        81.623,
        83.347,
        79.580,
        80.181
      ];
      for (var i = 0; i < smoothed.length; i++) {
        expect(smoothed[i], closeTo(expected[i], 1e-2));
      }
    });
  });
}
