import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

final data = <(EpochUTC, double)>[
  (EpochUTC.fromDateTimeString('2012-01-26T05:54:00.000Z'), 6781.681),
  (EpochUTC.fromDateTimeString('2012-01-26T10:31:00.000Z'), 6781.804),
  (EpochUTC.fromDateTimeString('2012-01-26T13:36:00.000Z'), 6781.788),
  (EpochUTC.fromDateTimeString('2012-01-26T15:08:00.000Z'), 6781.788),
  (EpochUTC.fromDateTimeString('2012-01-26T19:45:00.000Z'), 6781.788),
  (EpochUTC.fromDateTimeString('2012-01-27T03:27:00.000Z'), 6781.955),
  (EpochUTC.fromDateTimeString('2012-01-27T08:04:00.000Z'), 6782.060),
  (EpochUTC.fromDateTimeString('2012-01-27T11:09:00.000Z'), 6782.052),
  (EpochUTC.fromDateTimeString('2012-01-27T14:13:00.000Z'), 6782.255),
  (EpochUTC.fromDateTimeString('2012-01-27T20:23:00.000Z'), 6782.366),
  (EpochUTC.fromDateTimeString('2012-01-28T01:00:00.000Z'), 6782.316),
  (EpochUTC.fromDateTimeString('2012-01-28T05:37:00.000Z'), 6782.412),
  (EpochUTC.fromDateTimeString('2012-01-28T10:14:00.000Z'), 6782.398),
  (EpochUTC.fromDateTimeString('2012-01-28T14:51:00.000Z'), 6782.450),
  (EpochUTC.fromDateTimeString('2012-01-28T19:28:00.000Z'), 6782.542),
  (EpochUTC.fromDateTimeString('2012-01-29T00:05:00.000Z'), 6782.638),
  (EpochUTC.fromDateTimeString('2012-01-29T01:37:00.000Z'), 6782.815),
  (EpochUTC.fromDateTimeString('2012-01-29T06:16:00.000Z'), 6783.875),
  (EpochUTC.fromDateTimeString('2012-01-29T12:25:00.000Z'), 6783.991),
  (EpochUTC.fromDateTimeString('2012-01-29T13:58:00.000Z'), 6784.142),
  (EpochUTC.fromDateTimeString('2012-01-29T17:02:00.000Z'), 6783.975),
  (EpochUTC.fromDateTimeString('2012-01-29T20:07:00.000Z'), 6784.019),
  (EpochUTC.fromDateTimeString('2012-01-30T03:49:00.000Z'), 6784.227),
  (EpochUTC.fromDateTimeString('2012-01-30T08:26:00.000Z'), 6783.863),
  (EpochUTC.fromDateTimeString('2012-01-30T13:03:00.000Z'), 6784.016),
  (EpochUTC.fromDateTimeString('2012-01-30T17:41:00.000Z'), 6784.324),
  (EpochUTC.fromDateTimeString('2012-01-30T22:18:00.000Z'), 6784.374),
  (EpochUTC.fromDateTimeString('2012-01-31T04:27:00.000Z'), 6784.442),
  (EpochUTC.fromDateTimeString('2012-01-31T09:04:00.000Z'), 6784.290),
  (EpochUTC.fromDateTimeString('2012-01-31T13:42:00.000Z'), 6784.551),
  (EpochUTC.fromDateTimeString('2012-01-31T18:19:00.000Z'), 6784.606),
  (EpochUTC.fromDateTimeString('2012-02-01T02:01:00.000Z'), 6784.684),
  (EpochUTC.fromDateTimeString('2012-02-01T05:05:00.000Z'), 6784.735),
];

final epochs = data.map((final e) => e.$1).toList();
final values = data.map((final e) => e.$2).toList();
final expectedIndices = [0, 12, 24, 32];

void main() {
  group('Smoothing', () {
    test('Exponential Simple', () {
      final smoothed = ExponentialSmoothing.smooth(values, 0.8);
      final expected = [6781.681, 6782.397, 6783.998, 6784.720];
      for (var i = 0; i < expectedIndices.length; i++) {
        expect(smoothed[expectedIndices[i]], closeTo(expected[i], 1e-2));
      }
    });
  });

  test('Exponential Double', () {
    final smoothed = ExponentialSmoothing.smoothDouble(values, 0.9, 0.4);
    final expected = [6781.681, 6782.397, 6783.998, 6784.720];
    for (var i = 0; i < expectedIndices.length; i++) {
      expect(smoothed[expectedIndices[i]], closeTo(expected[i], 1e-2));
    }
  });

  test('Exponential Time', () {
    final smoothed = ExponentialSmoothing.smoothTime(epochs, values, 86400.0);
    final expected = [6781.681, 6782.397, 6783.998, 6784.720];
    for (var i = 0; i < expectedIndices.length; i++) {
      expect(smoothed[expectedIndices[i]], closeTo(expected[i], 1e-2));
    }
  });
}
