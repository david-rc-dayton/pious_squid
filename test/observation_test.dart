import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

final _state = J2000(
    EpochUTC.fromDateTimeString('1994-05-14T13:11:20.598Z'),
    Vector3D(5036.736529, -10806.660797, -4534.633784),
    Vector3D(2.6843855, -5.7595920, -2.4168093));

final _site = Geodetic.fromDegrees(39.007, -104.883, 2.19456)
    .toITRF(_state.epoch)
    .toJ2000();

void main() {
  group('Observation', () {
    test('RadecGeocentric', () {
      final radec = RadecGeocentric.fromStateVector(_state);
      expect(radec.rightAscensionDegrees, closeTo(294.9891458, 1e-3));
      expect(radec.declinationDegrees, closeTo(-20.8234944, 1e-3));
      expect(radec.range, closeTo(12756.000, 1e-3));
      expect(radec.rightAscensionRateDegrees, closeTo(-0.00000012244, 1e-6));
      expect(radec.declinationRateDegrees, closeTo(-0.00000001794, 1e-6));
      expect(radec.rangeRate, closeTo(6.7985140, 1e-6));
      expect(radec.position().distance(_state.position), lessThan(1e-3));
      expect(radec.velocity().distance(_state.velocity), lessThan(1e-6));
    });

    test('RadecTopocentric', () {
      final radec = RadecTopocentric.fromStateVectors(_state, _site);
      expect(radec.rightAscensionDegrees, closeTo(276.9337329, 1e-3));
      expect(radec.declinationDegrees, closeTo(-46.7583402, 1e-3));
      expect(radec.range, closeTo(11710.8120, 1e-3));
      expect(radec.rightAscensionRateDegrees, closeTo(0.01233970405, 1e-6));
      expect(radec.declinationRateDegrees, closeTo(0.01439246203, 1e-6));
      expect(radec.rangeRate, closeTo(6.0842826, 1e-6));
      expect(radec.position(_site).distance(_state.position), lessThan(1e-3));
      expect(radec.velocity(_site).distance(_state.velocity), lessThan(1e-6));
    });

    test('Razel', () {
      final razel = Razel.fromStateVectors(_state, _site);
      final state = razel.toStateVector(_site);
      expect(razel.azimuthDegrees, closeTo(210.8777747, 1e-3));
      expect(razel.elevationDegrees, closeTo(-5.9409535, 1e-3));
      expect(razel.range, closeTo(11710.8120, 1e-3));
      expect(razel.azimuthRateDegrees, closeTo(0.00384011466, 1e-6));
      expect(razel.elevationRateDegrees, closeTo(0.01495847759, 1e-6));
      expect(razel.rangeRate, closeTo(6.0842826, 1e-6));
      expect(state.position.distance(_state.position), lessThan(1e-3));
      expect(state.velocity.distance(_state.velocity), lessThan(1e-6));
    });
  });
}
