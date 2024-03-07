import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

final _state = J2000(
    EpochUTC.fromDate(2000, 2, 15, 14, 47, 39.570),
    Vector3D(-22864.644660, 35438.541539, -0.361127),
    Vector3D(-2.582939062, -1.666767504, -0.000795483));

final _covariance = StateCovariance.fromSigmas(
    [0.3943, 1.7769, 1.0018, 0.0001, 0.0000, 0.0002], CovarianceFrame.ric);

final _forceModel = ForceModel()
  ..setEarthGravity(8, 8)
  ..setThirdBodyGravity(sun: true, moon: true);

const _tolerance = 1e-4;

void main() {
  group('Covariance', () {
    test('CovarianceSample', () {
      final sample =
          CovarianceSample(_state, _covariance, originForceModel: _forceModel);
      var ricSig = sample.desampleRIC().sigmas();
      var j2kSig = sample.desampleJ2000().sigmas();

      expect(ricSig[0], closeTo(0.3943, _tolerance));
      expect(ricSig[1], closeTo(1.7769, _tolerance));
      expect(ricSig[2], closeTo(1.0018, _tolerance));
      expect(ricSig[3], closeTo(0.0002, _tolerance));
      expect(ricSig[4], closeTo(0.0001, _tolerance));
      expect(ricSig[5], closeTo(0.0002, _tolerance));

      expect(j2kSig[0], closeTo(1.5083, _tolerance));
      expect(j2kSig[1], closeTo(1.0187, _tolerance));
      expect(j2kSig[2], closeTo(1.0018, _tolerance));
      expect(j2kSig[3], closeTo(0.0001, _tolerance));
      expect(j2kSig[4], closeTo(0.0001, _tolerance));
      expect(j2kSig[5], closeTo(0.0002, _tolerance));

      sample.propagate(EpochUTC.fromDate(2000, 2, 22, 11, 47, 39.570));
      ricSig = sample.desampleRIC().sigmas();
      j2kSig = sample.desampleJ2000().sigmas();

      expect(ricSig[0], closeTo(2.2787, _tolerance));
      expect(ricSig[1], closeTo(199.6084, _tolerance));
      expect(ricSig[2], closeTo(1.8802, _tolerance));
      expect(ricSig[3], closeTo(0.0003, _tolerance));
      expect(ricSig[4], closeTo(0.0003, _tolerance));
      expect(ricSig[5], closeTo(0.0002, _tolerance));

      expect(j2kSig[0], closeTo(198.8608, _tolerance));
      expect(j2kSig[1], closeTo(17.4096, _tolerance));
      expect(j2kSig[2], closeTo(1.8801, _tolerance));
      expect(j2kSig[3], closeTo(0.0015, _tolerance));
      expect(j2kSig[4], closeTo(0.0142, _tolerance));
      expect(j2kSig[5], closeTo(0.0002, _tolerance));
    });
  });
}
