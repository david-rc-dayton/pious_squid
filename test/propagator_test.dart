import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

final startState = J2000(
    EpochUTC.fromDateTimeString('2017-01-07T05:31:00.243Z'),
    Vector3D(-5737.369776, -3423.651756, 364.099770),
    Vector3D(4.378112704, -6.646623519, 1.170571889));

final stopEpoch = EpochUTC.fromDateTimeString('2017-01-10T04:46:49.139Z');

final mass = 1400.0; // kilograms
final area = 16.0; // meters²
final forceModel = ForceModel()
  ..setEarthGravity(36, 36)
  ..setThirdBodyGravity(moon: true, sun: true)
  ..setSolarRadiationPressure(mass, area, coeff: 1.2)
  ..setAtmosphericDrag(mass, area, coeff: 2.2);

final tle = TLE(
    '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753',
    '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667');

void main() {
  group('Propagator', () {
    test('Kepler', () {
      final propagator = KeplerPropagator(startState.toClassicalElements());
      final expected = Vector3D(-251.600120, -6643.127745, 1031.665425);
      final actual = propagator.propagate(stopEpoch).position;
      expect(actual.distance(expected), lessThanOrEqualTo(0.01));
    });

    test('Runge-Kutta 4', () {
      final propagator = RungeKutta4(startState, forceModel);
      final expected = Vector3D(5059.691657, -4729.021976, 638.641366);
      final actual = propagator.propagate(stopEpoch).position;
      expect(actual.distance(expected), lessThanOrEqualTo(0.15));
    });

    test('Dormand-Prince 5(4)', () {
      final propagator = DormandPrince54Propagator(startState, forceModel);
      final expected = Vector3D(5059.691657, -4729.021976, 638.641366);
      final actual = propagator.propagate(stopEpoch).position;
      expect(actual.distance(expected), lessThanOrEqualTo(0.13));
    });

    test('Runge-Kutta 8(9)', () {
      final propagator = RungeKutta89Propagator(startState, forceModel);
      final expected = Vector3D(5059.691657, -4729.021976, 638.641366);
      final actual = propagator.propagate(stopEpoch).position;
      expect(actual.distance(expected), lessThanOrEqualTo(0.13));
    });

    test('SGP4', () {
      final propagator = Sgp4Propagator(tle);
      final expectedPos =
          Vector3D(-7154.03120202, -3783.17682504, -3536.19412294);
      final expectedVel = Vector3D(4.741887409, -4.151817765, -2.093935425);
      final actual = propagator.propagate(tle.epoch.roll(21600.0)).toTEME();
      expect(actual.position.distance(expectedPos), lessThanOrEqualTo(0.01));
      expect(actual.velocity.distance(expectedVel), lessThanOrEqualTo(0.01));
    });
  });
}
