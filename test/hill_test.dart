import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

void main() {
  group('Hill', () {
    test('solveManeuver', () {
      final target = J2000(
          EpochUTC.fromDateTimeString('2017-01-07T05:31:00.243Z'),
          Vector3D(-5737.369776, -3423.651756, 364.099770),
          Vector3D(4.378112704, -6.646623519, 1.170571889));

      final satellite = Hill.fromPerch(target, -10, 0.005, 0);

      final waypoint = Waypoint(target.epoch.roll(1200), Vector3D(0, 0, 0));
      final maneuver = satellite.solveManeuver(waypoint);
      final maneuverNoCt =
          satellite.solveManeuver(waypoint, ignoreCrosstrack: true);

      final expected = [-7.3287, 5.3057, -5.0000];
      expect(maneuver.center, equals(target.epoch));
      expect(maneuver.radial, closeTo(expected[0], 1e-3));
      expect(maneuver.intrack, closeTo(expected[1], 1e-3));
      expect(maneuver.crosstrack, closeTo(expected[2], 1e-3));
      expect(maneuver.duration, equals(0.0));
      expect(maneuverNoCt.center, equals(target.epoch));
      expect(maneuverNoCt.radial, closeTo(expected[0], 1e-3));
      expect(maneuverNoCt.intrack, closeTo(expected[1], 1e-3));
      expect(maneuverNoCt.crosstrack, equals(0.0));
      expect(maneuverNoCt.duration, equals(0.0));
    });
  });
}
