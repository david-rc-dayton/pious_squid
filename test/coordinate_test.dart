import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

void main() {
  group('ITRF', () {
    test('PolarHeight', () {
      final northPole =
          ITRF(EpochUTC(0), Vector3D(0, 0, Earth.radiusPolar), Vector3D.origin)
              .toGeodetic();
      final southPole =
          ITRF(EpochUTC(0), Vector3D(0, 0, -Earth.radiusPolar), Vector3D.origin)
              .toGeodetic();
      expect(northPole.altitude, equals(0));
      expect(southPole.altitude, equals(0));
    });
  });
  group('RelativeState', () {
    test('Transfer', () {
      final state = J2000.fromClassicalElements(
          ClassicalElements(EpochUTC(0), 42164, 1e-9, 1e-9, 0, 0, 0));
      final chief = RelativeState.fromPerch(state, 50, 0.01, 0,
          type: RelativeStateType.eqcm);

      final waypoint = Waypoint(state.epoch.roll(3600), Vector3D(0, 0, 0));
      final maneuver = chief.solveManeuver(waypoint);
      final maneuverNoCt =
          chief.solveManeuver(waypoint, ignoreCrosstrack: true);

      final relState = chief.maneuver(maneuver).propagate(waypoint.epoch);
      final relStateNoCt =
          chief.maneuver(maneuverNoCt).propagate(waypoint.epoch);

      expect(relState.position.magnitude(), lessThan(1e-6));
      expect(relStateNoCt.position.x.abs(), lessThan(1e-6));
      expect(relStateNoCt.position.y.abs(), lessThan(1e-6));
      expect(maneuverNoCt.crosstrack, equals(0));
    });
  });
}
