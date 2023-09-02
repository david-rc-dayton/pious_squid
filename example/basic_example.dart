import 'package:pious_squid/pious_squid.dart';

void main() {
  // Create a new J2000 inertial satellite state.
  final startState = J2000(
      EpochUTC.fromDateTimeString('2017-02-03T06:26:37.976Z'),
      Vector3D(-3134.15877, 7478.695162, 1568.694229),
      Vector3D(-5.227261462, -3.7717234, 2.643938099));

  // Define some spacecraft properties.
  final mass = 1400.0; // kilograms
  final area = 16.0; // meters²

  // Create a perturbation force model.
  final forceModel = ForceModel()
    // Model a 36x36 geopotential.
    ..setEarthGravity(36, 36)
    // Model Moon and Sun gravity.
    ..setThirdBodyGravity(moon: true, sun: true)
    // Model solar radiation pressure, with reflectivity coefficient 1.2.
    ..setSolarRadiationPressure(mass, area, coeff: 1.2)
    // Model atmospheric drag, with drag coefficient 2.2.
    ..setAtmosphericDrag(mass, area, coeff: 2.2);

  // Create a Runge-Kutta 8(9) propagator.
  final rk89Prop = RungeKutta89Propagator(startState, forceModel);

  // Propagate the start state to 1 day in the future.
  final oneDay = 86400.0; // seconds
  final finalState = rk89Prop.propagate(startState.epoch.roll(oneDay));

  print(finalState);
  // ->[J2000]
  // Epoch: 2017-02-04T06:26:37.976Z
  // Position: [5704.152590, -5470.867067, -3040.596164] km
  // Velocity: [4.554130436, 4.557924086, -2.152201166] km/s
}
