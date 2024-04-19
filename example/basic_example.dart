import 'dart:io';

import 'package:pious_squid/pious_squid.dart';

void main() {
  // Optionally, load Earth Orientation Parameter (EOP) data.
  DataHandler().updateEarthOrientationParametersFromCsv(
      File('external/EOP-All.csv').readAsStringSync());

  // Optionally, load Space Weather (SW) data.
  DataHandler().updateSpaceWeatherFromCsv(
      File('external/SW-All.csv').readAsStringSync());

  // Create a new J2000 inertial satellite state.
  final startState = J2000(
    EpochUTC.fromDateTimeString('2017-02-03T06:26:37.976Z'), // utc
    Vector3D(-3134.15877, 7478.695162, 1568.694229), // km
    Vector3D(-5.227261462, -3.7717234, 2.643938099), // km/s
  );

  // Define some spacecraft properties.
  final massArea = 87.5; // kg/m²

  // Create a perturbation force model.
  final forceModel = ForceModel()
    // Model a 36x36 geopotential.
    ..setEarthGravity(36, 36)
    // Model Moon and Sun gravity.
    ..setThirdBodyGravity(moon: true, sun: true)
    // Model solar radiation pressure, with reflectivity coefficient 1.2.
    ..setSolarRadiationPressure(massArea, reflectCoeff: 1.2)
    // Model atmospheric drag, with drag coefficient 2.2.
    ..setAtmosphericDrag(massArea, dragCoeff: 2.2);

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

  // Create a observer location.
  final observer = Geodetic.fromDegrees(-15, 80, 0.05);

  // Calculate look-angles from the observer to the satellite.
  print(finalState.toITRF().toGeodetic());
  final razel = Razel.fromStateVectors(
      finalState, observer.toITRF(finalState.epoch).toJ2000());
  print(razel);
  // => [RazEl]
  // Epoch:     2017-02-04T06:26:37.976Z
  // Azimuth:   141.6525°
  // Elevation: 60.3304°
  // Range:     2318.580 km
}
