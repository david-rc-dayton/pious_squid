import 'package:pious_squid/pious_squid.dart';

const secondsPerDay = 86400.0; // seconds

void main() {
  // Create a new J2000 state.
  final state = J2000(
    EpochUTC.fromDate(2023, 9, 2, 12, 19, 52.085), // utc
    Vector3D(-526.798998, 5726.586585, 3613.195304), // km
    Vector3D(-5.859872080, 2.247870322, -4.396729518), // km/s
  );
  print(state);
  // => [J2000]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [-526.798998, 5726.586585, 3613.195304] km
  // Velocity: [-5.859872080, 2.247870322, -4.396729518] km/s

  // Create a TLE object from an element set.
  final tle = TLE(
      '1 26388U 00034A   23245.65928248  .00000071  00000+0  00000+0 0  9990',
      '2 26388  11.0718  39.4320 0001820 217.0596 150.5171  1.00306462 84979');

  // Create a +1 m/s intrack maneuver, in RIC frame.
  final maneuver = Thrust(EpochUTC.fromDate(2023, 9, 2, 14, 0, 0), 0, 1, 0);

  // Define vehicle properties.
  final mass = 1400.0; // kg
  final area = 16.0; // m^2

  // Create a force model of satellite perturbation effects.
  final forceModel = ForceModel()
    ..setEarthGravity(36, 36)
    ..setThirdBodyGravity(sun: true, moon: true)
    ..setAtmosphericDrag(mass, area)
    ..setSolarRadiationPressure(mass, area);

  // Define a propagation epoch one day in the future.
  final propEpoch = state.epoch.roll(secondsPerDay);

  // --------------------------------------------------------------------------
  // Keplerian (Two-Body) Propagation
  // --------------------------------------------------------------------------

  // Create a new Kepler Propagator. This only models two-body gravitational
  // forces, and a force model cannot be used here.
  final keplerProp = KeplerPropagator(state.toClassicalElements());
  print(keplerProp.propagate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [437.289110, -5697.783123, -3684.370146] km
  // Velocity: [5.861151425, -2.338683120, 4.332153746] km/s

  // --------------------------------------------------------------------------
  // Fixed Step Propagation
  // --------------------------------------------------------------------------

  // Create a fixed step propagator, with an optional step size in seconds.
  final rk4Prop = RungeKutta4Propagator(state, forceModel, 15.0);
  print(rk4Prop.propagate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [380.649661, -5911.213183, -3332.559823] km
  // Velocity: [5.623275705, -2.267192474, 4.677707648] km/s

  // --------------------------------------------------------------------------
  // Adaptive Step Propagation
  // --------------------------------------------------------------------------

  // Create an adaptive step propagator, with an optional tolerance.
  final rk89Prop = RungeKutta89Propagator(state, forceModel, 1e-9);
  print(rk89Prop.propagate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [380.647375, -5911.212294, -3332.561750] km
  // Velocity: [5.623275881, -2.267195519, 4.677705900] km/s

  // --------------------------------------------------------------------------
  // SGP4 Propagation
  // --------------------------------------------------------------------------

  // Create an SGP4 propagator.
  final sgp4Prop = Sgp4Propagator(tle);
  print(sgp4Prop.propagate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [41658.245286, -3041.331455, -5702.914987] km
  // Velocity: [0.279968176, 3.032598531, 0.424305959] km/s

  // --------------------------------------------------------------------------
  // Propagator Methods
  // --------------------------------------------------------------------------

  // Propagate to apogee after the propagation epoch.
  print(rk89Prop.propagateApogee(propEpoch));
  // => [J2000]
  //   Epoch: 2023-09-03T13:30:27.907Z
  //   Position: [-4952.421703, 1625.602174, -4369.858629] km
  //   Velocity: [0.788540066, -6.796430678, -3.421958633] km/s

  // Propagate to perigee after the propagation epoch.
  print(rk89Prop.propagatePerigee(propEpoch));
  // => [J2000]
  //   Epoch: 2023-09-03T12:49:05.182Z
  //   Position: [4422.176988, 496.420118, 5121.699040] km
  //   Velocity: [-2.626183697, 7.026830728, 1.586424502] km/s

  // Propagate to the ascending node after the propagation epoch.
  print(rk89Prop.propagateAscendingNode(propEpoch));
  // => [J2000]
  //   Epoch: 2023-09-03T12:29:52.085Z
  //   Position: [3419.157113, -5868.411319, 0.000385] km
  //   Velocity: [4.114659551, 2.407441184, 6.003348144] km/s

  // Propagate to the descending node after the propagation epoch.
  print(rk89Prop.propagateDescendingNode(propEpoch));
  // => [J2000]
  //   Epoch: 2023-09-03T13:16:12.907Z
  //   Position: [-3405.923131, 5883.389223, -0.000025] km
  //   Velocity: [-4.126179601, -2.378569951, -5.997750061] km/s

  // Create a checkpoint at this state if you want to restore it later.
  final checkpoint = rk89Prop.checkpoint();
  print(checkpoint); // => 0;

  // Maneuver the propagator spacecraft, and propagate 5 minutes.
  rk89Prop.maneuver(maneuver);
  print(rk89Prop.propagate(rk89Prop.state.epoch.roll(300)));
  // => [J2000]
  // Epoch: 2023-09-02T14:05:00.000Z
  // Position: [-4167.329247, 5353.291789, -452.200230] km
  // Velocity: [-3.526422821, -3.239368097, -5.977094709] km/s

  // Restore our previous checkpoint.
  rk89Prop.restore(checkpoint);
  print(rk89Prop.state);
  // => [J2000]
  // Epoch: 2023-09-03T13:16:12.907Z
  // Position: [-3405.922639, 5883.389508, 0.000808] km
  // Velocity: [-4.126180407, -2.378569093, -5.997749849] km/s

  // Reset the propagator to its initial state in the constructor.
  rk89Prop.reset();
  print(rk89Prop.state);
  // => [J2000]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [-526.798998, 5726.586585, 3613.195304] km
  // Velocity: [-5.859872080, 2.247870322, -4.396729518] km/s

  // --------------------------------------------------------------------------
  // Interpolators
  // --------------------------------------------------------------------------

  // Create a week's worth freeflight ephemeris from a propagator, as a
  // Verlet-Blend interpolator.
  final ephemeris = rk89Prop.ephemeris(
      EpochUTC.fromDate(2023, 9, 2), EpochUTC.fromDate(2023, 9, 9));

  // Get the time window covered by this ephemeris.
  final (ephemStart, ephemStop) = ephemeris.window();
  print(ephemStart); // => 2023-09-02T00:00:00.000Z
  print(ephemStop); // => 2023-09-09T00:01:00.000Z

  // Check if an epoch is inside this interpolator's window.
  print(ephemeris.inWindow(EpochUTC.fromDate(2023, 9, 7))); // => true

  // Interpolate a state at an epoch.
  print(ephemeris.interpolate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [380.647219, -5911.212166, -3332.562213] km
  // Velocity: [5.623277599, -2.267213194, 4.677791562] km/s

  // Interpolation returns null if the epoch is not in the epoch window.
  print(ephemeris.interpolate(EpochUTC.fromDate(2023, 10, 4))); // => null

  // If using Verlet-Blend, you can get the closest non-interpolated state to
  // the provided epoch.
  print(ephemeris.getCachedState(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:20:00.000Z
  // Position: [425.139911, -5928.921910, -3295.405304] km
  // Velocity: [5.619228876, -2.207721799, 4.711094142] km/s

  // To make interpolation faster and more accurate for dense ephemeris, you
  // can convert to a Cubic-Spline.
  final cubicSpline = ephemeris.toCubicSpline();
  print(cubicSpline.interpolate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [380.647402, -5911.212239, -3332.561685] km
  // Velocity: [5.623276391, -2.267209727, 4.677697390] km/s

  // To get even more accurate interpolation but with reduced speed, convert
  // to Lagrange, and provide an optional interpolation order.
  final lagrange = ephemeris.toLagrange(10);
  print(lagrange.interpolate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [380.647404, -5911.212305, -3332.561725] km
  // Velocity: [5.623275894, -2.267195491, 4.677705927] km/s

  // For the fastest interpolation speed and smallest memory requirements you
  // can sacrifice accuracy and use Chebyshev with an optional coefficients
  // per revolution argument.
  final chebyshev = ChebyshevCompressor(cubicSpline).compress(21);
  print(chebyshev.interpolate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [380.647827, -5911.212855, -3332.561828] km
  // Velocity: [5.623279572, -2.267198003, 4.677705494] km/s

  // Interpolators can also contain maneuvers.
  rk89Prop.reset();
  final ephemerisManeuver = rk89Prop.ephemerisManeuver(
      EpochUTC.fromDate(2023, 9, 2), EpochUTC.fromDate(2023, 9, 4), [maneuver]);

  // The same interpolator logic applies.
  print(ephemerisManeuver.interpolate(propEpoch));
  // => [J2000]
  // Epoch: 2023-09-03T12:19:52.085Z
  // Position: [204.520172, -5839.315553, -3478.780291] km
  // Velocity: [5.633027814, -2.500614424, 4.539719953] km/s

  // Epoch window overlap between two interpolaters can be computed, or null is
  // returned if there is no overlap.
  final (overlapStart, overlapStop) = ephemerisManeuver.overlap(ephemeris)!;
  print(overlapStart); // => 2023-09-02T00:00:00.000Z
  print(overlapStop); // => 2023-09-04T00:00:00.000Z
}
