import 'dart:math';

import 'package:pious_squid/pious_squid.dart';

/// Convert radians to degrees.
const rad2deg = 180 / pi;

void main() {
  // ---- J2000 ---------------------------------------------------------------

  // Create a new J2000 state.
  final j2000 = J2000(
      EpochUTC.fromDate(2023, 9, 2, 12, 19, 52.085),
      Vector3D(-526.798998, 5726.586585, 3613.195304),
      Vector3D(-5.859872080, 2.247870322, -4.396729518));

  // Get state epoch (UTC).
  print(j2000.epoch); // => 2023-09-02T12:19:52.085Z

  // Get state position/velocity vectors.
  print(j2000.position);
  // => [-526.798998, 5726.586585, 3613.195304] (km)
  print(j2000.velocity);
  // => [-5.85987208, 2.247870322, -4.396729518] (km/s)

  // Get orbit period.
  print(j2000.period()); // => 5574.86307633312 (seconds)

  // Get orbit semimajor axis (km).
  print(j2000.semimajorAxis()); // => 6795.407153349586 (km)

  // Convert to classical elements.
  print(j2000.toClassicalElements());
  // => [ClassicalElements]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Semimajor Axis (a):       6795.4072 km
  // Eccentricity (e):         0.0015141
  // Inclination (i):          51.5361°
  // Right Ascension (Ω):      305.1981°
  // Argument of Perigee (ω):  68.5444°
  // True Anomaly (ν):         68.6553°

  // Convert to TEME inertial frame.
  print(j2000.toTEME());
  // => [TEME]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [-565.362599, 5723.543990, 3612.187805] km
  // Velocity: [-5.861620421, 2.217034320, -4.410035952] km/s

  // Convert to Earth fixed frame.
  print(j2000.toITRF());
  // => [ITRF]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [1898.887296, -5428.887336, 3612.187805] km
  // Velocity: [5.823207529, -0.910982178, -4.410035952] km/s

  // ---- Geodetic ------------------------------------------------------------

  // Create a new Geodetic object
  final geodetic = Geodetic.fromDegrees(
    38.4, // latitude (deg)
    -104.2, // longitude (deg)
    1.7, // altitude (km)
  );
  print(geodetic);
  // => [Geodetic]
  // Latitude:  38.4000°
  // Longitude: -104.2000°
  // Altitude:  1.700 km

  // Get geodetic fields.
  print(geodetic.latitude); // => 0.6702064327658225 (rad)
  print(geodetic.latitudeDegrees); // => 38.4 (deg)
  print(geodetic.longitude); // => -1.8186330805780915 (rad)
  print(geodetic.longitudeDegrees); // => -104.2 (deg)
  print(geodetic.altitude); // => 1.7 (km)

  // Calculate the angle between two geodetic points.
  final otherGeodetic = Geodetic.fromDegrees(40.1, -72.3, 0.006);
  print(geodetic.angle(otherGeodetic)); // => 0.42986063178253586 (rad)
  print(geodetic.angleDegrees(otherGeodetic)); // => 24.629199979966444 (deg)

  // Compute the surface distance between two points.
  print(geodetic.distance(otherGeodetic)); // => 2738.645555006403 (km)

  // Compute the surface field-of-view half angle at a geodetic point.
  final orbitalGeodetic = j2000.toITRF().toGeodetic();
  final fov = orbitalGeodetic.fieldOfView();
  print(fov); // => 0.3533714256661327 (rad)
  print(fov * rad2deg); // => 20.2466912911903 (deg)

  // Check if line of sight exists between two points.
  print(geodetic.sight(orbitalGeodetic)); // => false

  // Convert to ITRF Earth fixed frame at the provided epoch.
  final geodeticEpoch = EpochUTC.fromDate(2021, 3, 16, 22, 45, 12.321);
  print(geodetic.toITRF(geodeticEpoch));
  // => [ITRF]
  // Epoch: 2021-03-16T22:45:12.321Z
  // Position: [-1228.083275, -4853.337849, 3941.391545] km
  // Velocity: [0.000000000, 0.000000000, 0.000000000] km/s

  // ---- Hill ----------------------------------------------------------------

  // Create Hill state as a linear drift relative to an inertial state.
  final origin = j2000;
  final hill = Hill.fromLinearDrift(
    origin, // relative state origin
    -5, // radial (km)
    -15, // intrack (km)
    0, // node velocity (km/s)
    0, // node offset (seconds)
  );
  print(hill);
  // => [Hill]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [-5.000000, -15.000000, -0.000000] km
  // Velocity: [0.000000000, 0.008452923, 0.000000000] km/s

  // Propagate state using Clohessy-Wiltshire equations.
  final hillPropEpoch = origin.epoch.roll(
    300, // seconds
  );
  final hillPropState = hill.propagate(hillPropEpoch);
  print(hillPropState);
  // => [Hill]
  // Epoch: 2023-09-02T12:24:52.085Z
  // Position: [-5.000000, -12.464123, 0.000000] km
  // Velocity: [0.000000000, 0.008452923, 0.000000000] km/s

  // Convert back to a J2000 state, by using the origin state propagated
  // to the Hill state epoch.
  final originProp = RungeKutta4(origin);
  final originPropState = originProp.propagate(hillPropEpoch);
  final j2000FromHill = hillPropState.toJ2000(originPropState);
  print(j2000FromHill);
  // => [J2000]
  // Epoch: 2023-09-02T12:24:52.085Z
  // Position: [-2211.262831, 6059.131003, 2121.765016] km
  // Velocity: [-5.337273038, -0.009814813, -5.497286630] km/s

  // ---- ITRF ----------------------------------------------------------------

  // Create a new ITRF Earth fixed state.
  final itrf = ITRF(
      EpochUTC.fromDateTimeString('2023-09-02T12:19:52.085Z'),
      Vector3D(1898.887296, -5428.887336, 3612.187805),
      Vector3D(5.823207529, -0.910982178, -4.410035952));
  print(itrf);
  // => [ITRF]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [1898.887296, -5428.887336, 3612.187805] km
  // Velocity: [5.823207529, -0.910982178, -4.410035952] km/s

  // Get height above Earth's reference ellipsoid.
  print(itrf.getHeight()); // => 419.58336347002205 (km)

  // Convert to a J2000 inertial state.
  print(itrf.toJ2000());
  // => [J2000]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [-526.798998, 5726.586585, 3613.195304] km
  // Velocity: [-5.859872080, 2.247870322, -4.396729518] km/s

  // Convert to a geodetic location.
  print(itrf.toGeodetic());
  // => [Geodetic]
  // Latitude:  32.2939°
  // Longitude: -70.7215°
  // Altitude:  419.582 km

  // ---- RIC -----------------------------------------------------------------

  // Create a RIC relative state from two inertial states.
  final ric = RIC.fromJ2000(j2000FromHill, originPropState);
  print(ric);
  // => [RIC]
  // Position: [-5.011423, -12.454945, -0.000000] km
  // Velocity: [0.014054402, 0.002776669, 0.000000000] km/s

  // Calculate relative range and range rate.
  print(ric.range()); // => 13.425349328666746 (km)
  print(ric.rangeRate()); // => -0.007822203634006441 (km/s)

  // Convert back to a J2000 inertial state.
  print(ric.toJ2000(originPropState));
  // => [J2000]
  // Epoch: 2023-09-02T12:24:52.085Z
  // Position: [-2211.262831, 6059.131003, 2121.765016] km
  // Velocity: [-5.337273038, -0.009814813, -5.497286630] km/s

  // ---- TEME ----------------------------------------------------------------

  // Create a new TEME inertial state.
  final teme = TEME(
      EpochUTC.fromDateTimeString('2023-09-02T12:19:52.085Z'),
      Vector3D(-565.362599, 5723.543990, 3612.187804),
      Vector3D(-5.861620421, 2.217034319, -4.410035953));
  print(teme);
  // => [TEME]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [-565.362599, 5723.543990, 3612.187804] km
  // Velocity: [-5.861620421, 2.217034319, -4.410035953] km/s

  // Convert TEME to classical elements.
  print(teme.toClassicalElements());
  // => [ClassicalElements]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Semimajor Axis (a):       6795.4072 km
  // Eccentricity (e):         0.0015141
  // Inclination (i):          51.6444°
  // Right Ascension (Ω):      305.4428°
  // Argument of Perigee (ω):  68.6385°
  // True Anomaly (ν):         68.6553°

  // Convert to a J2000 inertial state.
  print(teme.toJ2000());
  // => [J2000]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [-526.798998, 5726.586585, 3613.195303] km
  // Velocity: [-5.859872080, 2.247870321, -4.396729519] km/s

  // ---- TLE -----------------------------------------------------------------

  // Create a new TLE object.
  final tle = TLE(
      '1 25544U 98067A   23245.51379729  .00012235  00000+0  22079-3 0  9994',
      '2 25544  51.6428 305.4682 0005253  24.3497 112.8045 15.50153461413785');

  // Get TLE epoch.
  print(tle.epoch); // => 2023-09-02T12:19:52.085Z

  // Convert to TEME state vector.
  print(tle.state);
  // => [TEME]
  // Epoch: 2023-09-02T12:19:52.085Z
  // Position: [-565.362599, 5723.543990, 3612.187804] km
  // Velocity: [-5.861620421, 2.217034319, -4.410035953] km/s

  // Propagate a TLE to a new epoch in TEME frame using SGP4.
  print(tle.propagate(tle.epoch.roll(300)));
  // => [TEME]
  // Epoch: 2023-09-02T12:24:52.085Z
  // Position: [-2258.454325, 6051.639512, 2109.027204] km
  // Velocity: [-5.317821061, -0.050598970, -5.515002994] km/s
}
