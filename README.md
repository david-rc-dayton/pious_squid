# pious_squid

An astrodynamics library for the Dart ecosystem, covering orbital mechanics
and satellite mission analysis logic.

## Features

- Celestial Bodies
  - Earth
  - Moon
  - Sun
- Coordinates
  - Classical Orbital Elements
  - Equinoctial Elements
  - Geodetic Coordinates
  - Hill Modified Equidistant Cylindrical Frame _(EQCM)_
  - International Terrestrial Reference Frame _(ITRF)_
  - J2000 Inertial Frame _(J2000)_
  - Radial-Intrack-Crosstrack Relative Frame _(RIC)_
  - True Equator Mean Equinox Inertial Frame _(TEME)_
  - Two-Line Element Set _(TLE)_
- Covariance
  - Cartesian Covariance
  - Classical Element Covariance
  - Equinoctial Element Covariance
  - Relative Covariance
- Perturbation Forces
  - Atmospheric Drag _(Harris-Priester)_
  - Earth Gravity _(up to 36x36 geopotential)_
  - Solar Radiation Pressure
  - Spacecraft Thrust
  - Spherical Body Gravity
  - Third Body Gravity _(Sun and Moon)_
- Interpolators
  - Chebyshev Ephemeris Interpolator
  - Cubic-Spline Ephemeris Interpolator
  - Lagrange Ephemeris Interpolator
  - Verlet-Blend Ephemeris Interpolator
- Maneuvers
  - Relative Waypoint Targeting
  - Two-Burn Transfer _(Hohmann)_
- Metric Observations
  - Geocentric Right-Ascension and Declination _(RaDec)_
  - Optical Observation
  - Radar Observation
  - Range-Azimuth-Elevation _(RAzEl)_
  - Topocentric Right-Ascension and Declination _(RaDec)_
- Math
  - Matrix Operations
  - Quaternion Operations
  - Vector Operations
- Optimization
  - Chebyshev Ephemeris Compression
  - Downhill Simplex _(Nelder-Mead)_
  - Golden Section
  - Polynomial Regression _(WIP)_
  - Simple Linear Regression
- Orbit Determination
  - Batch Least Squares Orbit Determination _(OD)_
  - Gibbs Initial Orbit Determination _(IOD)_
  - Gooding Initial Orbit Determination _(IOD)_
  - Herrik-Gibbs Initial Orbit Determination _(IOD)_
  - Lambert Initial Orbit Determination _(IOD)_
  - Modified Gooding Initial Orbit Determination _(IOD)_
- State Propagation
  - Dormand-Prince 5(4) Adaptive Numerical Propagator
  - Kepler Two-Body Analytical Propagator
  - Runge-Kutta 4 - Fixed Numerical Propagator
  - Runge-Kutta 8(9) Adaptive Numerical Propagator
  - Simplified Perturbations Model 4 _(SGP4)_
- Data Smoothing
  - Exponential Smoothing _(single/double/time)_
- Time
  - Barycentric Dynamical Time _(TDB)_
  - Global Positioning System Time _(GPS)_
  - International Atomic Time _(TAI)_
  - Terrestrial Time _(TT)_
  - Universal Coordinated Time _(UTC)_

## Usage

```dart
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

  // Create a observer location.
  final observer = Geodetic.fromDegrees(-15, 80, 0.05);

  // Calculate look-angles from the observer to the satellite .
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
```
More examples can be found in the
[example](https://github.com/david-rc-dayton/pious_squid/tree/master/example)
directory.

## License

```text
Copyright (c) 2023 David RC Dayton

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
