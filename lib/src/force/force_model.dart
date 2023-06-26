import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/atmospheric_drag.dart';
import 'package:pious_squid/src/force/earth_gravity.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/force/gravity.dart';
import 'package:pious_squid/src/force/solar_radiation_pressure.dart';
import 'package:pious_squid/src/force/third_body_gravity.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Base class for perturbation forces.
abstract class Force {
  /// Calculate the acceleration due to the perturbing force on a given
  /// [state] vector.
  Vector3D acceleration(final J2000 state);
}

/// Force model for spacecraft propagation.
class ForceModel {
  /// Create a new [ForceModel] object.
  ForceModel();

  /// Central body gravity model.
  Force? _centralGravity;

  /// Third body gravity model.
  Force? _thirdBodyGravity;

  /// Solar radiation pressure model.
  Force? _solarRadiationPressure;

  /// Atmospheric drag model.
  Force? _atmosphericDrag;

  /// Thrust model.
  Force? _maneuverThrust;

  /// Enable simple central-body gravity model.
  ///
  /// If gravitational parameter [mu] _(km²/s³)_ is not provided, Earth's
  /// parameter will be used.
  void setGravity([final double mu = Earth.mu]) {
    _centralGravity = Gravity(mu);
  }

  /// Enable complex Earth gravity model using the provided geopotential
  /// [degree] and [order] _(max: 36)_.
  void setEarthGravity(final int degree, final int order) {
    _centralGravity = EarthGravity(degree, order);
  }

  /// Enable third-body gravity for the provided bodies.
  void setThirdBodyGravity({final bool moon = false, final bool sun = false}) {
    _thirdBodyGravity = ThirdBodyGravity(moon: moon, sun: sun);
  }

  /// Enable solar radiation pressure using the spacecraft mass _(kg)_,
  /// cross-sectional area _(m²)_, and optional reflectivity coefficient.
  void setSolarRadiationPressure(final double mass, final double area,
      {final double coeff = 1.2}) {
    _solarRadiationPressure = SolarRadiationPressure(mass, area, coeff);
  }

  /// Enable atmospheric drag using the spacecraft mass _(kg)_,
  /// cross-sectional area _(m²)_, optional drag coefficient, and optional
  /// cosine exponent.
  ///
  /// The cosine exponent should be a number between `2` for low inclination
  /// orbits and `6` for polar orbits.
  ///
  /// **Note:** The atmospheric drag model assumes mean solar flux.
  void setAtmosphericDrag(final double mass, final double area,
      {final double coeff = 2.2, final int cosine = 4}) {
    _atmosphericDrag = AtmosphericDrag(mass, area, coeff, cosine);
  }

  /// Load a [maneuver] into the acceleration model.
  void loadManeuver(final Thrust maneuver) {
    _maneuverThrust = maneuver;
  }

  /// Clear a [Thrust] maneuver from the model if one is loaded.
  void clearManeuver() {
    _maneuverThrust = null;
  }

  /// Calculate inertial acceleration on a [state] using the forces in
  /// this model.
  Vector3D acceleration(final J2000 state) {
    var accVec = Vector3D.origin;
    if (_centralGravity != null) {
      accVec = accVec.add(_centralGravity!.acceleration(state));
    }
    if (_thirdBodyGravity != null) {
      accVec = accVec.add(_thirdBodyGravity!.acceleration(state));
    }
    if (_solarRadiationPressure != null) {
      accVec = accVec.add(_solarRadiationPressure!.acceleration(state));
    }
    if (_atmosphericDrag != null) {
      accVec = accVec.add(_atmosphericDrag!.acceleration(state));
    }
    if (_maneuverThrust != null) {
      accVec = accVec.add(_maneuverThrust!.acceleration(state));
    }
    return accVec;
  }

  /// Calculate inertial [state] derivative using the forces in this model.
  Vector derivative(final J2000 state) =>
      state.velocity.join(acceleration(state));
}
