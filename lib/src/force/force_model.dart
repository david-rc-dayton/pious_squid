import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/atmospheric_drag.dart';
import 'package:pious_squid/src/force/earth_gravity.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/force/gravity.dart';
import 'package:pious_squid/src/force/solar_radiation_pressure.dart';
import 'package:pious_squid/src/force/third_body_gravity.dart';
import 'package:pious_squid/src/operations/functions.dart';
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
  SolarRadiationPressure? _solarRadiationPressure;

  /// Atmospheric drag model.
  AtmosphericDrag? _atmosphericDrag;

  /// Thrust model.
  Force? _maneuverThrust;

  /// Ballistic coefficient _(kg/m²)_.
  double _bcoeffInv = 0.0;

  /// Solar radiation pressure coefficient _(kg/m²)_.
  double _srpcoeffInv = 0.0;

  /// Get this force model's ballistic coefficient _(kg/m²)_.
  double getBCoeff() => invertZero(_bcoeffInv);

  /// Get this force model's solar radiation pressure coefficient _(kg/m²)_.
  double getSrpCoeff() => invertZero(_srpcoeffInv);

  /// Set this force model's ballistic coefficient _(kg/m²)_.
  void setBCoeff(final double bcoeff) {
    _bcoeffInv = invertZero(bcoeff);
  }

  /// Set this force model's solar radiation pressure coefficient _(kg/m²)_.
  void setSrpCoeff(final double srpcoeff) {
    _srpcoeffInv = invertZero(srpcoeff);
  }

  /// Enable simple central-body gravity model.
  ///
  /// If gravitational parameter [mu] _(km²/s³)_ is not provided, Earth's
  /// parameter will be used.
  void setGravity([final double mu = Earth.mu]) {
    _centralGravity = Gravity(mu);
  }

  /// Enable complex Earth gravity model using the provided geopotential
  /// [degree] and [order] _(max: 70)_.
  void setEarthGravity(final int degree, final int order) {
    _centralGravity = EarthGravity(degree, order);
  }

  /// Enable third-body gravity for the provided bodies.
  void setThirdBodyGravity({final bool moon = false, final bool sun = false}) {
    _thirdBodyGravity = ThirdBodyGravity(moon: moon, sun: sun);
  }

  /// Enable solar radiation pressure using the spacecraft solar radiation
  /// pressure coefficient _()_, and optional reflectivity coefficient.
  void setSolarRadiationPressure(final double srpcoeff,
      {final double reflectCoeff = 1.2}) {
    _srpcoeffInv = invertZero(srpcoeff);
    _solarRadiationPressure = SolarRadiationPressure(reflectCoeff);
  }

  /// Enable atmospheric drag using the spacecraft ballistic coefficient
  /// _(kg/m²)_, optional drag coefficient, and optional cosine exponent.
  ///
  /// The cosine exponent should be a number between `2` for low inclination
  /// orbits and `6` for polar orbits.
  ///
  /// **Note:** The atmospheric drag model assumes mean solar flux if space
  /// weather data is not provided.
  void setAtmosphericDrag(final double bcoeff,
      {final double dragCoeff = 2.2, final int cosine = 4}) {
    _bcoeffInv = invertZero(bcoeff);
    _atmosphericDrag = AtmosphericDrag(dragCoeff, cosine);
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
    if (_atmosphericDrag != null) {
      accVec = accVec.add(_atmosphericDrag!.acceleration(state, _bcoeffInv));
    }
    if (_solarRadiationPressure != null) {
      accVec = accVec
          .add(_solarRadiationPressure!.acceleration(state, _srpcoeffInv));
    }
    if (_maneuverThrust != null) {
      accVec = accVec.add(_maneuverThrust!.acceleration(state));
    }
    return accVec;
  }

  /// Calculate inertial [state] derivative using the forces in this model.
  Vector derivative(final J2000 state) =>
      state.velocity.join(acceleration(state));

  /// Clone this force model into a new object.
  ForceModel clone() {
    final output = ForceModel();

    // gravity
    if (_centralGravity != null) {
      if (_centralGravity is EarthGravity) {
        final force = _centralGravity as EarthGravity;
        output.setEarthGravity(force.degree, force.order);
      } else if (_centralGravity is Gravity) {
        final force = _centralGravity as Gravity;
        output.setGravity(force.mu);
      }
    }
    // third body
    if (_thirdBodyGravity != null) {
      final force = _thirdBodyGravity as ThirdBodyGravity;
      output.setThirdBodyGravity(sun: force.sun, moon: force.moon);
    }
    // solar radiation pressure
    if (_solarRadiationPressure != null) {
      final force = _solarRadiationPressure as SolarRadiationPressure;
      output.setSolarRadiationPressure(getSrpCoeff(),
          reflectCoeff: force.reflectCoeff);
    }
    // atmospheric drag
    if (_atmosphericDrag != null) {
      final force = _atmosphericDrag as AtmosphericDrag;
      output.setAtmosphericDrag(getBCoeff(),
          dragCoeff: force.dragCoeff, cosine: force.cosine);
    }
    // thrust
    if (_maneuverThrust != null) {
      final force = _maneuverThrust as Thrust;
      output.loadManeuver(force);
    }

    // coefficients
    output.setBCoeff(invertZero(_bcoeffInv));
    output.setSrpCoeff(invertZero(_bcoeffInv));

    return output;
  }
}
