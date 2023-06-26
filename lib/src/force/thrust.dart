import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Thrust force model.
class Thrust implements Force {
  /// Create a new [Thrust] object given the maneuver [center] time and
  /// relative maneuver components _(m/s)_.
  ///
  /// An optional thrust [durationRate] _(s/m/s)_ can be provided to model
  /// finite maneuvers.
  Thrust(this.center, final double radial, final double intrack,
      final double crosstrack,
      {this.durationRate = 0.0}) {
    deltaV = Vector3D(radial * 1e-3, intrack * 1e-3, crosstrack * 1e-3);
  }

  /// Maneuver center time.
  final EpochUTC center;

  /// Delta-V vector _(km/s)_.
  late final Vector3D deltaV;

  /// Maneuver duration rate _(s/m/s)_.
  final double durationRate;

  /// Radial thrust component _(m/s)_.
  double get radial => deltaV.x * 1000.0;

  /// In-track thrust component _(m/s)_.
  double get intrack => deltaV.y * 1000.0;

  /// Cross-track thrust component _(m/s)_.
  double get crosstrack => deltaV.z * 1000.0;

  /// Thrust magnitude _(m/s)_.
  double get magnitude => deltaV.magnitude() * 1000.0;

  /// Maneuver duration _(seconds)_.
  double get duration => magnitude * durationRate;

  /// Maneuver start epoch.
  EpochUTC get start => center.roll(-0.5 * duration);

  /// Maneuver end epoch.
  EpochUTC get stop => center.roll(0.5 * duration);

  @override
  Vector3D acceleration(final J2000 state) {
    final relative = RIC(Vector3D.origin, deltaV.scale(1.0 / duration));
    return relative.toJ2000(state).velocity.subtract(state.velocity);
  }

  /// Return a copy of the provided [state] with this maneuver applied.
  J2000 apply(final J2000 state) {
    final relative = RIC(Vector3D.origin, deltaV);
    return relative.toJ2000(state);
  }

  /// Return `true` if this maneuver is impulsive, otherwise this maneuver
  /// is finite.
  bool get isImpulsive => duration <= 0;
}
