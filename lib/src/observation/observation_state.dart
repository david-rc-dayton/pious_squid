import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator.dart';
import 'package:pious_squid/src/time/epoch_utc.dart';

/// Class for ITRF state observations.
class ObservationState extends Observation {
  /// Create a new [ObservationState] from an ITRF state and noise matrix.
  ObservationState(this.observation, [final Matrix? noise])
      : _noise = noise ?? defaultNoise;

  /// ITRF state observation _(km,km/s)_
  final ITRF observation;

  /// ITRF state noise matrix _(km²,km/s²)_
  final Matrix _noise;

  /// Default state noise matrix _(km²,km/s²)_
  static final defaultNoise =
      observationNoiseFromSigmas([0.01, 0.01, 0.1, 0.001, 0.001, 0.001]);

  @override
  double clos(final Propagator propagator) => propagator
      .propagate(epoch)
      .toITRF()
      .position
      .distance(observation.position);

  @override
  EpochUTC get epoch => observation.epoch;

  @override
  J2000 get site => observation.toJ2000();

  @override
  Matrix get noise => _noise;

  /// Create a noise matrix from the position/velocity standard
  /// deviantions _(km,km/s)_.
  static Matrix noiseFromSigmas(final double rx, final double ry,
          final double rz, final double vx, final double vy, final double vz) =>
      observationNoiseFromSigmas([rx, ry, rz, vx, vy, vz]);

  @override
  Matrix jacobian(final PropagatorPairs propPairs) {
    final result = Matrix(6, 6);
    for (var i = 0; i < 6; i++) {
      final step = propPairs.step(i);
      final (high, low) = propPairs.get(i);
      final sl = low.propagate(epoch).toITRF();
      final sh = high.propagate(epoch).toITRF();
      result.set(
          0, i, observationDerivative(sh.position.x, sl.position.x, step));
      result.set(
          1, i, observationDerivative(sh.position.y, sl.position.y, step));
      result.set(
          2, i, observationDerivative(sh.position.z, sl.position.z, step));
      result.set(
          3, i, observationDerivative(sh.velocity.x, sl.velocity.x, step));
      result.set(
          4, i, observationDerivative(sh.velocity.y, sl.velocity.y, step));
      result.set(
          5, i, observationDerivative(sh.velocity.z, sl.velocity.z, step));
    }
    return result;
  }

  @override
  Matrix residual(final Propagator propagator) {
    final result = Matrix(6, 1);
    final state = propagator.propagate(epoch).toITRF();
    result.set(0, 0, observation.position.x - state.position.x);
    result.set(1, 0, observation.position.y - state.position.y);
    result.set(2, 0, observation.position.z - state.position.z);
    result.set(3, 0, observation.velocity.x - state.velocity.x);
    result.set(4, 0, observation.velocity.y - state.velocity.y);
    result.set(5, 0, observation.velocity.z - state.velocity.z);
    return result;
  }

  @override
  Vector3D ricDiff(final Propagator propagator) =>
      RelativeState.fromJ2000(site, propagator.propagate(epoch)).position;

  @override
  Observation sample(final RandomGaussianSource random,
      [final double sigma = 1.0]) {
    final result = sampleVector(random, sigma);
    return ObservationState(
        ITRF(epoch, result.toVector3D(0), result.toVector3D(3)), noise);
  }

  @override
  Vector toVector() => Vector.fromList([
        observation.position.x,
        observation.position.y,
        observation.position.z,
        observation.velocity.x,
        observation.velocity.y,
        observation.velocity.z
      ]);
}
