import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Radar observation data.
class ObservationRadar extends Observation {
  /// Create a new [ObservationRadar] object.
  ObservationRadar(this._site, this.observation, [final Matrix? noise])
      : _noise = noise ?? defaultNoise;

  /// Range-Azimuth-Elevation observation.
  final Razel observation;

  /// Site location.
  final J2000 _site;

  /// Noise matrix.
  final Matrix _noise;

  /// Default noise matrix.
  static final Matrix defaultNoise =
      noiseFromSigmas(0.32, 0.015 * deg2rad, 0.015 * deg2rad);

  @override
  EpochUTC get epoch => observation.epoch;

  @override
  J2000 get site => _site;

  @override
  Matrix get noise => _noise;

  @override
  Vector toVector() => Vector.fromList(
      [observation.range, observation.azimuth, observation.elevation]);

  @override
  double clos(final Propagator propagator) {
    final ri = propagator.propagate(epoch).position.subtract(site.position);
    return (observation.range - ri.magnitude()).abs();
  }

  @override
  Vector3D ricDiff(final Propagator propagator) {
    final r0 = propagator.propagate(epoch);
    final r1 = Razel.fromStateVectors(r0, site);
    final r2 = observation.position(site, r1.azimuth, r1.elevation);
    return RelativeState.fromJ2000(J2000(epoch, r2, Vector3D.origin), r0)
        .position;
  }

  @override
  Observation sample(final RandomGaussianSource random,
      [final double sigma = 1.0]) {
    final result = sampleVector(random, sigma);
    return ObservationRadar(
        site, Razel(observation.epoch, result[0], result[1], result[2]), noise);
  }

  @override
  Matrix jacobian(final PropagatorPairs propPairs) {
    final result = Matrix(3, 6);
    for (var i = 0; i < 6; i++) {
      final step = propPairs.step(i);
      final (high, low) = propPairs.get(i);
      final sl = low.propagate(epoch);
      final sh = high.propagate(epoch);
      final ol = Razel.fromStateVectors(sl, site);
      final oh = Razel.fromStateVectors(sh, site);
      result.set(0, i, observationDerivative(oh.range, ol.range, step));
      result.set(1, i,
          observationDerivative(oh.azimuth, ol.azimuth, step, isAngle: true));
      result.set(
          2,
          i,
          observationDerivative(oh.elevation, ol.elevation, step,
              isAngle: true));
    }
    return result;
  }

  @override
  Matrix residual(final Propagator propagator) {
    final result = Matrix(3, 1);
    final state = propagator.propagate(epoch);
    final razel = Razel.fromStateVectors(state, site);
    result.set(0, 0, observation.range - razel.range);
    result.set(1, 0, normalizeAngle(observation.azimuth, razel.azimuth));
    result.set(2, 0, normalizeAngle(observation.elevation, razel.elevation));
    return result;
  }

  /// Create a noise matrix from the range, azimuth, and elevation standard
  /// deviantions _(kilometers/radians)_.
  static Matrix noiseFromSigmas(
          final double rngSigma, final double azSigma, final double elSigma) =>
      observationNoiseFromSigmas([rngSigma, azSigma, elSigma]);
}
