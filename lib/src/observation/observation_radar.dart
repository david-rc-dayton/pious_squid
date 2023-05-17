import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/observation/observation_utils.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Radar observation data.
class ObservationRadar extends Observation {
  /// Create a new [ObservationRadar] object.
  ObservationRadar(this._site, this.observation, [final Matrix? noise]) {
    _noise = noise ?? defaultNoise;
  }

  /// Range-Azimuth-Elevation observation.
  final Razel observation;

  /// Site location.
  final Geodetic _site;

  /// Noise matrix.
  late final Matrix _noise;

  /// Default noise matrix.
  static final Matrix defaultNoise =
      noiseFromSigmas([0.32, 0.015 * deg2rad, 0.015 * deg2rad]);

  @override
  EpochUTC get epoch => observation.epoch;

  @override
  Geodetic get site => _site;

  @override
  Matrix get noise => _noise;

  @override
  Vector toVector() => Vector.fromList(
      [observation.range, observation.azimuth, observation.elevation]);

  @override
  double clos(final Propagator propagator) {
    final siteEci = site.toITRF(epoch).toJ2000().position;
    final ri = propagator.propagate(epoch).position.subtract(siteEci);
    return (observation.range - ri.magnitude()).abs();
  }

  @override
  Vector ricDiff(final Propagator propagator) {
    final r0 = propagator.propagate(epoch);
    final r1 = Razel.fromStateVectors(r0.toITRF(), site.toITRF(epoch));
    final r2 = observation.position(site, r1.azimuth, r1.elevation);
    return RIC.fromJ2000(J2000(epoch, r2, Vector.origin3), r0).position;
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
    final result = array2d(3, 6);
    final siteEcef = site.toITRF(epoch);
    for (var i = 0; i < 6; i++) {
      final step = propPairs.step(i);
      final (high, low) = propPairs.get(i);
      final sl = low.propagate(epoch).toITRF();
      final sh = high.propagate(epoch).toITRF();
      final ol = Razel.fromStateVectors(sl, siteEcef);
      final oh = Razel.fromStateVectors(sh, siteEcef);
      result[0][i] = observationDerivative(oh.range, ol.range, step);
      result[1][i] =
          observationDerivative(oh.azimuth, ol.azimuth, step, isAngle: true);
      result[2][i] = observationDerivative(oh.elevation, ol.elevation, step,
          isAngle: true);
    }
    return Matrix(result);
  }

  @override
  Matrix residual(final Propagator propagator) {
    final result = array2d(3, 1);
    final siteEcef = site.toITRF(epoch);
    final stateEcef = propagator.propagate(epoch).toITRF();
    final razel = Razel.fromStateVectors(stateEcef, siteEcef);
    result[0][0] = observation.range - razel.range;
    result[1][0] = normalizeAngle(observation.azimuth, razel.azimuth);
    result[2][0] = normalizeAngle(observation.elevation, razel.elevation);
    return Matrix(result);
  }
}
