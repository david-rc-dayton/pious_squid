import 'dart:math';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Optical observation data.
class ObservationOptical extends Observation {
  /// Create a new [ObservationOptical] object.
  ObservationOptical(this._site, this.observation, [final Matrix? noise])
      : _noise = noise ?? defaultNoise;

  /// Topocentric radec observation.
  final RadecTopocentric observation;

  /// Inertial site location.
  final J2000 _site;

  /// Noise matrix.
  final Matrix _noise;

  /// Default noise matrix _(right-ascension, declination)_.
  ///
  /// Based on the Maui Optical Site noise model.
  static final Matrix defaultNoise =
      noiseFromSigmas(0.0037 * deg2rad, 0.0030 * deg2rad);

  @override
  EpochUTC get epoch => observation.epoch;

  @override
  J2000 get site => _site;

  @override
  Matrix get noise => _noise;

  @override
  Vector toVector() =>
      Vector.fromList([observation.rightAscension, observation.declination]);

  @override
  double clos(final Propagator propagator) {
    final position = propagator.propagate(epoch).position;
    final offset = position.subtract(site.position);
    final actual = observation.lineOfSight().normalize();
    final expected = offset.normalize();
    final slantRange = offset.magnitude();
    final theta = actual.angle(expected);
    if (theta.isNaN) {
      return 0.0;
    }
    return 2.0 * slantRange * sin(theta * 0.5);
  }

  @override
  Vector3D ricDiff(final Propagator propagator) {
    final r0 = site;
    final r1 = propagator.propagate(epoch);
    final r2 = observation.position(site, r1.position.distance(r0.position));
    return RelativeState.fromJ2000(J2000(epoch, r2, Vector3D.origin), r1)
        .position;
  }

  @override
  Observation sample(final RandomGaussianSource random,
      [final double sigma = 1.0]) {
    final result = sampleVector(random, sigma);
    return ObservationOptical(
        site, RadecTopocentric(epoch, result[0], result[1]), noise);
  }

  @override
  Matrix jacobian(final PropagatorPairs propPairs) {
    final result = Matrix(2, 6);
    for (var i = 0; i < 6; i++) {
      final step = propPairs.step(i);
      final (high, low) = propPairs.get(i);
      final sl = low.propagate(epoch);
      final sh = high.propagate(epoch);
      final ol = RadecTopocentric.fromStateVectors(sl, site);
      final oh = RadecTopocentric.fromStateVectors(sh, site);
      result.set(
          0,
          i,
          observationDerivative(oh.rightAscension, ol.rightAscension, step,
              isAngle: true));
      result.set(
          1,
          i,
          observationDerivative(oh.declination, ol.declination, step,
              isAngle: true));
    }
    return result;
  }

  @override
  Matrix residual(final Propagator propagator) {
    final result = Matrix(2, 1);
    final state = propagator.propagate(epoch);
    final radec = RadecTopocentric.fromStateVectors(state, site);
    result.set(
        0, 0, normalizeAngle(observation.rightAscension, radec.rightAscension));
    result.set(
        1, 0, normalizeAngle(observation.declination, radec.declination));
    return result;
  }

  /// Create a noise matrix from right ascension and declination standard
  /// deviantions _(radians)_.
  static Matrix noiseFromSigmas(final double raSigma, final double decSigma) =>
      observationNoiseFromSigmas([raSigma, decSigma]);
}
