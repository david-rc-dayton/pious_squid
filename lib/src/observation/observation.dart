import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/observation/observation_utils.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Base class for observation types.
abstract class Observation {
  /// Observation epoch.
  EpochUTC get epoch;

  /// Inertial observer location.
  J2000 get site;

  /// Observation noise matrix.
  Matrix get noise;

  /// Return range-normalized cross line-of-sight residual for the observation
  /// when compared against a nominal state propagator.
  double clos(final Propagator propagator);

  /// Return relative state residual for the observation when compared against
  /// a nominal state propagator.
  Vector3D ricDiff(final Propagator propagator);

  /// Convert this observation to vector form.
  Vector toVector();

  /// Compute the rate derivative matrix for this observation.
  Matrix jacobian(final PropagatorPairs propPairs);

  /// Compute the state residual matrix for this observation.
  Matrix residual(final Propagator propagator);

  /// Convert this observation's noise matrix into a covariance matrix.
  Matrix noiseCovariance() => noise.reciprocal();

  Matrix _noiseSample(final double sigma) =>
      noiseCovariance().scale(sigma).cholesky();

  /// Randomly sample this observation in vector form within the
  /// observation noise.
  Vector sampleVector(final RandomGaussianSource random, final double sigma) {
    final chol = _noiseSample(sigma);
    final gauss = random.gaussVector(noise.columns);
    final meas = toVector();
    return meas.add(chol.multiplyVector(gauss));
  }

  /// Randomly sample this observation within the observation noise, scaled to
  /// the provided [sigma] value.
  Observation sample(final RandomGaussianSource random,
      [final double sigma = 1.0]);
}
