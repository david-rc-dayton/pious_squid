import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/operations_base.dart';

/// Covaraiance coordinate frame.
enum CovarianceFrame {
  /// J2000 inertial frame
  j2000,

  /// ITRF earth-fixed frame
  itrf,

  /// Radial-Intrack-Crosstrack frame
  ric,
}

/// State covariance.
class StateCovariance {
  /// Create a new [StateCovariance] object given its covariance [matrix] and
  /// [CovarianceFrame].
  StateCovariance(this.matrix, this.frame);

  /// Create a new [StateCovariance] object from 1-sigma values.
  factory StateCovariance.fromSigmas(
      final List<double> sigmas, final CovarianceFrame frame) {
    final n = sigmas.length;
    final output = Matrix.zero(n, n);
    for (var i = 0; i < n; i++) {
      output[i][i] = max(sigmas[i] * sigmas[i], 1e-32);
    }
    return StateCovariance(output, frame);
  }

  /// Covariance matrix.
  final Matrix matrix;

  /// Covariance frame.
  final CovarianceFrame frame;

  /// Calculate the standard deviations along the covariance
  /// diagonal variances.
  Vector sigmas([final double stddev = 1.0]) {
    final c = matrix.columns;
    final result = Float64List(c);
    for (var i = 0; i < c; i++) {
      final variance = matrix[i][i];
      result[i] = sqrt(variance) * stddev;
    }
    return Vector(result);
  }
}
