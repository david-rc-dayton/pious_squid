import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/operations_base.dart';

/// Covariance base class.
abstract class Covariance {
  /// Create a new [Covariance] object given its covariance [matrix].
  Covariance(this.matrix);

  /// Covariance matrix.
  final Matrix matrix;

  /// Calculate the standard deviations along the covariance
  /// diagonal variances.
  Vector sigmas() {
    final c = matrix.columns;
    final result = Float64List(c);
    for (var i = 0; i < c; i++) {
      final variance = matrix[i][i];
      result[i] = sqrt(variance);
    }
    return Vector(result);
  }
}
