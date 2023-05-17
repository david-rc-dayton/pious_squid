import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/covariance/covariance_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Relative state covariance.
class CovarianceRelative extends Covariance {
  /// Create a new [CovarianceRelative] given an inertial [state] and
  /// relative covariance [matrix].
  CovarianceRelative(this.state, final Matrix matrix) : super(matrix);

  /// Inertial state.
  final J2000 state;

  Matrix _relativeToCartesianTransform() {
    final m =
        RelativeState.createMatrix(state.position, state.velocity).transpose();
    return Matrix([
      [m[0][0], m[0][1], m[0][2], 0, 0, 0],
      [m[1][0], m[1][1], m[1][2], 0, 0, 0],
      [m[2][0], m[2][1], m[2][2], 0, 0, 0],
      [0, 0, 0, m[0][0], m[0][1], m[0][2]],
      [0, 0, 0, m[1][0], m[1][1], m[1][2]],
      [0, 0, 0, m[2][0], m[2][1], m[2][2]],
    ]).transpose();
  }

  /// Convert this to a [CovarianceCartesian] object.
  CovarianceCartesian toCartesian() {
    final j = _relativeToCartesianTransform();
    final p = j.multiply(matrix).multiply(j.transpose());
    return CovarianceCartesian(state, p);
  }
}
