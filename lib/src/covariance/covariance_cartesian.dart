import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/covariance/covariance_base.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Inertial cartesian state covariance.
class CovarianceCartesian extends Covariance {
  /// Create a new [CovarianceCartesian] given an inertial [state] and
  /// inertial covariance [matrix].
  CovarianceCartesian(this.state, final Matrix matrix) : super(matrix);

  /// Inertial state.
  final J2000 state;

  Matrix _cartesianToEquinoctialTransform() {
    Float64List f(final Float64List x) {
      final ct = J2000(
          state.epoch, Vector3D(x[0], x[1], x[2]), Vector3D(x[3], x[4], x[5]));
      final eq = ct.toClassicalElements().toEquinoctialElements();
      return Float64List.fromList([eq.af, eq.ag, eq.l, eq.n, eq.chi, eq.psi]);
    }

    return jacobian(f, 6, state.position.join(state.velocity).toArray());
  }

  Matrix _cartesianToClassicalTransform() {
    Float64List f(final Float64List x) {
      final ct = J2000(
          state.epoch, Vector3D(x[0], x[1], x[2]), Vector3D(x[3], x[4], x[5]));
      final coe = ct.toClassicalElements();
      return Float64List.fromList([
        coe.semimajorAxis,
        coe.eccentricity,
        coe.inclination,
        coe.rightAscension,
        coe.argPerigee,
        coe.trueAnomaly
      ]);
    }

    return jacobian(f, 6, state.position.join(state.velocity).toArray());
  }

  Matrix _cartesianToRelativeTransform() {
    final m =
        RelativeState.createMatrix(state.position, state.velocity).transpose();
    return Matrix([
      [m[0][0], m[0][1], m[0][2], 0, 0, 0],
      [m[1][0], m[1][1], m[1][2], 0, 0, 0],
      [m[2][0], m[2][1], m[2][2], 0, 0, 0],
      [0, 0, 0, m[0][0], m[0][1], m[0][2]],
      [0, 0, 0, m[1][0], m[1][1], m[1][2]],
      [0, 0, 0, m[2][0], m[2][1], m[2][2]],
    ]);
  }

  /// Convert this to a [CovarianceEquinoctial] object.
  CovarianceEquinoctial toEquinoctial() {
    final j = _cartesianToEquinoctialTransform();
    final p = j.multiply(matrix).multiply(j.transpose());
    return CovarianceEquinoctial(
        state.toClassicalElements().toEquinoctialElements(), p);
  }

  /// Convert this to a [CovarianceClassical] object.
  CovarianceClassical toClassical() {
    final j = _cartesianToClassicalTransform();
    final p = j.multiply(matrix).multiply(j.transpose());
    return CovarianceClassical(state.toClassicalElements(), p);
  }

  /// Convert this to a [CovarianceRelative] object.
  CovarianceRelative toRelative() {
    final j = _cartesianToRelativeTransform();
    final p = j.multiply(matrix).multiply(j.transpose());
    return CovarianceRelative(state, p);
  }
}
