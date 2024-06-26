import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Covaraiance coordinate frame.
enum CovarianceFrame {
  /// J2000 inertial frame
  j2000,

  /// ITRF earth-fixed frame
  itrf,

  /// Radial-Intrack-Crosstrack frame
  ric,

  /// Equinoctial elements
  equinoctial,
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
    final output = Matrix(n, n);
    for (var i = 0; i < n; i++) {
      output.set(i, i, max(sigmas[i] * sigmas[i], machineEpsilon));
    }
    return StateCovariance(output, frame);
  }

  /// Create a new [StateCovariance] object from lower-triangular values.
  factory StateCovariance.fromLowerTriangle(
      final List<double> values, final CovarianceFrame frame) {
    final n = (0.5 * (sqrt(8 * values.length + 1) - 1)).toInt();
    final output = Matrix(n, n);
    var dex = 0;
    for (var i = 0; i < n; i++) {
      for (var j = 0; j <= i; j++) {
        output.set(i, j, values[dex]);
        if (i != j) {
          output.set(j, i, values[dex]);
        }
        dex++;
      }
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
      final variance = matrix.get(i, i);
      result[i] = sqrt(variance) * stddev;
    }
    return Vector(result);
  }

  /// Convert RIC covariance to J2000 frame.
  ///
  /// Will throw an error if this covariance isn't in RIC or J2000 frame.
  StateCovariance toJ2000(final J2000 state) {
    if (frame == CovarianceFrame.j2000) {
      return StateCovariance(matrix, frame);
    }
    if (frame != CovarianceFrame.ric) {
      throw 'Conversion to RIC only supported for J2000 covariance.';
    }
    final transform =
        RelativeState.createMatrix6(state.position, state.velocity);
    return StateCovariance(
        transform.transpose().multiply(matrix).multiply(transform),
        CovarianceFrame.j2000);
  }

  /// Convert J2000 covariance to RIC frame.
  ///
  /// Will throw an error if this covariance isn't in RIC or J2000 frame.
  StateCovariance toRIC(final J2000 state) {
    if (frame == CovarianceFrame.ric) {
      return StateCovariance(matrix, frame);
    }
    if (frame != CovarianceFrame.j2000) {
      throw 'Conversion to RIC only supported for J2000 covariance.';
    }
    final transform =
        RelativeState.createMatrix6(state.position, state.velocity);
    return StateCovariance(
        transform.multiply(matrix).multiply(transform.transpose()),
        CovarianceFrame.ric);
  }

  /// Convert covariance to right ascension, declination, and range covariance.
  ///
  /// Will throw an error if this covariance isn't in J2000 frame.
  Matrix toRadec(final J2000 state, final J2000 site) {
    if (frame != CovarianceFrame.j2000) {
      throw 'Conversion to Radec only supported for J2000 covariance.';
    }
    final x = state.position.x - site.position.x;
    final y = state.position.y - site.position.y;
    final z = state.position.z - site.position.z;
    final x2y2 = x * x + y * y;
    final r = sqrt(x * x + y * y + z * z);
    final r3 = r * r * r;
    final rxy = sqrt(x * x + y * y);
    final hArr = Float64List(6);
    hArr[0] = -y / x2y2;
    hArr[1] = x / x2y2;
    hArr[2] = 0;
    hArr[3] = -(x * z) / r3;
    hArr[4] = -(y * z) / r3;
    hArr[5] = rxy / r3;
    final h = Matrix(2, 3, hArr);
    return h.multiply(matrix.getBlock(0, 0, 3, 3)).multiply(h.transpose());
  }

  /// Convert to range, azimuth, and elevation covariance.
  ///
  /// Will throw an error if this covariance isn't in J2000 frame.
  Matrix toRazel(final J2000 state, final J2000 site) {
    if (frame != CovarianceFrame.j2000) {
      throw 'Conversion to Razel only supported for J2000 covariance.';
    }
    final x = state.position.x - site.position.x;
    final y = state.position.y - site.position.y;
    final z = state.position.z - site.position.z;
    final x2y2 = x * x + y * y;
    final r = sqrt(x * x + y * y + z * z);
    final r3 = r * r * r;
    final rxy = sqrt(x * x + y * y);
    final r3sz = r3 * sqrt(1 - pow(z / r, 2));
    final hArr = Float64List(9);
    hArr[0] = x / r;
    hArr[1] = y / r;
    hArr[2] = z / r;
    hArr[3] = -y / x2y2;
    hArr[4] = -x / x2y2;
    hArr[6] = -(x * z) / r3sz;
    hArr[7] = -(y * z) / r3sz;
    hArr[8] = rxy / r3;
    final h = Matrix(3, 3, hArr);
    return h.multiply(matrix.getBlock(0, 0, 3, 3)).multiply(h.transpose());
  }
}
