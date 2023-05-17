import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/covariance/covariance_base.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Classical element set covariance.
class CovarianceClassical extends Covariance {
  /// Create a new [CovarianceClassical] object given a set of [elements] and
  /// covariance [matrix].
  CovarianceClassical(this.elements, final Matrix matrix) : super(matrix);

  /// Classical element set.
  ClassicalElements elements;

  Matrix _classicalToCartesianTransform() {
    Float64List f(final Float64List x) {
      final coe =
          ClassicalElements(elements.epoch, x[0], x[1], x[2], x[3], x[4], x[5]);
      final ct = J2000.fromClassicalElements(coe);
      return ct.position.join(ct.velocity).toArray();
    }

    return jacobian(
        f,
        6,
        Float64List.fromList([
          elements.semimajorAxis,
          elements.eccentricity,
          elements.inclination,
          elements.rightAscension,
          elements.argPerigee,
          elements.trueAnomaly
        ]));
  }

  /// Conver this to a [CovarianceCartesian] object.
  CovarianceCartesian toCartesian() {
    final j = _classicalToCartesianTransform();
    final p = j.multiply(matrix).multiply(j.transpose());
    return CovarianceCartesian(J2000.fromClassicalElements(elements), p);
  }
}
