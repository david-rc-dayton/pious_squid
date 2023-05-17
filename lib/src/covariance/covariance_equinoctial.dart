import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/covariance/covariance_base.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Equinoctial element set covariance.
class CovarianceEquinoctial extends Covariance {
  /// Create a new [CovarianceEquinoctial] object given a set of [elements] and
  /// covariance [matrix].
  CovarianceEquinoctial(this.elements, final Matrix matrix) : super(matrix);

  /// Equinoctial elements.
  final EquinoctialElements elements;

  Matrix _equinoctialToCartesianTransform() {
    Float64List f(final Float64List x) {
      final eq = EquinoctialElements(
          elements.epoch, x[0], x[1], x[2], x[3], x[4], x[5]);
      final ct = J2000.fromClassicalElements(eq.toClassicalElements());
      return ct.position.join(ct.velocity).toArray();
    }

    return jacobian(
        f,
        6,
        Float64List.fromList([
          elements.af,
          elements.ag,
          elements.l,
          elements.n,
          elements.chi,
          elements.psi
        ]));
  }

  /// Convert this to a [CovarianceCartesian] object.
  CovarianceCartesian toCartesian() {
    final j = _equinoctialToCartesianTransform();
    final p = j.multiply(matrix).multiply(j.transpose());
    return CovarianceCartesian(
        J2000.fromClassicalElements(elements.toClassicalElements()), p);
  }
}
