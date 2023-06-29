import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Interpolator base class.
abstract class Interpolator {
  /// Return the start and end epoch covered by this interpolator.
  EpochWindow window();

  /// Return `true` if the provided [epoch] is within this interpolator's
  /// cached value range.
  bool inWindow(final EpochUTC epoch) {
    final (start, stop) = window();
    return start <= epoch && epoch <= stop;
  }

  /// Calculate the start/stop epoch between this and another [Interpolator].
  ///
  /// Returns `null` if there is no overlap between interpolators.
  EpochWindow? overlap(final Interpolator interpolator) {
    final (x1, x2) = window();
    final (y1, y2) = interpolator.window();
    if ((x1 <= y2) && (y1 <= x2)) {
      final e1 = EpochUTC(max(x1.posix, y1.posix));
      final e2 = EpochUTC(min(x2.posix, y2.posix));
      return (e1, e2);
    }
    return null;
  }
}

/// Base class for state vector interpolators.
abstract class StateInterpolator extends Interpolator {
  /// Interpolate an intertial state vector at the provided [epoch].
  ///
  /// Returns `null` if the epoch is outside the cached ephemeris range.
  J2000? interpolate(final EpochUTC epoch);

  /// Return the size _(bytes)_ of this interpolator's cached data.
  int get sizeBytes;
}

/// Base class for attitude interpolators.
abstract class AttitudeInterpolator extends Interpolator {
  /// Interpolate an attitude [Quaternion] at the provided [epoch].
  ///
  /// Returns `null` if the epoch is outside the cached attitude range.
  Quaternion? interpolate(final EpochUTC epoch);
}

/// Base class for arbitrary field interpolators.
abstract class FieldInterpolator extends Interpolator {
  /// Interpolate field values at the provided [epoch].
  ///
  /// Returns `null` if the epoch is outside the cached field range.
  Float64List? interpolate(final EpochUTC epoch);
}
