import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';

/// Propagator pairs used for numerically finding the observation Jacobian.
class PropagatorPairs {
  /// Create a new [PropagatorPairs] object given the position _(km)_ and
  /// velocity _(km/s)_ step.
  PropagatorPairs(this._posStep, this._velStep);
  final List<Propagator?> _high = List.filled(6, null, growable: false);
  final List<Propagator?> _low = List.filled(6, null, growable: false);
  final double _posStep;
  final double _velStep;

  /// Set the [high] and [low] propagator at the provided state [index].
  void set(final int index, final Propagator high, final Propagator low) {
    _high[index] = high;
    _low[index] = low;
  }

  /// Get the high and low propagator at the provided [index].
  (Propagator, Propagator) get(final int index) =>
      (_high[index]!, _low[index]!);

  /// Get the step size at the provided index.
  double step(final int index) => index < 3 ? _posStep : _velStep;
}

/// Convert right-ascension [ra] _(rad)_, declination [dec] _(rad)_, and slant
/// range [r] _(km)_ into a position vector relative to the observer.
Vector radecToPosition(final double ra, final double dec, final double r) {
  final ca = cos(ra);
  final sa = sin(ra);
  final cd = cos(dec);
  final sd = sin(dec);
  final v = Float64List(3);
  v[0] = r * cd * ca;
  v[1] = r * cd * sa;
  v[2] = r * sd;
  return Vector(v);
}

/// Convert right-ascension [ra] _(rad)_, declination [dec] _(rad)_, slant
/// range [r] _(km)_, right-ascension rate [raDot] _(rad/s)_, declination rate
/// [decDot] _(rad/s)_, and slant range rate _(km/s)_ into a velocity vector
/// relative to the observer.
Vector radecToVelocity(final double ra, final double dec, final double r,
    final double raDot, final double decDot, final double rDot) {
  final ca = cos(ra);
  final sa = sin(ra);
  final cd = cos(dec);
  final sd = sin(dec);
  final v = Float64List(3);
  v[0] = rDot * cd * ca - r * sd * ca * decDot - r * cd * sa * raDot;
  v[1] = rDot * cd * sa - r * sd * sa * decDot + r * cd * ca * raDot;
  v[2] = rDot * sd + r * cd * decDot;
  return Vector(v);
}

/// Calculate the difference between two angles.
double normalizeAngle(final double a, final double b) {
  final x = a - b;
  return atan2(sin(x), cos(x));
}

/// Calculate the component derivative for the high [xh] and low [xl]
/// observation values and step size.
///
/// Set [isAngle] to true, if the component is an angular value.
double observationDerivative(
        final double xh, final double xl, final double step,
        {final bool isAngle = false}) =>
    (isAngle ? normalizeAngle(xh, xl) : (xh - xl)) / step;

/// Create a noise matrix from the provided list of standard deviation
/// values [sigmas].
Matrix noiseFromSigmas(final List<double> sigmas) {
  final n = sigmas.length;
  final result = array2d(n, n);
  for (var i = 0; i < n; i++) {
    final s = sigmas[i];
    result[i][i] = 1 / (s * s);
  }
  return Matrix(result);
}
