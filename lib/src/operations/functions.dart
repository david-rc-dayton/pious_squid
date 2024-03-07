import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Differentiable function that takes a [double] and returns a [double].
typedef DifferentiableFunction = double Function(double x);

/// Jacobian function that takes an array of [double] values and returns an
/// array of [double] values.
typedef JacobianFunction = Float64List Function(Float64List xs);

/// Method for calculating angular distance.
enum AngularDistanceMethod {
  /// Law of cosines, fast but inaccurate.
  cosine,

  /// Haversine formula, slow but accurate.
  haversine
}

/// Method for calculating angular diameter.
enum AngularDiameterMethod {
  /// Circle formula, for flat objects.
  circle,

  /// Sphere formula, for spherical objects.
  sphere
}

/// Calculate the factorial of input [n].
double factorial(final int n) {
  final nAbs = n.abs();
  var result = 1.0;
  for (var i = 2; i <= nAbs; i++) {
    result *= i;
  }
  return result;
}

/// Calculate the log10 of [x].
double log10(final double x) => log(x) / ln10;

/// Calculate the cube root of [n].
double cbrt(final double n) => pow(n, 1 / 3) as double;

/// Calculate the hyperbolic sine of [x].
double sinh(final double x) => 0.5 * (exp(x) - exp(-x));

/// Calculate the hyperbolic cosine of [x].
double cosh(final double x) => 0.5 * (exp(x) + exp(-x));

/// Calculate the hyperbolic tangent of [x].
double tanh(final double x) => sinh(x) / cosh(x);

/// Calculate the hyperbolic cotangent of [x].
double coth(final double x) => cosh(x) / sinh(x);

/// Calculate the hyperbolic secant of [x].
double sech(final double x) => 1 / cosh(x);

/// Calculate the hyperbolic cosecant of [x].
double csch(final double x) => 1 / sinh(x);

/// Calculate the hyperbolic arcsine of [x].
double asinh(final double x) => log(x + sqrt(x * x + 1));

/// Calculate the hyperbolic arccosine of [x].
double acosh(final double x) => log(x + sqrt(x * x - 1));

/// Calculate the hyperbolic arctangent of [x].
double atanh(final double x) => 0.5 * log((1 + x) / (1 - x));

/// Calculate the hyperbolic arccosecant of [x].
double acsch(final double x) => log((1 / x) + sqrt((1 / (x * x)) + 1));

/// Calculate the hyperbolic arcsecant of [x].
double asech(final double x) => log((1 / x) + sqrt((1 / (x * x)) - 1));

/// Calculate the hyperbolic arccotangent of [x].
double acoth(final double x) => 0.5 * log((x + 1) / (x - 1));

/// Copy the sign of [sgn] to the magnitude of [mag].
double copySign(final double mag, final double sgn) => mag.abs() * sgn.sign;

/// Evaluate polynomial coefficients for the provied value.
///
/// Polynomial coefficients must be ordered from largest exponent to smallest
/// exponent; the last coefficient must relate to `x⁰ ≡ 1`.
double evalPoly(final double x, final Float64List coeffs) {
  var result = coeffs[0];
  for (var i = 1; i < coeffs.length; i++) {
    result = result * x + coeffs[i];
  }
  return result;
}

/// Translate the provided [angle] _(rad)_ to be in the same half-plane as the
/// [match] _(rad)_ argument.
double matchHalfPlane(final double angle, final double match) {
  final a1 = angle;
  final a2 = twoPi - angle;
  final d1 = atan2(sin(a1 - match), cos(a1 - match));
  final d2 = atan2(sin(a2 - match), cos(a2 - match));
  return d1.abs() < d2.abs() ? a1 : a2;
}

/// Wrap a (0, 2π) angle to a (-π, π) range.
double wrapAngle(final double theta) {
  final result = ((theta + pi) % twoPi) - pi;
  if (result == -pi) {
    return pi;
  }
  return result;
}

/// Calculate the angular distance _(rad)_ between two phi and lambda angle
/// sets using the law of cosines.
double _angularDistanceCosine(final double lam1, final double phi1,
    final double lam2, final double phi2) {
  final a = sin(phi1) * sin(phi2);
  final b = cos(phi1) * cos(phi2) * cos(lam2 - lam1);
  return acos(a + b);
}

/// Calculate the angular distance _(rad)_ between two phi and lambda angle
/// sets using the Haversine formula.
double _angularDistanceHaversine(final double lam1, final double phi1,
    final double lam2, final double phi2) {
  final dlam = lam2 - lam1;
  final dphi = phi2 - phi1;
  final sdlam = sin(0.5 * dlam);
  final sdphi = sin(0.5 * dphi);
  final a = sdphi * sdphi + cos(phi1) * cos(phi2) * sdlam * sdlam;
  return 2.0 * asin(min(1.0, sqrt(a)));
}

/// Calculate the angular distance _(rad)_ between two phi and lambda angle
/// sets using the selected [method].
double angularDistance(
    final double lam1, final double phi1, final double lam2, final double phi2,
    {final AngularDistanceMethod method = AngularDistanceMethod.cosine}) {
  switch (method) {
    case AngularDistanceMethod.cosine:
      return _angularDistanceCosine(lam1, phi1, lam2, phi2);
    case AngularDistanceMethod.haversine:
      return _angularDistanceHaversine(lam1, phi1, lam2, phi2);
  }
}

/// Calculate the angular diameter _(rad)_ given the actual [diameter],
/// [distance], and calculation [method].
double angularDiameter(final double diameter, final double distance,
    [final AngularDiameterMethod method = AngularDiameterMethod.sphere]) {
  switch (method) {
    case AngularDiameterMethod.circle:
      return 2 * atan(diameter / (2 * distance));
    case AngularDiameterMethod.sphere:
      return 2 * asin(diameter / (2 * distance));
  }
}

/// Linearly interpolate the output value for [x] given two xy pairs.
double linearInterpolate(final double x, final double x0, final double y0,
        final double x1, final double y1) =>
    (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);

/// Calculate the mean from a list of [values].
double mean(final List<double> values) {
  final n = values.length;
  var sum = 0.0;
  for (final v in values) {
    sum += v;
  }
  return sum / n;
}

/// Calculate the standard deviation from a list of [values].
///
/// Set [isSample] to `true` if not using the entire value population.
double standardDeviation(final List<double> values,
    {final bool isSample = false}) {
  final mu = mean(values);
  final n = values.length;
  var sum = 0.0;
  for (final v in values) {
    final sub = v - mu;
    sum += sub * sub;
  }
  final m = isSample ? 1 : 0;
  return sqrt((1.0 / (n - m)) * sum);
}

/// Calculate the covariance between the values of list [a] and [b].
///
/// Set [isSample] to `true` if not using the entire value population.
double covariance(final List<double> a, final List<double> b,
    {final bool isSample = false}) {
  final n = a.length;
  final am = mean(a);
  final bm = mean(b);
  var result = 0.0;
  for (var i = 0; i < n; i++) {
    result += (a[i] - am) * (b[i] - bm);
  }
  final m = isSample ? 1 : 0;
  return result / (n - m);
}

/// Create a new function that performs central finite differencing to return
/// the derivative of the provided function [f].
///
/// The finite differencing step [h] can be passed as an optional parameter.
DifferentiableFunction derivative(final DifferentiableFunction f,
    {final double h = 1e-3}) {
  double df(final double x) {
    final hh = h * 0.5;
    return (f(x + hh) - f(x - hh)) / h;
  }

  return df;
}

/// Calculate the gamma result of [n].
double gamma(final int n) => factorial(n - 1);

/// Return the inverted value [x], or zero if the value is zero.
double invertZero(final double x) => x == 0 ? 0.0 : 1.0 / x;

/// Convert eccentricity [ecc] and mean anomaly [m] _(rad)_ into eccentric
/// anomaly and true anomaly _(rad)_;
({double e0, double nu}) newtonM(final double ecc, final double m) {
  final numiter = 50;
  final small = 1e-8;
  double e0;
  double nu;

  // elliptical
  if (ecc > small) {
    if (((m < 0.0) && (m > -pi)) || (m > pi)) {
      e0 = m - ecc;
    } else {
      e0 = m + ecc;
    }
    var ktr = 1;
    var e1 = e0 + (m - e0 + ecc * sin(e0)) / (1.0 - ecc * cos(e0));

    while (((e1 - e0).abs() > small) && (ktr <= numiter)) {
      ktr++;
      e0 = e1;
      e1 = e0 + (m - e0 + ecc * sin(e0)) / (1.0 - ecc * cos(e0));
    }
    final sinv = sqrt(1.0 - ecc * ecc) * sin(e1) / (1.0 - ecc * cos(e1));
    final cosv = (cos(e1) - ecc) / (1.0 - ecc * cos(e1));
    nu = atan2(sinv, cosv);
  } else {
    // circular
    nu = m;
    e0 = m;
  }
  return (e0: e0, nu: nu);
}

/// Convert eccentricity [ecc] and true anomaly [nu] _(rad)_ into eccentric
/// anomaly and mean anomaly _(rad)_;
({double e0, double m}) newtonNu(final double ecc, final double nu) {
  final small = 1e-8;
  var e0 = 0.0;
  var m = 0.0;

  // circular
  if (ecc.abs() < small) {
    m = nu;
    e0 = nu;
  } else {
    // elliptical
    if (ecc < 1.0 - small) {
      final sine = (sqrt(1.0 - ecc * ecc) * sin(nu)) / (1.0 + ecc * cos(nu));
      final cose = (ecc + cos(nu)) / (1.0 + ecc * cos(nu));
      e0 = atan2(sine, cose);
      m = e0 - ecc * sin(e0);
    }
  }
  if (ecc < 1.0) {
    m -= (m / twoPi).floor() * twoPi;
    if (m < 0.0) {
      m += 2.0 * pi;
    }
    e0 -= (e0 / twoPi).floor() * twoPi;
  }

  return (e0: e0, m: m);
}

/// Create a 2D array with the provided [rows] and [columns] and
/// default [value].
List<List<T>> array2d<T>(final int rows, final int columns, final T value,
    {final bool growable = false}) {
  final output = <List<T>>[];
  for (var i = 0; i < rows; i++) {
    output.add(List.filled(columns, value, growable: growable));
  }
  return output;
}

/// Copy [n] values from a [sourceArray] starting at [sourcePosition] into
/// the [destinationArray] starting at [destinationPositon].
void arraycopy<T>(final List<T> sourceArray, final int sourcePosition,
    final List<T> destinationArray, final int destinationPositon, final int n) {
  var sourceDex = sourcePosition;
  var destDex = destinationPositon;
  var iter = 0;
  while (iter < n) {
    destinationArray[destDex] = sourceArray[sourceDex];
    sourceDex++;
    destDex++;
    iter++;
  }
}

/// Use central finite differencing to generate a Jacobian matrix for
/// function [f], output length [m], and input value [x0].
Matrix jacobian(final JacobianFunction f, final int m, final Float64List x0,
    {final double step = 1e-5}) {
  final n = x0.length;
  final j = array2d(m, n, 0.0);
  final h = 0.5 * step;

  for (var k = 0; k < n; k++) {
    // forward
    final xp = x0.sublist(0);
    xp[k] += h;
    final fp = Vector(f(xp));

    // backward
    final xm = x0.sublist(0);
    xm[k] -= h;
    final fm = Vector(f(xm));

    // central difference
    final cd = fp.subtract(fm).scale(1.0 / step);

    // update matrix
    for (var i = 0; i < m; i++) {
      j[i][k] = cd[i];
    }
  }
  return Matrix(j);
}

/// Compute the Mahalanobis distance between an [actual] and [expected] vector
/// given [covariance].
double mahalanobisDistance(
    final Vector actual, final Vector expected, final Matrix covariance) {
  final z = actual.subtract(expected).toColumnMatrix();
  final t = z.transpose().multiply(covariance.inverse()).multiply(z);
  return sqrt(t[0][0]);
}
