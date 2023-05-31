import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Attitude Quaternion.
class Quaternion {
  /// Create a new [Quaternion] object, real component last.
  Quaternion(this.x, this.y, this.z, this.w);

  /// Create a new [Quaternion] object from the provided [axis] and angle
  /// [theta] _(deg)_.
  factory Quaternion.fromAxisAngle(final Vector axis, final double theta) {
    final w = cos(theta * 0.5);
    final hSin = sin(theta * 0.5);
    final x = axis.x * hSin;
    final y = axis.y * hSin;
    final z = axis.z * hSin;
    return Quaternion(x, y, z, w);
  }

  /// Create a new [Quaternion] that will point from an observer position to
  /// a target position given a forward unit vector.
  factory Quaternion.lookAt(
      final Vector observer, final Vector target, final Vector forward) {
    final forwardVector = target.subtract(observer).normalize();
    final rotAxis = forward.cross(forwardVector);
    final dot = forward.dot(forwardVector);
    return Quaternion(rotAxis.x, rotAxis.y, rotAxis.z, dot + 1.0).normalize();
  }

  /// X-axis.
  final double x;

  /// Y-axis.
  final double y;

  /// Z-axis.
  final double z;

  /// Real component.
  final double w;

  /// Zero quaternion.
  static final Quaternion zero = Quaternion(0, 0, 0, 0);

  /// Identity quaternion.
  static final Quaternion one = Quaternion(0, 0, 0, 1);

  /// X-axis unit quaternion.
  static final Quaternion xAxis = Quaternion(1, 0, 0, 0);

  /// Y-axis unit quaternion.
  static final Quaternion yAxis = Quaternion(0, 1, 0, 0);

  /// Z-axis unit quaternion.
  static final Quaternion zAxis = Quaternion(0, 0, 1, 0);

  /// Convert this [Quaternion] to its positive polar form.
  Quaternion positivePolar() => (w >= 0) ? normalize() : negate().normalize();

  /// Calculate the magnitude squared of this [Quaternion].
  double magnitudeSquared() => w * w + x * x + y * y + z * z;

  /// Calculate the magnitude of this [Quaternion].
  double magnitude() => sqrt(magnitudeSquared());

  /// Multiply this [Quaternion] by a scalar value.
  Quaternion scale(final double n) => Quaternion(n * x, n * y, n * z, n * w);

  /// Negate the components of this [Quaternion].
  Quaternion negate() => scale(-1);

  /// Convert this to a unit [Quaternion].
  Quaternion normalize() {
    final m = magnitude();
    if (m == 0) {
      return zero;
    }
    return scale(1 / m);
  }

  /// Return the conjugate of this [Quaternion].
  Quaternion conjugate() => Quaternion(-x, -y, -z, w);

  /// Return the inverse of this [Quaternion].
  Quaternion reciprocal() => conjugate().scale(1 / magnitudeSquared());

  /// Add this and another [Quaternion].
  Quaternion add(final Quaternion q) =>
      Quaternion(x + q.x, y + q.y, z + q.z, w + q.w);

  /// Subtract this and another quaternion.
  Quaternion subtract(final Quaternion q) =>
      Quaternion(x - q.x, y - q.y, z - q.z, w - q.w);

  /// Add the provided number [n] to the real component of this [Quaternion].
  Quaternion addReal(final double n) => Quaternion(x, y, z, w + n);

  /// Multiply this by another [Quaternion].
  Quaternion multiply(final Quaternion q) {
    final mx = w * q.x + x * q.w + y * q.z - z * q.y;
    final my = w * q.y - x * q.z + y * q.w + z * q.x;
    final mz = w * q.z + x * q.y - y * q.x + z * q.w;
    final mw = w * q.w - x * q.x - y * q.y - z * q.z;
    return Quaternion(mx, my, mz, mw);
  }

  /// Calculate the dot product of this and another [Quaternion].
  double dot(final Quaternion q) => x * q.x + y * q.y + z * q.z + w * q.w;

  /// Apply this rotation to a [Vector].
  Vector rotateVector(final Vector v) {
    final q = multiply(Quaternion(v.x, v.y, v.z, 0)).multiply(conjugate());
    final result = Float64List(3);
    result[0] = q.x;
    result[1] = q.y;
    result[2] = q.z;
    return Vector(result);
  }

  /// Linearly interpolate between this and another quaternion by the given
  /// ratio [t] _(0.0, 1.0)_.
  ///
  /// Returns the positive polar form of the resulting quaternion.
  Quaternion lerp(final Quaternion q, final double t) {
    final f = 1.0 - t;
    return Quaternion(
            f * x + t * q.x, f * y + t * q.y, f * z + t * q.z, f * w + t * q.w)
        .positivePolar();
  }

  ///  Spherically interpolate between this and another quaternion by the given
  ///  ratio [t] _(0.0, 1.0)_.
  ///
  ///  Returns the positive polar form of the resulting quaternion.
  Quaternion slerp(final Quaternion q, final double t) {
    var qp = q;
    var dotP = dot(qp);
    if (dotP < 0) {
      dotP = -dotP;
      qp = qp.negate();
    }
    if (dotP > 0.9995) {
      return lerp(qp, t);
    }
    final theta = acos(dotP);
    final sinTheta = sin(theta);
    final f1 = sin((1.0 - t) * theta) / sinTheta;
    final f2 = sin(t * theta) / sinTheta;
    return Quaternion(f1 * x + f2 * qp.x, f1 * y + f2 * qp.y,
            f1 * z + f2 * qp.z, f1 * w + f2 * qp.w)
        .positivePolar();
  }

  /// Return the imaginary components of this [Quaternion] as a [Vector].
  Vector toVector() {
    final v = Float64List(3);
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return Vector(v);
  }

  /// Calculate the angle _(rad)_ between this and another [Quaternion].
  double angle(final Quaternion q) {
    final c = multiply(q.conjugate()).normalize();
    return wrapAngle(2 * atan2(c.toVector().magnitude(), c.w));
  }

  /// Calculate the geodesic angle _(rad)_ between this and
  /// another [Quaternion].
  double geodesicAngle(final Quaternion q) {
    final p = dot(q);
    return wrapAngle(acos((2 * p * p) - 1.0));
  }

  /// Calculate the distance between this and another [Quaternion].
  double distance(final Quaternion q) {
    final m01 = subtract(q).magnitude();
    final p01 = add(q).magnitude();
    return m01 < p01 ? m01 : p01;
  }

  /// Return a delta between this and another [Quaternion].
  Quaternion delta(final Quaternion qTo) => reciprocal().multiply(qTo);

  /// Convert this [Quaternion] to a direction cosine [Matrix].
  Matrix toDirectionCosineMatrix() {
    final w2 = w * w;
    final x2 = x * x;
    final y2 = y * y;
    final z2 = z * z;
    final m = [
      [w2 + x2 - y2 - z2, 2 * (x * y + w * z), 2.0 * (x * z - w * y)],
      [2.0 * (x * y - w * z), w2 - x2 + y2 - z2, 2.0 * (y * z + w * x)],
      [2.0 * (x * z + w * y), 2 * (y * z - w * x), w2 - x2 - y2 + z2],
    ];
    return Matrix(m);
  }

  /// Calculate the angle between this and the [Quaternion] pointing between
  /// the [observer] and [target] positions, given the [forward] vector.
  double vectorAngle(
      final Vector observer, final Vector target, final Vector forward) {
    final delta = target.subtract(observer);
    final transform = toDirectionCosineMatrix().multiplyVector(delta);
    return forward.angle(transform);
  }
}
