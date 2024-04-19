import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/operations_base.dart';

/// Attitude Quaternion.
class Quaternion {
  /// Create a new [Quaternion] object, real component last.
  const Quaternion(this.x, this.y, this.z, this.w);

  /// Create a new [Quaternion] object from the provided [axis] and angle
  /// [theta] _(deg)_.
  factory Quaternion.fromAxisAngle(final Vector3D axis, final double theta) {
    final w = cos(theta * 0.5);
    final hSin = sin(theta * 0.5);
    final x = axis.x * hSin;
    final y = axis.y * hSin;
    final z = axis.z * hSin;
    return Quaternion(x, y, z, w);
  }

  /// Create a new [Quaternion] object from the provided direction cosine
  /// matrix [c].
  factory Quaternion.fromDirectionCosineMatrix(final Matrix c) {
    final c00 = c.get(0, 0);
    final c11 = c.get(1, 1);
    final c22 = c.get(2, 2);
    final q2 = Float64List(4);
    q2[0] = 0.25 * (1 + c00 - c11 - c22); // qx
    q2[1] = 0.25 * (1 - c00 + c11 - c22); // qy
    q2[2] = 0.25 * (1 - c00 - c11 + c22); // qz
    q2[3] = 0.25 * (1 + c00 + c11 + c22); // qw
    var n = 0;
    var qn2 = q2[0];
    for (var i = 1; i < 4; i++) {
      if (q2[i] > qn2) {
        n = i;
        qn2 = q2[i];
      }
    }
    final s = 1.0 / (4.0 * sqrt(qn2));
    var output = Quaternion.one;
    switch (n) {
      case 0: // qx
        output = Quaternion(4.0 * qn2, c.get(0, 1) + c.get(1, 0),
            c.get(2, 0) + c.get(0, 2), c.get(1, 2) - c.get(2, 1));
        break;
      case 1: // qy
        output = Quaternion(c.get(0, 1) + c.get(1, 0), 4.0 * qn2,
            c.get(1, 2) + c.get(2, 1), c.get(2, 0) - c.get(0, 2));
        break;
      case 2: // qz
        output = Quaternion(c.get(2, 0) + c.get(0, 2),
            c.get(1, 2) + c.get(2, 1), 4.0 * qn2, c.get(0, 1) - c.get(1, 0));
        break;
      case 3: // qw
        output = Quaternion(c.get(1, 2) - c.get(2, 1),
            c.get(2, 0) - c.get(0, 2), c.get(0, 1) - c.get(1, 0), 4.0 * qn2);
        break;
    }
    return output.scale(s).positivePolar();
  }

  /// Create a new [Quaternion] object from a 3-2-1 ordered [EulerAngles].
  factory Quaternion.fromEulerAngles321(final EulerAngles ea) =>
      Quaternion.fromDirectionCosineMatrix(ea.dcm321());

  /// Create a new [Quaternion] that will point from an observer position to
  /// a target position given a forward unit vector.
  factory Quaternion.lookAt(
      final Vector3D observer, final Vector3D target, final Vector3D forward) {
    final forwardVector = target.subtract(observer).normalize();
    final rotAxis = forward.cross(forwardVector);
    final dot = forward.dot(forwardVector);
    return Quaternion(rotAxis.x, rotAxis.y, rotAxis.z, dot + 1.0).normalize();
  }

  /// Perform Tri-Axial Attitude Determination (TRIAD) to find the rotation
  /// quaternion that transforms reference unit vectors [v1] and [v2] to
  /// the new frame defined by [w1] and [w2].
  factory Quaternion.triad(final Vector3D v1, final Vector3D v2,
      final Vector3D w1, final Vector3D w2) {
    final rr = v1.cross(v2).normalize();
    final sr = v1.cross(rr);
    final mrArr = Float64List(9);
    mrArr[0] = v1.x;
    mrArr[1] = v1.y;
    mrArr[2] = v1.z;
    mrArr[3] = rr.x;
    mrArr[4] = rr.y;
    mrArr[5] = rr.z;
    mrArr[6] = sr.x;
    mrArr[7] = sr.y;
    mrArr[8] = sr.z;
    final mr = Matrix(3, 3, mrArr);
    final rb = w1.cross(w2).normalize();
    final sb = w1.cross(rb);
    final mbArr = Float64List(9);
    mbArr[0] = w1.x;
    mbArr[1] = w1.y;
    mbArr[2] = w1.z;
    mbArr[3] = rb.x;
    mbArr[4] = rb.y;
    mbArr[5] = rb.z;
    mbArr[6] = sb.x;
    mbArr[7] = sb.y;
    mbArr[8] = sb.z;
    final mb = Matrix(3, 3, mbArr);
    final a = mr.transpose().multiply(mb);
    return Quaternion.fromDirectionCosineMatrix(a);
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
  static const Quaternion zero = Quaternion(0, 0, 0, 0);

  /// Identity quaternion.
  static const Quaternion one = Quaternion(0, 0, 0, 1);

  /// X-axis unit quaternion.
  static const Quaternion xAxis = Quaternion(1, 0, 0, 0);

  /// Y-axis unit quaternion.
  static const Quaternion yAxis = Quaternion(0, 1, 0, 0);

  /// Z-axis unit quaternion.
  static const Quaternion zAxis = Quaternion(0, 0, 1, 0);

  @override
  String toString({final int precision = 8}) {
    final xStr = x.toStringAsFixed(precision);
    final yStr = y.toStringAsFixed(precision);
    final zStr = z.toStringAsFixed(precision);
    final wStr = w.toStringAsFixed(precision);
    return 'Q(x: $xStr, y: $yStr, z: $zStr, w: $wStr)';
  }

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
  Quaternion inverse() => conjugate().scale(1 / magnitudeSquared());

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

  /// Apply this rotation to a [Vector3D].
  Vector3D rotateVector3D(final Vector3D v) {
    final q = multiply(Quaternion(v.x, v.y, v.z, 0)).multiply(conjugate());
    return Vector3D(q.x, q.y, q.z);
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

  /// Return the imaginary components of this [Quaternion] as a [Vector3D].
  Vector3D toVector3D() => Vector3D(x, y, z);

  /// Calculate the angle _(rad)_ between this and another [Quaternion].
  double angle(final Quaternion q) {
    final c = multiply(q.conjugate()).normalize();
    return 2 * atan2(c.toVector3D().magnitude(), c.w);
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
  Quaternion delta(final Quaternion qTo) => inverse().multiply(qTo);

  /// Convert this [Quaternion] to a direction cosine [Matrix].
  Matrix toDirectionCosineMatrix() {
    final w2 = w * w;
    final x2 = x * x;
    final y2 = y * y;
    final z2 = z * z;
    final mArr = Float64List(9);
    mArr[0] = w2 + x2 - y2 - z2;
    mArr[1] = 2.0 * (x * y + z * w);
    mArr[2] = 2.0 * (x * z - y * w);
    mArr[3] = 2.0 * (x * y - z * w);
    mArr[4] = w2 - x2 + y2 - z2;
    mArr[5] = 2.0 * (y * z + x * w);
    mArr[6] = 2.0 * (x * z + y * w);
    mArr[7] = 2.0 * (y * z - x * w);
    mArr[8] = w2 - x2 - y2 + z2;
    return Matrix(3, 3, mArr);
  }

  /// Convert this [Quaternion] to a rotation [Matrix].
  Matrix toRotationMatrix() => toDirectionCosineMatrix().transpose();

  /// Calculate the angle between this and the [Quaternion] pointing between
  /// the [observer] and [target] positions, given the [forward] vector.
  double vectorAngle(
      final Vector3D observer, final Vector3D target, final Vector3D forward) {
    final delta = target.subtract(observer);
    final transform = toDirectionCosineMatrix().multiplyVector3D(delta);
    return forward.angle(transform);
  }

  /// Calculate the rate of change for this [Quaternion] given an
  /// [angularVelocity] vector _(rad/s)_.
  Quaternion kinematics(final Vector3D angularVelocity) {
    final wPrime = Vector.fromList(
        [0, angularVelocity.x, angularVelocity.y, angularVelocity.z]);
    final qArr = Float64List(16);
    qArr[0] = x;
    qArr[1] = w;
    qArr[2] = -z;
    qArr[3] = y;
    qArr[4] = y;
    qArr[5] = z;
    qArr[6] = w;
    qArr[7] = -x;
    qArr[8] = z;
    qArr[9] = -y;
    qArr[10] = x;
    qArr[11] = w;
    qArr[12] = w;
    qArr[13] = -x;
    qArr[14] = -y;
    qArr[15] = -z;
    final qMat = Matrix(4, 4, qArr);
    final result = qMat.multiplyVector(wPrime).scale(0.5);
    return Quaternion(result[0], result[1], result[2], result[3]);
  }
}
