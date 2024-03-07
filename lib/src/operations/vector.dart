import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/matrix.dart';

/// 3-dimensional vector.
class Vector3D {
  //// Create a new [Vector3D] object.
  const Vector3D(this.x, this.y, this.z);

  /// Create a new [Vector3D] object from the first three elements of a
  /// [Vector] object.
  factory Vector3D.fromVector(final Vector v) => Vector3D(v[0], v[1], v[2]);

  /// X-axis component.
  final double x;

  /// Y-axis component.
  final double y;

  /// Z-axis component.
  final double z;

  /// Origin vector.
  static const origin = Vector3D(0, 0, 0);

  /// X-axis unit vector.
  static const xAxis = Vector3D(1, 0, 0);

  /// Y-axis unit vector.
  static const yAxis = Vector3D(0, 1, 0);

  /// Z-axis unit vector.
  static const zAxis = Vector3D(0, 0, 1);

  /// Negative x-axis unit vector.
  static const xAxisNeg = Vector3D(-1, 0, 0);

  /// Negative y-axis unit vector.
  static const yAxisNeg = Vector3D(0, -1, 0);

  /// Negative z-axis unit vector.
  static const zAxisNeg = Vector3D(0, 0, -1);

  /// Convert this to a [List] of doubles.
  List<double> toList() => [x, y, z];

  /// Convert this to a [Float64List] object.
  Float64List toArray() {
    final output = Float64List(3);
    output[0] = x;
    output[1] = y;
    output[2] = z;
    return output;
  }

  /// Return the [Vector3D] element at the provided [index].
  double operator [](final int index) {
    switch (index) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        throw 'Index $index outside 3D vector bounds.';
    }
  }

  /// Convert this to a [Vector] object.
  Vector toVector() => Vector(toArray());

  @override
  String toString([final int fixed = -1]) {
    if (fixed < 0) {
      return '[${toList().join(', ')}]';
    }
    final output = toList().map((final e) => e.toStringAsFixed(fixed));
    return '[${output.join(", ")}]';
  }

  /// Return the magnitude of this vector.
  double magnitude() => sqrt(x * x + y * y + z * z);

  /// Return the result of adding this to another [Vector3D].
  Vector3D add(final Vector3D v) => Vector3D(x + v.x, y + v.y, z + v.z);

  /// Return the result of subtracting this and another [Vector3D].
  Vector3D subtract(final Vector3D v) => Vector3D(x - v.x, y - v.y, z - v.z);

  /// Return a copy of this [Vector3D] scaled by [n];
  Vector3D scale(final double n) => Vector3D(x * n, y * n, z * n);

  /// Return a copy of this [Vector3D] with the elements negated.
  Vector3D negate() => scale(-1);

  /// Return the Euclidean distance between this and another [Vector3D].
  double distance(final Vector3D v) {
    final dx = x - v.x;
    final dy = y - v.y;
    final dz = z - v.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
  }

  /// Convert this to a unit [Vector3D].
  Vector3D normalize() {
    final m = magnitude();
    if (m == 0) {
      return Vector3D.origin;
    }
    return Vector3D(x / m, y / m, z / m);
  }

  /// Calculate the dot product of this and another [Vector3D].
  double dot(final Vector3D v) => x * v.x + y * v.y + z * v.z;

  /// Calculate the outer product between this and another [Vector3D].
  Matrix outer(final Vector3D v) => Matrix([
        [x * v.x, x * v.y, x * v.z],
        [y * v.x, y * v.y, y * v.z],
        [z * v.x, z * v.y, z * v.z],
      ]);

  /// Calculate the cross product of this and another [Vector3D].
  Vector3D cross(final Vector3D v) =>
      Vector3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);

  /// Calculate the skew-symmetric matrix for this [Vector3D].
  Matrix skewSymmetric() => Matrix([
        [0, -z, y],
        [z, 0, -x],
        [-y, x, 0]
      ]);

  /// Create a copy of this [Vector3D] rotated in the x-axis by angle
  /// [theta] _(rad)_.
  Vector3D rotX(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    return Vector3D(x, cosT * y + sinT * z, -sinT * y + cosT * z);
  }

  /// Create a copy of this [Vector3D] rotated in the y-axis by angle
  /// [theta] _(rad)_.
  Vector3D rotY(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    return Vector3D(cosT * x + -sinT * z, y, sinT * x + cosT * z);
  }

  /// Create a copy of this [Vector3D] rotated in the z-axis by angle
  /// [theta] _(rad)_.
  Vector3D rotZ(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    return Vector3D(cosT * x + sinT * y, -sinT * x + cosT * y, z);
  }

  /// Calculate the angle _(rad)_ between this and another [Vector3D].
  double angle(final Vector3D v) {
    final theta = atan2(cross(v).magnitude(), dot(v));
    return theta.isNaN ? 0.0 : theta;
  }

  /// Calculate the angle _(°)_ between this and another [Vector3D].
  double angleDegrees(final Vector3D v) => angle(v) * rad2deg;

  /// Return `true` if line-of-sight exists between this and another [Vector3D]
  /// with a central body of the given [radius].
  bool sight(final Vector3D v, final double radius) {
    final r1Mag2 = pow(magnitude(), 2);
    final r2Mag2 = pow(v.magnitude(), 2);
    final rDot = dot(v);
    final tMin = (r1Mag2 - rDot) / (r1Mag2 + r2Mag2 - 2.0 * rDot);
    var los = false;
    if (tMin < 0 || tMin > 1) {
      los = true;
    } else {
      final c = (1.0 - tMin) * r1Mag2 + rDot * tMin;
      if (c >= radius * radius) {
        los = true;
      }
    }
    return los;
  }

  /// Return the unit vector that bisects this and another [Vector3D].
  Vector3D bisect(final Vector3D v) =>
      scale(v.magnitude()).add(v.scale(magnitude())).normalize();

  /// Convert this [Vector3D] into a row [Matrix].
  Matrix row() => Matrix([
        [x, y, z]
      ]);

  /// Convert this [Vector3D] into a column [Matrix].
  Matrix column() => Matrix([
        [x],
        [y],
        [z]
      ]);

  /// Join this and another [Vector3D] into a new [Vector] object.
  Vector join(final Vector3D v) {
    final output = Float64List(6);
    output[0] = x;
    output[1] = y;
    output[2] = z;
    output[3] = v.x;
    output[4] = v.y;
    output[5] = v.z;
    return Vector(output);
  }
}

/// Vector operations.
class Vector {
  /// Create a new [Vector] object from an array of [_elements].
  const Vector(this._elements) : length = _elements.length;

  /// Create a zero-filled [Vector] of the provided [length];
  factory Vector.zero(final int length) => Vector(Float64List(length));

  /// Create a [Vector] of the provided [length], filled with the
  /// provided [value].
  factory Vector.filled(final int length, final double value) {
    final output = Vector.zero(length);
    for (var i = 0; i < length; i++) {
      output._elements[i] = value;
    }
    return output;
  }

  /// Create a [Vector] from the provided [elements] list.
  factory Vector.fromList(final List<double> elements) =>
      Vector(Float64List.fromList(elements));

  /// Vector elements.
  final Float64List _elements;

  /// Vector length.
  final int length;

  /// 3-dimensional origin.
  static final origin3 = Vector(Float64List.fromList([0, 0, 0]));

  /// 6-dimensional origin.
  static final origin6 = Vector(Float64List.fromList([0, 0, 0, 0, 0, 0]));

  /// X-axis unit vector.
  static final xAxis = Vector(Float64List.fromList([1, 0, 0]));

  /// Y-axis unit vector.
  static final yAxis = Vector(Float64List.fromList([0, 1, 0]));

  /// Z-axis unit vector.
  static final zAxis = Vector(Float64List.fromList([0, 0, 1]));

  /// Negative x-axis unit vector.
  static final xAxisNeg = Vector.fromList([-1, 0, 0]);

  /// Negative y-axis unit vector.
  static final yAxisNeg = Vector.fromList([0, -1, 0]);

  /// Negative z-axis unit vector.
  static final zAxisNeg = Vector.fromList([0, 0, -1]);

  @override
  String toString([final int fixed = -1]) {
    if (fixed < 0) {
      return '[${_elements.join(', ')}]';
    }
    final output =
        _elements.toList().map((final e) => e.toStringAsFixed(fixed));
    return '[${output.join(", ")}]';
  }

  /// X-axis component.
  double get x => _elements[0];

  /// Y-axis component.
  double get y => _elements[1];

  /// Z-axis component.
  double get z => _elements[2];

  /// Return the [Vector] element at the provided [index].
  double operator [](final int index) => _elements[index];

  /// Set [Vector] element [value] at the provided [index].
  void operator []=(final int index, final double value) =>
      _elements[index] = value;

  /// Convert the elements of this [Vector] to a list object.
  List<double> toList({final bool growable = false}) =>
      _elements.toList(growable: true);

  /// Copy the elements of this [Vector] to a new array.
  Float64List toArray() => _elements.sublist(0);

  /// Return the magnitude of this vector.
  double magnitude() {
    var total = 0.0;
    for (final x in _elements) {
      total += x * x;
    }
    return sqrt(total);
  }

  /// Return the result of adding this to another [Vector].
  Vector add(final Vector v) {
    final output = Float64List(length);
    for (var i = 0; i < length; i++) {
      output[i] = _elements[i] + v._elements[i];
    }
    return Vector(output);
  }

  /// Return the result of subtracting this and another [Vector].
  Vector subtract(final Vector v) {
    final output = Float64List(length);
    for (var i = 0; i < length; i++) {
      output[i] = _elements[i] - v._elements[i];
    }
    return Vector(output);
  }

  /// Return a copy of this [Vector] scaled by [n];
  Vector scale(final double n) {
    final output = Float64List(length);
    for (var i = 0; i < length; i++) {
      output[i] = _elements[i] * n;
    }
    return Vector(output);
  }

  /// Return a copy of this [Vector] with the elements negated.
  Vector negate() => scale(-1);

  /// Return the Euclidean distance between this and another [Vector].
  double distance(final Vector v) => subtract(v).magnitude();

  /// Convert this to a unit [Vector].
  Vector normalize() {
    final m = magnitude();
    if (m == 0) {
      return Vector.zero(length);
    }
    return scale(1.0 / m);
  }

  /// Calculate the dot product of this and another [Vector];
  double dot(final Vector v) {
    var total = 0.0;
    for (var i = 0; i < length; i++) {
      total += _elements[i] * v._elements[i];
    }
    return total;
  }

  /// Calculate the outer product between this and another [Vector].
  Matrix outer(final Vector v) {
    final result = array2d(length, v.length, 0.0);
    for (var i = 0; i < length; i++) {
      for (var j = 0; j < v.length; j++) {
        result[i][j] = _elements[i] * v._elements[j];
      }
    }
    return Matrix(result);
  }

  /// Calculate the cross product of this and another [Vector];
  Vector cross(final Vector v) {
    final output = Float64List(length);
    for (var i = 0; i < length; i++) {
      output[i] = _elements[(i + 1) % length] * v._elements[(i + 2) % length] -
          _elements[(i + 2) % length] * v._elements[(i + 1) % length];
    }
    return Vector(output);
  }

  /// Calculate the skew-symmetric matrix for this [Vector].
  ///
  /// An error will be thrown if the vector is not length 3.
  Matrix skewSymmetric() {
    if (length != 3) {
      throw 'Skew-symmetric matrix requires a vector of length 3.';
    }
    return Matrix([
      [0, -_elements[2], _elements[1]],
      [_elements[2], 0, -_elements[0]],
      [-_elements[1], _elements[0], 0]
    ]);
  }

  /// Create a copy of this [Vector] rotated in the x-axis by angle
  /// [theta] _(rad)_.
  Vector rotX(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final output = Float64List(3);
    output[0] = _elements[0];
    output[1] = cosT * _elements[1] + sinT * _elements[2];
    output[2] = -sinT * _elements[1] + cosT * _elements[2];
    return Vector(output);
  }

  /// Create a copy of this [Vector] rotated in the y-axis by angle
  /// [theta] _(rad)_.
  Vector rotY(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final output = Float64List(3);
    output[0] = cosT * _elements[0] + -sinT * _elements[2];
    output[1] = _elements[1];
    output[2] = sinT * _elements[0] + cosT * _elements[2];
    return Vector(output);
  }

  /// Create a copy of this [Vector] rotated in the z-axis by angle
  /// [theta] _(rad)_.
  Vector rotZ(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final output = Float64List(3);
    output[0] = cosT * _elements[0] + sinT * _elements[1];
    output[1] = -sinT * _elements[0] + cosT * _elements[1];
    output[2] = _elements[2];
    return Vector(output);
  }

  /// Calculate the angle _(rad)_ between this and another [Vector].
  double angle(final Vector v) {
    // better than acos for small angles
    final theta = atan2(cross(v).magnitude(), dot(v));
    if (theta.isNaN) {
      return 0.0;
    }
    return theta;
  }

  /// Calculate the angle _(°)_ between this and another [Vector].
  double angleDegrees(final Vector v) => angle(v) * rad2deg;

  /// Return `true` if line-of-sight exists between this and another [Vector]
  /// with a central body of the given [radius].
  bool sight(final Vector v, final double radius) {
    final r1Mag2 = pow(magnitude(), 2);
    final r2Mag2 = pow(v.magnitude(), 2);
    final rDot = dot(v);
    final tMin = (r1Mag2 - rDot) / (r1Mag2 + r2Mag2 - 2.0 * rDot);
    var los = false;
    if (tMin < 0 || tMin > 1) {
      los = true;
    } else {
      final c = (1.0 - tMin) * r1Mag2 + rDot * tMin;
      if (c >= radius * radius) {
        los = true;
      }
    }
    return los;
  }

  /// Return the unit vector that bisects this and another [Vector].
  Vector bisect(final Vector v) =>
      scale(v.magnitude()).add(v.scale(magnitude())).normalize();

  /// Create a new [Vector] combining the elements of this and
  /// another [Vector].
  Vector join(final Vector v) =>
      Vector(Float64List.fromList(toList() + v.toList()));

  /// Create a new [Vector] containing a subset of this object's elements,
  /// from index [start] to [end] _(exclusive)_.
  Vector slice(final int start, final int end) =>
      Vector(_elements.sublist(start, end));

  /// Convert this [Vector] into a row [Matrix].
  Matrix row() => Matrix([_elements.toList()]);

  /// Convert this [Vector] into a column [Matrix].
  Matrix column() => Matrix(_elements.map((final e) => [e]).toList());

  /// Convert elements from this to a [Vector3D], starting at the
  /// provided [index].
  Vector3D toVector3D(final int index) =>
      Vector3D(_elements[index], _elements[index + 1], _elements[index + 2]);

  /// Convert this [Vector] into a row matrix.
  Matrix toRowMatrix() {
    final output = Matrix.zero(1, length);
    for (var i = 0; i < length; i++) {
      output[0][i] = _elements[i];
    }
    return output;
  }

  /// Convert this [Vector] into a column matrix.
  Matrix toColumnMatrix() {
    final output = Matrix.zero(length, 1);
    for (var i = 0; i < length; i++) {
      output[i][0] = _elements[i];
    }
    return output;
  }
}
