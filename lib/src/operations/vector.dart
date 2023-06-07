import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/matrix.dart';

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

  /// Calculate the cross product of this and another [Vector];
  Vector cross(final Vector v) {
    final output = Float64List(length);
    for (var i = 0; i < length; i++) {
      output[i] = _elements[(i + 1) % length] * v._elements[(i + 2) % length] -
          _elements[(i + 2) % length] * v._elements[(i + 1) % length];
    }
    return Vector(output);
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

  /// Convert this [Vector] into a column [Matrix].
  Matrix column() => Matrix(_elements.map((final e) => [e]).toList());
}
