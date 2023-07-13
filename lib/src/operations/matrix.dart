import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Matrix operations.
class Matrix {
  /// Create a matrix from a nested array.
  Matrix(this._elements)
      : rows = _elements.length,
        columns = _elements.first.length;

  /// Create a 3x3 x-axis rotation matrix, rotated by angle [theta] _(rad)_.
  factory Matrix.rotX(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final result = array2d(3, 3);
    result[0][0] = 1.0;
    result[1][1] = cosT;
    result[1][2] = sinT;
    result[2][1] = -sinT;
    result[2][2] = cosT;
    return Matrix(result);
  }

  /// Create a 3x3 y-axis rotation matrix, rotated by angle [theta] _(rad)_.
  factory Matrix.rotY(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final result = array2d(3, 3);
    result[0][0] = cosT;
    result[0][2] = -sinT;
    result[1][1] = 1.0;
    result[2][0] = sinT;
    result[2][2] = cosT;
    return Matrix(result);
  }

  /// Create a 3x3 y-axis rotation matrix, rotated by angle [theta] _(rad)_.
  factory Matrix.rotZ(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final result = array2d(3, 3);
    result[0][0] = cosT;
    result[0][1] = sinT;
    result[1][0] = -sinT;
    result[1][1] = cosT;
    result[2][2] = 1.0;
    return Matrix(result);
  }

  /// Create a new zero-filled matrix of the provided [rows] and
  /// [columns] dimensions.
  factory Matrix.zero(final int rows, final int columns) =>
      Matrix(array2d(rows, columns));

  /// Create a new square identity matrix of the provided [dimension].
  factory Matrix.identity(final int dimension) {
    final output = array2d(dimension, dimension);
    for (var i = 0; i < dimension; i++) {
      output[i][i] = 1.0;
    }
    return Matrix(output);
  }

  /// Create a new square [Matrix] using elements [d] for the diagonal
  /// components.
  factory Matrix.diagonal(final List<double> d) {
    final output = array2d(d.length, d.length);
    for (var i = 0; i < d.length; i++) {
      output[i][i] = d[i];
    }
    return Matrix(output);
  }

  /// Matrix elements.
  final List<List<double>> _elements;

  /// Number of rows in this matrix.
  final int rows;

  /// Number of columns in this matrix.
  final int columns;

  /// Return the matrix row at the provided [index].
  List<double> operator [](final int index) => _elements[index];

  @override
  String toString() => _elements.toString();

  /// Return the result of adding this and another [Matrix].
  Matrix add(final Matrix m) {
    final result = array2d(rows, columns);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        result[i][j] = _elements[i][j] + m._elements[i][j];
      }
    }
    return Matrix(result);
  }

  /// Return the result of subtracting this and another [Matrix].
  Matrix subtract(final Matrix m) {
    final result = array2d(rows, columns);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        result[i][j] = _elements[i][j] - m._elements[i][j];
      }
    }
    return Matrix(result);
  }

  /// Return a copy of this [Matrix], scaled by [n].
  Matrix scale(final double n) {
    final result = array2d(rows, columns);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        result[i][j] = _elements[i][j] * n;
      }
    }
    return Matrix(result);
  }

  /// Return a copy of this [Matrix] with all elements negated.
  Matrix negate() => scale(-1);

  /// Return the result of multiplying this by another [Matrix];
  Matrix multiply(final Matrix m) {
    final result = array2d(rows, m.columns);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < m.columns; j++) {
        for (var k = 0; k < columns; k++) {
          result[i][j] += _elements[i][k] * m._elements[k][j];
        }
      }
    }
    return Matrix(result);
  }

  /// Return the result of element-wise multiplying this and another [Matrix].
  Matrix outerProduct(final Matrix m) {
    final result = array2d(rows, columns);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        result[i][j] = _elements[i][j] * m._elements[i][j];
      }
    }
    return Matrix(result);
  }

  /// Return the result of multiplying this with the [Vector] argument.
  Vector multiplyVector(final Vector v) {
    final result = Float64List(rows);
    for (var i = 0; i < rows; i++) {
      var total = 0.0;
      for (var j = 0; j < columns; j++) {
        total += _elements[i][j] * v[j];
      }
      result[i] = total;
    }
    return Vector(result);
  }

  /// Return the result of multiplying this with the [Vector3D] argument.
  Vector3D multiplyVector3D(final Vector3D v) {
    final result = Float64List(3);
    for (var i = 0; i < rows; i++) {
      var total = 0.0;
      for (var j = 0; j < columns; j++) {
        total += _elements[i][j] * v[j];
      }
      result[i] = total;
    }
    return Vector3D(result[0], result[1], result[2]);
  }

  /// Return a copy of this [Matrix] with all elements inverted.
  Matrix reciprocal() {
    final output = array2d(rows, columns);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        if (_elements[i][j] != 0) {
          output[i][j] = 1 / _elements[i][j];
        }
      }
    }
    return Matrix(output);
  }

  /// Return the transpose of this [Matrix].
  Matrix transpose() {
    final result = array2d(columns, rows);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        result[j][i] = _elements[i][j];
      }
    }
    return Matrix(result);
  }

  /// Return the lower-triangular Cholesky decomposed form of this [Matrix];
  Matrix cholesky() {
    final result = array2d(rows, rows);
    for (var i = 0; i < rows; i++) {
      for (var k = 0; k < (i + 1); k++) {
        var total = 0.0;
        for (var j = 0; j < k; j++) {
          total += result[i][j] * result[k][j];
        }
        result[i][k] = (i == k)
            ? sqrt(_elements[i][i] - total)
            : (1 / result[k][k] * (_elements[i][k] - total));
      }
    }
    return Matrix(result);
  }

  /// Swap row [i] and [j] in-place.
  void _swapRows(final int i, final int j) {
    if (i == j) {
      return;
    }
    final tmp = _elements[i];
    _elements[i] = _elements[j];
    _elements[j] = tmp;
  }

  /// Convert this [Matrix] to reduced row echelon form in-place.
  void _toReducedRowEchelonForm() {
    for (var row = 0, lead = 0; row < rows && lead < columns; ++row, ++lead) {
      var i = row;
      while (_elements[i][lead] == 0) {
        if (++i == rows) {
          i = row;
          if (++lead == columns) {
            return;
          }
        }
      }
      _swapRows(i, row);
      if (_elements[row][lead] != 0) {
        final f = _elements[row][lead];
        for (var column = 0; column < columns; ++column) {
          _elements[row][column] /= f;
        }
      }
      for (var j = 0; j < rows; ++j) {
        if (j == row) {
          continue;
        }
        final f = _elements[j][lead];
        for (var column = 0; column < columns; ++column) {
          _elements[j][column] -= f * _elements[row][column];
        }
      }
    }
  }

  /// Return the inverse of this [Matrix] using the Gauss-Jackson method.
  Matrix inverse() {
    final tmp = Matrix.zero(rows, columns * 2);
    for (var row = 0; row < rows; ++row) {
      for (var column = 0; column < columns; ++column) {
        tmp._elements[row][column] = _elements[row][column];
      }
      tmp._elements[row][row + columns] = 1.0;
    }
    tmp._toReducedRowEchelonForm();
    final inv = Matrix.zero(rows, columns);
    for (var row = 0; row < rows; ++row) {
      for (var column = 0; column < columns; ++column) {
        inv._elements[row][column] = tmp._elements[row][column + columns];
      }
    }
    return inv;
  }
}
