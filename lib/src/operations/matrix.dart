import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/constants.dart';
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
    final result = array2d(3, 3, 0.0);
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
    final result = array2d(3, 3, 0.0);
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
    final result = array2d(3, 3, 0.0);
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
      Matrix(array2d(rows, columns, 0.0));

  /// Create a new square identity matrix of the provided [dimension].
  factory Matrix.identity(final int dimension) {
    final output = array2d(dimension, dimension, 0.0);
    for (var i = 0; i < dimension; i++) {
      output[i][i] = 1.0;
    }
    return Matrix(output);
  }

  /// Create a new square [Matrix] using elements [d] for the diagonal
  /// components.
  factory Matrix.diagonal(final List<double> d) {
    final output = array2d(d.length, d.length, 0.0);
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
  String toString() => _elements.join('\n');

  /// Copy the elements of this [Matrix] into a [List].
  List<List<double>> toList() {
    final output = array2d(rows, columns, 0.0);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        output[i][j] = this[i][j];
      }
    }
    return output;
  }

  /// Clone this [Matrix] into a new object.
  Matrix clone() => Matrix(toList());

  /// Return the column at the provided [index] as a [Vector].
  Vector column(final int index) {
    final output = Float64List(rows);
    for (var i = 0; i < rows; i++) {
      output[i] = _elements[i][index];
    }
    return Vector(output);
  }

  /// Get a block of this [Matrix], given the start row and column index and
  /// the row and column length.
  Matrix getBlock(
      final int rDex, final int cDex, final int rLen, final int cLen) {
    final output = Matrix.zero(rLen, cLen);
    for (var i = 0; i < rLen; i++) {
      for (var j = 0; j < cLen; j++) {
        output[i][j] = _elements[rDex + i][cDex + j];
      }
    }
    return output;
  }

  /// Set a block of this [Matrix] from input [m], at the provided start row
  /// and column index.
  void setBlock(final Matrix m, final int rDex, final int cDex) {
    for (var i = 0; i < m.rows; i++) {
      for (var j = 0; j < m.columns; j++) {
        _elements[rDex + i][cDex + j] = m[i][j];
      }
    }
  }

  /// Return the result of adding this and another [Matrix].
  Matrix add(final Matrix m) {
    final result = array2d(rows, columns, 0.0);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        result[i][j] = _elements[i][j] + m._elements[i][j];
      }
    }
    return Matrix(result);
  }

  /// Return the result of subtracting this and another [Matrix].
  Matrix subtract(final Matrix m) {
    final result = array2d(rows, columns, 0.0);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        result[i][j] = _elements[i][j] - m._elements[i][j];
      }
    }
    return Matrix(result);
  }

  /// Return a copy of this [Matrix], scaled by [n].
  Matrix scale(final double n) {
    final result = array2d(rows, columns, 0.0);
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
  ///
  /// Throws an error if the column length of this matrix is not equal to
  /// the row length of the argument matrix.
  Matrix multiply(final Matrix m) {
    final rowsA = rows;
    final colsA = columns;
    final rowsB = m.rows;
    final colsB = m.columns;
    if (colsA != rowsB) {
      throw 'Invalid dimensions for matrix multiplication';
    }
    final result = array2d(rowsA, colsB, 0.0);
    for (var i = 0; i < rowsA; i++) {
      for (var j = 0; j < colsB; j++) {
        var sum = 0.0;
        for (var k = 0; k < colsA; k++) {
          sum += _elements[i][k] * m._elements[k][j];
        }
        result[i][j] = sum;
      }
    }
    return Matrix(result);
  }

  /// Return the result of element-wise multiplying this and another [Matrix].
  Matrix outerProduct(final Matrix m) {
    final result = array2d(rows, columns, 0.0);
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
    final output = array2d(rows, columns, 0.0);
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
    final result = array2d(columns, rows, 0.0);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        result[j][i] = _elements[i][j];
      }
    }
    return Matrix(result);
  }

  /// Return the lower-triangular Cholesky decomposed form of this [Matrix];
  Matrix cholesky() {
    final result = array2d(rows, rows, 0.0);
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

  /// Return the Gauss-Jackson inverse of this [Matrix].
  ///
  /// Throws an error if this matrix is non-square.
  Matrix inverse() {
    final n = rows;

    // Check if the matrix is square
    if (n != columns) {
      throw ArgumentError('Matrix must be square for inversion.');
    }

    // Initialize the inverted matrix
    final invertedMatrix = array2d(n, n, 0.0);

    // Augment the matrix with the identity matrix
    final augmentedMatrix = array2d(n, 2 * n, 0.0);
    for (var i = 0; i < n; i++) {
      arraycopy(_elements[i], 0, augmentedMatrix[i], 0, n);
      augmentedMatrix[i][n + i] = 1;
    }

    // Perform Gaussian elimination with partial pivoting
    for (var i = 0; i < n; i++) {
      // Find pivot row
      var maxRowIndex = i;
      for (var k = i + 1; k < n; k++) {
        if (augmentedMatrix[k][i].abs() >
            augmentedMatrix[maxRowIndex][i].abs()) {
          maxRowIndex = k;
        }
      }
      // Swap current row with pivot row
      if (maxRowIndex != i) {
        for (var j = 0; j < 2 * n; j++) {
          final temp = augmentedMatrix[i][j];
          augmentedMatrix[i][j] = augmentedMatrix[maxRowIndex][j];
          augmentedMatrix[maxRowIndex][j] = temp;
        }
      }
      // Make the diagonal element 1
      final divisor = augmentedMatrix[i][i];
      if (divisor == 0) {
        throw ArgumentError('Matrix is singular.');
      }
      for (var j = 0; j < 2 * n; j++) {
        augmentedMatrix[i][j] /= divisor;
      }
      // Make the rest of the column zero
      for (var k = 0; k < n; k++) {
        if (k != i) {
          final factor = augmentedMatrix[k][i];
          for (var j = 0; j < 2 * n; j++) {
            augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
          }
        }
      }
    }

    // Extract the inverted matrix from the augmented matrix
    for (var i = 0; i < n; i++) {
      arraycopy(augmentedMatrix[i], n, invertedMatrix[i], 0, n);
    }

    return Matrix(invertedMatrix);
  }

  /// Attempt to compute the inverse of this matrix using matrix augmentation
  /// to "fix" a singular matrix.
  Matrix inverseSingular() {
    final m = clone();
    for (var i = 0; i < rows; i++) {
      m[i][i] = m[i][i].abs() < machineEpsilon ? machineEpsilon : m[i][i];
    }
    return m.inverse();
  }

  /// Compute the Moore-Penrose pseudoinverse of this matrix.
  Matrix pseudoinverse() {
    final numRows = rows;
    final numCols = columns;
    final at = transpose();

    // calculate A^T * A for wide matrices and A * A^T for tall matrices
    final product = numRows >= numCols ? at.multiply(this) : multiply(at);

    // augment the product matrix to make in non-singular, and invert
    final inverted = product.inverseSingular();

    // compute A^+ = (A^T * A)^-1 * A^T for wide matrices and A^T * (A * A^T)^-1
    // for tall matrices
    return numRows >= numCols ? inverted.multiply(at) : at.multiply(inverted);
  }

  /// Perform LU decomposition on this [Matrix].
  (Matrix l, Matrix u) luDecomposition() {
    if (rows != columns) {
      throw ArgumentError('Matrix must be square for LU decomposition.');
    }
    final n = rows;
    final l = Matrix.zero(rows, columns);
    final u = Matrix.zero(rows, columns);

    for (var i = 0; i < n; i++) {
      l[i][i] = 1.0;
      for (var j = i; j < n; j++) {
        var sum1 = 0.0;
        for (var k = 0; k < i; k++) {
          sum1 += l[i][k] * u[k][j];
        }
        u[i][j] = this[i][j] - sum1;
      }
      for (var j = i + 1; j < n; j++) {
        var sum2 = 0.0;
        for (var k = 0; k < i; k++) {
          sum2 += l[j][k] * u[k][i];
        }
        if (u[i][i] == 0.0) {
          throw ArgumentError(
              'Division by zero encountered in LU decomposition.');
        }
        l[j][i] = (this[j][i] - sum2) / u[i][i];
      }
    }
    return (l, u);
  }

  Vector _forwardSubstitution(final Vector b) {
    if (rows != b.length) {
      throw ArgumentError('The matrix and vector dimensions do not match.');
    }
    final n = rows;
    final y = Vector.zero(b.length);
    for (var i = 0; i < n; i++) {
      var sum = 0.0;
      for (var j = 0; j < i; j++) {
        sum += this[i][j] * y[j];
      }
      if (this[i][i] == 0.0) {
        throw ArgumentError(
            'Division by zero encountered in forward substitution.');
      }
      y[i] = (b[i] - sum) / this[i][i];
    }
    return y;
  }

  Vector _backSubstitution(final Vector y) {
    if (rows != y.length) {
      throw ArgumentError('The matrix and vector dimensions do not match.');
    }
    final n = rows;
    final x = Vector.zero(n);
    for (var i = n - 1; i >= 0; i--) {
      var sum = 0.0;
      for (var j = i + 1; j < n; j++) {
        sum += this[i][j] * x[j];
      }
      if (this[i][i] == 0.0) {
        throw ArgumentError(
            'Division by zero encountered in back substitution.');
      }
      x[i] = (y[i] - sum) / this[i][i];
    }
    return x;
  }

  /// Solve the set of linear equations represented by this [Matrix] for input
  /// vector [b] using LU decomposition.
  Vector solve(final Vector b) {
    if (rows != b.length) {
      throw ArgumentError(
          'The matrix rows and vector length must match to solve.');
    }
    final (l, u) = luDecomposition();
    final y = l._forwardSubstitution(b);
    return u._backSubstitution(y);
  }
}
