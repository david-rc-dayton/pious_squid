import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Matrix operations.
class Matrix {
  /// Create a matrix from a nested array.
  Matrix(this.rows, this.columns, [final Float64List? elements])
      : _elements = elements ?? Float64List(rows * columns);

  /// Create a 3x3 x-axis rotation matrix, rotated by angle [theta] _(rad)_.
  factory Matrix.rotX(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final result = Matrix(3, 3);
    result._elements[0] = 1.0;
    result._elements[4] = cosT;
    result._elements[5] = sinT;
    result._elements[7] = -sinT;
    result._elements[8] = cosT;
    return result;
  }

  /// Create a 3x3 y-axis rotation matrix, rotated by angle [theta] _(rad)_.
  factory Matrix.rotY(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final result = Matrix(3, 3);
    result._elements[0] = cosT;
    result._elements[2] = -sinT;
    result._elements[4] = 1.0;
    result._elements[6] = sinT;
    result._elements[8] = cosT;
    return result;
  }

  /// Create a 3x3 y-axis rotation matrix, rotated by angle [theta] _(rad)_.
  factory Matrix.rotZ(final double theta) {
    final cosT = cos(theta);
    final sinT = sin(theta);
    final result = Matrix(3, 3);
    result._elements[0] = cosT;
    result._elements[1] = sinT;
    result._elements[3] = -sinT;
    result._elements[4] = cosT;
    result._elements[8] = 1.0;
    return result;
  }

  /// Create a new square identity matrix of the provided [dimension].
  factory Matrix.identity(final int size) {
    final output = Matrix(size, size);
    for (var i = 0; i < size; i++) {
      output.set(i, i, 1.0);
    }
    return output;
  }

  /// Create a new square [Matrix] using elements [d] for the diagonal
  /// components.
  factory Matrix.diagonal(final List<double> d) {
    final output = Matrix(d.length, d.length);
    for (var i = 0; i < d.length; i++) {
      output.set(i, i, d[i]);
    }
    return output;
  }

  /// Create a new matrix from a 2D array of elements.
  factory Matrix.fromList(final List<List<double>> list) {
    final rows = list.length;
    final cols = list.first.length;
    final elements =
        Float64List.fromList(list.expand((final element) => element).toList());
    return Matrix(rows, cols, elements);
  }

  /// Matrix elements.
  final Float64List _elements;

  /// Number of rows in this matrix.
  final int rows;

  /// Number of columns in this matrix.
  final int columns;

  /// Get the element value at the provided [row] and [column];
  double get(final int row, final int column) {
    if (row >= rows || column >= columns || row < 0 || column < 0) {
      throw 'Index out of bounds';
    }
    return _elements[row * columns + column];
  }

  /// Set the element [value] at the provided [row] and [column].
  void set(final int row, final int column, final double value) {
    if (row >= rows || column >= columns || row < 0 || column < 0) {
      throw 'Index out of bounds';
    }
    _elements[row * columns + column] = value;
  }

  @override
  String toString() => toList().join('\n');

  /// Return `true` if this [Matrix] is square.
  bool isSquare() => rows == columns;

  /// Copy the elements of this [Matrix] into a [List].
  List<List<double>> toList() {
    final arr = <List<double>>[];
    for (var i = 0; i < rows; i++) {
      arr.add(List.from(_elements.sublist(i * columns, (i + 1) * columns)));
    }
    return arr;
  }

  /// Clone this [Matrix] into a new object.
  Matrix clone() => Matrix(rows, columns, _elements.sublist(0));

  /// Convert the row at the provided [index] into a [Vector].
  Vector row(final int index) {
    final output = Float64List(columns);
    for (var i = 0; i < columns; i++) {
      output[i] = get(index, i);
    }
    return Vector(output);
  }

  /// Return the column at the provided [index] as a [Vector].
  Vector column(final int index) {
    final output = Float64List(rows);
    for (var i = 0; i < rows; i++) {
      output[i] = get(i, index);
    }
    return Vector(output);
  }

  /// Get a block of this [Matrix], given the start row and column index and
  /// the row and column length.
  Matrix getBlock(
      final int rDex, final int cDex, final int rLen, final int cLen) {
    final output = Matrix(rLen, cLen);
    for (var i = 0; i < rLen; i++) {
      for (var j = 0; j < cLen; j++) {
        output.set(i, j, get(rDex + i, cDex + j));
      }
    }
    return output;
  }

  /// Set a block of this [Matrix] from input [m], at the provided start row
  /// and column index.
  void setBlock(final Matrix m, final int rDex, final int cDex) {
    for (var i = 0; i < m.rows; i++) {
      for (var j = 0; j < m.columns; j++) {
        set(rDex + i, cDex + j, m.get(i, j));
      }
    }
  }

  /// Return the result of adding this and another [Matrix].
  Matrix add(final Matrix m) {
    if (rows != m.rows || columns != m.columns) {
      throw 'Matrix dimensions must match.';
    }
    final result = Matrix(rows, columns);
    for (var i = 0; i < _elements.length; i++) {
      result._elements[i] = _elements[i] + m._elements[i];
    }
    return result;
  }

  /// Return the result of subtracting this and another [Matrix].
  Matrix subtract(final Matrix m) {
    if (rows != m.rows || columns != m.columns) {
      throw 'Matrix dimensions must match.';
    }
    final result = Matrix(rows, columns);
    for (var i = 0; i < _elements.length; i++) {
      result._elements[i] = _elements[i] - m._elements[i];
    }
    return result;
  }

  /// Return a copy of this [Matrix], scaled by [n].
  Matrix scale(final double n) {
    final result = Matrix(rows, columns);
    for (var i = 0; i < _elements.length; i++) {
      result._elements[i] = _elements[i] * n;
    }
    return result;
  }

  /// Return a copy of this [Matrix] with all elements negated.
  Matrix negate() => scale(-1);

  /// Return the result of multiplying this by another [Matrix];
  ///
  /// Throws an error if the column length of this matrix is not equal to
  /// the row length of the argument matrix.
  Matrix multiply(final Matrix m) {
    if (columns != m.rows) {
      throw 'Matrix multiplication dimensions do not match.';
    }
    final result = Matrix(rows, m.columns, Float64List(rows * m.columns));
    for (var i = 0; i < rows; i++) {
      for (var k = 0; k < columns; k++) {
        for (var j = 0; j < m.columns; j++) {
          result._elements[i * m.columns + j] += get(i, k) * m.get(k, j);
        }
      }
    }
    return result;
  }

  /// Return the result of element-wise multiplying this and another [Matrix].
  Matrix outerProduct(final Matrix m) {
    final result = Matrix(rows, columns);
    for (var i = 0; i < _elements.length; i++) {
      result._elements[i] = _elements[i] * m._elements[i];
    }
    return result;
  }

  /// Return the result of multiplying this with the [Vector] argument.
  Vector multiplyVector(final Vector v) {
    final result = Float64List(rows);
    for (var i = 0; i < rows; i++) {
      var total = 0.0;
      for (var j = 0; j < columns; j++) {
        total += get(i, j) * v[j];
      }
      result[i] = total;
    }
    return Vector(result);
  }

  /// Return the result of multiplying this with the [Vector3D] argument.
  Vector3D multiplyVector3D(final Vector3D v) {
    if (rows != 3 || columns != 3) {
      throw 'Invalid dimensions for multiplication';
    }
    final x = _elements[0] * v.x + _elements[1] * v.y + _elements[2] * v.z;
    final y = _elements[3] * v.x + _elements[4] * v.y + _elements[5] * v.z;
    final z = _elements[6] * v.x + _elements[7] * v.y + _elements[8] * v.z;
    return Vector3D(x, y, z);
  }

  /// Return a copy of this [Matrix] with all elements inverted.
  Matrix reciprocal() {
    final n = rows * columns;
    final output = Float64List(n);
    for (var i = 0; i < n; i++) {
      if (_elements[i] != 0.0) {
        output[i] = 1.0 / _elements[i];
      }
    }
    return Matrix(rows, columns, output);
  }

  /// Return the transpose of this [Matrix].
  Matrix transpose() {
    final output = Float64List(rows * columns);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < columns; j++) {
        output[j * rows + i] = get(i, j);
      }
    }
    return Matrix(columns, rows, output);
  }

  /// Return the lower-triangular Cholesky decomposed form of this [Matrix];
  Matrix cholesky() {
    if (!isSquare()) {
      throw 'Matrix must be square';
    }

    final n = rows;
    final L = Matrix(n, n, Float64List(n * n));

    for (var i = 0; i < n; i++) {
      for (var j = 0; j <= i; j++) {
        var sum = 0.0;
        if (j == i) {
          for (var k = 0; k < j; k++) {
            sum += L._elements[j * n + k] * L._elements[j * n + k];
          }
          L._elements[j * n + j] = sqrt(_elements[j * n + j] - sum);
        } else {
          for (var k = 0; k < j; k++) {
            sum += L._elements[i * n + k] * L._elements[j * n + k];
          }
          if (L._elements[j * n + j] == 0) {
            throw 'Matrix is not positive definite';
          }
          L._elements[i * n + j] =
              (_elements[i * n + j] - sum) / L._elements[j * n + j];
        }
      }
    }

    return L;
  }

  /// Return the Gauss-Jackson inverse of this [Matrix].
  ///
  /// Throws an error if this matrix is non-square.
  Matrix inverse() {
    final n = rows;

    // Check if the matrix is square
    if (!isSquare()) {
      throw 'Matrix must be square for inversion.';
    }

    // Augment the matrix with the identity matrix
    final augmentedMatrix = Matrix(n, 2 * n);
    augmentedMatrix.setBlock(this, 0, 0);
    for (var i = 0; i < n; i++) {
      augmentedMatrix._elements[i * augmentedMatrix.columns + n + i] = 1.0;
    }

    // Perform Gaussian elimination with partial pivoting
    for (var i = 0; i < n; i++) {
      // Find pivot row
      var maxRowIndex = i;
      for (var k = i + 1; k < n; k++) {
        if (augmentedMatrix.get(k, i).abs() >
            augmentedMatrix.get(maxRowIndex, i).abs()) {
          maxRowIndex = k;
        }
      }
      // Swap current row with pivot row
      if (maxRowIndex != i) {
        for (var j = 0; j < 2 * n; j++) {
          final temp = augmentedMatrix.get(i, j);
          augmentedMatrix.set(i, j, augmentedMatrix.get(maxRowIndex, j));
          augmentedMatrix.set(maxRowIndex, j, temp);
        }
      }
      // Make the diagonal element 1
      final divisor = augmentedMatrix.get(i, i);
      if (divisor == 0) {
        throw 'Matrix is singular.';
      }
      for (var j = 0; j < 2 * n; j++) {
        augmentedMatrix._elements[i * augmentedMatrix.columns + j] /= divisor;
      }
      // Make the rest of the column zero
      for (var k = 0; k < n; k++) {
        if (k != i) {
          final factor = augmentedMatrix.get(k, i);
          for (var j = 0; j < 2 * n; j++) {
            augmentedMatrix._elements[k * augmentedMatrix.columns + j] -=
                factor * augmentedMatrix.get(i, j);
          }
        }
      }
    }

    return augmentedMatrix.getBlock(0, n, n, n);
  }

  /// Augment this matrix along th diagonal to "fix" a sigular matrix.
  ///
  /// Throws an error if the matrix is not square.
  Matrix augmented() {
    if (!isSquare()) {
      throw 'Only square matrices can be augmented';
    }
    final m = clone();
    for (var i = 0; i < rows; i++) {
      final value = get(i, i);
      m.set(i, i, value.abs() < machineEpsilon ? machineEpsilon : value);
    }
    return m;
  }

  /// Attempt to compute the inverse of this matrix using matrix augmentation
  /// to "fix" a singular matrix.
  Matrix inverseSingular() => augmented().inverse();

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
    if (!isSquare()) {
      throw 'Matrix must be square for LU decomposition.';
    }
    final n = rows;
    final l = Matrix(rows, columns);
    final u = Matrix(rows, columns);

    for (var i = 0; i < n; i++) {
      l.set(i, i, 1.0);
      for (var j = i; j < n; j++) {
        var sum1 = 0.0;
        for (var k = 0; k < i; k++) {
          sum1 += l.get(i, k) * u.get(k, j);
        }
        u.set(i, j, get(i, j) - sum1);
      }
      for (var j = i + 1; j < n; j++) {
        var sum2 = 0.0;
        for (var k = 0; k < i; k++) {
          sum2 += l.get(j, k) * u.get(k, i);
        }
        if (u.get(i, i) == 0.0) {
          throw 'Division by zero encountered in LU decomposition.';
        }
        l.set(j, i, (get(j, i) - sum2) / u.get(i, i));
      }
    }
    return (l, u);
  }

  Vector _forwardSubstitution(final Vector b) {
    if (rows != b.length) {
      throw 'The matrix and vector dimensions do not match.';
    }
    final n = rows;
    final y = Vector.zero(b.length);
    for (var i = 0; i < n; i++) {
      var sum = 0.0;
      for (var j = 0; j < i; j++) {
        sum += get(i, j) * y[j];
      }
      if (get(i, i) == 0) {
        throw 'Division by zero encountered in forward substitution.';
      }
      y[i] = (b[i] - sum) / get(i, i);
    }
    return y;
  }

  Vector _backSubstitution(final Vector y) {
    if (rows != y.length) {
      throw 'The matrix and vector dimensions do not match.';
    }
    final n = rows;
    final x = Vector.zero(n);
    for (var i = n - 1; i >= 0; i--) {
      var sum = 0.0;
      for (var j = i + 1; j < n; j++) {
        sum += get(i, j) * x[j];
      }
      if (get(i, i) == 0.0) {
        throw 'Division by zero encountered in back substitution.';
      }
      x[i] = (y[i] - sum) / get(i, i);
    }
    return x;
  }

  /// Compute the dot product of a [column] of this [Matrix] and vector [a].
  double dotColumn(final Vector a, final int column) {
    var dot = 0.0;
    for (var i = 0; i < a.length; i++) {
      dot += get(i, column) * a[i];
    }
    return dot;
  }

  /// Compute the QR decomposition of this [Matrix].
  (Matrix q, Matrix r) qrDecomposition() {
    final m = rows;
    final n = columns;

    final q = Matrix(m, n);
    final r = Matrix(n, n);

    for (var j = 0; j < n; j++) {
      final a = Vector.zero(m);
      for (var i = 0; i < m; i++) {
        a[i] = get(i, j);
      }

      for (var i = 0; i < j; i++) {
        final dot = q.dotColumn(a, i);
        for (var k = 0; k < m; k++) {
          a[k] -= dot * q.get(k, i);
        }
        r.set(i, j, dot);
      }

      final norm = a.magnitude();
      for (var i = 0; i < m; i++) {
        q.set(i, j, a[i] / norm);
      }
      r.set(j, j, norm);
    }

    return (q, r);
  }

  /// Solve the set of linear equations represented by this [Matrix] for input
  /// vector [b] using LU decomposition.
  Vector solve(final Vector b) {
    if (rows != b.length) {
      throw 'The matrix rows and vector length must match to solve.';
    }
    // solve with lu decomposition if square
    if (isSquare()) {
      final (l, u) = luDecomposition();
      final y = l._forwardSubstitution(b);
      return u._backSubstitution(y);
    }
    // solve with qr decomposition if rectangular
    final (q, r) = qrDecomposition();
    final m = rows;
    final n = columns;
    final y = Vector.zero(n);
    for (var i = 0; i < n; i++) {
      y[i] = 0.0;
      for (var j = 0; j < m; j++) {
        final qji = q.get(j, i);
        final qv = qji.isFinite && qji != 0 ? qji : machineEpsilon;
        y[i] += qv * b[j];
      }
    }
    final x = Vector.zero(n);
    for (var i = x.length - 1; i >= 0; i--) {
      x[i] = y[i];
      for (var j = i + 1; j < x.length; j++) {
        final rij = r.get(i, j);
        final rv = rij.isFinite && rij != 0 ? rij : machineEpsilon;
        x[i] -= rv * x[j];
      }
      x[i] /= r.get(i, i);
    }

    return x;
  }
}
