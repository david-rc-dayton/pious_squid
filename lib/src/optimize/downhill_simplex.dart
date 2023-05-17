import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/operations_base.dart';

/// Optimization cost function that takes in parameter array [xs] and returns
/// a score.
typedef CostFunction = double Function(Float64List xs);

/// Simplex coordinate and score data.
class SimplexEntry {
  /// Create a new [SimplexEntry] object, with score function [_f] and
  /// an array of simplex [points].
  SimplexEntry(this._f, final Float64List points) {
    _x = Vector(points);
    _score = _f(points);
  }

  /// Cost function.
  final CostFunction _f;

  /// Cost function result for the simplex point.
  late final double _score;

  /// Values for the simplex point.
  late final Vector _x;

  /// Simplex point element array.
  Float64List get points => _x.toArray();

  /// Cost function result.
  double get score => _score;

  /// Modify the simplex entry in accordance with the Nelder-Mead algorithm.
  SimplexEntry modify(
          final double n, final SimplexEntry xa, final SimplexEntry xb) =>
      SimplexEntry(_f, _x.add(xa._x.add(xb._x.scale(-1.0)).scale(n)).toArray());

  /// Calculate the Euclidean distance between this and another [SimplexEntry];
  double distance(final SimplexEntry se) => _x.distance(se._x);
}

/// Derivative-free Nelder-Mead simplex optimizer.
class DownhillSimplex {
  DownhillSimplex._(); // disable constructor

  /// Compute the centroid from a list of [SimplexEntry] objects, using cost
  /// function [f].
  static SimplexEntry centroid(
      final CostFunction f, final List<SimplexEntry> xss) {
    final n = xss[0].points.length;
    final m = xss.length - 1;
    final output = Float64List(n);
    for (var i = 0; i < m; i++) {
      for (var j = 0; j < n; j++) {
        output[j] += xss[i].points[j];
      }
    }
    for (var i = 0; i < n; i++) {
      output[i] /= m;
    }
    return SimplexEntry(f, output);
  }

  static void _shrink(final double s, final List<SimplexEntry> xss) {
    final x1 = xss[0];
    for (var i = 1; i < xss.length; i++) {
      final xi = xss[i];
      xss[i] = x1.modify(s, xi, x1);
    }
  }

  /// Generate a new simplex from initial guess [x0], and an optional
  /// simples [step] value.
  static List<Float64List> generateSimplex(final Float64List x0,
      {final double step = 0.01}) {
    final output = [x0.sublist(0)];
    for (var i = 0; i < x0.length; i++) {
      final tmp = x0.sublist(0);
      tmp[i] = tmp[i] + (tmp[i] * step);
      output.add(tmp);
    }
    return output;
  }

  /// Perform derivative-free Nelder-Mead simplex optimization to minimize the
  /// cost function [f] for the initial simplex [xs].
  ///
  /// Optional arguments:
  ///   - `xTolerance`: centroid delta termination criteria
  ///   - `fTolerance`: cost function delta termination criteria
  ///   - `maxIter`: maximum number of optimization iterations
  ///   - `adaptive`: use adaptive coefficients if possible
  ///   - `printIter`: print a debug statement after each iteration
  static Float64List solveSimplex(
      final CostFunction f, final List<Float64List> xs,
      {final double xTolerance = 1e-12,
      final double fTolerance = 1e-12,
      final int maxIter = 10000,
      final bool adaptive = false,
      final bool printIter = false}) {
    double a;
    double g;
    double p;
    double s;
    final n = xs.length - 1;
    if (adaptive && n >= 2) {
      a = 1.0;
      g = 1.0 + (2.0 / n);
      p = 0.75 - (1.0 / (2.0 * n));
      s = 1.0 - (1.0 / n);
    } else {
      a = 1.0;
      g = 2.0;
      p = 0.5;
      s = 0.5;
    }
    var iter = 0;
    var action = 'init';
    final ordered = <SimplexEntry>[];
    for (final x in xs) {
      ordered.add(SimplexEntry(f, x));
    }
    while (true) {
      ordered.sort((final x, final y) => x.score.compareTo(y.score));
      final x0 = centroid(f, ordered);
      // update exit criterea
      var xd = 0.0;
      var fd = 0.0;
      for (var i = 1; i < ordered.length; i++) {
        xd = max(xd, x0.distance(ordered[i]));
        fd = max(fd, (x0.score - ordered[i].score).abs());
      }
      if (printIter) {
        print('$iter: score=${x0.score} xd=$xd fd=$fd [$action]');
      }
      if (iter != 0 && (xd < xTolerance || fd < fTolerance)) {
        return ordered.first.points;
      }
      if (iter >= maxIter) {
        return ordered.first.points;
      }
      iter++;
      // reflection
      final xr = x0.modify(a, x0, ordered.last);
      if (ordered.first.score <= xr.score &&
          xr.score < ordered[ordered.length - 2].score) {
        ordered.last = xr;
        action = 'reflect';
        continue;
      }
      // expansion
      if (xr.score < ordered.first.score) {
        final xe = x0.modify(g, xr, x0);
        if (xe.score < xr.score) {
          ordered.last = xe;
        } else {
          ordered.last = xr;
        }
        action = 'expand';
        continue;
      }
      // contraction
      if (xr.score < ordered.last.score) {
        final xc = x0.modify(p, xr, x0);
        if (xc.score < xr.score) {
          ordered.last = xc;
          action = 'contract';
          continue;
        } else {
          _shrink(s, ordered);
          action = 'shrink';
          continue;
        }
      } else if (xr.score >= ordered.last.score) {
        final xc = x0.modify(p, ordered.last, x0);
        if (xc.score < ordered.last.score) {
          ordered.last = xc;
          action = 'contract';
          continue;
        } else {
          _shrink(s, ordered);
          action = 'shrink';
          continue;
        }
      }
    }
  }
}
