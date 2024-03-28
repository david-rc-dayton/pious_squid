import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/operations_base.dart';

/// Polynomial regression optimization result.
class PolynomicalRegressionResult {
  /// Create a new [PolynomicalRegressionResult] object, containing the
  /// polynomial [coefficients], root-sum-squared error [rss], and Bayes
  /// information criterea [bic].
  PolynomicalRegressionResult(this.coefficients, this.rss, this.bic);

  /// Polynomial coefficients.
  Float64List coefficients;

  /// Polynomial fit root-sum-squared.
  double rss;

  /// Polynomial fit Bays information criterea.
  double bic;

  /// Evaluate value [x] using this result's coefficients.
  double evaluate(final double x) => evalPoly(x, coefficients);
}

/// Polynomial regression optimizer.
class PolynomialRegression {
  PolynomialRegression._(); // disable constructor

  static double _bayesInformationCriterea(
          final int n, final int k, final double sse) =>
      n * log(sse) + k * log(n);

  /// Optimize polynomial coefficients to fit data series [xs] and [ys] for the
  /// provided polynomial [order].
  static PolynomicalRegressionResult solve(
      final Float64List xs, final Float64List ys, final int order) {
    final n = xs.length;
    final m = order;

    final xMat = Matrix(n, m + 1);
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < m + 1; j++) {
        xMat.set(i, j, pow(xs[i], j).toDouble());
      }
    }

    final bMat = xMat.solve(Vector.fromList(ys));
    final result = Float64List.fromList(bMat.toArray().reversed.toList());

    var sse = 0.0;
    for (var i = 0; i < xs.length; i++) {
      sse += pow(evalPoly(xs[i], result) - ys[i], 2).toDouble();
    }

    return PolynomicalRegressionResult(
        result, sqrt(sse), _bayesInformationCriterea(xs.length, order, sse));
  }

  /// Optimize polynomial coefficients to fit data series [xs] and [ys], and
  /// attempt to find an optimal order within the [minOrder] and
  /// [maxOrder] bounds.
  static PolynomicalRegressionResult solveOrder(
    final Float64List xs,
    final Float64List ys,
    final int minOrder,
    final int maxOrder,
  ) {
    final cache = <PolynomicalRegressionResult>[];
    for (var order = minOrder; order <= maxOrder; order++) {
      cache.add(solve(xs, ys, order));
    }
    cache.sort((final a, final b) => a.bic.compareTo(b.bic));
    return cache.first;
  }
}
