import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/optimize/optimize_base.dart';

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
      final Float64List xs, final Float64List ys, final int order,
      {final bool printIter = false}) {
    final simplex = DownhillSimplex.generateSimplex(
        Vector.filled(order + 1, 1.0).toArray());
    double f(final Float64List coeffs) {
      var sse = 0.0;
      for (var i = 0; i < xs.length; i++) {
        final diff = ys[i] - evalPoly(xs[i], coeffs);
        sse += diff * diff;
      }
      return sse;
    }

    final result = DownhillSimplex.solveSimplex(f, simplex,
        adaptive: true, printIter: printIter);
    final sse = f(result);
    return PolynomicalRegressionResult(
        result, sqrt(sse), _bayesInformationCriterea(xs.length, order, sse));
  }

  /// Optimize polynomial coefficients to fit data series [xs] and [ys], and
  /// attempt to find an optimal order within the [minOrder] and
  /// [maxOrder] bounds.
  static PolynomicalRegressionResult solveOrder(final Float64List xs,
      final Float64List ys, final int minOrder, final int maxOrder,
      {final bool printIter = false}) {
    final cache = <PolynomicalRegressionResult>[];
    for (var order = minOrder; order <= maxOrder; order++) {
      cache.add(solve(xs, ys, order, printIter: printIter));
    }
    cache.sort((final a, final b) => a.bic.compareTo(b.bic));
    return cache.first;
  }
}
