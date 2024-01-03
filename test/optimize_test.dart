import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

double _rosenbrock(final Float64List xs) =>
    pow(1 - xs[0], 2) + 100 * pow(xs[1] - xs[0] * xs[0], 2) as double;

double _himmelblau(final Float64List xs) =>
    pow(xs[0] * xs[0] + xs[1] - 11, 2) + pow(xs[0] + xs[1] * xs[1] - 7, 2)
        as double;

void main() {
  group('Optimize', () {
    test('DownhillSimplex', () {
      final rosenbrockSimplex =
          DownhillSimplex.generateSimplex(Float64List.fromList([-1, 1]));
      final rosenbrockSolve =
          DownhillSimplex.solveSimplex(_rosenbrock, rosenbrockSimplex);
      expect(rosenbrockSolve[0], closeTo(1, 1e-5));
      expect(rosenbrockSolve[1], closeTo(1, 1e-5));

      final himmelblauSimplex =
          DownhillSimplex.generateSimplex(Float64List.fromList([-1, 1]));
      final himmelblauSolve =
          DownhillSimplex.solveSimplex(_himmelblau, himmelblauSimplex);
      expect(himmelblauSolve[0], closeTo(-2.805118, 1e-5));
      expect(himmelblauSolve[1], closeTo(3.131312, 1e-5));
    });

    test('SimpleLinearRegression', () {
      final xs = <double>[17, 13, 12, 15, 16, 14, 16, 16, 18, 19];
      final ys = <double>[94, 73, 59, 80, 93, 85, 66, 79, 77, 91];
      final f = SimpleLinearRegression(xs, ys);
      expect(f.evaluate(15), closeTo(77.792, 1e-3));
      expect(f.error, closeTo(9.295, 1e-3));
      final f1 = f.filterOutliers(1.0);
      expect(f1.length, equals(5));
    });

    test('PolynomialRegression', () {
      final xs = <double>[-1, 0, 1, 2, 3, 5, 7, 9];
      final ys = <double>[-1, 3, 2, 5, 4, 2, 5, 4];
      final result = PolynomialRegression.solve(
          Float64List.fromList(xs), Float64List.fromList(ys), 4);
      expect(result.coefficients.length, equals(5));
      expect(result.coefficients[0], closeTo(-0.0108, 1e-3));
      expect(result.coefficients[1], closeTo(0.1984, 1e-3));
      expect(result.coefficients[2], closeTo(-1.1341, 1e-3));
      expect(result.coefficients[3], closeTo(2.193, 1e-3));
      expect(result.coefficients[4], closeTo(2.5004, 1e-3));
      expect(result.rss, closeTo(2.6417, 1e-3));
      expect(result.bic, closeTo(23.8608, 1e-3));
      expect(result.evaluate(0), closeTo(2.5003, 1e-3));
    });
  });
}
