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
    });
  });
}
