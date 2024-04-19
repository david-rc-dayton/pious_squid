import 'package:pious_squid/pious_squid.dart';
import 'package:test/test.dart';

void main() {
  group('Functions', () {
    test('mahalanobisDistance', () {
      final z = Vector.fromList([5, 3]);
      final zHat = Vector.fromList([4, 4]);
      final s = Matrix.fromList([
        [2, 0.5],
        [0.5, 1]
      ]);
      final d = mahalanobisDistance(z, zHat, s);
      expect(d, closeTo(1.5118, 1e-3));
    });
  });
}
