import 'package:pious_squid/pious_squid.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:test/test.dart';

void main() {
  group('Operations', () {
    test('Quaternion', () {
      final p = Vector.fromList([1, 0, 0]);
      final r = Quaternion(0.0, 0.707, 0.0, 0.707);
      final c = r.toDirectionCosineMatrix();
      final pc = Quaternion.fromDirectionCosineMatrix(c);

      final pRotR = r.rotateVector(p);
      expect(pRotR.x, closeTo(0.0, 1e-5));
      expect(pRotR.y, closeTo(0.0, 1e-5));
      expect(pRotR.z, closeTo(-0.999697, 1e-5));

      final pRotC = c.transpose().multiplyVector(p);
      expect(pRotC.x, closeTo(0.0, 1e-5));
      expect(pRotC.y, closeTo(0.0, 1e-5));
      expect(pRotC.z, closeTo(-0.999697, 1e-5));

      final pRotPc = pc.negate().rotateVector(p);
      expect(pRotPc.x, closeTo(0.0, 1e-3));
      expect(pRotPc.y, closeTo(0.0, 1e-3));
      expect(pRotPc.z, closeTo(-0.999697, 1e-3));

      final q = Quaternion(1, 2, 3, 4).normalize();
      final w = Vector3D(deg2rad, deg2rad, deg2rad);
      final k = q.kinematics(w);
      expect(k.x, closeTo(0.00477978, 1e-7));
      expect(k.y, closeTo(0.00955956, 1e-7));
      expect(k.z, closeTo(0.00477978, 1e-7));
      expect(k.w, closeTo(-0.00955956, 1e-7));
    });
  });
}
