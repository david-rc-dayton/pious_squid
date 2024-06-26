import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/matrix.dart';
import 'package:pious_squid/src/operations/vector.dart';

/// Class containing Euler angles.
class EulerAngles {
  /// Create a new [EulerAngles] object from [roll], [pitch], and [yaw]
  /// angles _(rad)_.
  const EulerAngles(this.roll, this.pitch, this.yaw);

  /// Create a new [EulerAngles] object from roll, pitch, and yaw
  /// angles _(deg)_.
  factory EulerAngles.fromDegrees(
          final double rollDeg, final double pitchDeg, final double yawDeg) =>
      EulerAngles(rollDeg * deg2rad, pitchDeg * deg2rad, yawDeg * deg2rad);

  /// Create a new [EulerAngles] object from 3-2-1 ordered direction cosine
  /// matrix [c].
  factory EulerAngles.fromDcm321(final Matrix c) => EulerAngles(
      atan(c.get(1, 2) / c.get(2, 2)),
      -asin(c.get(0, 2)),
      atan(c.get(0, 1) / c.get(0, 0)));

  /// Roll component _(rad)_.
  final double roll;

  /// Pitch component _(rad)_.
  final double pitch;

  /// Yaw component _(rad)_.
  final double yaw;

  /// Roll component _(deg)_.
  double get rollDegrees => roll * rad2deg;

  /// Pitch component _(deg)_.
  double get pitchDegrees => pitch * rad2deg;

  /// Yaw component _(deg)_.
  double get yawDegrees => yaw * rad2deg;

  /// Roll component alias _(rad)_.
  double get phi => roll;

  /// Pitch component alias _(rad)_.
  double get theta => pitch;

  /// Yaw component alias _(rad)_.
  double get psi => yaw;

  /// Roll component alias _(deg)_.
  double get phiDegrees => phi * rad2deg;

  /// Pitch component alias _(deg)_.
  double get thetaDegrees => theta * rad2deg;

  /// Yaw component alias _(deg)_.
  double get psiDegrees => psi * rad2deg;

  @override
  String toString({final int precision = 6}) {
    final rollStr = rollDegrees.toStringAsFixed(precision);
    final pitchStr = pitchDegrees.toStringAsFixed(precision);
    final yawStr = yawDegrees.toStringAsFixed(precision);
    return 'Euler(roll: $rollStr°, pitch: $pitchStr°, yaw: $yawStr°)';
  }

  /// Convert this to a 3-2-1 ordered direction cosine [Matrix].
  Matrix dcm321() {
    final sPhi = sin(phi);
    final cPhi = cos(phi);
    final sTheta = sin(theta);
    final cTheta = cos(theta);
    final sPsi = sin(psi);
    final cPsi = cos(psi);
    final output = Float64List(9);
    output[0] = cTheta * cPsi;
    output[1] = cTheta * sPsi;
    output[2] = -sTheta;
    output[3] = sPhi * sTheta * cPsi - cPhi * sPsi;
    output[4] = sPhi * sTheta * sPsi + cPhi * cPsi;
    output[5] = sPhi * cTheta;
    output[6] = cPhi * sTheta * cPsi + sPhi * sPsi;
    output[7] = cPhi * sTheta * sPsi - sPhi * cPsi;
    output[8] = cPhi * cTheta;
    return Matrix(3, 3, output);
  }

  /// Perform a 3-2-1 ordered rotation on provided vector [v].
  Vector3D rotateVector321(final Vector3D v) => dcm321().multiplyVector3D(v);
}
