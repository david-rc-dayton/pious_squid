import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Equinoctial element set.
class EquinoctialElements {
  /// Create a new [EquinoctialElements] object from orbital elements.
  EquinoctialElements(
      this.epoch, this.af, this.ag, this.l, this.n, this.chi, this.psi,
      {this.mu = Earth.mu, this.fr = 1});

  /// UTC epoch.
  final EpochUTC epoch;

  /// Af parameter.
  final double af;

  /// Ag parameter.
  final double ag;

  /// L parameter.
  final double l;

  /// N parameter.
  final double n;

  /// Chi parameter.
  final double chi;

  /// Psi parameter.
  final double psi;

  /// Gravitational parameter _(km³/s²)_.
  final double mu;

  /// Posigrade parameter _(-1, 1)_.
  final int fr;

  @override
  String toString() => [
        '[EquinoctialElements]',
        '  Epoch: $epoch',
        '  Af:   $af',
        '  Ag:   $ag',
        '  L:    $l rad',
        '  N:    $n rad/s',
        '  Chi:  $chi',
        '  Psi:  $psi'
      ].join('\n');

  /// Return the orbit semimajor-axis _(km)_.
  double semimajorAxis() => cbrt(mu / (n * n));

  /// Compute the mean motion _(rad/s)_ of this orbit.
  double meanMotion() {
    final a = semimajorAxis();
    return sqrt(mu / (a * a * a));
  }

  /// Compute the period _(seconds)_ of this orbit.
  double period() => twoPi * sqrt(pow(semimajorAxis(), 3) / mu);

  /// Compute the number of revolutions completed per day for this orbit.
  double revsPerDay() => secondsPerDay / period();

  /// Convert this to [ClassicalElements].
  ClassicalElements toClassicalElements() {
    final a = semimajorAxis();
    final e = sqrt(af * af + ag * ag);
    final i =
        pi * ((1.0 - fr) * 0.5) + 2.0 * fr * atan(sqrt(chi * chi + psi * psi));
    final o = atan2(chi, psi);
    final w = atan2(ag, af) - fr * atan2(chi, psi);
    final m = l - fr * o - w;
    final v = newtonM(e, m).nu;
    return ClassicalElements(epoch, a, e, i, o, w, v, mu: mu);
  }

  /// Convert this to inertial position and velocity vectors.
  ({Vector position, Vector velocity}) toPositionVelocity() =>
      toClassicalElements().toPositionVelocity();
}
