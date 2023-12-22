import 'dart:math';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/covariance/covariance_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/epoch.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Sigma point covariance sample.
class CovarianceSample {
  /// Create a new [CovarianceSample] object from an inertial state, covariance
  /// and optional force models for the origin state and samples.
  ///
  /// Two-body physics will be used if a force model is not provided.
  CovarianceSample(final J2000 state, final StateCovariance covariance,
      {ForceModel? originForceModel, ForceModel? sampleForceModel}) {
    originForceModel ??= (ForceModel()..setGravity());
    sampleForceModel ??= (ForceModel()..setGravity());
    _origin = RungeKutta89Propagator(state, originForceModel);

    // build sigma points
    final s = covariance.matrix.cholesky();
    final sqrt6 = sqrt(6.0);
    for (var i = 0; i < 6; i++) {
      for (var j = 0; j < 6; j++) {
        s[i][j] *= sqrt6;
      }
    }

    final sigmapts = Matrix.zero(6, 12);
    for (var i = 0; i < 6; i++) {
      final jj = (i - 1) * 2 + 2;
      for (var j = 0; j < 3; j++) {
        sigmapts[j][jj] = s[j][i];
        sigmapts[j + 3][jj] = s[j + 3][i];

        sigmapts[j][jj + 1] = -s[j][i];
        sigmapts[j + 3][jj + 1] = -s[j + 3][i];
      }
    }

    // build propagators from points
    for (var i = 0; i < 12; i++) {
      final sampleR = Vector3D(sigmapts[0][i], sigmapts[1][i], sigmapts[2][i]);
      final sampleV = Vector3D(sigmapts[3][i], sigmapts[4][i], sigmapts[5][i]);
      if (covariance.frame == CovarianceFrame.eci) {
        final sample = J2000(state.epoch, state.position.add(sampleR),
            state.velocity.add(sampleV));
        _samples.add(RungeKutta89Propagator(sample, sampleForceModel));
      } else if (covariance.frame == CovarianceFrame.ecef) {
        final itrf = state.toITRF();
        final sample = ITRF(state.epoch, itrf.position.add(sampleR),
                itrf.velocity.add(sampleV))
            .toJ2000();
        _samples.add(RungeKutta89Propagator(sample, sampleForceModel));
      } else if (covariance.frame == CovarianceFrame.ric) {
        final sample = RIC(sampleR, sampleV).toJ2000(state);
        _samples.add(RungeKutta89Propagator(sample, sampleForceModel));
      }
    }
  }

  late RungeKutta89Propagator _origin;
  final List<RungeKutta89Propagator> _samples = [];
  final Matrix _pts = Matrix.zero(6, 12);

  /// Current covariance sample epoch.
  Epoch get epoch => _origin.state.epoch;

  /// Current covariance sample origin state.
  J2000 get state => _origin.state;

  /// Rebuild covariance from sigma points.
  Matrix _rebuildCovariance(final Matrix pts) {
    const c = 1.0 / 12.0;
    final yu = Vector.zero(6);
    final y = Matrix.zero(6, 12);

    for (var i = 0; i < 12; i++) {
      for (var j = 0; j < 6; j++) {
        yu[j] = yu[j] + pts[j][i];
      }
    }

    for (var j = 0; j < 6; j++) {
      yu[j] = yu[j] * c;
    }

    for (var i = 0; i < 12; i++) {
      for (var j = 0; j < 6; j++) {
        y[j][i] = pts[j][i] - yu[j];
      }
    }

    final yt = y.transpose();
    final tmp = y.multiply(yt);
    return tmp.scale(c);
  }

  /// Propagate covariance to a new epoch.
  void propagate(final EpochUTC epoch) {
    _origin.propagate(epoch);
    for (final sample in _samples) {
      sample.propagate(epoch);
    }
  }

  /// Apply a maneuver to this covariance.
  void maneuver(final Thrust maneuver) {
    _origin.maneuver(maneuver);
    for (final sample in _samples) {
      sample.maneuver(maneuver);
    }
  }

  /// Desample covariance in J2000 frame.
  StateCovariance desampleJ2000() {
    for (var i = 0; i < 12; i++) {
      final state = _samples[i].state;
      _pts[0][i] = state.position[0];
      _pts[1][i] = state.position[1];
      _pts[2][i] = state.position[2];
      _pts[3][i] = state.velocity[0];
      _pts[4][i] = state.velocity[1];
      _pts[5][i] = state.velocity[2];
    }
    final matrix = _rebuildCovariance(_pts);
    return StateCovariance(matrix, CovarianceFrame.eci);
  }

  /// Desample covariance in RIC frame.
  StateCovariance desampleRIC() {
    final rot = RelativeState.createMatrix(
        _origin.state.position, _origin.state.velocity);
    for (var i = 0; i < 12; i++) {
      final state = RIC.fromJ2000Matrix(_samples[i].state, _origin.state, rot);
      _pts[0][i] = state.position[0];
      _pts[1][i] = state.position[1];
      _pts[2][i] = state.position[2];
      _pts[3][i] = state.velocity[0];
      _pts[4][i] = state.velocity[1];
      _pts[5][i] = state.velocity[2];
    }
    final matrix = _rebuildCovariance(_pts);
    return StateCovariance(matrix, CovarianceFrame.ric);
  }

  /// Desample covariance in ECEF frame.
  StateCovariance desampleECEF() {
    for (var i = 0; i < 12; i++) {
      final state = _samples[i].state.toITRF();
      _pts[0][i] = state.position[0];
      _pts[1][i] = state.position[1];
      _pts[2][i] = state.position[2];
      _pts[3][i] = state.velocity[0];
      _pts[4][i] = state.velocity[1];
      _pts[5][i] = state.velocity[2];
    }
    final matrix = _rebuildCovariance(_pts);
    return StateCovariance(matrix, CovarianceFrame.ecef);
  }
}
