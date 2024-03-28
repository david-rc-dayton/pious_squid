import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/covariance/covariance_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Batch least squares orbit determination result.
class BatchLeastSquaresResult {
  /// Create a new [BatchLeastSquaresResult] object, containing the solved
  /// [state], [covariance], and root-mean-squared error [rms].
  BatchLeastSquaresResult(this.state, this.covariance, this.rms);

  /// Solved inertial state vector.
  final J2000 state;

  /// Solved covariance matrix.
  final StateCovariance covariance;

  /// Root-mean-squared error.
  final double rms;
}

/// Batch least squares orbit determination.
class BatchLeastSquaresOD {
  /// Create a new [BatchLeastSquaresOD] object from a list of [Observation]
  /// objects, an [apriori] state estimate, and an optional
  /// spacecraft [forceModel].
  BatchLeastSquaresOD(this._observations, final J2000 apriori,
      {final ForceModel? forceModel,
      final double posStep = 1e-5,
      final double velStep = 1e-5,
      final bool fastDerivatives = false}) {
    _observations
        .sort((final a, final b) => b.epoch.posix.compareTo(a.epoch.posix));
    _start = _observations.first.epoch;
    _propPairs = PropagatorPairs(posStep, velStep);
    _forceModel = forceModel ?? (ForceModel()..setGravity());
    _propagator = RungeKutta89Propagator(apriori, _forceModel);
    _nominal = _propagator.propagate(_start);
    _fastDerivatives = fastDerivatives;
  }

  /// Propagator pair cache, for generating observation Jacobians.
  late final PropagatorPairs _propPairs;

  /// Observations to use in the solve.
  final List<Observation> _observations;

  /// Nominal state propagator.
  late Propagator _propagator;

  /// State estimate during solve.
  late final J2000 _nominal;

  /// Spacecraft force model.
  late final ForceModel _forceModel;

  /// Solve start epoch.
  late final EpochUTC _start;

  /// Use Keplerian logic for derivatives if `true`.
  late final bool _fastDerivatives;

  Propagator _buildPropagator(final Float64List x0, final bool simple) {
    final state = J2000(_nominal.epoch, Vector3D(x0[0], x0[1], x0[2]),
        Vector3D(x0[3], x0[4], x0[5]));
    if (simple) {
      return KeplerPropagator(state.toClassicalElements());
    }
    return RungeKutta89Propagator(state, _forceModel);
  }

  static Float64List _stateToX0(final J2000 state) =>
      state.position.join(state.velocity).toArray();

  void _setPropagatorPairs(final Float64List x0) {
    final pl = _buildPropagator(x0, _fastDerivatives);
    for (var i = 0; i < 6; i++) {
      final step = _propPairs.step(i);
      final xh = x0.sublist(0);
      xh[i] += step;
      final ph = _buildPropagator(xh, _fastDerivatives);
      _propPairs.set(i, ph, pl);
    }
  }

  /// Attempt to solve a state estimate with the given root-mean-squared
  /// delta [tolerance].
  BatchLeastSquaresResult solve(
      {final double tolerance = 1e-6,
      final int maxIter = 250,
      final bool printIter = false}) {
    var breakFlag = false;
    final xNom = _stateToX0(_nominal);
    var weightedRms = double.infinity;
    final atwaMatInit = Matrix(6, 6);
    final atwbMatInit = Matrix(6, 1);
    var atwaMat = atwaMatInit;
    for (var iter = 0; iter < maxIter; iter++) {
      atwaMat = atwaMatInit;
      var atwbMat = atwbMatInit;
      _propagator = _buildPropagator(xNom, false);
      _setPropagatorPairs(xNom);
      var rmsTotal = 0.0;
      var measCount = 0;
      for (final ob in _observations) {
        final noise = ob.noise;
        final aMat = ob.jacobian(_propPairs);
        final aMatTN = aMat.transpose().multiply(noise);
        final bMat = ob.residual(_propagator);
        atwaMat = atwaMat.add(aMatTN.multiply(aMat));
        atwbMat = atwbMat.add(aMatTN.multiply(bMat));
        rmsTotal += bMat.transpose().multiply(noise).multiply(bMat).get(0, 0);
        measCount += noise.rows;
      }
      final newWeightedRms = sqrt(rmsTotal / measCount);
      if (printIter) {
        print('${iter + 1}: rms=$newWeightedRms x=${Vector(xNom)}');
      }
      if (((weightedRms - newWeightedRms) / weightedRms).abs() <= tolerance) {
        breakFlag = true;
      }
      weightedRms = newWeightedRms;
      final atwbVec = atwbMat.transpose().row(0);
      final dX = atwaMat.solve(atwbVec);
      for (var i = 0; i < 6; i++) {
        xNom[i] += dX[i];
      }
      if (breakFlag) {
        break;
      }
    }
    final p = atwaMat.pseudoinverse();
    final covariance = StateCovariance(p, CovarianceFrame.j2000);
    return BatchLeastSquaresResult(
        _buildPropagator(xNom, false).propagate(_start),
        covariance,
        weightedRms);
  }
}
