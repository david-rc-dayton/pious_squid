import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/covariance/covariance_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/optimize/optimize_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Gauss-Newton orbit determination result.
class GaussNewtonODResult {
  /// Create a new [GaussNewtonODResult] from a solved inertial [state]
  /// _(km,km/s)_, ballistic coefficient _(kg/m²)_, and solar radiation
  /// pressure coefficient _(kg/m²)_.
  GaussNewtonODResult(
      this.state, this.covariance, this.rms, this.bcoeff, this.srpcoeff);

  /// Solved inertial state vector _(km,km/s)_.
  final J2000 state;

  /// Solved inertial state covariance.
  final StateCovariance covariance;

  /// Solve inertial state Root-Mean-Squared error.
  final double rms;

  /// Solved ballistic coefficient _(kg/m²)_.
  final double bcoeff;

  /// Solved solar radiation pressure coefficient _(kg/m²)_.
  final double srpcoeff;

  @override
  String toString() => '$state\n'
      'bcoeff:     ${bcoeff.toStringAsFixed(3)} kg/m²\n'
      'srpcoeff:   ${srpcoeff.toStringAsFixed(3)} kg/m²\n'
      'rms:        ${rms.toStringAsFixed(6)}\n'
      'sigmas:\n'
      '${covariance.sigmas()}\n'
      'covariance:\n'
      '${covariance.matrix}';
}

/// Gauss-Newton differential corrector orbit determination.
class GaussNewtonOD {
  GaussNewtonOD._(); // disable constructor

  /// Convert a parameter vector to an inertial state object.
  static J2000 _vecToState(final Vector x, final EpochUTC epoch) =>
      J2000(epoch, Vector3D(x[0], x[1], x[2]), Vector3D(x[3], x[4], x[5]));

  /// Convert a parameter vector and force model to a propagator object.
  static Propagator _vecToProp(
      final Vector x, final EpochUTC epoch, final ForceModel forceModel) {
    final state = _vecToState(x, epoch);

    /// note: parameter coefficients are inverted during the solve
    final fm = forceModel.clone()
      ..setBCoeff(invertZero(x[6]))
      ..setSrpCoeff(invertZero(x[7]));
    return RungeKutta4Propagator(state, fm);
  }

  /// Convert an inertial state and force model into a parameter vector.
  static Vector _stateToVec(final J2000 state, final ForceModel forceModel) {
    final output = Float64List(8);
    output[0] = state.position.x;
    output[1] = state.position.y;
    output[2] = state.position.z;
    output[3] = state.velocity.x;
    output[4] = state.velocity.y;
    output[5] = state.velocity.z;
    // solve for the coefficient inverse
    output[6] = invertZero(forceModel.getBCoeff());
    output[7] = invertZero(forceModel.getSrpCoeff());
    return Vector(output);
  }

  /// Estimae covariance from a solve vector, and objective/jacobian functions.
  static Matrix _solveToCov(final Vector residuals, final Matrix jacobian) {
    var sigsq = 0.0;
    for (var i = 0; i < residuals.length; i++) {
      sigsq += residuals[i] * residuals[i];
    }
    sigsq /= residuals.length - 8;
    final cov =
        jacobian.transpose().multiply(jacobian).inverseSingular().scale(sigsq);
    return cov.getBlock(0, 0, 6, 6);
  }

  static double _solveToRms(final Vector residuals) {
    final n = residuals.length;
    var rms = 0.0;
    for (var i = 0; i < n; i++) {
      rms += residuals[i] * residuals[i];
    }
    rms /= n;
    return sqrt(rms);
  }

  /// Create a solve objective function for [ObservationState] objects.
  static Vector Function(Vector) _stateObjective(
      final List<ObservationState> obs, final ForceModel forceModel) {
    const m = 6;
    final n = obs.length;
    final epoch = obs.first.epoch;
    return (final Vector x) {
      final propagator = _vecToProp(x, epoch, forceModel);
      final elements = Vector.zero(n * m);
      for (var i = 0; i < n; i++) {
        final expected = obs[i].observation;
        final actual = propagator.propagate(expected.epoch).toITRF();
        elements[m * i + 0] = expected.position.x - actual.position.x;
        elements[m * i + 1] = expected.position.y - actual.position.y;
        elements[m * i + 2] = expected.position.z - actual.position.z;
        elements[m * i + 3] = expected.velocity.x - actual.velocity.x;
        elements[m * i + 4] = expected.velocity.y - actual.velocity.y;
        elements[m * i + 5] = expected.velocity.z - actual.velocity.z;
      }
      return elements;
    };
  }

  /// Create a solve jacobian function for [ObservationState] objects.
  static Matrix Function(Vector) _stateJacobian(
      final List<ObservationState> obs, final ForceModel forceModel) {
    const delta = 1e-5;
    const m = 6;
    final n = obs.length;
    final epoch = obs.first.epoch;
    return (final Vector x) {
      final low = _vecToProp(x, epoch, forceModel);
      final high = <Propagator>[];
      for (var i = 0; i < 8; i++) {
        final xh = x.toArray();
        xh[i] += delta;
        high.add(_vecToProp(Vector(xh), epoch, forceModel));
      }
      final matrix = Matrix(n * m, 8);
      for (var i = 0; i < n; i++) {
        final ob = obs[i];
        final t = ob.epoch;
        final check = ob.observation.posvel();
        final lowItrf = low.propagate(t).toITRF().posvel();
        final dl = check.subtract(lowItrf);
        for (var j = 0; j < 8; j++) {
          final highItrf = high[j].propagate(t).toITRF().posvel();
          final dh = check.subtract(highItrf);
          matrix.set(m * i + 0, j, (dh[0] - dl[0]) / delta);
          matrix.set(m * i + 1, j, (dh[1] - dl[1]) / delta);
          matrix.set(m * i + 2, j, (dh[2] - dl[2]) / delta);
          matrix.set(m * i + 3, j, (dh[3] - dl[3]) / delta);
          matrix.set(m * i + 4, j, (dh[4] - dl[4]) / delta);
          matrix.set(m * i + 5, j, (dh[5] - dl[5]) / delta);
        }
      }
      return matrix;
    };
  }

  /// Create a solve objective function for [ObservationOptical] objects.
  static Vector Function(Vector) _opticalObjective(
      final List<ObservationOptical> obs, final ForceModel forceModel) {
    const m = 2;
    final n = obs.length;
    final epoch = obs.first.epoch;
    return (final Vector x) {
      final propagator = _vecToProp(x, epoch, forceModel);
      final elements = Vector.zero(n * m);
      for (var i = 0; i < n; i++) {
        final expected = obs[i];
        final actual = propagator.propagate(expected.epoch);
        final losExp = expected.observation;
        final losAct = RadecTopocentric.fromStateVectors(actual, expected.site);
        elements[m * i + 0] =
            normalizeAngle(losExp.rightAscension, losAct.rightAscension);
        elements[m * i + 1] =
            normalizeAngle(losExp.declination, losAct.declination);
      }
      return elements;
    };
  }

  /// Create a solve jacobian function for [ObservationOptical] objects.
  static Matrix Function(Vector) _opticalJacobian(
      final List<ObservationOptical> obs, final ForceModel forceModel) {
    const delta = 1e-5;
    const m = 2;
    final n = obs.length;
    final epoch = obs.first.epoch;
    return (final Vector x) {
      final low = _vecToProp(x, epoch, forceModel);
      final high = <Propagator>[];
      for (var i = 0; i < 8; i++) {
        final xh = x.toArray();
        xh[i] += delta;
        high.add(_vecToProp(Vector(xh), epoch, forceModel));
      }
      final matrix = Matrix(n * m, 8);
      for (var i = 0; i < n; i++) {
        final ob = obs[i];
        final expected = ob.observation;
        final lowState = low.propagate(ob.epoch);
        final actualLow = RadecTopocentric.fromStateVectors(lowState, ob.site);
        final dlRa =
            normalizeAngle(expected.rightAscension, actualLow.rightAscension);
        final dlDec =
            normalizeAngle(expected.declination, actualLow.declination);
        for (var j = 0; j < 8; j++) {
          final highState = high[j].propagate(ob.epoch);
          final actualHigh =
              RadecTopocentric.fromStateVectors(highState, ob.site);
          final dhRa = normalizeAngle(
              expected.rightAscension, actualHigh.rightAscension);
          final dhDec =
              normalizeAngle(expected.declination, actualHigh.declination);
          matrix.set(m * i + 0, j, (dhRa - dlRa) / delta);
          matrix.set(m * i + 1, j, (dhDec - dlDec) / delta);
        }
      }
      return matrix;
    };
  }

  /// Create a solve objective function for [ObservationRadar] objects.
  static Vector Function(Vector) _radarObjective(
      final List<ObservationRadar> obs, final ForceModel forceModel) {
    const m = 3;
    final n = obs.length;
    final epoch = obs.first.epoch;
    return (final Vector x) {
      final propagator = _vecToProp(x, epoch, forceModel);
      final elements = Vector.zero(n * m);
      for (var i = 0; i < n; i++) {
        final expected = obs[i];
        final actual = propagator.propagate(expected.epoch);
        final losExp = expected.observation;
        final losAct = Razel.fromStateVectors(actual, expected.site);
        elements[m * i + 0] = losExp.range - losAct.range;
        elements[m * i + 1] = normalizeAngle(losExp.azimuth, losAct.azimuth);
        elements[m * i + 2] =
            normalizeAngle(losExp.elevation, losAct.elevation);
      }
      return elements;
    };
  }

  /// Create a solve jacobian function for [ObservationRadar] objects.
  static Matrix Function(Vector) _radarJacobian(
      final List<ObservationRadar> obs, final ForceModel forceModel) {
    const delta = 1e-5;
    const m = 3;
    final n = obs.length;
    final epoch = obs.first.epoch;
    return (final Vector x) {
      final low = _vecToProp(x, epoch, forceModel);
      final high = <Propagator>[];
      for (var i = 0; i < 8; i++) {
        final xh = x.toArray();
        xh[i] += delta;
        high.add(_vecToProp(Vector(xh), epoch, forceModel));
      }
      final matrix = Matrix(n * m, 8);
      for (var i = 0; i < n; i++) {
        final ob = obs[i];
        final expected = ob.observation;
        final lowState = low.propagate(ob.epoch);
        final actualLow = Razel.fromStateVectors(lowState, ob.site);
        final dlR = expected.range - actualLow.range;
        final dlAz = normalizeAngle(expected.azimuth, actualLow.azimuth);
        final dlEl = normalizeAngle(expected.elevation, actualLow.elevation);
        for (var j = 0; j < 8; j++) {
          final highState = high[j].propagate(ob.epoch);
          final actualHigh = Razel.fromStateVectors(highState, ob.site);
          final dhR = expected.range - actualHigh.range;
          final dhAz = normalizeAngle(expected.azimuth, actualHigh.azimuth);
          final dhEl = normalizeAngle(expected.elevation, actualHigh.elevation);
          matrix.set(m * i + 0, j, (dhR - dlR) / delta);
          matrix.set(m * i + 1, j, (dhAz - dlAz) / delta);
          matrix.set(m * i + 2, j, (dhEl - dlEl) / delta);
        }
      }
      return matrix;
    };
  }

  /// Solve state and parameters for a list of [ObservationState] objects, for
  /// the provided force model.
  static GaussNewtonODResult solveState(
      final List<ObservationState> obs, final ForceModel forceModel,
      {final bool printIter = false}) {
    obs.sort((final a, final b) => b.epoch.compareTo(a.epoch));
    final init = obs.first.observation.toJ2000();
    final guess = _stateToVec(init, forceModel);
    final objective = _stateObjective(obs, forceModel);
    final jacobian = _stateJacobian(obs, forceModel);
    final solve =
        GaussNewton.solve(guess, objective, jacobian, printIter: printIter);
    final solveState = _vecToState(solve, init.epoch);
    final residuals = objective(solve);
    final derivatives = jacobian(solve);
    final covariance = _solveToCov(residuals, derivatives);
    final rms = _solveToRms(residuals);
    return GaussNewtonODResult(
        solveState,
        StateCovariance(covariance, CovarianceFrame.j2000),
        rms,
        invertZero(solve[6]),
        invertZero(solve[7]));
  }

  /// Solve state and parameters for a list of [ObservationOptical] objects,
  /// for the provided force model.
  static GaussNewtonODResult solveOptical(final List<ObservationOptical> obs,
      final J2000 apriori, final ForceModel forceModel,
      {final bool printIter = false}) {
    obs.sort((final a, final b) => b.epoch.compareTo(a.epoch));
    final aprop = RungeKutta89Propagator(apriori, forceModel);
    final init = aprop.propagate(obs.first.epoch);
    final guess = _stateToVec(init, forceModel);
    final objective = _opticalObjective(obs, forceModel);
    final jacobian = _opticalJacobian(obs, forceModel);
    final solve =
        GaussNewton.solve(guess, objective, jacobian, printIter: printIter);
    final solveState = _vecToState(solve, init.epoch);
    final residuals = objective(solve);
    final derivatives = jacobian(solve);
    final covariance = _solveToCov(residuals, derivatives);
    final rms = _solveToRms(residuals);
    return GaussNewtonODResult(
        solveState,
        StateCovariance(covariance, CovarianceFrame.j2000),
        rms,
        invertZero(solve[6]),
        invertZero(solve[7]));
  }

  /// Solve state and parameters for a list of [ObservationRadar] objects,
  /// for the provided force model.
  static GaussNewtonODResult solveRadar(final List<ObservationRadar> obs,
      final J2000 apriori, final ForceModel forceModel,
      {final bool printIter = false}) {
    obs.sort((final a, final b) => b.epoch.compareTo(a.epoch));
    final aprop = RungeKutta89Propagator(apriori, forceModel);
    final init = aprop.propagate(obs.first.epoch);
    final guess = _stateToVec(init, forceModel);
    final objective = _radarObjective(obs, forceModel);
    final jacobian = _radarJacobian(obs, forceModel);
    final solve =
        GaussNewton.solve(guess, objective, jacobian, printIter: printIter);
    final solveState = _vecToState(solve, init.epoch);
    final residuals = objective(solve);
    final derivatives = jacobian(solve);
    final covariance = _solveToCov(residuals, derivatives);
    final rms = _solveToRms(residuals);
    return GaussNewtonODResult(
        solveState,
        StateCovariance(covariance, CovarianceFrame.j2000),
        rms,
        invertZero(solve[6]),
        invertZero(solve[7]));
  }
}
