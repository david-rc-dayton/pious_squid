import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/observation/observation_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/orbit_determination/orbit_determination_base.dart';
import 'package:pious_squid/src/propagator/propagator_base.dart';

/// Gooding angles-only initial orbit determination.
class GoodingIOD {
  /// Create a new [GoodingIOD] object from three optical observations and
  /// optional gravitational parameter [_mu].
  GoodingIOD(this._o1, this._o2, this._o3, [this._mu = Earth.mu]);
  final double _mu;
  final ObservationOptical _o1;
  final ObservationOptical _o2;
  final ObservationOptical _o3;
  Vector _vObserverPosition1 = Vector.origin3;
  Vector _vObserverPosition2 = Vector.origin3;
  Vector _vObserverPosition3 = Vector.origin3;
  double _r = 0.0;
  double _v = 0.0;
  double _t = 0.0;
  double _r1 = 0.0;
  double _r2 = 0.0;
  double _r3 = 0.0;
  double _rho1 = 0.0;
  double _rho2 = 0.0;
  double _rho3 = 0.0;
  double _d1 = 0.0;
  double _d3 = 0.0;
  double _facFiniteDiff = 0.0;
  final ForceModel _forceModel = ForceModel()..setGravity(1.0);

  Vector? _getPositionOnLoS2(
      final Vector e1,
      final double r01,
      final Vector e3,
      final double r03,
      final double t13,
      final double t12,
      final int nRev,
      final bool posigrade) {
    final p1 = _vObserverPosition1.add(e1.scale(r01));
    _r1 = p1.magnitude();

    final p3 = _vObserverPosition3.add(e3.scale(r03));
    _r3 = p3.magnitude();

    final p13 = p1.cross(p3);

    var th = atan2(p13.magnitude(), p1.dot(p3));

    if (!posigrade) {
      th = twoPi - th;
    }

    final v1 = Float64List(2);
    final exitflag = LambertIOD.solve(_r1, _r3, th, t13, nRev, v1);

    if (exitflag) {
      final pn = p1.cross(p3);
      final pt = pn.cross(p1);

      var rt = pt.magnitude();
      if (!posigrade) {
        rt = -rt;
      }

      final vel1 = p1.scale(v1[0] / _r1).add(pt.scale(v1[1] / rt));

      final p2 = RungeKutta89Propagator(J2000(_o1.epoch, p1, vel1), _forceModel)
          .propagate(_o1.epoch.roll(t12))
          .position;

      return p2;
    }

    return null;
  }

  void _modifyIterate(final Vector lineOfSight1, final Vector lineOfSight3) {
    final r13 = _vObserverPosition3.add(_vObserverPosition1.negate());
    _d1 = r13.dot(lineOfSight1);
    _d3 = r13.dot(lineOfSight3);
    final d2 = lineOfSight1.dot(lineOfSight3);
    final d4 = 1.0 - d2 * d2;
    _rho1 = max((_d1 - _d3 * d2) / d4, 0.0);
    _rho3 = max((_d1 * d2 - _d3) / d4, 0.0);
  }

  void _computeDerivatives(
      final double x,
      final double y,
      final Vector lineOfSight1,
      final Vector lineOfSight3,
      final Vector pin,
      final Vector ein,
      final double t13,
      final double t12,
      final int nrev,
      final bool direction,
      final Float64List fd,
      final Float64List gd) {
    final p = pin.normalize();
    final en = ein.normalize();

    final dx = _facFiniteDiff * x;
    final dy = _facFiniteDiff * y;

    final cm1 = (_getPositionOnLoS2(
            lineOfSight1, x - dx, lineOfSight3, y, t13, t12, nrev, direction)!)
        .subtract(_vObserverPosition2);

    final fm1 = p.dot(cm1);
    final gm1 = en.dot(cm1);

    final cp1 = (_getPositionOnLoS2(
            lineOfSight1, x + dx, lineOfSight3, y, t13, t12, nrev, direction)!)
        .subtract(_vObserverPosition2);

    final fp1 = p.dot(cp1);
    final gp1 = en.dot(cp1);

    final fx = (fp1 - fm1) / (2.0 * dx);
    final gx = (gp1 - gm1) / (2.0 * dx);

    final cm3 = (_getPositionOnLoS2(
            lineOfSight1, x, lineOfSight3, y - dy, t13, t12, nrev, direction)!)
        .subtract(_vObserverPosition2);

    final fm3 = p.dot(cm3);
    final gm3 = en.dot(cm3);

    final cp3 = (_getPositionOnLoS2(
            lineOfSight1, x, lineOfSight3, y + dy, t13, t12, nrev, direction)!)
        .subtract(_vObserverPosition2);

    final fp3 = p.dot(cp3);
    final gp3 = en.dot(cp3);

    final fy = (fp3 - fm3) / (2.0 * dy);
    final gy = (gp3 - gm3) / (2.0 * dy);

    fd[0] = fx;
    fd[1] = fy;
    gd[0] = gx;
    gd[1] = gy;
  }

  /// Attempt to solve a state estimate given a range guess for the first
  /// and last optical observation.
  ///
  /// Adjust the number of revolutions [nRev] and [direction] if the solve
  /// takes more than one revolution or is retrograde.
  J2000 solve(final double r1Init, final double r3Init,
      {final int nRev = 0, final bool direction = true}) {
    final lineOfSight1 = _o1.observation.lineOfSight();
    final lineOfSight2 = _o2.observation.lineOfSight();
    final lineOfSight3 = _o3.observation.lineOfSight();

    _r = max(r1Init, r3Init);
    _v = sqrt(_mu / _r);
    _t = _r / _v;

    _vObserverPosition1 =
        _o1.site.toITRF(_o1.epoch).toJ2000().position.scale(1.0 / _r);
    _vObserverPosition2 =
        _o2.site.toITRF(_o2.epoch).toJ2000().position.scale(1.0 / _r);
    _vObserverPosition3 =
        _o3.site.toITRF(_o3.epoch).toJ2000().position.scale(1.0 / _r);

    final maxiter = 100;
    _solveRangeProblem(
        r1Init / _r,
        r3Init / _r,
        _o3.epoch.difference(_o1.epoch) / _t,
        _o2.epoch.difference(_o1.epoch) / _t,
        nRev,
        direction,
        lineOfSight1,
        lineOfSight2,
        lineOfSight3,
        maxiter);
    final gibbs = GibbsIOD(_mu);
    final p1 = _vObserverPosition1.add(lineOfSight1.scale(_rho1)).scale(_r);
    final p2 = _vObserverPosition2.add(lineOfSight2.scale(_rho2)).scale(_r);
    final p3 = _vObserverPosition3.add(lineOfSight3.scale(_rho3)).scale(_r);
    return gibbs.solve(p1, p2, p3, _o2.epoch, _o3.epoch);
  }

  void _solveRangeProblem(
      final double rho1init,
      final double rho3init,
      final double t13,
      final double t12,
      final int nrev,
      final bool direction,
      final Vector lineOfSight1,
      final Vector lineOfSight2,
      final Vector lineOfSight3,
      final int maxIterations) {
    final arbf = 1e-6;
    final cvtol = 1e-14;

    _rho1 = rho1init;
    _rho3 = rho3init;

    var iter = 0;
    var stoppingCriterion = 10.0 * cvtol;
    while (iter < maxIterations && stoppingCriterion.abs() > cvtol) {
      _facFiniteDiff = arbf;

      final p2 = _getPositionOnLoS2(
          lineOfSight1, _rho1, lineOfSight3, _rho3, t13, t12, nrev, direction);

      if (p2 == null) {
        _modifyIterate(lineOfSight1, lineOfSight3);
      } else {
        _r2 = p2.magnitude();
        final c = p2.subtract(_vObserverPosition2);
        _rho2 = c.magnitude();
        final cr = lineOfSight2.dot(c);

        final u = lineOfSight2.cross(c);
        final p = u.cross(lineOfSight2).normalize();
        final ent = lineOfSight2.cross(p);

        final enr = ent.magnitude();
        if (enr == 0.0) {
          return;
        }

        final en = ent.normalize();

        final fc = p.dot(c);

        final fd = Float64List(2);
        final gd = Float64List(2);
        _computeDerivatives(_rho1, _rho3, lineOfSight1, lineOfSight3, p, en,
            t13, t12, nrev, direction, fd, gd);

        final fr1 = fd[0];
        final fr3 = fd[1];
        final gr1 = gd[0];
        final gr3 = gd[1];
        final detj = fr1 * gr3 - fr3 * gr1;

        _d3 = -gr3 * fc / detj;
        _d1 = gr1 * fc / detj;

        _rho1 += _d3;
        _rho3 += _d1;

        final den = max(cr, _r2);
        stoppingCriterion = fc / den;
      }

      ++iter;
    }

    return;
  }
}
