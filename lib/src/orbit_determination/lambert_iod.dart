import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Lambert two-position and time initial orbit determination.
class LambertIOD {
  /// Create a new [LambertIOD] object with optional gravitational
  /// parameter [mu].
  LambertIOD([this.mu = Earth.mu]);

  /// Gravitational parameter _(km²/s³)_.
  final double mu;

  /// Try to guess the short path argument given an [interceptor] and
  /// [target] state.
  static bool useShortPath(final J2000 interceptor, final J2000 target) {
    final transN = interceptor.position.cross(target.position);
    final h = interceptor.position.cross(interceptor.velocity);
    return h.dot(transN) >= 0;
  }

  static double _timeOfFlight(final double x, final int longway, final int mrev,
      final double minSma, final double speri, final double chord) {
    final a = minSma / (1.0 - x * x);
    double tof;
    if (x.abs() < 1) {
      final beta = longway * 2.0 * asin(sqrt((speri - chord) / (2.0 * a)));
      final alpha = 2.0 * acos(x);
      tof = a *
          sqrt(a) *
          (alpha - sin(alpha) - (beta - sin(beta)) + 2.0 * pi * mrev);
    } else {
      final alpha = 2.0 * acosh(x);
      final beta = longway * 2.0 * asinh(sqrt((speri - chord) / (-2.0 * a)));
      tof = -a * sqrt(-a) * (sinh(alpha) - alpha - (sinh(beta) - beta));
    }
    return tof;
  }

  /// Attempt to solve output velocity [v1] _(km/s)_ given radii [r1] and
  /// [r2] _(canonical)_, sweep angle [dth] _(rad)_, time of flight [tau]
  /// _(canonical)_, and number of revolutions _(mRev)_.
  static bool solve(final double r1, final double r2, final double dth,
      final double tau, final int mRev, final Float64List v1) {
    final leftBranch = dth < pi;
    var longway = 1;
    if (dth > pi) {
      longway = -1;
    }

    final m = mRev.abs();
    final rtof = tau.abs();
    final theta = dth;

    final chord = sqrt(r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * cos(theta));
    final speri = 0.5 * (r1 + r2 + chord);

    final minSma = 0.5 * speri;
    final lambda = longway * sqrt(1.0 - chord / speri);
    final logt = log(rtof);

    double in1;
    double in2;
    double x1;
    double x2;
    if (m == 0) {
      in1 = -0.6523333;
      in2 = 0.6523333;
      x1 = log(1.0 + in1);
      x2 = log(1.0 + in2);
    } else {
      if (!leftBranch) {
        in1 = -0.523334;
        in2 = -0.223334;
      } else {
        in1 = 0.723334;
        in2 = 0.523334;
      }
      x1 = tan(in1 * halfPi);
      x2 = tan(in2 * halfPi);
    }
    final tof1 = _timeOfFlight(in1, longway, m, minSma, speri, chord);
    final tof2 = _timeOfFlight(in2, longway, m, minSma, speri, chord);

    double y1;
    double y2;
    if (m == 0) {
      y1 = log(tof1) - logt;
      y2 = log(tof2) - logt;
    } else {
      y1 = tof1 - rtof;
      y2 = tof2 - rtof;
    }
    var err = 1e20;
    var iterations = 0;
    final tol = 1e-13;
    final maxiter = 50;
    var xnew = 0.0;
    while (err > tol && iterations < maxiter) {
      xnew = (x1 * y2 - y1 * x2) / (y2 - y1);
      double xt;
      if (m == 0) {
        xt = exp(xnew) - 1.0;
      } else {
        xt = atan(xnew) * 2.0 / pi;
      }

      final tof = _timeOfFlight(xt, longway, m, minSma, speri, chord);

      double ynew;
      if (m == 0) {
        ynew = log(tof) - logt;
      } else {
        ynew = tof - rtof;
      }

      x1 = x2;
      x2 = xnew;
      y1 = y2;
      y2 = ynew;

      err = (x1 - xnew).abs();
      ++iterations;
    }

    if (err > tol) {
      return false;
    }

    double x;
    if (m == 0) {
      x = exp(xnew) - 1.0;
    } else {
      x = atan(xnew) * 2.0 / pi;
    }

    final sma = minSma / (1.0 - x * x);

    double eta;
    if (x < 1) {
      final alfa = 2.0 * acos(x);
      final beta = longway * 2.0 * asin(sqrt((speri - chord) / (2.0 * sma)));
      final psi = (alfa - beta) / 2.0;
      final sinPsi = sin(psi);
      final etaSq = 2.0 * sma * sinPsi * sinPsi / speri;
      eta = sqrt(etaSq);
    } else {
      final gamma = 2.0 * acosh(x);
      final delta = longway * 2.0 * asinh(sqrt((chord - speri) / (2.0 * sma)));
      final psi = (gamma - delta) / 2.0;
      final sinhPsi = sinh(psi);
      final etaSq = -2.0 * sma * sinhPsi * sinhPsi / speri;
      eta = sqrt(etaSq);
    }

    final vr1 = 1.0 /
        eta *
        sqrt(1.0 / minSma) *
        (2.0 * lambda * minSma / r1 - (lambda + x * eta));
    final vt1 = 1.0 / eta * sqrt(1.0 / minSma) * sqrt(r2 / r1) * sin(dth / 2.0);
    v1[0] = vr1;
    v1[1] = vt1;

    return true;
  }

  /// Estimate a state vector for inertial position [p1] _(km)_ given the
  /// two epoch and positions.
  J2000? estimate(
      final Vector p1, final Vector p2, final EpochUTC t1, final EpochUTC t2,
      {final bool posigrade = true, final int nRev = 0}) {
    final r1 = p1.magnitude();
    final r2 = p2.magnitude();
    final tof = t2.difference(t1);

    final r = max(r1, r2);
    final v = sqrt(mu / r);
    final t = r / v;

    var dth = p1.angle(p2);
    if (!posigrade) {
      dth = twoPi - dth;
    }

    final vDep = Float64List(2);
    final exitFlag = solve(r1 / r, r2 / r, dth, tof / t, nRev, vDep);

    if (exitFlag) {
      final pn = p1.cross(p2);
      final pt = pn.cross(p1);
      var rt = pt.magnitude();
      if (!posigrade) {
        rt = -rt;
      }
      final vel1 = p1.scale(v * vDep[0] / r1).add(pt.scale(v * vDep[1] / rt));
      return J2000(t1, p1, vel1);
    }
    return null;
  }
}
