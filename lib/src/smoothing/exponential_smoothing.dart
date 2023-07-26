import 'dart:math';

import 'package:pious_squid/src/time/time_base.dart';

/// Exponential data smoothing methods.
class ExponentialSmoothing {
  ExponentialSmoothing._(); // disable constructor

  /// Perform exponential smoothing on data set [xs] using the provided
  /// data smoothing factor [alpha] _(0.0 <= alpha <= 1.0)_.
  static List<double> smooth(final List<double> xs, final double alpha) {
    final ss = <double>[];
    for (var i = 0; i < xs.length; i++) {
      if (i == 0) {
        ss.add(xs[0]);
        continue;
      }
      ss.add(alpha * xs[i] + (1 - alpha) * ss[i - 1]);
    }
    return ss;
  }

  /// Perform exponential smoothing on data set [xs] using the provided data
  /// smoothing factor [alpha] _(0.0 <= alpha <= 1.0)_ and trend smoothing
  /// factor [beta] _(0.0 <= beta <= 1.0)_.
  static List<double> smoothDouble(
      final List<double> xs, final double alpha, final double beta) {
    final bs = <double>[];
    final ss = <double>[];
    for (var i = 0; i < xs.length; i++) {
      if (i == 0) {
        ss.add(xs[0]);
        bs.add(xs[1] - xs[0]);
        continue;
      }
      ss.add(alpha * xs[i] + (1 - alpha) * (ss[i - 1] + bs[i - 1]));
      bs.add(beta * (ss[i] - ss[i - 1]) + (1 - beta) * bs[i - 1]);
    }
    return ss;
  }

  /// Perform exponential smoothing on time series data set [xs] correlated
  /// with the [epochs] array using the provided [timeConstant]
  /// value _(seconds)_.
  ///
  /// Note: the [timeConstant] is the amount of time for the smoothed response
  /// of a unit step function to reach ~63.2% of the original signal.
  static List<double> smoothTime(final List<EpochUTC> epochs,
      final List<double> xs, final double timeConstant) {
    final ts = epochs.map((final e) => e.posix).toList();
    final ss = <double>[];
    for (var i = 0; i < xs.length; i++) {
      if (i == 0) {
        ss.add(xs[0]);
        continue;
      }
      final a = 1 - exp(-(ts[i] - ts[i - 1]) / timeConstant);
      ss.add(a * xs[i] + (1 - a) * ss[i - 1]);
    }
    return ss;
  }
}
