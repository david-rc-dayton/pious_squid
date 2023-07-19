import 'dart:math';

/// Simple linear regression _(y = mx + b)_.
class SimpleLinearRegression {
  /// Create a new [SimpleLinearRegression] object from lists of x and y
  /// values.
  SimpleLinearRegression(this.xs, this.ys) {
    update();
  }

  /// X-axis values
  final List<double> xs;

  /// Y-axis values
  final List<double> ys;

  /// Line slope
  double _slope = 0.0;

  /// Y-axis intercept
  double _intercept = 0.0;

  /// Line slope
  double get slope => _slope;

  /// Y-axis intercept
  double get intercept => _intercept;

  /// Update the linear fit with this object's current [xs] and [ys] values.
  void update() {
    final n = min(xs.length, ys.length);
    var xMu = 0.0;
    var yMu = 0.0;
    for (var i = 0; i < n; i++) {
      xMu += xs[i];
      yMu += ys[i];
    }
    xMu /= n;
    yMu /= n;
    var pa = 0.0;
    var xSig = 0.0;
    var ySig = 0.0;
    for (var i = 0; i < n; i++) {
      final xd = xs[i] - xMu;
      final yd = ys[i] - yMu;
      pa += xd * yd;
      xSig += xd * xd;
      ySig += yd * yd;
    }
    final p = pa / (sqrt(xSig) * sqrt(ySig));
    xSig = sqrt(xSig / (n - 1));
    ySig = sqrt(ySig / (n - 1));
    _slope = p * (ySig / xSig);
    _intercept = yMu - (_slope * xMu);
  }

  /// Evaluate this linear fit for y, given an [x] value.
  double evaluate(final double x) => _slope * x + _intercept;
}
