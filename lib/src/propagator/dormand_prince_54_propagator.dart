import 'dart:typed_data';

import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/force/force_base.dart';
import 'package:pious_squid/src/propagator/runge_kutta_adaptive.dart';

final Float64List _a = Float64List.fromList(
    [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0]);

final List<Float64List> _b = [
  Float64List(0),
  Float64List.fromList([1.0 / 5.0]),
  Float64List.fromList([3.0 / 40.0, 9.0 / 40.0]),
  Float64List.fromList([44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0]),
  Float64List.fromList(
      [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0]),
  Float64List.fromList([
    9017.0 / 3168.0,
    -355.0 / 33.0,
    46732.0 / 5247.0,
    49.0 / 176.0,
    -5103.0 / 18656.0
  ]),
  Float64List.fromList([
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0
  ]),
];

final Float64List _ch = Float64List.fromList([
  35.0 / 384.0,
  0.0,
  500.0 / 1113.0,
  125.0 / 192.0,
  -2187.0 / 6784.0,
  11.0 / 84.0,
  0.0
]);

final Float64List _c = Float64List.fromList([
  5179.0 / 57600.0,
  0.0,
  7571.0 / 16695.0,
  393.0 / 640.0,
  -92097.0 / 339200.0,
  187.0 / 2100.0,
  1.0 / 40.0
]);

/// Dormand-Prince 5(4) adaptive numerical propagator.
class DormandPrince54Propagator extends RungeKuttaAdaptive {
  /// Create a new [DormandPrince54Propagator] from an initial [state] and
  /// an optional [ForceModel].
  DormandPrince54Propagator(final J2000 state, [final ForceModel? forceModel])
      : super(state, forceModel);

  @override
  Float64List get a => _a;

  @override
  List<Float64List> get b => _b;

  @override
  Float64List get ch => _ch;

  @override
  Float64List get c => _c;

  @override
  int get order => 5;
}
