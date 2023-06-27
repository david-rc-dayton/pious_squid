import 'dart:math';
import 'dart:typed_data';

import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/operations_base.dart';

/// Box-Muller random Gaussian number generator.
class BoxMuller {
  /// Create a new [BoxMuller] object with mean [mu], standard deviation
  /// [sigma], and [seed] number.
  BoxMuller(this.mu, this.sigma, [final int seed = 0]) : rand = Random(seed) {
    _generate();
  }

  /// Mean value.
  final double mu;

  /// Standard deviation.
  final double sigma;

  /// Random number cache.
  final Float64List _cache = Float64List(2);

  /// Cache index.
  int _index = 0;

  /// Uniform random number generator.
  final Random rand;

  /// Refill the cache with random Gaussian numbers.
  void _generate() {
    _index = 0;
    final u1 = rand.nextDouble();
    final u2 = rand.nextDouble();
    final mag = sigma * sqrt(-2.0 * log(u1));
    _cache[0] = mag * cos(twoPi * u2) + mu;
    _cache[1] = mag * sin(twoPi * u2) + mu;
  }

  /// Generate a gaussian number, with mean [mu] and standard
  /// deviation [sigma].
  double nextGauss() {
    if (_index > 1) {
      _generate();
    }
    final result = _cache[_index];
    _index++;
    return result;
  }

  /// Generate a [Vector] of gaussian numbers, with mean [mu] and standard
  /// deviation [sigma].
  Vector gaussVector(final int n) {
    final result = Float64List(n);
    for (var i = 0; i < n; i++) {
      result[i] = nextGauss();
    }
    return Vector(result);
  }
}

/// Random Gaussian number generation with a fixed mean of 0 and standard
/// deviation of 1.
class RandomGaussianSource {
  /// Create a new [RandomGaussianSource] with an optional [seed] value.
  RandomGaussianSource([final int seed = 0])
      : _boxMuller = BoxMuller(0, 1, seed);

  /// Box-Muller number generator.
  final BoxMuller _boxMuller;

  /// Generate a gaussian number, with mean 0 and standard
  /// deviation 1.
  double nextGauss() => _boxMuller.nextGauss();

  /// Generate a [Vector] of gaussian numbers, with mean 0, standard
  /// deviation 1, and length [n].
  Vector gaussVector(final int n) {
    final result = Float64List(n);
    for (var i = 0; i < n; i++) {
      result[i] = nextGauss();
    }
    return Vector(result);
  }

  /// Generate a uniformly distributed random vector, with points on the
  /// surface of a sphere with the provided [radius].
  Vector3D gaussSphere([final double radius = 1.0]) =>
      gaussVector(3).toVector3D(0).normalize().scale(radius);
}
