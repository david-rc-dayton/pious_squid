import 'dart:math';
import 'dart:typed_data';

/// Particle swarm optimizer.
class ParticleSwarm {
  ParticleSwarm._(); // disable constructor

  /// Attempt to find the minima of an [objective] function, within the provided
  /// [bounds].
  static Float64List optimize(
      final double Function(Float64List) objective,
      final List<List<double>> bounds,
      final int particles,
      final int iterations,
      {final double wMax = 0.9, // max inertia weight
      final double wMin = 0.4, // min inertia weight
      final double c1 = 2.0, // cognitive coefficient
      final double c2 = 2.0, // social coefficient
      final int seed = 0,
      final bool debug = false}) {
    final dim = bounds.first.length;
    final random = Random(0);
    final positions = <Float64List>[];
    final velocities = <Float64List>[];
    for (var i = 0; i < particles; i++) {
      final p = Float64List(dim);
      final v = Float64List(dim);
      for (var j = 0; j < dim; j++) {
        p[j] =
            random.nextDouble() * (bounds[j][1] - bounds[j][0]) + bounds[j][0];
        v[j] = 0.0;
      }
      positions.add(p);
      velocities.add(v);
    }
    final pBestPositions = List<Float64List>.from(positions);
    final pBestScores = positions.map(objective).toList();
    var gBestScore = pBestScores.reduce(min);
    var gBestPosition =
        List<double>.from(pBestPositions[pBestScores.indexOf(gBestScore)]);

    final step = (wMax - wMin) / max(iterations, 1);
    var w = wMax; // inertia weight

    for (var iteration = 0; iteration < iterations; iteration++) {
      for (var i = 0; i < particles; i++) {
        for (var j = 0; j < dim; j++) {
          // update velocity
          velocities[i][j] = (w * velocities[i][j]) +
              (c1 *
                  random.nextDouble() *
                  (pBestPositions[i][j] - positions[i][j])) +
              (c2 * random.nextDouble() * (gBestPosition[j] - positions[i][j]));

          // update position
          positions[i][j] += velocities[i][j];

          // apply bounds
          positions[i][j] = positions[i][j].clamp(bounds[j][0], bounds[j][1]);
        }

        // update personal best
        final score = objective(positions[i]);
        if (score < pBestScores[i]) {
          pBestScores[i] = score;
          pBestPositions[i] = positions[i].sublist(0);
        }
      }

      // update global best
      final bestScoreIdx = pBestScores.indexOf(pBestScores.reduce(min));
      if (pBestScores[bestScoreIdx] < gBestScore) {
        gBestScore = pBestScores[bestScoreIdx];
        gBestPosition = List.from(pBestPositions[bestScoreIdx]);
      }

      // reduce inertia weight
      w -= step;

      if (debug) {
        print('${iteration + 1}/$iterations, x=$gBestPosition fx=$gBestScore');
      }
    }

    return Float64List.fromList(gBestPosition);
  }
}
